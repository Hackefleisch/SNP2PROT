"""A grid that keeps its numbers and throws its models away cannot be asked anything later."""

from __future__ import annotations

import numpy as np
import pytest

from snp2prot.config import checkpoint_file

torch = pytest.importorskip("torch")

from snp2prot import training  # noqa: E402
from snp2prot.data.matrix import KmerMatrix  # noqa: E402
from snp2prot.embeddings import DomainEmbeddings  # noqa: E402

CONFIG = {
    "model": {
        "width": 8,
        "dna": {"channels": 4, "layers": 1},
        "protein": {"hidden": 0, "dropout": 0.0},
        "temperature": 0.07,
        "learn_temperature": True,
        "max_logit_scale": 100.0,
    },
    "training": {
        "steps": 6,
        "batch_domains": 2,
        "learning_rate": 1e-3,
        "weight_decay": 0.01,
        "warmup_steps": 1,
        "eval_every": 2,
        "patience": 0,
    },
}


@pytest.fixture
def trainer():
    rng = np.random.default_rng(0)
    kmers = np.array(["ACGTACGT", "TTTTAAAA", "GGGGCCCC", "ACACACAC"])
    label = np.array([[1, 0, 0, -1], [0, 1, 0, 0], [1, 1, 0, -1], [0, 0, 0, 0]], dtype=np.int8)
    matrix = KmerMatrix(
        domains=np.array(["AAAA", "CCCC", "DDDD", "EEEE"]),
        kmers=kmers,
        escore=rng.random((4, 4)).astype(np.float32),
        label=label,
    )
    vectors = DomainEmbeddings("A1", "test", matrix.domains, rng.random((4, 6)).astype(np.float32))
    return training.Trainer(matrix, vectors, CONFIG, device=torch.device("cpu"), seed=0)


def test_a_saved_checkpoint_reloads_into_a_model_that_scores_identically(trainer, tmp_path):
    result = trainer.train(np.array([0, 1, 2]), np.array([3]), seed=0)
    before = trainer.predict(np.array([0, 1, 2, 3]))

    path = trainer.save(tmp_path / "fold.pt", result, {"arm": "A1", "fold": "all"})
    model, meta = training.load_checkpoint(path)

    table = model.dna_table(trainer.tokens.cpu())
    after = model.score(trainer.proteins.cpu(), table)[:, :-1].detach().numpy()
    np.testing.assert_allclose(before, after, rtol=1e-5, atol=1e-6)
    assert meta["arm"] == "A1" and meta["fold"] == "all"
    assert "state_dict" not in meta and "final_state_dict" not in meta


def test_both_the_selected_and_the_budget_model_are_kept(trainer, tmp_path):
    """A run produces two models and which is better is a per-fold empirical question, so the
    checkpoint carries both and `run_fold` scores both."""
    result = trainer.train(np.array([0, 1, 2]), np.array([3]), seed=0)
    assert result.best_state and result.final_state
    path = trainer.save(tmp_path / "f.pt", result, {})

    best, meta = training.load_checkpoint(path, which="best")
    final, _ = training.load_checkpoint(path, which="final")
    assert meta["best_step"] <= meta["final_step"] == result.steps_run

    # Each loaded model must match the state it was saved from, whether or not the two coincide.
    for model, state in ((best, result.best_state), (final, result.final_state)):
        for name, tensor in model.state_dict().items():
            torch.testing.assert_close(tensor, state[name], rtol=0, atol=0)


def test_a_validation_slice_with_no_positives_reports_the_budget_not_step_zero(trainer):
    """Every eval returns nan, so nothing is ever selected. The model at the budget is then the
    selection, and `best_step` must say so instead of claiming step 0."""
    result = trainer.train(np.array([0, 1, 2]), np.array([3]), seed=0)
    assert all(np.isnan(h["validation_aupr"]) for h in result.history)
    assert result.best_step == result.steps_run == CONFIG["training"]["steps"]


def test_an_unknown_checkpoint_selector_is_an_error(trainer, tmp_path):
    result = trainer.train(np.array([0, 1, 2]), np.array([3]), seed=0)
    path = trainer.save(tmp_path / "f.pt", result, {})
    with pytest.raises(ValueError, match="which must be one of"):
        training.load_checkpoint(path, which="latest")


def test_the_checkpoint_carries_the_split_and_the_code_it_came_from(trainer, tmp_path):
    """The same rule `tracking.start_run` enforces: weights alone cannot say what they never saw."""
    meta = {"digests": {"test": "abc123"}, "code": {"commit": "9390e41", "dirty": "true"}}
    result = trainer.train(np.array([0, 1, 2]), np.array([3]), seed=0)
    _, read = training.load_checkpoint(trainer.save(tmp_path / "f.pt", result, meta))
    assert read["digests"]["test"] == "abc123"
    assert read["code"]["commit"] == "9390e41"
    assert read["model_config"]["width"] == 8
    assert read["parameter_counts"]["total"] > 0


def test_a_fold_name_with_a_colon_becomes_a_usable_filename():
    assert checkpoint_file("A1", "P3", "half:draw-0").name == "A1_P3_half-draw-0.pt"
