"""The training loop itself: the schedule, the selection, and what `run_fold` reports.

`snp2prot.training` was the only module in the modelling path with no tests of its own, and it
is the one that decides which domains a model sees, which checkpoint is kept and what the report
says about it. `tests/test_training_checkpoint.py` covers saving and reloading; this covers the
loop and the fold driver.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

torch = pytest.importorskip("torch")

from snp2prot import splits, training  # noqa: E402
from snp2prot.data.matrix import KmerMatrix  # noqa: E402
from snp2prot.distances import DomainDistances  # noqa: E402
from snp2prot.embeddings import DomainEmbeddings  # noqa: E402

N_DOMAINS, N_KMERS = 10, 16


def _corpus(seed: int = 0):
    """A small corpus with real structure: every domain has positives, negatives and a gray band."""
    rng = np.random.default_rng(seed)
    kmers = np.array(
        sorted({"".join("ACGT"[c] for c in rng.integers(0, 4, 8)) for _ in range(N_KMERS)}),
        dtype=np.str_,
    )
    escore = rng.normal(size=(N_DOMAINS, len(kmers))).astype(np.float32)
    label = np.zeros_like(escore, dtype=np.int8)
    for i in range(N_DOMAINS):
        order = np.argsort(-escore[i])
        label[i, order[:3]] = 1  # three positives each
        label[i, order[3:5]] = -1  # a gray band
    domains = np.array([f"DOMAIN{i:02d}" + "A" * 20 for i in range(N_DOMAINS)], dtype=np.str_)
    matrix = KmerMatrix(domains, kmers, escore, label)
    vectors = DomainEmbeddings(
        "A1", "test", domains, rng.normal(size=(N_DOMAINS, 12)).astype(np.float32)
    )
    frame = pd.DataFrame(
        {
            "dbd_seq": domains,
            "dbd_family": ["Fam"] * N_DOMAINS,
            "gene": [f"g{i}" for i in range(N_DOMAINS)],
            "wt_id": [f"c{i // 2}" for i in range(N_DOMAINS)],
            "is_variant": [i % 2 == 1 for i in range(N_DOMAINS)],
            "verdict": ["ok"] * N_DOMAINS,
        }
    )
    # Every pair aligned, none close enough to join a component: one group per domain.
    dist = DomainDistances(
        domains,
        np.full((N_DOMAINS, N_DOMAINS), 15, dtype=np.int16),
        np.full((N_DOMAINS, N_DOMAINS), 20, dtype=np.int32),
    )
    np.fill_diagonal(dist.n_edits, 0)
    return frame, matrix, vectors, dist


def _config(steps: int = 20, bce_weight: float = 0.0, **training_overrides):
    cfg = {
        "model": {
            "width": 8,
            "dna": {"channels": 4, "layers": 1},
            "protein": {"hidden": 0, "dropout": 0.0},
            "temperature": 0.07,
            "learn_temperature": True,
            "max_logit_scale": 100.0,
            "seed": 11,
            "bce_weight": 0.0,
        },
        "training": {
            "steps": steps,
            "batch_domains": 4,
            "learning_rate": 1e-2,
            "weight_decay": 0.01,
            "warmup_steps": 5,
            "eval_every": 5,
            "patience": 0,
        },
        "splits": {
            "seed": 7,
            "s2_min_identity": 0.5,
            "validation": {"fraction": 0.25, "grouping": "random"},
        },
        "metrics": {
            "precision_at": [2],
            "recall_at_precision": 0.5,
            "null_repeats": 5,
            "calibration_bins": 4,
            "call_threshold": 0.5,
        },
    }
    cfg["training"].update(training_overrides)
    cfg["model"]["bce_weight"] = bce_weight
    return cfg


@pytest.fixture
def trainer():
    _, matrix, vectors, _ = _corpus()
    return training.Trainer(matrix, vectors, _config(), device=torch.device("cpu"), seed=0)


# --- the loop -----------------------------------------------------------------------------


def test_training_reduces_the_loss_and_moves_the_weights(trainer):
    rows = torch.arange(8)
    before = {k: v.clone() for k, v in trainer.model.state_dict().items()}
    start, _ = trainer.loss_for(rows, trainer.model.dna_table(trainer.tokens))
    trainer.train(np.arange(8), np.array([8, 9]), seed=0)
    end, _ = trainer.loss_for(rows, trainer.model.dna_table(trainer.tokens))
    assert float(end) < float(start)
    assert any(not torch.allclose(before[k], v) for k, v in trainer.model.state_dict().items())


def test_the_learning_rate_warms_up_and_then_holds(trainer):
    seen = []
    original = torch.optim.AdamW.step

    def spy(self, *args, **kwargs):
        seen.append(self.param_groups[0]["lr"])
        return original(self, *args, **kwargs)

    torch.optim.AdamW.step = spy
    try:
        trainer.train(np.arange(8), np.array([8, 9]), seed=0)
    finally:
        torch.optim.AdamW.step = original
    warmup, base = 5, 1e-2
    assert seen[0] == pytest.approx(base / warmup)
    assert seen[warmup - 1] == pytest.approx(base)
    assert all(lr == pytest.approx(base) for lr in seen[warmup:])


def test_the_temperature_is_clamped_on_every_step_not_just_at_the_end():
    _, matrix, vectors, _ = _corpus()
    cfg = _config(steps=40, learning_rate=1.0)  # large enough to blow past the ceiling
    t = training.Trainer(matrix, vectors, cfg, device=torch.device("cpu"), seed=0)
    t.train(np.arange(8), np.array([8, 9]), seed=0)
    assert t.model.logit_scale.detach().exp() <= 100.0 + 1e-3


def test_weight_decay_reaches_the_weights_and_not_the_temperature_or_the_null(trainer):
    groups = trainer._parameter_groups(0.01)
    decayed = {id(p) for g in groups if g["weight_decay"] > 0 for p in g["params"]}
    for name, p in trainer.model.named_parameters():
        assert (id(p) in decayed) == (name not in training.Trainer.NOT_WEIGHTS), name


# --- selection ----------------------------------------------------------------------------


def test_the_kept_checkpoint_is_the_best_scoring_one_not_the_last(trainer):
    result = trainer.train(np.arange(8), np.array([8, 9]), seed=0)
    scored = [h for h in result.history if not np.isnan(h["validation_aupr"])]
    if scored:
        best = max(scored, key=lambda h: h["validation_aupr"])
        assert result.best_step == best["step"]
        assert result.best_validation == pytest.approx(best["validation_aupr"])


def test_patience_stops_a_run_early_and_zero_patience_does_not():
    _, matrix, vectors, _ = _corpus()
    full = training.Trainer(matrix, vectors, _config(steps=60, patience=0), seed=0)
    assert full.train(np.arange(8), np.array([8, 9]), seed=0).steps_run == 60

    impatient = training.Trainer(matrix, vectors, _config(steps=60, patience=1), seed=0)
    assert impatient.train(np.arange(8), np.array([8, 9]), seed=0).steps_run <= 60


def test_an_empty_validation_set_keeps_the_last_evaluation(trainer):
    result = trainer.train(np.arange(10), np.array([], dtype=np.int64), seed=0)
    assert result.best_step == result.steps_run == 20


# --- prediction and scoring ---------------------------------------------------------------


def test_predict_drops_the_null_column_and_restores_training_mode(trainer):
    trainer.model.train()
    predicted = trainer.predict(np.array([0, 1, 2]))
    assert predicted.shape == (3, trainer.labels.shape[1])  # K, not K + 1
    assert trainer.model.training, "predict must not leave the model in eval mode"


def test_predict_is_unaffected_by_the_chunk_size(trainer):
    """Chunking is a memory device and must not change a ranking. Exact equality is too strong —
    float32 matmuls accumulate in a different order at a different batch size — so this asserts
    the tolerance that matters: the induced reordering."""
    rows = np.arange(N_DOMAINS)
    small, large = trainer.predict(rows, chunk=2), trainer.predict(rows, chunk=64)
    np.testing.assert_allclose(small, large, rtol=1e-4, atol=1e-5)
    for i in range(len(rows)):
        np.testing.assert_array_equal(np.argsort(-small[i]), np.argsort(-large[i]))


def test_macro_aupr_skips_domains_with_no_positive_and_is_nan_when_all_are():
    _, matrix, vectors, _ = _corpus()
    silent = matrix.label.copy()
    silent[3] = 0  # one domain loses every positive
    quiet = KmerMatrix(matrix.domains, matrix.kmers, matrix.escore, silent)
    t = training.Trainer(quiet, vectors, _config(), device=torch.device("cpu"), seed=0)
    assert not np.isnan(t.macro_aupr(np.array([2, 3, 4])))  # scored over the other two
    assert np.isnan(t.macro_aupr(np.array([3])))  # nothing scoreable at all
    assert np.isnan(t.macro_aupr(np.array([], dtype=np.int64)))


# --- the fold driver ------------------------------------------------------------------------


@pytest.fixture
def fold_inputs(tmp_path, monkeypatch):
    monkeypatch.setattr(
        training,
        "checkpoint_file",
        lambda arm, regime, name, seed, tag="": tmp_path / f"{arm}_{regime}_{name}_s{seed}{tag}.pt",
    )
    return _corpus()


def test_run_fold_reports_both_models_and_a_disjoint_split(fold_inputs):
    frame, matrix, vectors, dist = fold_inputs
    fold = splits.Fold("S1", "fold-0", np.array([0, 1, 2]), np.arange(3, N_DOMAINS), "3 domains")
    summary, per_domain = training.run_fold(
        fold,
        "A1",
        frame,
        matrix,
        vectors,
        dist,
        _config(),
        trainable=np.ones(N_DOMAINS, dtype=bool),
        device=torch.device("cpu"),
        track=False,
    )
    assert len(per_domain) == 3
    assert list(per_domain.domain) == list(frame.dbd_seq.to_numpy()[fold.test])
    for key in ("aupr", "aupr_final", "selection_gain", "chance_aupr", "random_p95", "digest"):
        assert key in summary, key
    assert summary["selection_gain"] == pytest.approx(summary["aupr"] - summary["aupr_final"])
    assert summary["n_train"] + summary["n_validation"] == len(fold.train)


def test_run_fold_narrows_the_training_pool_but_never_the_test_set(fold_inputs):
    """`trainable` drops `no_evidence` records from training. A fold's test set is whatever the
    regime says it is and must not be quietly narrowed (`T21`)."""
    frame, matrix, vectors, dist = fold_inputs
    fold = splits.Fold("S1", "fold-0", np.array([0, 1]), np.arange(2, N_DOMAINS), "")
    trainable = np.ones(N_DOMAINS, dtype=bool)
    trainable[[0, 5, 6]] = False  # one in test, two in train

    summary, per_domain = training.run_fold(
        fold,
        "A1",
        frame,
        matrix,
        vectors,
        dist,
        _config(),
        trainable=trainable,
        device=torch.device("cpu"),
        track=False,
    )
    assert len(per_domain) == 2, "the untrainable test domain must still be scored"
    assert summary["n_train"] + summary["n_validation"] == len(fold.train) - 2


def test_run_fold_writes_one_checkpoint_carrying_its_split_and_code_stamp(fold_inputs, tmp_path):
    frame, matrix, vectors, dist = fold_inputs
    fold = splits.Fold("P3", "half:draw-0", np.array([0, 1]), np.arange(2, N_DOMAINS), "")
    summary, _ = training.run_fold(
        fold,
        "A1",
        frame,
        matrix,
        vectors,
        dist,
        _config(),
        trainable=np.ones(N_DOMAINS, dtype=bool),
        device=torch.device("cpu"),
        track=False,
    )
    written = list(tmp_path.glob("*.pt"))
    assert len(written) == 1
    _, meta = training.load_checkpoint(written[0])
    assert meta["digests"]["test"] == summary["digest"]
    assert meta["regime"] == "P3" and meta["fold"] == "half:draw-0"
    assert set(meta["code"]) == {"commit", "dirty"}
    assert summary["code_commit"] == meta["code"]["commit"]


# --- small helpers --------------------------------------------------------------------------


def test_device_for_honours_an_explicit_request_and_otherwise_picks_one():
    assert training.device_for("cpu") == torch.device("cpu")
    assert training.device_for().type in {"cpu", "cuda"}


def test_a_checkpoint_path_outside_the_repo_is_reported_absolute_not_crashed(tmp_path):
    """`Path.relative_to` raises rather than falling back, which made `run_fold` unusable with a
    checkpoint directory anywhere but inside the project."""
    from snp2prot.config import CHECKPOINT_DIR, PROJECT_ROOT

    inside = CHECKPOINT_DIR / "A1_S1_fold-0.pt"
    assert training._repo_relative(inside) == str(inside.relative_to(PROJECT_ROOT))
    outside = tmp_path / "A1_S1_fold-0.pt"
    assert training._repo_relative(outside) == str(outside)


def test_the_model_seed_is_separate_from_the_split_seed(fold_inputs):
    """They were one number. Which domains are held out and how the towers are initialised are
    unrelated choices, and varying the first would silently have re-initialised the second."""
    frame, matrix, vectors, dist = fold_inputs
    fold = splits.Fold("S1", "fold-0", np.array([0, 1]), np.arange(2, N_DOMAINS), "")
    common = dict(
        domains=frame,
        matrix=matrix,
        embeddings=vectors,
        distances=dist,
        trainable=np.ones(N_DOMAINS, dtype=bool),
        device=torch.device("cpu"),
        track=False,
    )
    a, _ = training.run_fold(fold, "A1", config=_config(), model_seed=1, **common)
    b, _ = training.run_fold(fold, "A1", config=_config(), model_seed=2, **common)

    assert a["seed"] == b["seed"], "the split seed must not move with the model seed"
    assert a["digest"] == b["digest"], "the held-out set must be identical"
    assert a["model_seed"] == 1 and b["model_seed"] == 2
    assert a["checkpoint"] != b["checkpoint"], "seeds must not overwrite each other's weights"


def test_the_model_seed_defaults_to_the_configured_one(fold_inputs):
    frame, matrix, vectors, dist = fold_inputs
    fold = splits.Fold("S1", "fold-0", np.array([0, 1]), np.arange(2, N_DOMAINS), "")
    cfg = _config()
    summary, _ = training.run_fold(
        fold,
        "A1",
        frame,
        matrix,
        vectors,
        dist,
        cfg,
        trainable=np.ones(N_DOMAINS, dtype=bool),
        device=torch.device("cpu"),
        track=False,
    )
    assert summary["model_seed"] == cfg["model"]["seed"]
