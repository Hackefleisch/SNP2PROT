"""The two towers, and the properties the design depends on.

What matters is not the numbers a randomly initialised model produces but the invariants: that
the DNA encoder is strand-symmetric by construction, that the null anchor lives in the same
space as the 8-mers, and that both towers reach the same width.
"""

from __future__ import annotations

import numpy as np
import pytest

torch = pytest.importorskip("torch")

from snp2prot.evaluation import metrics  # noqa: E402
from snp2prot.models import DNAEncoder, ProteinTower, TwoTower, tokenise  # noqa: E402
from snp2prot.models.encoders import reverse_complement  # noqa: E402

KMERS = np.array(["AAAAAAAC", "ACGTACGT", "GGGGCCCC", "TTTTTTTT"])


def test_tokenise_maps_bases_so_that_the_complement_is_three_minus_the_index():
    """What lets the reverse complement be a flip and a subtraction rather than a lookup."""
    tokens = tokenise(np.array(["ACGT"]))
    assert list(tokens[0]) == [0, 1, 2, 3]
    assert list(reverse_complement(torch.tensor(tokens))[0]) == [0, 1, 2, 3]  # ACGT is a palindrome


def test_reverse_complement_of_a_known_pair():
    tokens = torch.tensor(tokenise(np.array(["AAAAAAAC"])))
    assert list(reverse_complement(tokens)[0].numpy()) == list(tokenise(np.array(["GTTTTTTT"]))[0])


def test_the_dna_encoder_gives_a_stored_8mer_and_its_reverse_complement_one_vector():
    """The array measured both strands together and the vocabulary keeps only one of each pair,
    so the embedding must be strand-symmetric *by construction* rather than learned
    (`ML_PLAN.md` §4.1). Architectural pooling makes the difference exactly zero."""
    encoder = DNAEncoder(width=16, channels=8).eval()
    tokens = torch.tensor(tokenise(KMERS))
    with torch.no_grad():
        assert torch.allclose(encoder(tokens), encoder(reverse_complement(tokens)), atol=1e-6)


def test_both_towers_produce_unit_vectors_of_the_shared_width():
    model = TwoTower(protein_features=32, width=16, dna_channels=8).eval()
    with torch.no_grad():
        dna = model.dna(torch.tensor(tokenise(KMERS)))
        protein = model.protein(torch.randn(3, 32))
    assert dna.shape == (len(KMERS), 16)
    assert protein.shape == (3, 16)
    assert torch.allclose(dna.norm(dim=1), torch.ones(len(KMERS)), atol=1e-5)
    assert torch.allclose(protein.norm(dim=1), torch.ones(3), atol=1e-5)


def test_the_dna_table_carries_exactly_one_extra_row_for_the_null():
    model = TwoTower(protein_features=32, width=16, dna_channels=8).eval()
    with torch.no_grad():
        table = model.dna_table(torch.tensor(tokenise(KMERS)))
    assert table.shape == (len(KMERS) + 1, 16)
    # The null is normalised alongside the real 8-mers, so it is a point on the same sphere and
    # every score in the model is a cosine.
    assert float(table[-1].norm()) == pytest.approx(1.0, abs=1e-5)


def test_the_temperature_is_clamped_so_a_learned_scale_cannot_saturate_the_softmax():
    model = TwoTower(protein_features=8, width=4, dna_channels=4)
    with torch.no_grad():
        model.logit_scale.fill_(50.0)  # far past any sane value
    model.clamp_temperature()
    with torch.no_grad():
        logits = model(torch.randn(2, 8), torch.tensor(tokenise(KMERS)))
    assert torch.isfinite(logits).all()
    assert logits.abs().max() <= 100.0 + 1e-3
    assert model.temperature_is_clamped


def test_the_scale_keeps_its_gradient_at_the_ceiling():
    """Clamping the forward pass instead left a dead zone: past the ceiling the gradient was
    exactly 0, so the parameter oscillated on weight decay alone, or froze solid without it
    (measured on `P3/all`, 2026-08-26). Clamping the value keeps the loss able to see it."""
    model = TwoTower(protein_features=8, width=4, dna_channels=4)
    with torch.no_grad():
        model.logit_scale.fill_(50.0)
    model.clamp_temperature()
    model(torch.randn(2, 8), torch.tensor(tokenise(KMERS))).sum().backward()
    assert model.logit_scale.grad is not None
    assert model.logit_scale.grad.abs() > 0


def test_clamping_never_raises_the_scale():
    """It is a ceiling, not a target: a model below it must be left alone."""
    model = TwoTower(protein_features=8, width=4, dna_channels=4, temperature=0.5)
    before = float(model.logit_scale)
    model.clamp_temperature()
    assert float(model.logit_scale) == before
    assert not model.temperature_is_clamped


def test_every_reported_metric_is_invariant_to_the_temperature():
    """`score` is `scale x cosine` and every metric ranks within one domain, so the temperature
    cannot reorder anything. This is why the clamp value is not a hyperparameter that needs
    sweeping — measured end to end at 0.8558 clamped to 100 against 0.8554 running free to 148."""
    rng = np.random.default_rng(0)
    labels = np.zeros(500, dtype=np.int64)
    labels[rng.choice(500, 20, replace=False)] = 1
    cosine = rng.normal(size=500) * 0.3
    reference = metrics.average_precision(labels, cosine)
    for scale in (1.0, 100.0, 10_000.0):
        assert metrics.average_precision(labels, cosine * scale) == reference
        assert metrics.auroc(labels, cosine * scale) == metrics.auroc(labels, cosine)


def test_parameter_counts_are_reported_per_tower():
    """Equal `D` controls the shared space but not the capacity feeding it, so a cross-arm
    comparison is only readable next to these (`ML_PLAN.md` §4.2)."""
    counts = TwoTower(protein_features=1280, width=256, dna_channels=64).parameter_counts()
    # projection + bias, plus the input LayerNorm's gain and shift.
    assert counts["protein_tower"] == 1280 * 256 + 256 + 2 * 1280
    assert (
        counts["total"] == sum(v for k, v in counts.items() if k != "total") + 1
    )  # the learned temperature


def test_a_protein_tower_with_a_hidden_layer_is_opt_in():
    """With 1,338 training proteins a hidden layer is an experiment, not a default
    (`docs/TRAINING.md` §5)."""
    linear = ProteinTower(1280, 256)
    deep = ProteinTower(1280, 256, hidden=512)
    assert sum(p.numel() for p in deep.parameters()) > sum(p.numel() for p in linear.parameters())


def test_tokenise_rejects_a_non_acgt_base():
    with pytest.raises(ValueError, match="non-ACGT"):
        tokenise(np.array(["AAAANAAA"]))


def test_the_protein_tower_is_invariant_to_a_rescale_of_its_input():
    """`A1`'s pooled norms are 4.8-9.9 and `A4`'s are 0.8-1.2, and the delta between the two arms
    is supposed to isolate the pretraining corpus. A tower that answers differently to the same
    directions at a different scale cannot do that.

    Invariance is exact up to `LayerNorm`'s `eps`, which is a fixed absolute term added to a
    variance that scales with the input. On the real arms the residual is 8e-7 for `A1` and 3e-5
    for `A4` — the arm whose per-dimension variance is 5.4e-4, so `eps` is 1.85% of it — against
    unit-norm outputs, which is four orders of magnitude below anything that reorders a ranking.
    """
    torch.manual_seed(0)
    tower = ProteinTower(64, 16)
    vectors = torch.randn(12, 64)
    for scale in (0.11, 0.5, 7.0, 100.0):
        torch.testing.assert_close(tower(vectors), tower(vectors * scale), rtol=0, atol=1e-3)


def test_a_better_separated_embedding_stays_better_separated_through_the_tower():
    """Without the input LayerNorm this failed: at ESM-DBP's scale the projection bias outweighed
    the signal, and both arms came out of an untrained tower at 0.8930 whatever went in."""

    def spread(out):
        gram = out @ out.T
        return float(gram[torch.triu_indices(len(out), len(out), offset=1).unbind()].mean())

    torch.manual_seed(0)
    common = torch.randn(1, 64)
    tight = torch.nn.functional.normalize(common + 0.2 * torch.randn(40, 64), dim=-1)
    loose = torch.nn.functional.normalize(common + 2.0 * torch.randn(40, 64), dim=-1)
    assert spread(tight) > spread(loose)

    tower = ProteinTower(64, 16)
    # …and at a tenth the norm, which is the case that used to collapse.
    assert spread(tower(tight)) > spread(tower(loose * 0.1))


def test_both_protein_tower_shapes_normalise_their_input():
    assert isinstance(ProteinTower(32, 8).net[0], torch.nn.LayerNorm)
    assert isinstance(ProteinTower(32, 8, hidden=16).net[0], torch.nn.LayerNorm)
