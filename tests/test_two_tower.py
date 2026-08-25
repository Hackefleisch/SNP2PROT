"""The two towers, and the properties the design depends on.

What matters is not the numbers a randomly initialised model produces but the invariants: that
the DNA encoder is strand-symmetric by construction, that the null anchor lives in the same
space as the 8-mers, and that both towers reach the same width.
"""

from __future__ import annotations

import numpy as np
import pytest

torch = pytest.importorskip("torch")

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
        logits = model(torch.randn(2, 8), torch.tensor(tokenise(KMERS)))
    assert torch.isfinite(logits).all()
    assert logits.abs().max() <= 100.0 + 1e-3


def test_parameter_counts_are_reported_per_tower():
    """Equal `D` controls the shared space but not the capacity feeding it, so a cross-arm
    comparison is only readable next to these (`ML_PLAN.md` §4.2)."""
    counts = TwoTower(protein_features=1280, width=256, dna_channels=64).parameter_counts()
    assert counts["protein_tower"] == 1280 * 256 + 256
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
