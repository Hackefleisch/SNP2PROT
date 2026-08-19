"""Pooled protein embeddings, and the batching that produces them.

The embeddings themselves are a model's output and not something to pin in a test. What is
pinned here is everything around them that could silently corrupt an arm: which tokens the
pooling averages, that the cached artifact survives a round trip, and that a mismatched domain
set is refused rather than quietly misaligned.
"""

from __future__ import annotations

import numpy as np
import pytest

from snp2prot import embeddings


def test_pooling_drops_the_bos_and_eos_tokens():
    """They carry sequence-level summary information every domain has equally, so averaging
    them in dilutes exactly the per-residue signal the pooling exists to capture."""
    from scripts.build_embeddings import pool

    torch = pytest.importorskip("torch")
    # One sequence of 3 residues: [BOS, r1, r2, r3, EOS], each token a distinct constant.
    representation = torch.tensor([[[0.0], [1.0], [2.0], [3.0], [99.0]]])
    assert pool(representation, [3])[0, 0] == pytest.approx(2.0)  # mean(1, 2, 3), not mean(all)


def test_batches_cover_every_sequence_exactly_once():
    from scripts.build_embeddings import batches

    sequences = ["A" * n for n in (10, 300, 45, 45, 200, 7)]
    seen = [i for batch in batches(sequences, batch_tokens=600) for i in batch]
    assert sorted(seen) == list(range(len(sequences)))


def test_a_batch_stays_within_its_token_budget_once_padded():
    """Padding is to the longest member, so the cost of a batch is longest x count — which is
    what has to be bounded, not the sum of the lengths."""
    from scripts.build_embeddings import batches

    sequences = ["A" * n for n in (400, 380, 100, 90, 80, 20)]
    for batch in batches(sequences, batch_tokens=800):
        longest = max(len(sequences[i]) for i in batch)
        assert longest * len(batch) <= 800 or len(batch) == 1


def test_the_table_round_trips_without_pickle(tmp_path):
    """An object array in an `.npz` needs `allow_pickle` on load, and a cached artifact should
    never require that to be read."""
    table = embeddings.DomainEmbeddings(
        arm="A1",
        model="esm2_t33_650M_UR50D",
        domains=np.array(["ACDE", "ACDF"]),
        vectors=np.array([[1.0, 0.0], [0.0, 1.0]], dtype=np.float32),
    )
    path = table.save(tmp_path / "A1.npz")
    back = embeddings.DomainEmbeddings.load("A1", path)
    assert list(back.domains) == ["ACDE", "ACDF"]
    assert back.model == "esm2_t33_650M_UR50D"
    assert back.width == 2


def test_a_domain_without_an_embedding_is_an_error_not_a_silent_gap():
    table = embeddings.DomainEmbeddings("A1", "m", np.array(["ACDE"]), np.zeros((1, 2), np.float32))
    with pytest.raises(KeyError, match="no embedding"):
        table.rows_for(["ACDE", "WWWW"])


def test_cosine_distance_ignores_vector_length():
    """PLM embedding norms vary with sequence length, and a variant is by construction almost
    the same length as its reference — a norm difference would report as a distance."""
    a = np.array([[1.0, 2.0, 3.0]])
    assert embeddings.cosine_distance(a, a * 7.0)[0] == pytest.approx(0.0, abs=1e-12)
    assert embeddings.cosine_distance(a, -a)[0] == pytest.approx(2.0)


def test_pairwise_distance_agrees_with_the_rowwise_one():
    rng = np.random.default_rng(0)
    vectors = rng.normal(size=(5, 8)).astype(np.float32)
    full = embeddings.pairwise_cosine_distance(vectors)
    for i in range(5):
        for j in range(5):
            expected = embeddings.cosine_distance(vectors[i : i + 1], vectors[j : j + 1])[0]
            assert full[i, j] == pytest.approx(expected, abs=1e-9)
