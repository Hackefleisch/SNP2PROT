"""The one probe that survived, and the property that makes it readable.

`rank_ceiling` is a *lower* bound on what rank `r` can express — an oracle fit, and Frobenius-
optimal rather than AUPR-optimal. Both errors point the same way, which is what lets its height
rule the shared space out as the constraint. See `snp2prot.evaluation.ceilings`.
"""

from __future__ import annotations

import numpy as np

from snp2prot.data.matrix import KmerMatrix
from snp2prot.evaluation import ceilings


def _matrix(n_domains: int = 12, n_kmers: int = 60, rank: int = 3, seed: int = 0) -> KmerMatrix:
    """A LABEL matrix that is genuinely low rank. `rank_ceiling` factorises the calls, not the
    E-scores, so the low-rank structure has to live in the labels for the probe to find it."""
    rng = np.random.default_rng(seed)
    escore = (rng.normal(size=(n_domains, rank)) @ rng.normal(size=(rank, n_kmers))).astype(
        np.float32
    )
    # Rank-`rank` labels: threshold a rank-`rank` score matrix per domain, so each row's call
    # pattern is a function of the same few factors.
    label = (escore > np.quantile(escore, 0.9, axis=1, keepdims=True)).astype(np.int8)
    return KmerMatrix(
        domains=np.array([f"D{i}" for i in range(n_domains)]),
        kmers=np.array([f"K{j}" for j in range(n_kmers)]),
        escore=escore,
        label=label,
    )


def test_more_rank_never_scores_worse():
    scores = ceilings.rank_ceiling(_matrix(), ranks=(1, 2, 4, 8))
    values = [scores[r] for r in (1, 2, 4, 8)]
    assert values == sorted(values)


def test_low_rank_label_structure_is_recovered_at_full_rank():
    """The floor has to be tight where the data really is low rank, or its height means nothing.
    A binary matrix is not itself low rank even when the scores behind it are, so the tight case
    is full rank rather than the rank of the generating factors."""
    matrix = _matrix(rank=3)
    scores = ceilings.rank_ceiling(matrix, ranks=(3, len(matrix.domains)))
    assert scores[len(matrix.domains)] > 0.99
    assert scores[len(matrix.domains)] >= scores[3]


def test_the_reported_number_is_an_in_sample_oracle():
    """It scores the rows it was fitted on, which is why the report says so in as many words."""
    matrix = _matrix()
    full = ceilings.rank_ceiling(matrix, ranks=(len(matrix.domains),))
    assert full[len(matrix.domains)] > 0.99


def test_the_removed_probes_are_gone():
    """They were lower bounds read as upper ones; `T36`'s ablation replaces them."""
    assert not hasattr(ceilings, "ridge_probe")
    assert not hasattr(ceilings, "kernel_probe")
