"""The nearest-neighbour lookup baseline.

It is the bar every model is read against, so what it does has to be exactly what it claims:
copy the most identical *eligible* training domain's E-score profile, and never copy across an
alignment the overlap guard rejects.
"""

from __future__ import annotations

import numpy as np
import pytest

from snp2prot import distances
from snp2prot.baselines import nn_lookup
from snp2prot.data.matrix import KmerMatrix

BASE = "ACDEFGHIKLMNPQRSTVWY" * 4  # 80 aa
NEAR = BASE[:3] + "P" + BASE[4:]  # one edit
FAR = "".join("P" if i % 4 == 0 and i < 40 else c for i, c in enumerate(BASE))  # ten edits
SLIVER = "WYWYWYWY" * 10  # unrelated, and 80 aa so length alone does not exclude it


@pytest.fixture
def corpus():
    """The query, a one-edit relative, a ten-edit relative, and an unrelated domain."""
    sequences = np.array([BASE, NEAR, FAR, SLIVER])
    d = distances.build(sequences, processes=1)
    escore = np.array(
        [
            [0.0, 0.0, 0.0, 0.0, 0.0, 0.0],  # the query's own, never used as a prediction
            [0.5, 0.4, 0.3, 0.2, 0.1, 0.0],  # the one-edit relative's
            [0.4, 0.5, 0.2, 0.3, 0.0, 0.1],  # the ten-edit relative's
            [0.0, 0.1, 0.2, 0.3, 0.4, 0.5],  # the unrelated one's, deliberately reversed
        ],
        dtype=np.float32,
    )
    label = np.zeros((4, 6), dtype=np.int8)
    return KmerMatrix(sequences, np.array(list("ABCDEF")), escore, label), d


def test_it_copies_the_most_identical_training_domain(corpus):
    matrix, d = corpus
    got = nn_lookup.fit_predict(matrix, d, np.array([0]), np.array([1, 2, 3]), min_overlap=0.6)
    assert list(got.profile[0]) == list(matrix.escore[1])
    assert got.neighbours.neighbour[0] == NEAR
    assert got.neighbours.n_edits[0] == 1


def test_an_unrelated_domain_reports_perfect_identity_over_a_sliver(corpus):
    """The `T25` degeneracy, stated as a fact about the data rather than about the code: with
    free terminal gaps the aligner keeps two residues of each sequence and discards the rest
    for nothing, so `identity` is 1.0 and `n_edits` is 0 between two unrelated domains."""
    _, d = corpus
    assert d.n_edits[0, 3] == 0
    assert d.identity()[0, 3] == 1.0
    assert d.overlap()[0, 3] == pytest.approx(2 / 80)


def test_the_guard_is_what_stops_that_sliver_being_copied_from(corpus):
    matrix, d = corpus
    unguarded = nn_lookup.fit_predict(matrix, d, np.array([0]), np.array([1, 3]), min_overlap=0.0)
    assert unguarded.neighbours.neighbour[0] == SLIVER  # 1.0 identity beats the relative's 0.99

    guarded = nn_lookup.fit_predict(matrix, d, np.array([0]), np.array([1, 3]), min_overlap=0.6)
    assert guarded.neighbours.neighbour[0] == NEAR


def test_a_query_with_no_eligible_neighbour_gets_the_mean_profile_and_is_counted(corpus):
    """Never the best *ineligible* candidate: that would copy whichever unrelated domain the
    aligner mangled most favourably. Measured on the real corpus this never fires."""
    matrix, d = corpus
    got = nn_lookup.fit_predict(matrix, d, np.array([3]), np.array([0, 1]), min_overlap=0.6)
    assert got.n_without_neighbour == 1
    assert list(got.profile[0]) == list(matrix.escore[[0, 1]].mean(axis=0))


def test_top_k_averages_the_neighbours_but_still_reports_the_closest(corpus):
    matrix, d = corpus
    got = nn_lookup.fit_predict(matrix, d, np.array([0]), np.array([1, 2]), 0.6, k=2)
    assert got.neighbours.neighbour[0] == NEAR
    between = (matrix.escore[1] + matrix.escore[2]) / 2
    assert got.profile[0] == pytest.approx(between, abs=0.01)  # identity-weighted, near equal


def test_an_empty_training_pool_is_an_error_rather_than_a_silent_mean(corpus):
    matrix, d = corpus
    with pytest.raises(ValueError, match="nothing to copy"):
        nn_lookup.fit_predict(matrix, d, np.array([0]), np.array([], dtype=int), 0.6)
