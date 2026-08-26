"""Retrieval metrics.

The dataset is 466:1 negative to positive, so the whole point of these is that they cannot be
satisfied by getting the negatives right. Each test states the property that matters at that
imbalance rather than pinning a number for its own sake.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from snp2prot.evaluation import metrics


def test_a_perfect_ranking_scores_one_and_a_reversed_one_scores_near_the_positive_rate():
    labels = np.zeros(1000, dtype=int)
    labels[:5] = 1
    perfect = -np.arange(1000.0)  # positives first
    assert metrics.average_precision(labels, perfect) == 1.0
    assert metrics.average_precision(labels, -perfect) < 0.02


def test_average_precision_is_the_mean_precision_at_each_positive():
    """Hand-computable: hits at ranks 1, 3 and 4 give (1/1 + 2/3 + 3/4) / 3."""
    labels = np.array([1, 0, 1, 1, 0])
    scores = np.array([5.0, 4.0, 3.0, 2.0, 1.0])
    assert metrics.average_precision(labels, scores) == pytest.approx((1 + 2 / 3 + 3 / 4) / 3)


def test_a_constant_prediction_scores_the_positive_rate_and_nothing_more():
    """The failure mode §5.3 names, in ranking form.

    A prediction that expresses no ordering at all must score chance. It only does so if tied
    scores are collapsed into one threshold: sorting them by array position would put this
    corpus's positives first and return AUPR 1.0.
    """
    labels = np.zeros(10_000, dtype=int)
    labels[:21] = 1  # the corpus's 0.21% positive rate
    constant = np.zeros(10_000)
    assert metrics.average_precision(labels, constant) == pytest.approx(0.0021)
    assert metrics.precision_at_k(labels, constant, 10) == pytest.approx(0.0021)


def test_a_domain_with_no_positives_has_no_aupr_rather_than_a_zero():
    """The 54 zero-positive records: `nan`, so a macro-average skips them and says so."""
    labels = np.zeros(100, dtype=int)
    assert math.isnan(metrics.average_precision(labels, np.random.default_rng(0).normal(size=100)))
    assert math.isnan(metrics.auroc(labels, np.arange(100.0)))
    assert math.isnan(metrics.recall_at_precision(labels, np.arange(100.0), 0.5))


def test_recall_at_precision_is_nan_when_the_precision_is_never_reached():
    """Different from recall 0, and reported differently."""
    labels = np.zeros(1000, dtype=int)
    labels[-1] = 1  # the single positive ranks last
    scores = -np.arange(1000.0)
    assert math.isnan(metrics.recall_at_precision(labels, scores, 0.5))


def test_the_no_call_band_is_excluded_from_the_ranking_metrics():
    """`label == -1` is absent evidence, not a negative (§3.1: two masks, not filters)."""
    labels = np.array([1, -1, -1, 0, 0])
    escore = np.array([0.5, 0.4, 0.4, 0.1, 0.0])
    predicted = np.array([0.9, 0.8, 0.7, 0.1, 0.0])  # the no-calls rank second and third
    scored = metrics.score_domain(labels, escore, predicted)
    assert scored.n_pos == 1
    assert scored.aupr == 1.0  # unaffected by the two no-calls ranked above the negatives


def test_ties_share_their_mean_rank():
    assert list(metrics._average_ranks(np.array([1.0, 1.0, 2.0, 3.0, 3.0, 3.0]))) == [
        1.5,
        1.5,
        3.0,
        5.0,
        5.0,
        5.0,
    ]


def test_spearman_is_one_against_itself_and_minus_one_reversed():
    x = np.array([0.1, 0.4, 0.2, 0.9, 0.3])
    assert metrics.spearman(x, x) == 1.0
    assert metrics.spearman(x, -x) == -1.0


def test_the_macro_average_reports_how_many_domains_it_skipped():
    scored = [
        metrics.score_domain(np.array([1, 0, 0, 0]), np.zeros(4), np.array([4.0, 3, 2, 1])),
        metrics.score_domain(np.array([0, 0, 0, 0]), np.zeros(4), np.array([4.0, 3, 2, 1])),
    ]
    summary = metrics.macro_average(scored)
    assert summary["n_domains"] == 2
    assert summary["n_scored"] == 1
    assert summary["n_undefined"] == 1
    assert summary["aupr"] == 1.0


def test_chance_is_the_mean_positive_rate_of_the_domains_the_aupr_mean_covers():
    """A domain with no positives has no AUPR, so it must not be averaged into the null either.

    `run_nn_baseline.positive_rate` used to include them at rate 0, which pushed the null down
    and every reported `x chance` ratio up — by 15% on the P3 folds, where 19 of 173 held-out
    variants have no positive 8-mer.
    """
    labels = np.array(
        [
            [1, 1, 0, 0],  # rate 0.5
            [1, 0, 0, 0],  # rate 0.25
            [0, 0, 0, 0],  # no positives — skipped, not counted as 0.0
        ]
    )
    assert metrics.chance_aupr(labels) == pytest.approx(0.375)


def test_chance_ignores_the_no_call_band_like_every_other_ranking_metric():
    labels = np.array([[1, 0, -1, -1]])  # 1 of 2 scored cells, not 1 of 4
    assert metrics.chance_aupr(labels) == pytest.approx(0.5)


def test_chance_is_nan_when_no_domain_has_a_positive():
    assert np.isnan(metrics.chance_aupr(np.array([[0, 0, -1], [0, -1, 0]])))


def test_the_null_band_brackets_the_chance_level_and_a_random_ranking_falls_inside_it():
    rng = np.random.default_rng(0)
    labels = np.zeros((6, 400), dtype=np.int8)
    for i in range(6):
        labels[i, rng.choice(400, 8, replace=False)] = 1

    band = metrics.null_aupr(labels, repeats=200, seed=1)
    assert band["null_mean"] < band["null_p95"]

    drawn = [
        np.mean([metrics.average_precision(labels[i], rng.random(400)) for i in range(6)])
        for _ in range(40)
    ]
    # A genuinely random ranking should clear the 95th percentile about 5% of the time.
    assert np.mean(np.array(drawn) > band["null_p95"]) < 0.25


def test_a_perfect_ranking_sits_far_above_the_null_and_the_null_is_near_chance():
    labels = np.zeros((4, 500), dtype=np.int8)
    labels[:, :10] = 1
    band = metrics.null_aupr(labels, repeats=100, seed=0)
    assert band["null_p95"] < 0.15
    assert metrics.chance_aupr(labels) == pytest.approx(0.02)
    assert band["null_mean"] > metrics.chance_aupr(labels)  # AP is biased up at small n_pos


def test_the_null_is_nan_rather_than_zero_when_there_is_nothing_to_rank():
    band = metrics.null_aupr(np.array([[0, 0, -1]]), repeats=5)
    assert np.isnan(band["null_mean"]) and np.isnan(band["null_p95"])
