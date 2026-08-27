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


def test_recall_at_precision_is_zero_when_the_precision_is_never_reached():
    """A failure with a value, not an absence. Returning nan here and letting `macro_average`
    drop it reported the NN baseline's `P1` recall as 0.130 over 24 of 412 domains, where over
    all of them it is 0.008."""
    labels = np.zeros(1000, dtype=int)
    labels[-1] = 1  # the single positive ranks last: no cut ever reaches precision 0.5
    scores = -np.arange(1000.0)
    assert metrics.recall_at_precision(labels, scores, 0.5) == 0.0


def test_recall_at_precision_is_nan_only_when_there_is_no_positive_to_recall():
    labels = np.zeros(100, dtype=int)
    assert math.isnan(
        metrics.recall_at_precision(labels, np.random.default_rng(0).random(100), 0.5)
    )


def test_a_perfect_ranking_still_reaches_full_recall_at_the_target():
    labels = np.zeros(100, dtype=int)
    labels[:10] = 1
    assert metrics.recall_at_precision(labels, -np.arange(100.0), 0.5) == pytest.approx(1.0)


def test_failures_stay_in_the_macro_average_instead_of_vanishing_from_it():
    """One domain reaches the target and one cannot. The mean must be over both."""
    good = np.zeros(100, dtype=int)
    good[:5] = 1
    bad = np.zeros(100, dtype=int)
    bad[-1] = 1
    scores = -np.arange(100.0)
    rows = [
        metrics.score_domain(np.where(v == 1, 1, 0), scores, scores, precision_at=(10,))
        for v in (good, bad)
    ]
    summary = metrics.macro_average(rows)
    assert summary["n_scored_recall_at_precision"] == 2
    assert summary["recall_at_precision"] == pytest.approx(0.5)  # (1.0 + 0.0) / 2


def test_every_metric_reports_how_many_domains_its_mean_is_over():
    """`n_scored` counts AUPR only, and was read as covering everything."""
    labels = np.zeros(50, dtype=int)
    labels[:4] = 1
    scores = -np.arange(50.0)
    rows = [
        metrics.score_domain(labels, scores, scores, precision_at=(5,)),
        metrics.score_domain(np.zeros(50, dtype=int), scores, scores, precision_at=(5,)),
    ]
    summary = metrics.macro_average(rows)
    for name in ("aupr", "auroc", "recall_at_precision", "precision_at_5"):
        assert f"n_scored_{name}" in summary, name
    assert summary["n_scored_aupr"] == 1  # the second domain has no positives
    assert summary["n_scored_precision_at_5"] == 2  # precision@k is defined for it regardless


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


def test_spearman_survives_as_a_utility_but_is_no_longer_a_model_metric():
    """`check_pooling.py` needs a rank correlation for embedding displacement vs edit count.
    Nothing scores a prediction against the raw E-score any more."""
    assert not hasattr(metrics.DomainScores("d", 0, 0.0, 0.0, 0.0), "spearman")


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

    band = metrics.random_baseline(labels, repeats=200, seed=1)
    assert band["random_mean"] < band["random_p95"]

    drawn = [
        np.mean([metrics.average_precision(labels[i], rng.random(400)) for i in range(6)])
        for _ in range(40)
    ]
    # A genuinely random ranking should clear the 95th percentile about 5% of the time.
    assert np.mean(np.array(drawn) > band["random_p95"]) < 0.25


def test_a_perfect_ranking_sits_far_above_the_null_and_the_null_is_near_chance():
    labels = np.zeros((4, 500), dtype=np.int8)
    labels[:, :10] = 1
    band = metrics.random_baseline(labels, repeats=100, seed=0)
    assert band["random_p95"] < 0.15
    assert metrics.chance_aupr(labels) == pytest.approx(0.02)
    assert band["random_mean"] > metrics.chance_aupr(labels)  # AP is biased up at small n_pos


def test_the_null_is_nan_rather_than_zero_when_there_is_nothing_to_rank():
    band = metrics.random_baseline(np.array([[0, 0, -1]]), repeats=5)
    assert np.isnan(band["random_mean"]) and np.isnan(band["random_p95"])


# --- suppression: the metric for the domains AUPR cannot reach --------------------------


def _dead_case(n=200, n_pos=20):
    labels = np.zeros(n, dtype=np.int64)
    labels[:n_pos] = 1
    wild_type = labels.astype(np.float64)  # the wild type's own calls
    return labels, wild_type


def test_copying_the_wild_type_scores_exactly_zero():
    """What the nearest-neighbour baseline does, so it is C1's null hypothesis literally."""
    labels, wild_type = _dead_case()
    assert metrics.suppression(wild_type, wild_type, labels) == 0.0


def test_ranking_every_wild_type_site_lower_scores_one():
    labels, wild_type = _dead_case()
    variant = -wild_type  # every site the wild type binds is now bottom-ranked
    assert metrics.suppression(variant, wild_type, labels) == 1.0


def test_a_global_downward_shift_is_not_mistaken_for_sensitivity():
    """Comparing raw scores would hand 1.0 to a model that merely scores the variant lower
    everywhere. Ranks are invariant to that, so it scores 0 — no site actually moved."""
    labels, wild_type = _dead_case()
    assert metrics.suppression(wild_type - 100.0, wild_type, labels) == 0.0
    assert metrics.suppression(wild_type * 0.001, wild_type, labels) == 0.0


def test_random_rearrangement_lands_near_a_half():
    rng = np.random.default_rng(0)
    labels, wild_type = _dead_case(n=2000, n_pos=200)
    draws = [
        metrics.suppression(rng.normal(size=2000), rng.normal(size=2000), labels) for _ in range(40)
    ]
    assert 0.4 < float(np.mean(draws)) < 0.6


def test_suppression_is_undefined_when_the_reference_binds_nothing_either():
    labels = np.zeros(50, dtype=np.int64)
    assert np.isnan(metrics.suppression(np.zeros(50), np.zeros(50), labels))
