"""The decision-rule metrics: `snp2prot.evaluation.calibration`.

Every test pins a property the module docstring gives a reason for. The ones that matter most
are the two `nan` conventions — "the model called nothing" and "the protein binds nothing" are
different statements and must not collapse into one number — and the direction of
`detection_auroc`, which is the metric `T38` was raised to obtain and would be silently
backwards if nothing checked it.
"""

from __future__ import annotations

import numpy as np
import pytest

from snp2prot.evaluation import calibration

# --- calibration error --------------------------------------------------------------------


def test_a_perfectly_calibrated_prediction_scores_zero():
    labels = np.array([1, 1, 0, 0, 0, 0])
    probabilities = np.array([1.0, 1.0, 0.0, 0.0, 0.0, 0.0])
    assert calibration.expected_calibration_error(labels, probabilities, bins=2) == pytest.approx(
        0.0, abs=1e-9
    )


def test_a_confidently_wrong_prediction_scores_one():
    labels = np.array([0, 0, 1, 1])
    probabilities = np.array([1.0, 1.0, 0.0, 0.0])
    assert calibration.expected_calibration_error(labels, probabilities, bins=2) == pytest.approx(
        1.0, abs=1e-9
    )


def test_the_bins_are_equal_mass_so_the_negative_cloud_is_not_one_bin():
    """The reason for quantile binning: at this base rate an equal-width binning would put every
    cell below 0.067 into one bucket and report a single number for the part being asked about."""
    rng = np.random.default_rng(0)
    probabilities = np.concatenate([rng.uniform(0.0, 0.01, 990), rng.uniform(0.9, 1.0, 10)])
    labels = np.concatenate([np.zeros(990), np.ones(10)])
    # Well calibrated, and an equal-width binning would still see only two occupied bins.
    assert calibration.expected_calibration_error(labels, probabilities, bins=10) < 0.02


# --- the call ------------------------------------------------------------------------------


def test_calling_nothing_leaves_precision_undefined_but_recall_zero():
    """Two different absences. Calling nothing is a decision the model made — its recall is 0,
    a failure with a value — while its precision is over an empty set and does not exist."""
    scores = calibration.score_calls(np.array([1, 0, 0]), np.array([0.1, 0.1, 0.1]), 0.5)
    assert np.isnan(scores.precision)
    assert scores.recall == 0.0
    assert scores.f1 == 0.0
    assert scores.n_called == 0


def test_a_protein_that_binds_nothing_leaves_recall_undefined():
    """The complement: no positive 8-mer is a property of the assay, not a model failure."""
    scores = calibration.score_calls(np.array([0, 0, 0]), np.array([0.9, 0.1, 0.1]), 0.5)
    assert np.isnan(scores.recall)
    assert scores.precision == 0.0
    assert scores.f1 == 0.0


def test_both_sets_empty_is_the_only_case_where_every_set_metric_is_nan():
    scores = calibration.score_calls(np.array([0, 0]), np.array([0.1, 0.1]), 0.5)
    assert np.isnan(scores.precision) and np.isnan(scores.recall)
    assert np.isnan(scores.f1) and np.isnan(scores.jaccard)


def test_a_perfect_call_scores_one_everywhere():
    scores = calibration.score_calls(np.array([1, 1, 0, 0]), np.array([0.9, 0.8, 0.1, 0.0]), 0.5)
    assert scores.precision == 1.0 and scores.recall == 1.0
    assert scores.f1 == 1.0 and scores.jaccard == 1.0


def test_the_gray_band_never_reaches_score_calls():
    """It is masked by the caller, so a `-1` must not be counted as a positive here if one slips
    through — `labels == 1` is the truth set, not `labels != 0`."""
    scores = calibration.score_calls(np.array([1, -1, 0]), np.array([0.9, 0.9, 0.1]), 0.5)
    assert scores.n_true == 1
    assert scores.precision == pytest.approx(0.5)


# --- the parameter-free rule ---------------------------------------------------------------


def test_the_expected_count_rule_calls_exactly_the_expected_number():
    """`Σ p` is the expected number of binders when the probabilities are calibrated, so the cut
    that keeps that many is the rule with no free parameter."""
    probabilities = np.array([0.9, 0.8, 0.2, 0.05, 0.05])  # sums to 2.0
    cut = calibration.expected_count_rule(probabilities)
    assert int((probabilities >= cut).sum()) == 2


def test_a_model_that_believes_nothing_binds_calls_nothing():
    """The behaviour the dead-variant case needs: rounding to zero must produce an empty set
    rather than a top-1 fallback that always names something."""
    probabilities = np.array([0.001, 0.002, 0.001])
    cut = calibration.expected_count_rule(probabilities)
    assert int((probabilities >= cut).sum()) == 0


def test_calling_nothing_survives_a_float32_profile():
    """The regression for a bug that inverted the rule exactly where it matters most. NumPy's
    weak scalar promotion casts a Python float down to a `float32` array's dtype, so a cut
    expressed as "one ulp above the maximum" — a denormal double — became `0.0` and called
    **every** 8-mer instead of none. Measured on the baseline: 32,894 calls where 0 was meant.
    Both the model's probabilities and the baseline's profiles are `float32`."""
    probabilities = np.zeros(2048, dtype=np.float32)
    cut = calibration.expected_count_rule(probabilities)
    assert int((probabilities >= cut).sum()) == 0
    assert calibration.score_calls(np.zeros(2048, dtype=np.int8), probabilities, cut).n_called == 0


def test_the_rule_never_asks_for_more_kmers_than_exist():
    probabilities = np.full(4, 0.99)
    cut = calibration.expected_count_rule(probabilities)
    assert int((probabilities >= cut).sum()) == 4


# --- the protein-level statistic -----------------------------------------------------------


def test_interaction_power_is_the_expected_count():
    assert calibration.interaction_power(np.array([0.5, 0.25, 0.25])) == pytest.approx(1.0)


def test_detection_auroc_is_one_when_low_power_means_dead():
    """The direction check. Dead is the positive class and *low* interaction power is the
    evidence for it, so a perfect separation must score 1.0 and not 0.0."""
    is_dead = np.array([1, 1, 0, 0])
    power = np.array([0.1, 0.5, 40.0, 90.0])
    assert calibration.detection_auroc(is_dead, power) == pytest.approx(1.0)


def test_detection_auroc_is_zero_when_the_signal_is_backwards():
    is_dead = np.array([1, 1, 0, 0])
    power = np.array([90.0, 40.0, 0.5, 0.1])
    assert calibration.detection_auroc(is_dead, power) == pytest.approx(0.0)


def test_detection_auroc_is_undefined_without_both_classes():
    """A fold whose held-out domains all bind, or all do not, has no detection to measure — and
    reporting 0.5 for it would read as "no information" rather than "not asked"."""
    assert np.isnan(calibration.detection_auroc(np.array([0, 0, 0]), np.array([1.0, 2.0, 3.0])))
    assert np.isnan(calibration.detection_auroc(np.array([1, 1]), np.array([1.0, 2.0])))


def test_power_ratio_is_one_when_the_mutation_is_predicted_to_do_nothing():
    """What copying the wild type scores, which is what the nearest-neighbour baseline does by
    construction under `P3/all`."""
    assert calibration.power_ratio(50.0, 50.0) == pytest.approx(1.0)
    assert calibration.power_ratio(0.0, 50.0) == pytest.approx(0.0)


def test_power_ratio_is_undefined_when_the_wild_type_is_predicted_dead_too():
    """The ratio would then be a statement about the model's opinion of the wild type rather
    than about the mutation."""
    assert np.isnan(calibration.power_ratio(1.0, 0.0))
