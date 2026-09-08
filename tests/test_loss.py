"""Multi-positive InfoNCE with a null anchor.

The loss is the part of the model most able to be quietly wrong: it would still train, and the
numbers would still look plausible. Every test here pins a property `docs/TRAINING.md` §2 or §3
gives a reason for, rather than a value.
"""

from __future__ import annotations

import math

import pytest

torch = pytest.importorskip("torch")

from snp2prot.models.loss import multi_positive_infonce  # noqa: E402


def logits_for(rows, null):
    """`(B, K + 1)` logits from per-row 8-mer scores plus a null score."""
    return torch.tensor([list(r) + [n] for r, n in zip(rows, null, strict=True)], dtype=torch.float)


def test_ranking_every_positive_above_every_negative_drives_the_loss_to_zero():
    """The property that makes the objective agree with the metric: this configuration is
    AUPR = 1, and it must also be the loss's minimum."""
    logits = logits_for([[20.0, 20.0, -20.0, -20.0]], [-20.0])
    labels = torch.tensor([[1, 1, 0, 0]], dtype=torch.int8)
    loss, _ = multi_positive_infonce(logits, labels)
    assert float(loss) == pytest.approx(0.0, abs=1e-6)


def test_a_reversed_ranking_is_penalised_heavily():
    logits = logits_for([[-20.0, -20.0, 20.0, 20.0]], [-20.0])
    labels = torch.tensor([[1, 1, 0, 0]], dtype=torch.int8)
    loss, _ = multi_positive_infonce(logits, labels)
    assert float(loss) > 20.0


def test_sibling_positives_do_not_compete_with_each_other():
    """The departure from SupCon's `L_out` (§2.3b). Under `L_out` the loss of a domain with
    many positives floors at `log|P|`; here a domain with two positives ranked top scores the
    same as a domain with one, because its positives are not in each other's denominators."""
    one = multi_positive_infonce(
        logits_for([[20.0, -20.0, -20.0]], [-20.0]), torch.tensor([[1, 0, 0]], dtype=torch.int8)
    )[0]
    many = multi_positive_infonce(
        logits_for([[20.0, 20.0, -20.0]], [-20.0]), torch.tensor([[1, 1, 0]], dtype=torch.int8)
    )[0]
    assert float(many) == pytest.approx(float(one), abs=1e-6)
    assert float(many) < math.log(2)  # what L_out would floor at


def test_the_gray_band_enters_neither_the_numerator_nor_the_denominator():
    """`label == -1` is absent evidence, not a negative (`ML_PLAN.md` §3.1)."""
    labels = torch.tensor([[1, 0, -1]], dtype=torch.int8)
    scored = multi_positive_infonce(logits_for([[2.0, 0.0, 0.0]], [0.0]), labels)[0]
    # Moving the gray 8-mer's score anywhere must not change the loss.
    moved = multi_positive_infonce(logits_for([[2.0, 0.0, 50.0]], [0.0]), labels)[0]
    assert float(scored) == pytest.approx(float(moved), abs=1e-6)


def test_a_domain_with_no_positives_still_produces_a_gradient():
    """The defect the null anchor exists to fix: plain multi-positive InfoNCE averages over an
    empty set and the row contributes nothing at all (§3.2)."""
    logits = logits_for([[0.0, 0.0, 0.0]], [0.0])
    logits.requires_grad_(True)
    loss, per_domain = multi_positive_infonce(logits, torch.tensor([[0, 0, 0]], dtype=torch.int8))
    assert math.isfinite(float(per_domain[0]))
    loss.backward()
    assert logits.grad.abs().sum() > 0


def test_for_an_all_negative_row_the_null_is_the_target():
    """Scoring the null above every 8-mer is the correct answer for a domain that binds
    nothing, and must be the low-loss configuration."""
    labels = torch.tensor([[0, 0, 0]], dtype=torch.int8)
    good = multi_positive_infonce(logits_for([[-20.0, -20.0, -20.0]], [20.0]), labels)[0]
    bad = multi_positive_infonce(logits_for([[20.0, 20.0, 20.0]], [-20.0]), labels)[0]
    assert float(good) == pytest.approx(0.0, abs=1e-6)
    assert float(bad) > 20.0


def test_the_null_competes_with_the_positives_of_an_ordinary_row():
    """Its second job (§3.3): a real positive must be ranked above the null, so raising the
    null's score on a row that has positives must increase the loss."""
    labels = torch.tensor([[1, 0]], dtype=torch.int8)
    low = multi_positive_infonce(logits_for([[5.0, 0.0]], [-5.0]), labels)[0]
    high = multi_positive_infonce(logits_for([[5.0, 0.0]], [4.9]), labels)[0]
    assert float(high) > float(low)


def test_every_domain_counts_once_regardless_of_its_positive_count():
    """Positives per domain run 1 to 236, and without the `1/|P|` a 236-positive domain would
    outweigh 236 one-positive domains (§2.3c).

    The two rows here are built to differ in **positive count alone**: both have two negatives
    and one null, both score every positive at 1.0 and every negative at 0.0, and the second
    row's spare column is gray so that it does not enlarge the denominator.
    """
    labels = torch.tensor([[1, 1, 0, 0], [1, -1, 0, 0]], dtype=torch.int8)
    logits = logits_for([[1.0, 1.0, 0.0, 0.0], [1.0, 0.0, 0.0, 0.0]], [0.0, 0.0])
    _, per_domain = multi_positive_infonce(logits, labels)
    assert float(per_domain[0]) == pytest.approx(float(per_domain[1]), abs=1e-6)


def test_more_negatives_make_a_row_harder_even_at_the_same_positive_count():
    """The complement of the test above, so that it cannot pass by ignoring the denominator."""
    labels = torch.tensor([[1, 0, -1, -1], [1, 0, 0, 0]], dtype=torch.int8)
    logits = logits_for([[1.0, 0.0, 0.0, 0.0], [1.0, 0.0, 0.0, 0.0]], [0.0, 0.0])
    _, per_domain = multi_positive_infonce(logits, labels)
    assert float(per_domain[1]) > float(per_domain[0])


def test_a_labels_and_logits_shape_disagreement_is_an_error():
    """The null column is easy to forget, and forgetting it would silently treat the last
    8-mer as the anchor."""
    with pytest.raises(ValueError, match="null anchor"):
        multi_positive_infonce(torch.zeros(1, 3), torch.zeros(1, 3, dtype=torch.int8))


# --- the calibration term (T38) ---------------------------------------------------------------


def test_infonce_is_exactly_invariant_to_a_per_row_shift():
    """**The diagnosis the whole hybrid loss rests on.** A per-row softmax depends only on
    differences within its row, so adding a constant to a whole protein's scores changes nothing
    — which is why no threshold can be recovered from a model trained on this term alone, and
    why the null anchor never became one (`T38`)."""
    labels = torch.tensor([[1, 0, 0, -1]], dtype=torch.int8)
    base = logits_for([[3.0, 1.0, 0.5, 9.0]], [0.0])
    shifted = base + 7.5
    assert float(multi_positive_infonce(base, labels)[0]) == pytest.approx(
        float(multi_positive_infonce(shifted, labels)[0]), abs=1e-5
    )


def test_the_bce_term_is_not_invariant_to_a_per_row_shift():
    """The complement, and the reason binary cross-entropy is the term that was added: it is the
    cheapest objective that *does* pin where a row sits."""
    from snp2prot.models.loss import calibration_bce

    labels = torch.tensor([[1, 0, 0]], dtype=torch.int8)
    base = torch.tensor([[3.0, -1.0, -2.0]])
    assert float(calibration_bce(base, labels)[0]) != pytest.approx(
        float(calibration_bce(base + 5.0, labels)[0]), abs=1e-3
    )


def test_the_gray_band_is_masked_in_the_bce_as_well():
    """`label == -1` is absent evidence in both terms, or the calibration head would be fitted
    against cells the assay never called."""
    from snp2prot.models.loss import calibration_bce

    labels = torch.tensor([[1, 0, -1]], dtype=torch.int8)
    here = calibration_bce(torch.tensor([[2.0, -2.0, 0.0]]), labels)[0]
    moved = calibration_bce(torch.tensor([[2.0, -2.0, 40.0]]), labels)[0]
    assert float(here) == pytest.approx(float(moved), abs=1e-6)


def test_the_bce_weights_every_domain_once_regardless_of_its_width():
    """Same normalisation as the InfoNCE term and as the reported metric (§2.3c): a domain with
    more scored cells must not count for more."""
    from snp2prot.models.loss import calibration_bce

    labels = torch.tensor([[1, 0, -1, -1], [1, 0, 0, 0]], dtype=torch.int8)
    logits = torch.tensor([[2.0, -2.0, 0.0, 0.0], [2.0, -2.0, -2.0, -2.0]])
    _, per_domain = calibration_bce(logits, labels)
    assert float(per_domain[0]) == pytest.approx(float(per_domain[1]), abs=1e-6)


def test_lambda_zero_is_the_pure_ranking_loss_to_the_last_bit():
    """What makes this a one-knob ablation rather than a redesign: at `λ = 0` the objective is
    the one that existed before the calibration term, to the last bit.

    This is a statement about the *loss*, not about a training run. GPU training is not
    reproducible across processes — measured 2026-08-28, identical code and seeds gave losses
    8.96644592285 and 8.96644687653 at step 50 — so a `λ = 0` run does not reproduce an older
    run's numbers and the sweep re-runs its own control."""
    from snp2prot.models.loss import hybrid_loss

    labels = torch.tensor([[1, 0, 0]], dtype=torch.int8)
    logits = logits_for([[2.0, 0.0, -1.0]], [-3.0])
    calibrated = torch.tensor([[50.0, 50.0, 50.0]])  # would be a huge BCE if it counted
    total, parts = hybrid_loss(logits, calibrated, labels, 0.0)
    assert float(total) == float(multi_positive_infonce(logits, labels)[0])
    assert float(parts["bce"]) == 0.0


def test_a_positive_lambda_adds_exactly_lambda_times_the_bce():
    from snp2prot.models.loss import calibration_bce, hybrid_loss

    labels = torch.tensor([[1, 0, 0]], dtype=torch.int8)
    logits = logits_for([[2.0, 0.0, -1.0]], [-3.0])
    calibrated = torch.tensor([[1.0, -4.0, -5.0]])
    total, parts = hybrid_loss(logits, calibrated, labels, 30.0)
    expected = (
        multi_positive_infonce(logits, labels)[0] + 30.0 * calibration_bce(calibrated, labels)[0]
    )
    assert float(total) == pytest.approx(float(expected), abs=1e-6)
    assert float(parts["infonce"]) == pytest.approx(
        float(multi_positive_infonce(logits, labels)[0]), abs=1e-6
    )


def test_the_null_anchor_column_must_be_dropped_before_the_bce():
    """The anchor carries no label, so giving it a target would invent evidence. The shape check
    is the only thing standing between that and a silent off-by-one across the whole vocabulary."""
    from snp2prot.models.loss import calibration_bce

    with pytest.raises(ValueError, match="null anchor"):
        calibration_bce(torch.zeros(1, 4), torch.zeros(1, 3, dtype=torch.int8))
