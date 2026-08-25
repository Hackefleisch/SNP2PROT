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
