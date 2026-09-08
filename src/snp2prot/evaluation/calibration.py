"""Decision-rule metrics: does a calibrated probability mean the same thing across proteins?

**Everything in `snp2prot.evaluation.metrics` is rank-based, and that is the gap this module
exists to close** (`T38`, raised 2026-08-28). AUPR, AUROC, precision@k and R@P all score an
*ordering*, so none of them ever requires the model to commit to a call — while the deliverable
is a call, and the capability that matters most is noticing that a variant has lost binding
entirely, which is a statement about a protein and not about its internal ranking.

Three questions, three families of metric here:

1. **Is the probability honest?** `expected_calibration_error` — of the cells the model calls
   `p ≈ 0.3`, do about 30% of them bind.
2. **Is the call good?** `score_calls` — precision, recall, F1 and Jaccard of the set the model
   actually names, against the set the assay names. This is the first metric in the project a
   nearest-neighbour lookup and a neural model can be compared on *as decisions*: `k = 1` emits a
   set natively, which is exactly why its AUPR was never a like-for-like number.
3. **Did the interaction survive?** `interaction_power` and `detection_auroc` — a per-protein
   scalar, and whether it separates the variants that bind nothing from the ones that still bind.

**The threshold is not chosen here.** `score_calls` takes one; `expected_count_rule` derives one
per protein from the model's own probabilities, which is the rule with no free parameter at all.
Which of them to report is a decision for `docs/DECISIONS.md`, not for a metrics module.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from snp2prot.evaluation.metrics import auroc

#: Quantile bins for the calibration curve, not equal-width ones. At a 0.0021 base rate an
#: equal-width binning puts 99.9% of the cells in `[0, 0.067)` and reports one number for the
#: whole negative cloud, which is the part being asked about. Equal-mass bins spend their
#: resolution where the data is.
CALIBRATION_BINS = 15

#: The naive fixed probability cut, reported beside `expected_count_rule` so that what choosing
#: a threshold costs is visible rather than argued about.
DEFAULT_CALL_THRESHOLD = 0.5


def expected_calibration_error(
    labels: np.ndarray, probabilities: np.ndarray, bins: int = CALIBRATION_BINS
) -> float:
    """Mean |predicted - observed| over equal-mass bins of the predicted probability.

    The standard ECE with quantile binning. 0 is perfect; the value is on the probability scale,
    so 0.01 means the average cell's stated probability is off by a hundredth.

    **Read it beside `score_calls`, never alone.** At 466:1 a model that predicts the base rate
    for every single cell — the constant predictor, which is useless — has an ECE near zero. ECE
    says the numbers are not lies; it does not say they are informative.
    """
    labels = np.asarray(labels, dtype=np.float64)
    probabilities = np.asarray(probabilities, dtype=np.float64)
    if not len(labels):
        return float("nan")
    order = np.argsort(probabilities, kind="stable")
    total = 0.0
    for chunk in np.array_split(order, min(bins, len(order))):
        if not len(chunk):
            continue
        total += len(chunk) * abs(probabilities[chunk].mean() - labels[chunk].mean())
    return float(total / len(labels))


@dataclass(frozen=True)
class CallScores:
    """One protein's predicted binding set, scored against the measured one."""

    n_called: int
    n_true: int
    precision: float
    recall: float
    f1: float
    jaccard: float
    threshold: float

    def as_dict(self) -> dict[str, float]:
        return {
            "n_called": float(self.n_called),
            "n_true": float(self.n_true),
            "call_precision": self.precision,
            "call_recall": self.recall,
            "call_f1": self.f1,
            "call_jaccard": self.jaccard,
            "call_threshold": self.threshold,
        }


def score_calls(labels: np.ndarray, probabilities: np.ndarray, threshold: float) -> CallScores:
    """Precision, recall, F1 and Jaccard of the set `probabilities >= threshold`.

    **Precision is `nan` when nothing is called and recall is `nan` when nothing binds**, and the
    two absences are different: calling nothing is a decision the model made, while a protein
    with no positive 8-mer is a property of the assay. F1 and Jaccard are 0 when the model called
    something and hit nothing, which is a failure with a value — the distinction
    `metrics.recall_at_precision` had to learn (see its docstring). Both are `nan` only when the
    truth set and the called set are *both* empty, where every set metric is genuinely undefined
    and the right question is `interaction_power` instead.
    """
    called = np.asarray(probabilities) >= threshold
    truth = np.asarray(labels) == 1
    hits = int((called & truth).sum())
    n_called, n_true = int(called.sum()), int(truth.sum())
    union = n_called + n_true - hits
    return CallScores(
        n_called=n_called,
        n_true=n_true,
        precision=float(hits / n_called) if n_called else float("nan"),
        recall=float(hits / n_true) if n_true else float("nan"),
        f1=float(2 * hits / (n_called + n_true)) if (n_called + n_true) else float("nan"),
        jaccard=float(hits / union) if union else float("nan"),
        threshold=float(threshold),
    )


def interaction_power(probabilities: np.ndarray) -> float:
    """The expected number of 8-mers this protein binds: `Σ p`.

    **The protein-level statistic, and the one the pathogenic-variant question needs.** A
    calibrated probability makes this a count with units rather than a score: a domain whose
    every cell is near the base rate has an interaction power near `K · rate ≈ 69`, and one the
    model believes binds nothing has one near 0.

    It is the sum and not the maximum on purpose. `max p` answers "is there any site at all",
    which a single confident false positive destroys; the sum is the expectation of a count, is
    what the calibration term is actually fitting, and degrades gracefully. `max p` is reported
    beside it by the callers as a second opinion, not as the headline.
    """
    return float(np.asarray(probabilities, dtype=np.float64).sum())


def expected_count_rule(probabilities: np.ndarray) -> float:
    """The self-consistent per-protein threshold: keep the top `round(Σ p)` 8-mers.

    **A decision rule with no free parameter.** If the probabilities are calibrated then `Σ p` is
    the expected number of binders, so a model that both ranks well and calibrates well already
    contains its own answer to "how many" — and taking exactly that many converts the ranking
    into a set without anyone choosing a cut. It also cannot produce the failure `T38` recorded
    for a global cut, where one threshold called 56 8-mers for one protein and 408 for another
    whose true counts were comparable: here the count *is* the per-protein quantity.

    Returns the probability at the cut, so it can be passed straight to `score_calls`. A protein
    whose expected count rounds to 0 gets `inf` — it calls nothing, which is the correct output
    for a domain the model believes binds nothing.

    **`inf` rather than "one ulp above the maximum", which was wrong and silently so.** Under
    NumPy's weak scalar promotion a Python float compared against a `float32` array is cast to
    `float32`, and `nextafter(0.0, inf)` — a denormal double — becomes exactly `0.0` there. The
    cut meant to call *nothing* then called **everything**: measured on the nearest-neighbour
    baseline, a held-out domain whose neighbour binds no 8-mer was reported as 32,894 calls
    rather than 0. `inf` survives any downcast and cannot invert.
    """
    probabilities = np.asarray(probabilities, dtype=np.float64)
    k = int(round(interaction_power(probabilities)))
    if k <= 0:
        return float("inf")
    if k >= len(probabilities):
        return float(probabilities.min())
    return float(np.partition(probabilities, len(probabilities) - k)[len(probabilities) - k])


def detection_auroc(is_dead: np.ndarray, power: np.ndarray) -> float:
    """AUROC for calling a protein dead from its interaction power alone.

    1.0 means every variant that binds nothing has a lower predicted interaction power than every
    variant that still binds. 0.5 is no information.

    **This is the metric `T38` asked for and the one `metrics.suppression` could not be.**
    Suppression compares a variant against its own wild type down the protein axis, which needs
    no threshold but also needs the wild type — so it is undefined for any protein presented on
    its own, which is every protein a user would actually bring. This asks the question the user
    asks: given this domain and nothing else, has the interaction gone.

    Reported with its `n`. On this corpus that `n` is small — 18 dead variants pass
    `label_health` — so a value here is a direction and not a result, and the callers print the
    count beside it for exactly that reason.
    """
    dead = np.asarray(is_dead).astype(np.int64)
    power = np.asarray(power, dtype=np.float64)
    keep = np.isfinite(power)
    if keep.sum() < 2 or not (0 < dead[keep].sum() < keep.sum()):
        return float("nan")
    # Dead is the positive class and low power is the evidence for it, so the score is negated.
    return auroc(dead[keep], -power[keep])


def power_ratio(variant_power: float, reference_power: float) -> float:
    """The variant's predicted interaction power as a fraction of its wild type's.

    **The paired form of `detection_auroc`, and the stronger read at this sample size.** 0.0 says
    the model predicts the mutation abolished binding; 1.0 says it predicts no effect, which is
    what copying the wild type does and therefore what the nearest-neighbour baseline scores by
    construction under `P3/all`.

    `nan` when the reference is predicted to bind nothing either — the ratio is then a
    statement about the model's opinion of the wild type, not about the mutation.
    """
    if not np.isfinite(reference_power) or reference_power <= 0:
        return float("nan")
    return float(variant_power / reference_power)
