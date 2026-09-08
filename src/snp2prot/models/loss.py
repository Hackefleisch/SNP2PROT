"""Multi-positive InfoNCE over the complete 8-mer axis, with a null anchor.

The full derivation and the reasons are `docs/TRAINING.md` §2 and §3. In brief, for a domain
`p` with positives `P_p`, negatives `N_p`, and the null anchor `ν`:

```
             1                    exp s(p, k⁺)
L_p  =  − ───────  ·   Σ    log ──────────────────────────────────────
            |P_p|    k⁺ ∈ P_p    exp s(p, k⁺) + Σ_{k ∈ N_p ∪ {ν}} exp s(p, k)
```

Three properties, each load-bearing:

**The denominator is exact.** Every InfoNCE you will read approximates the negative set because
in the usual setting it is unbounded. Here the complete label matrix is resident — 1,338 ×
32,896, every cell present — so sampling a negative set we already know buys nothing and costs
bias.

**Sibling positives are not in the denominator.** This departs from SupCon's `L_out`, and the
reason is concrete: with them in, a 236-positive domain has a loss that floors at `log|P_p|` and
a gradient that pushes its own positives apart. That is right for representation learning and
wrong for retrieval. Excluded, the loss reaches zero exactly when every positive outranks every
negative — which is AUPR = 1, so the objective and the metric agree at their optimum.

**The gray band is masked, not filtered.** `label == -1` is absent evidence, not a negative, so
it enters neither numerator nor denominator (`ML_PLAN.md` §3.1: two masks, not filters).

## The null anchor, and how it unifies the empty case

A domain with no positives should push *every* 8-mer away from itself, but the loss above
averages over `P_p` and an empty average is vacuous — the row silently contributes no gradient.
The fix is to make `ν` the target when `P_p` is empty, which drives every real 8-mer's score down
relative to a learned reference point. Writing `P̃_p = P_p` if non-empty else `{ν}` covers both
cases in one expression, so there is no branch on "does this domain have positives".

For a domain that *does* have positives, `ν` is one more competitor in the denominator.

**It was expected to settle into a decision threshold, and it does not.** Measured on a trained
checkpoint (2026-08-26): the anchor sits at −24.3 where the negatives' median is −26.4 and the
positives' median is +9.6, and a median of **15,159 of 32,460 8-mers score above it**. It lands
inside the negative cloud, not between the two classes. The reason is structural — for the ~1,318
domains that have positives it is simply another negative being pushed down, and the only force
raising it comes from the 20 domains with none, where it is the target. Nothing downstream may
treat it as an operating point; `snp2prot.evaluation.metrics.suppression` exists because of this.

What `ν` does do is the paragraph above: it gives a domain with no positives a target, so such a
row contributes a gradient instead of nothing at all. That job is real and is why it stays.

## The decision rule lives in a second term, not in the anchor

`calibration_bce` below, and `hybrid_loss` which adds it — `L = L_infonce + λ · L_bce`. The
diagnosis that motivates it is one line: **the loss above is exactly invariant to adding a
constant to a whole row**, being a softmax over that row, so nothing in it ever says where a
protein's scores should sit relative to another protein's. A per-protein threshold cannot be
recovered afterwards from a quantity training never constrained, which is why the anchor failed
and why global Platt scaling on the validation slice failed after it (`T38`). Binary
cross-entropy is not shift-invariant; `λ = 0` recovers this file's original behaviour exactly.

## Computed in O(K), not O(|P| · K)

The negative part of the denominator is shared across all of a domain's positives, so with
`logC_p = logsumexp` over `N_p ∪ {ν}` the whole thing collapses to

```
L_p  =  mean over k⁺ ∈ P̃_p  of  softplus( logC_p − s(p, k⁺) )
```

One masked `logsumexp` per row, one `softplus` over the positives.
"""

from __future__ import annotations

import torch
from torch import nn

#: Label values, matching `snp2prot.schema`.
BIND, NONBIND, GRAY = 1, 0, -1


def multi_positive_infonce(
    logits: torch.Tensor, labels: torch.Tensor
) -> tuple[torch.Tensor, torch.Tensor]:
    """The loss, and the per-domain values behind it.

    `logits` is `(B, K + 1)` — every 8-mer plus the null anchor in the final column — and
    `labels` is `(B, K)` in `{1, 0, -1}`. Returns `(mean, per_domain)`; the per-domain vector is
    what a caller logs a distribution of, since a mean over proteins hides which proteins fail.
    """
    if logits.shape[1] != labels.shape[1] + 1:
        raise ValueError(
            f"logits has {logits.shape[1]} columns and labels {labels.shape[1]}; "
            "the logits must carry exactly one extra column for the null anchor"
        )

    positive = labels == BIND
    negative = labels == NONBIND
    has_positive = positive.any(dim=1, keepdim=True)

    # Targets: a domain's own positives, or the null anchor when it has none.
    targets = torch.cat([positive, ~has_positive], dim=1)
    # Denominator: the negatives, plus the null whenever the null is not itself the target.
    pool = torch.cat([negative, has_positive], dim=1)

    masked = logits.masked_fill(~pool, float("-inf"))
    log_c = torch.logsumexp(masked, dim=1, keepdim=True)

    # softplus(logC - s) = -log( exp(s) / (exp(s) + exp(logC)) ), the per-positive term.
    per_positive = nn.functional.softplus(log_c - logits)
    n_targets = targets.sum(dim=1).clamp(min=1)
    per_domain = (per_positive * targets).sum(dim=1) / n_targets
    return per_domain.mean(), per_domain


def calibration_bce(
    calibrated: torch.Tensor, labels: torch.Tensor
) -> tuple[torch.Tensor, torch.Tensor]:
    """Masked binary cross-entropy over the real 8-mers: the term that pins the row offsets.

    **This exists because `multi_positive_infonce` is exactly invariant to a per-domain score
    offset, and that is why the model has no decision rule** (`T38`). Its per-domain loss is a
    softmax over one row, so it depends only on `logC_p - s(p, k⁺)` — differences *within* a
    row. Add a constant to every entry of a row and the loss does not change by a float. The
    objective therefore never says where a row should sit, only how it should be ordered, and
    no post-hoc calibration can recover an offset the training signal never constrained. Measured
    consequence: a globally-calibrated cut called between 56 and 408 8-mers per domain against a
    true median of 44, and the null anchor — designed to be the threshold — settled inside the
    negative cloud.

    Binary cross-entropy is not shift-invariant. It is the cheapest term with that property, and
    it produces the quantity the deliverable actually needs: `sigmoid(calibrated)` is a
    probability that *this* protein binds *this* 8-mer, comparable across proteins, so a protein
    that binds nothing is one whose probabilities are all low rather than one whose ranking has
    to be interpreted.

    **Unweighted, deliberately.** At 466:1 the negatives are 99.8% of every row and the usual
    reflex is to class-weight or to reach for focal loss. Both re-bias the output away from the
    empirical rate, which is precisely the thing being asked for here — a weighted BCE is
    calibrated to a class balance that does not exist. The ranking is already supplied by the
    InfoNCE term; this term's only job is to be honest about the base rate, so it is left alone
    and `λ` (`training.bce_weight`) sets how much it counts.

    **Per-domain mean, then mean over domains**, matching §2.3(c): every protein counts once
    whatever its positive count, exactly as the InfoNCE term and the reported metric do. Gray
    cells (`label == -1`) are masked out of both the sum and the count, so a domain with a wide
    no-call band is not quietly averaged over cells it has no evidence for.

    `calibrated` is `(B, K)` — **no null-anchor column**. The anchor is not a labelled 8-mer, so
    it has no target here and taking one would invent evidence.
    """
    if calibrated.shape != labels.shape:
        raise ValueError(
            f"calibrated logits are {tuple(calibrated.shape)} and labels {tuple(labels.shape)}; "
            "the BCE term is over the real 8-mers only, so the null anchor's column must be "
            "dropped before calling this"
        )
    scored = labels != GRAY
    target = (labels == BIND).to(calibrated.dtype)
    per_cell = nn.functional.binary_cross_entropy_with_logits(calibrated, target, reduction="none")
    counts = scored.sum(dim=1).clamp(min=1)
    per_domain = (per_cell * scored).sum(dim=1) / counts
    return per_domain.mean(), per_domain


def hybrid_loss(
    logits: torch.Tensor,
    calibrated: torch.Tensor,
    labels: torch.Tensor,
    bce_weight: float,
) -> tuple[torch.Tensor, dict[str, torch.Tensor]]:
    """`L = L_infonce + λ · L_bce`, and the two terms behind it.

    **`λ = 0` is the pure-ranking objective exactly** — the calibration head takes no gradient at
    all and the returned total *is* the InfoNCE value, to the last bit. That is what makes this a
    one-knob ablation rather than a redesign (`docs/DECISIONS.md` §13). It does not make an old
    run's *numbers* reproducible: GPU training is not deterministic across processes, so a `λ = 0`
    control has to be re-run alongside whatever it is controlling for.

    The two terms are on very different scales and that is expected rather than a defect. At
    initialisation the InfoNCE term is about `log K ≈ 10.4` while the BCE term, predicting the
    base rate, is about `0.014` — the entropy of a 0.0015 Bernoulli. So the useful `λ` grid is
    logarithmic and centred near their ratio, ~10³, not near 1. Both terms are returned
    separately and logged separately, so a run says which of them it was actually minimising.
    """
    infonce, _ = multi_positive_infonce(logits, labels)
    if not bce_weight:
        return infonce, {"infonce": infonce.detach(), "bce": torch.zeros((), device=logits.device)}
    bce, _ = calibration_bce(calibrated, labels)
    return infonce + bce_weight * bce, {"infonce": infonce.detach(), "bce": bce.detach()}
