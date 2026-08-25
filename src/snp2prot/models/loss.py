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

For a domain that *does* have positives, `ν` is one more competitor in the denominator, and it
learns a position meaning "good enough to call binding" — a global decision threshold, obtained
for free (`docs/TRAINING.md` §3.3).

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
