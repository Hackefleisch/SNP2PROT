"""The C1 evaluation set: the held-out variants where copying the wild type fails.

**Why the set has to exist** (`D6`, decided 2026-08-19). Claim C1 is that a model is sensitive
to single-residue change. The regime that tests it is `P3/all`, which holds out every variant
and leaves every wild type in training — so the nearest-neighbour baseline copies each held-out
variant's own wild type, and *is* the hypothesis "the mutation has no effect".

That baseline scores **mean AUPR 0.928 and median 1.000**. The median is the problem: for most
of the 173 variants the wild-type profile simply is the right answer, because most single
substitutions do not measurably change which 8-mers a domain binds. A model could be perfect
and gain two points on the mean. **The claim does not live in the mean, it lives in the
variants where the wild type is wrong** — and pooled with the rest they are invisible.

## What the set is

Held-out variants under `P3/all` in one of two states:

- **the wild-type copy scores badly** — baseline AUPR below `c1_set.max_baseline_aupr`;
- **the wild-type copy scores nothing at all** — the variant has no positive 8-mer, so its AUPR
  is undefined rather than zero. These are `label_health`'s `dead_variant` records: a variant
  that measurably lost binding, with a control in its own series (`T21`, `T30`). They are the
  extreme of the same tail and the sharpest C1 evidence in the corpus.

`no_evidence` records are excluded, as everywhere: with no control there is no way to tell a
true non-binder from an assay that failed.

## What it is not

**It is not a test the baseline can be said to fail.** The set is *defined* by the baseline
scoring badly on it, so the baseline scores badly on it by construction and "the model beats
the baseline here" is close to vacuous. Two things make it meaningful instead, and both belong
on the figure:

1. **the model's absolute number on the set** — can it recover these profiles at all;
2. **the model's number on the complement** — a model that learned to distrust wild types
   everywhere would gain here and lose there, which is not C1 but a shifted prior.

The selection uses only measured data — a variant's own labels against its wild type's measured
E-scores — and no model output, so it is fixed before any model runs and does not move between
arms.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from snp2prot import experiment, label_health

#: Why a variant is in the set. `dead` is the stronger case: not a poorly predicted profile but
#: no profile at all to predict.
LOST_BINDING = "dead"
POORLY_PREDICTED = "poorly_predicted"


def select(
    per_domain: pd.DataFrame,
    max_baseline_aupr: float | None = None,
    keep_no_evidence: bool | None = None,
) -> pd.DataFrame:
    """The C1 set, from `run_nn_baseline.py`'s per-domain frame.

    `per_domain` needs the `P3/all` rows at `k == 1`: one row per held-out variant with its
    `aupr`, `verdict` and `spearman`. Returns those rows plus a `c1_reason` column, ordered
    worst-predicted first.
    """
    cfg = experiment.section("c1_set")
    if max_baseline_aupr is None:
        max_baseline_aupr = float(cfg["max_baseline_aupr"])
    if keep_no_evidence is None:
        keep_no_evidence = bool(cfg["keep_no_evidence"])

    rows = per_domain[(per_domain.k == 1) & (per_domain.fold == "all")].copy()
    if not keep_no_evidence:
        rows = rows[rows.verdict != label_health.NO_EVIDENCE]

    undefined = ~np.isfinite(rows.aupr)
    poor = rows.aupr < max_baseline_aupr
    selected = rows[undefined | poor].copy()
    selected["c1_reason"] = np.where(~np.isfinite(selected.aupr), LOST_BINDING, POORLY_PREDICTED)
    # Worst first, and the ones with no AUPR at all before the ones that merely score low.
    return selected.sort_values(["aupr", "spearman"], na_position="first").reset_index(drop=True)


def summarise(selected: pd.DataFrame, all_variants: pd.DataFrame) -> pd.DataFrame:
    """Per family: how many variants there are, and how many of them the wild type fails on.

    Stratified because 84 of the 173 variants are homeodomain, so a pooled number is a
    homeodomain number wearing a general claim's clothes (`docs/ML_PLAN.md` §6.1).
    """
    total = all_variants.groupby("family").size().rename("n_variants")
    chosen = selected.groupby("family").size().rename("n_in_set")
    dead = selected[selected.c1_reason == LOST_BINDING].groupby("family").size().rename("n_dead")
    out = pd.concat([total, chosen, dead], axis=1).fillna(0).astype(int)
    out["fraction"] = out.n_in_set / out.n_variants
    return out.sort_values(["n_in_set", "n_variants"], ascending=False).reset_index()
