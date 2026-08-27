#!/usr/bin/env python
"""Train every arm on every fold, and write the phase deliverable.

    python scripts/run_grid.py [--arm A1 --arm A4] [--regime S2] [--out reports/training.md]

Phase 7, build step 1. **19 folds per arm** (`docs/TRAINING.md` §6.1): 5 `S1` + 5 `S2` + `P1` +
`P2` + 7 `P3`, which are exactly the folds `snp2prot.splits.all_regimes` produces and exactly
the ones `reports/nn_baseline.md` scored. With both sequence arms built that is **38 runs**, at
roughly a minute and a half each.

The grid is *enumerated*, not searched — which is why `ML_PLAN.md` §9.2 rejected Hydra and
Optuna and adopted a loop over configs.

**Every number here is reported against the nearest-neighbour baseline on the same fold.** A
degradation curve is not diagnostic on its own: a family-holdout AUPR of 0.22 means the model is
a lookup table if the baseline scores 0.21, and is the talk if the baseline scores 0.05
(`ML_PLAN.md` §8.1). The baseline's per-fold numbers are read from
`data/processed/nn_baseline_domains.parquet`, so the comparison is always on the same held-out
domains.

⚠️ **This takes about an hour for two arms and is the owner's to run** (`CLAUDE.md` rule 10).
`scripts/train.py --fold P3/all` is the fast single-fold check.
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np
import pandas as pd
from train import add_common_arguments, apply_overrides, load_everything

from snp2prot import experiment, splits, tracking, training
from snp2prot.config import PROCESSED_DIR, REPORTS_DIR

DEFAULT_OUT = REPORTS_DIR / "training.md"
SUMMARY_OUT = PROCESSED_DIR / "training_folds.parquet"
PER_DOMAIN_OUT = PROCESSED_DIR / "training_domains.parquet"
BASELINE = PROCESSED_DIR / "nn_baseline_domains.parquet"


def baseline_by_fold() -> pd.DataFrame:
    """The NN baseline's macro AUPR per fold, for the horizontal line on every figure."""
    if not BASELINE.exists():
        raise SystemExit(
            f"no baseline at {BASELINE}\nrun scripts/run_nn_baseline.py — every number in "
            "this report is read against it (ML_PLAN.md §8.1)"
        )
    frame = pd.read_parquet(BASELINE)
    frame = frame[frame.k == 1]
    grouped = frame.groupby(["regime", "fold"], sort=False)
    return pd.DataFrame(
        {
            "baseline_aupr": grouped.aupr.mean(),
            "baseline_median": grouped.aupr.median(),
        }
    ).reset_index()


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    add_common_arguments(ap)
    ap.add_argument("--arm", action="append", dest="arms", help="repeatable; default A1 and A4")
    ap.add_argument("--regime", action="append", help="restrict to these regimes")
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = ap.parse_args()
    arms = args.arms or ["A1", "A4"]

    config = apply_overrides(experiment.load(), args)
    device = training.device_for(args.device)
    started = time.time()

    summaries: list[dict] = []
    per_domain: list[pd.DataFrame] = []
    for arm in arms:
        domains, matrix, vectors, dist, trainable = load_everything(arm)
        folds = [f for f in splits.all_regimes(domains, dist)]
        if args.regime:
            folds = [f for f in folds if f.regime in set(args.regime)]
        print(f"{arm} ({vectors.model}): {len(folds)} folds")
        for fold in folds:
            summary, rows = training.run_fold(
                fold,
                arm,
                domains,
                matrix,
                vectors,
                dist,
                config,
                trainable,
                device=device,
                track=not args.no_track,
            )
            summaries.append(summary)
            per_domain.append(rows)
            print(
                f"  {fold.label:<20} AUPR {summary['aupr']:.4f}  "
                f"val {summary['validation_aupr']:.4f}  "
                f"{int(summary['steps_run'])} steps  ({time.time() - started:.0f}s)",
                flush=True,
            )

    frame = pd.DataFrame(summaries).merge(baseline_by_fold(), on=["regime", "fold"], how="left")
    frame["delta"] = frame.aupr - frame.baseline_aupr
    domains_frame = pd.concat(per_domain, ignore_index=True)

    SUMMARY_OUT.parent.mkdir(parents=True, exist_ok=True)
    frame.to_parquet(SUMMARY_OUT, index=False)
    domains_frame.to_parquet(PER_DOMAIN_OUT, index=False)
    write_report(args.out, frame, domains_frame, config)
    print(f"wrote {args.out}, {SUMMARY_OUT} and {PER_DOMAIN_OUT} in {time.time() - started:.0f}s")


def _fmt(value: float, places: int = 4) -> str:
    return "n/a" if not np.isfinite(value) else f"{value:.{places}f}"


def _lift(value: float, chance: float) -> str:
    """An AUPR as a multiple of the fold's random-ranking null.

    A decimal below 10x, because that is where the number decides something: `A4` on `P1` is
    1.4x and rendering it as `1x` reads as "exactly chance" when the null test says otherwise.
    Above 10x the decimal is noise.
    """
    if not (np.isfinite(value) and np.isfinite(chance)) or chance <= 0:
        return "n/a"
    lift = value / chance
    return f"{lift:.0f}x" if lift >= 10 else f"{lift:.1f}x"


def _selection_section(frame: pd.DataFrame) -> list[str]:
    """Did selecting on validation beat just taking the model at the step budget?

    Every run stores both models, so this is measured per fold rather than assumed. It is the
    only readout of whether the validation slice is a useful selector — on `S2` it may not be,
    since a slice carved from the training pool cannot imitate an unseen-component test task.
    """
    gain = frame.selection_gain
    lines = [
        "",
        "## Does selecting on validation beat the model at the budget?",
        "",
        "`model AUPR` is the validation-selected checkpoint and `at budget` the model at step",
        "`training.steps`. Both are stored with every run, so the question is measured rather",
        "than assumed. A negative mean would say the validation slice is selecting worse than",
        "not selecting at all.",
        "",
        "| regime | folds | selected | at budget | mean gain | folds where selection lost |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for regime, group in frame.groupby("regime", sort=False):
        lines.append(
            f"| {regime} | {len(group)} | {_fmt(group.aupr.mean())} | "
            f"{_fmt(group.aupr_final.mean())} | {group.selection_gain.mean():+.4f} | "
            f"{int((group.selection_gain < 0).sum())} |"
        )
    lines += ["", f"Overall mean gain from selecting: **{gain.mean():+.4f}**.", ""]
    return lines


def _null_section(frame: pd.DataFrame, repeats: int) -> list[str]:
    """Which results are, and are not, distinguishable from a random ranking.

    A multiple near 1 is an eyeball, not a test. `metrics.null_aupr` samples the macro statistic
    under a random ranking of the same held-out labels, and a result at or below the 95th
    percentile of that null is reported as **at chance** rather than as a small number.
    """
    at_chance = frame[frame.aupr <= frame.null_p95]
    lines = [
        "",
        "## Is it better than random?",
        "",
        "The null is the macro AUPR of a *random* ranking of the same held-out domains, sampled",
        f"{repeats} times (`snp2prot.evaluation.metrics.null_aupr`). A result at or below the",
        "95th percentile of that null is indistinguishable from guessing, whatever its delta",
        "against the baseline looks like.",
        "",
    ]
    if not len(at_chance):
        return lines + ["Every fold scored above its null. ✅", ""]

    lines += [
        f"⚠️ **{len(at_chance)} of {len(frame)} run{'s' if len(at_chance) != 1 else ''}"
        f"{' are' if len(at_chance) != 1 else ' is'} at chance.**",
        "",
        "| arm | regime | fold | model AUPR | null mean | null 95th | verdict |",
        "|---|---|---|---:|---:|---:|---|",
    ]
    for row in at_chance.itertuples():
        lines.append(
            f"| `{row.arm}` | {row.regime} | `{row.fold}` | {_fmt(row.aupr)} | "
            f"{_fmt(row.null_mean)} | {_fmt(row.null_p95)} | **at chance** |"
        )
    return lines + [""]


def _provenance_line(frame: pd.DataFrame) -> str:
    """Which commit produced these rows, and whether the tree was modified.

    A digest pins the held-out domains but not the procedure: `OVERSHOOT` landed in the same
    commit as the first version of this report, so one row of it came from a `validation_split`
    that no longer existed (`snp2prot.tracking.code_version`). Every row carries its own stamp,
    so a table assembled from more than one state of the tree says so rather than looking whole.
    """
    stamps = sorted({(r.code_commit, r.code_dirty) for r in frame.itertuples()})
    rendered = ", ".join(tracking.stamp(c, d) for c, d in stamps)
    if len(stamps) > 1:
        return (
            f"⚠️ **These rows were not all produced by the same code**: {rendered}. "
            "Re-run the grid before reading the table as one experiment."
        )
    return f"Produced from {rendered}."


def write_report(path: Path, frame: pd.DataFrame, per_domain: pd.DataFrame, config: dict) -> None:
    lines = [
        "# Two-tower contrastive model — results",
        "",
        "Generated by `scripts/run_grid.py`. The design is",
        "[`docs/TRAINING.md`](../docs/TRAINING.md):",
        "multi-positive InfoNCE over the complete 8-mer axis with a learned null anchor, batches",
        "drawn on the protein axis with every domain scored against all 32,896 8-mers.",
        "",
        "**Every number is against the nearest-neighbour baseline on the same held-out domains**",
        "([`nn_baseline.md`](nn_baseline.md)). That comparison is the whole point: a degradation",
        "curve is not diagnostic on its own, because every imperfectly-generalising model has one.",
        "What separates a model that learned biophysics from a lookup table is the *level*.",
        "",
        "## Per fold",
        "",
        "| arm | regime | fold | model AUPR | baseline | delta | median | AUROC "
        "| Spearman | steps |",
        "|---|---|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in frame.itertuples():
        lines.append(
            f"| `{row.arm}` | {row.regime} | `{row.fold}` | **{_fmt(row.aupr)}** | "
            f"{_fmt(row.baseline_aupr)} | {row.delta:+.4f} | {_fmt(row.aupr_median)} | "
            f"{_fmt(row.auroc, 3)} | {_fmt(row.spearman, 3)} | {int(row.steps_run)} |"
        )

    lines += [
        "",
        "## Per regime, averaged over folds",
        "",
        "| arm | regime | folds | model AUPR | baseline | delta |",
        "|---|---|---:|---:|---:|---:|",
    ]
    for (arm, regime), group in frame.groupby(["arm", "regime"], sort=False):
        lines.append(
            f"| `{arm}` | {regime} | {len(group)} | **{_fmt(group.aupr.mean())}** | "
            f"{_fmt(group.baseline_aupr.mean())} | {group.delta.mean():+.4f} |"
        )

    lines += [
        "",
        "## The C1 evaluation set",
        "",
        "The 29 variants whose binding measurably changed, so that their own wild type does not",
        "predict them (`D6`). They are never trained on. **The baseline scores badly on them by",
        'construction — the set is defined that way — so "the model beats the baseline here" is',
        "vacuous.** What is readable is the model's absolute number, and its number on the",
        "complement: a model that had merely learned to distrust wild types everywhere would gain",
        "here and lose there, which is a shifted prior and not C1.",
        "",
        "| arm | set | n | AUPR | median |",
        "|---|---|---:|---:|---:|",
    ]
    p3 = per_domain[(per_domain.regime == "P3") & (per_domain.fold == "all")]
    if len(p3):
        wanted = set(pd.read_parquet(PROCESSED_DIR / "c1_variants.parquet").domain)
        for arm, group in p3.groupby("arm", sort=False):
            inside = group[group.domain.isin(wanted)]
            outside = group[~group.domain.isin(wanted)]
            for name, part in (("C1 set", inside), ("the other variants", outside)):
                lines.append(
                    f"| `{arm}` | {name} | {len(part)} | {_fmt(np.nanmean(part.aupr))} | "
                    f"{_fmt(np.nanmedian(part.aupr))} |"
                )

    lines += [
        "",
        "## Configuration",
        "",
        "```yaml",
        f"width: {config['model']['width']}",
        f"dna: {config['model']['dna']}",
        f"protein: {config['model']['protein']}",
        f"training: {config['training']}",
        f"validation: {config['splits']['validation']}",
        "```",
        "",
        "Fold membership is hashed rather than described, because a regime name and a seed do not",
        "pin down which domains were held out once the corpus changes (`ML_PLAN.md` §9.2). The",
        "digests are in `data/processed/training_folds.parquet` and in each MLflow run.",
        "",
        _provenance_line(frame),
        "",
        "```bash",
        "python scripts/run_grid.py          # about an hour for two arms",
        "```",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
