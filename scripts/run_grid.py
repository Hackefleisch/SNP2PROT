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
    out = {}
    for k, prefix in ((1, "baseline"), (5, "baseline_k5")):
        sub = frame[frame.k == k]
        if not len(sub):
            continue
        grouped = sub.groupby(["regime", "fold"], sort=False)
        out[f"{prefix}_aupr"] = grouped.aupr.mean()
        out[f"{prefix}_median"] = grouped.aupr.median()
    return pd.DataFrame(out).reset_index()


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    add_common_arguments(ap)
    ap.add_argument("--arm", action="append", dest="arms", help="repeatable; default A1 and A4")
    ap.add_argument("--regime", action="append", help="restrict to these regimes")
    ap.add_argument(
        "--seeds",
        type=int,
        default=1,
        help="run each fold at this many model seeds, starting from model.seed (default 1). "
        "The only way to get an error bar: with one seed a delta between arms cannot be "
        "told from initialisation noise.",
    )
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = ap.parse_args()
    arms = args.arms or ["A1", "A4"]

    config = apply_overrides(experiment.load(), args)
    device = training.device_for(args.device)
    if args.seeds < 1:
        raise SystemExit("--seeds must be at least 1")
    seeds = [int(config["model"]["seed"]) + i for i in range(args.seeds)]
    started = time.time()

    summaries: list[dict] = []
    per_domain: list[pd.DataFrame] = []
    for arm in arms:
        domains, matrix, vectors, dist, trainable = load_everything(arm)
        folds = [f for f in splits.all_regimes(domains, dist)]
        if args.regime:
            folds = [f for f in folds if f.regime in set(args.regime)]
        seeded = f" x {len(seeds)} seeds" if len(seeds) > 1 else ""
        print(f"{arm} ({vectors.model}): {len(folds)} folds{seeded}")
        for fold in folds:
            for model_seed in seeds:
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
                    model_seed=model_seed,
                )
                summaries.append(summary)
                per_domain.append(rows)
                tag = f"  seed {model_seed}" if len(seeds) > 1 else ""
                print(
                    f"  {fold.label:<20} AUPR {summary['aupr']:.4f}  "
                    f"val {summary['validation_aupr']:.4f}  "
                    f"{int(summary['steps_run'])} steps{tag}  ({time.time() - started:.0f}s)",
                    flush=True,
                )

    frame = pd.DataFrame(summaries).merge(baseline_by_fold(), on=["regime", "fold"], how="left")
    frame["delta"] = frame.aupr - frame.baseline_aupr
    if "baseline_k5_aupr" in frame:
        frame["delta_k5"] = frame.aupr - frame.baseline_k5_aupr
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


def _seed_section(frame: pd.DataFrame) -> list[str]:
    """How much of a fold's number is the initialisation rather than the data.

    With one seed there is no error bar and a delta between arms cannot be told from noise, so
    the report says so rather than letting the reader assume otherwise.
    """
    per_fold = frame.groupby(["arm", "regime", "fold"], sort=False).aupr
    n = int(per_fold.count().max())
    if n < 2:
        return [
            "",
            "**One seed per fold.** `model.seed` fixes the initialisation and the batch order, so",
            "every number here is a single draw and none of them carry an error bar. A difference",
            "between arms or between folds smaller than the run-to-run spread cannot be told from",
            "initialisation noise, and that spread is unmeasured. `run_grid.py --seeds 3` measures",
            "it (`TODO.md` `T37`).",
            "",
        ]
    spread = per_fold.std(ddof=1)
    return [
        "",
        f"## Across seeds ({n} per fold)",
        "",
        "`model.seed` varies the initialisation and the batch order; the split is identical, so",
        "this isolates run-to-run noise. A delta smaller than the spread below is not a result.",
        "",
        "| regime | folds | mean spread (sd) | worst fold |",
        "|---|---:|---:|---:|",
        *[
            f"| {regime} | {len(g)} | {g.mean():.4f} | {g.max():.4f} |"
            for regime, g in spread.groupby(level="regime", sort=False)
        ],
        "",
    ]


def _suppression_section(frame: pd.DataFrame) -> list[str]:
    """The dead variants — the domains AUPR cannot reach, and the sharpest C1 evidence.

    Of the wild type's binding sites, the fraction the model ranks lower in the variant. The
    nearest-neighbour baseline scores exactly 0 wherever the wild type is in its training pool,
    because it predicts the variant by copying it — so this is C1's null hypothesis literally
    rather than by interpretation, and 0.5 is what an untargeted downward shift scores.
    """
    scored = frame[frame.n_scored_suppression > 0]
    lines = [
        "",
        "## Did it notice the mutation? — the dead variants",
        "",
        "20 held-out records have **no positive 8-mer at all**: variants whose binding measurably",
        "vanished (`T21`). AUPR, AUROC and R@P0.5 are undefined for every one of them, and they",
        "are the sharpest evidence the corpus holds for claim **C1**. They are scored instead by",
        "`suppression` — of the sites the wild type binds, the fraction the model ranks *lower*",
        "in the variant (`snp2prot.evaluation.metrics.suppression`).",
        "",
        "**1.0** the model saw the mutation abolish binding · **0.5** the sites moved at random ·",
        "**0.0** the variant is predicted exactly like its wild type, which is what the",
        "nearest-neighbour baseline does by construction.",
        "",
    ]
    if not len(scored):
        return lines + [
            "No fold scored one: the wild type was held out alongside every variant.",
            "",
        ]
    lines += [
        "| arm | regime | fold | dead variants | suppression | baseline |",
        "|---|---|---|---:|---:|---:|",
    ]
    for row in scored.itertuples():
        lines.append(
            f"| `{row.arm}` | {row.regime} | `{row.fold}` | {int(row.n_scored_suppression)} | "
            f"**{_fmt(row.suppression)}** | 0.0000 |"
        )
    return lines + [""]


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

    A multiple near 1 is an eyeball, not a test. `metrics.random_baseline` samples the macro
    statistic under a random ranking of the same held-out labels, and a result at or below the
    95th percentile of it is reported as **at chance** rather than as a small number.
    """
    at_chance = frame[frame.aupr <= frame.random_p95]
    lines = [
        "",
        "## Is it better than random?",
        "",
        "The reference is the macro AUPR of a *random* ranking of the same held-out domains,",
        f"sampled {repeats} times (`metrics.random_baseline`). A result at or below the 95th",
        "percentile of it is indistinguishable from guessing, whatever its delta against the",
        "baseline looks like.",
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
            f"{_fmt(row.random_mean)} | {_fmt(row.random_p95)} | **at chance** |"
        )
    return lines + [""]


def _calibration_section(frame: pd.DataFrame) -> list[str]:
    """The decision rule, in the fold table's own report.

    Blank at `λ = 0`, and that is a statement rather than a gap: the calibration head takes no
    gradient there, so its output is `sigmoid` of an initialisation and every number computed
    from it would be a constant dressed up as a measurement.
    """
    if "cal_ece" not in frame or not np.isfinite(frame.cal_ece).any():
        return [
            "",
            "## The decision rule",
            "",
            "**Not measured: every run in this table had `model.bce_weight = 0`.** The model",
            "emits a ranking and no call — `snp2prot.models.loss.multi_positive_infonce` is a",
            "per-row softmax and so is exactly invariant to a per-protein offset, which is why",
            "the null anchor never became a threshold (`T38`). `scripts/run_lambda_sweep.py`",
            "sweeps `λ` and writes [`calibration.md`](calibration.md).",
            "",
        ]
    scored = frame[np.isfinite(frame.cal_ece)]
    return [
        "",
        "## The decision rule",
        "",
        "Where `model.bce_weight > 0` the model emits a calibrated `P(binds)` and therefore a",
        "*call*, not only a ranking. `calls` is the median number of 8-mers named per protein",
        "under the parameter-free rule — keep the top `round(Σ p)` — against the median number",
        "that truly bind; `spread` is the gap between the most and fewest called in the fold,",
        "which is where `T38` recorded a global cut failing (56 to 408 for comparable proteins).",
        "`dead` is the protein-level AUROC for calling a variant that binds nothing from its",
        "predicted interaction power alone, over `n` such variants.",
        "",
        "| arm | regime | fold | λ | ECE | calls | true | spread | F1 | dead | n |",
        "|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|",
        *[
            f"| `{r.arm}` | {r.regime} | `{r.fold}` | {r.bce_weight:g} | {_fmt(r.cal_ece)} | "
            f"{r.cal_n_called_median:.0f} | {r.cal_n_true_median:.0f} | "
            f"{r.cal_count_spread:.0f} | {_fmt(r.cal_call_f1, 3)} | "
            f"**{_fmt(r.cal_dead_auroc, 3)}** | {r.cal_n_dead:.0f} |"
            for r in scored.itertuples()
        ],
        "",
    ]


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
        "`baseline` copies the neighbour's **binary calls** — the same information the model is",
        "trained on — so `delta` compares methods rather than inputs",
        "(`snp2prot.baselines.nn_lookup`). `chance` is what a random ranking scores on that",
        "fold's held-out domains, and the `x` columns are each AUPR as a multiple of it: it",
        "ranges 0.0015-0.0034 across the regimes, so two folds reporting the same AUPR are not",
        "reporting the same thing. `at budget` is the model at `training.steps` rather than the",
        "checkpoint the validation slice selected.",
        "",
        "## Per fold",
        "",
        "| arm | regime | fold | chance | model AUPR | x | `k=1` | delta | **`k=5`** "
        "| **delta k=5** | at budget | best step | AUROC |",
        "|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for keys, group in frame.groupby(["arm", "regime", "fold"], sort=False):
        arm, regime, fold = keys
        row = group.iloc[0]
        aupr = f"**{_fmt(group.aupr.mean())}**"
        if len(group) > 1:
            aupr += f" ±{group.aupr.std(ddof=1):.4f}"
        lines.append(
            f"| `{arm}` | {regime} | `{fold}` | {_fmt(row.chance_aupr)} | "
            f"{aupr} | {_lift(group.aupr.mean(), row.chance_aupr)} | "
            f"{_fmt(row.baseline_aupr)} | {group.delta.mean():+.4f} | "
            f"**{_fmt(getattr(row, 'baseline_k5_aupr', float('nan')))}** | "
            f"**{group.delta_k5.mean():+.4f}** | {_fmt(group.aupr_final.mean())} | "
            f"{int(group.best_step.mean())} | {_fmt(group.auroc.mean(), 3)} |"
        )

    lines += _seed_section(frame)
    lines += _null_section(frame, int(config["metrics"]["null_repeats"]))
    lines += _calibration_section(frame)
    lines += _suppression_section(frame)
    lines += _selection_section(frame)

    lines += [
        "",
        "## Per regime, averaged over folds",
        "",
        "| arm | regime | folds | model AUPR | `k=1` | delta | **`k=5`** | **delta k=5** |",
        "|---|---|---:|---:|---:|---:|---:|---:|",
    ]
    for (arm, regime), group in frame.groupby(["arm", "regime"], sort=False):
        lines.append(
            f"| `{arm}` | {regime} | {len(group)} | **{_fmt(group.aupr.mean())}** | "
            f"{_fmt(group.baseline_aupr.mean())} | {group.delta.mean():+.4f} | "
            f"**{_fmt(group.baseline_k5_aupr.mean())}** | "
            f"**{group.delta_k5.mean():+.4f}** |"
        )

    lines += [
        "",
        "## The C1 evaluation set",
        "",
        "The variants whose binding measurably changed, so that their own wild type does not",
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
        "python scripts/run_grid.py          # about 5 h for two arms at 15,000 steps",
        "```",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")


# --- the lambda sweep's report (`scripts/run_lambda_sweep.py`) -------------------------------


def _lambda_table(group: pd.DataFrame) -> list[str]:
    """One row per λ, aggregated over whatever seeds and folds the caller passed in."""
    lines = [
        "| λ | runs | AUPR | vs λ=0 | ECE | calls | true | spread | F1 | dead AUROC | ratio |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    control = group[group.bce_weight == 0].aupr.mean()
    for lam, part in group.groupby("bce_weight", sort=True):
        aupr = f"**{_fmt(part.aupr.mean())}**"
        if len(part) > 1:
            aupr += f" ±{part.aupr.std(ddof=1):.4f}"
        delta = "—" if lam == 0 else f"{part.aupr.mean() - control:+.4f}"
        lines.append(
            f"| {lam:g} | {len(part)} | {aupr} | {delta} | "
            f"{_fmt(part.cal_ece.mean())} | "
            f"{_fmt(part.cal_n_called_median.mean(), 0)} | "
            f"{_fmt(part.cal_n_true_median.mean(), 0)} | "
            f"{_fmt(part.cal_count_spread.mean(), 0)} | "
            f"{_fmt(part.cal_call_f1.mean(), 3)} | "
            f"**{_fmt(part.cal_dead_auroc.mean(), 3)}** | "
            f"{_fmt(part.cal_power_ratio.mean(), 3)} |"
        )
    return lines


def _reading(frame: pd.DataFrame) -> list[str]:
    """The arithmetic a reader would do first, done here so it is reproducible.

    **This is a computation, not a decision.** It names the λ values whose ranking is
    statistically indistinguishable from the `λ = 0` control and orders those by the capability
    the sweep exists to obtain. Which one to adopt — and whether any of them is good enough to
    adopt at all — belongs in `docs/DECISIONS.md` and is the owner's.
    """
    control = frame[frame.bce_weight == 0]
    if not len(control):
        return []
    # The bar the ranking must clear: the control's mean, less one seed-spread. With one seed
    # per cell there is no spread to subtract and the bar is the control itself, which is the
    # strict reading and is stated as such.
    spread = frame.groupby(["arm", "regime", "fold", "bce_weight"]).aupr.std(ddof=1).mean()
    measured = bool(np.isfinite(spread))
    spread = float(spread) if measured else 0.0
    bar = control.aupr.mean() - spread
    per_lambda = frame.groupby("bce_weight").agg(
        aupr=("aupr", "mean"), dead=("cal_dead_auroc", "mean"), ece=("cal_ece", "mean")
    )
    survivors = per_lambda[(per_lambda.index > 0) & (per_lambda.aupr >= bar)]
    lines = [
        "",
        "## Reading it",
        "",
        f"The `λ = 0` control ranks at **{control.aupr.mean():.4f}** over these folds, and the",
        (
            f"mean across-seed spread is **{spread:.4f}** — so a λ whose AUPR is at or above"
            if measured
            else "across-seed spread is **not yet measured** — no cell has a second seed, so the"
            " bar below is the control itself, which is the strict reading. A λ whose AUPR is"
            " at or above"
        ),
        f"**{bar:.4f}** has not measurably cost anything in ranking. Among those, the ones that",
        "buy the most dead-variant detection:",
        "",
    ]
    if not len(survivors):
        return lines + [
            "**None.** Every λ above 0 cost more ranking than the seed spread can explain, which",
            "is the outcome that says the two objectives genuinely conflict on this data. The",
            "honest reading is then `T38`'s third option: the output is a ranking, and the",
            "deliverable is top-`k` retrieval.",
            "",
        ]
    ordered = survivors.sort_values("dead", ascending=False)
    lines += [
        "| λ | AUPR | vs control | dead AUROC | ECE |",
        "|---:|---:|---:|---:|---:|",
        *[
            f"| {lam:g} | {_fmt(r.aupr)} | {r.aupr - control.aupr.mean():+.4f} | "
            f"**{_fmt(r.dead, 3)}** | {_fmt(r.ece)} |"
            for lam, r in ordered.iterrows()
        ],
        "",
        "A dead-variant AUROC at or below 0.5 means the model cannot tell a variant that binds",
        "nothing from one that still binds, whatever its AUPR — which is the capability `T38`",
        "was raised about, and a good ranking does not substitute for it.",
        "",
    ]
    return lines


def write_calibration_report(
    path: Path,
    frame: pd.DataFrame,
    per_domain: pd.DataFrame,
    config: dict,
    baseline: pd.DataFrame,
    grid_lambdas: tuple[float, ...] = (),
) -> None:
    """`reports/calibration.md` — the λ sweep, rewritten after every run so a partial sweep reads.

    Shares `write_report`'s helpers deliberately: the two reports must anchor to the same null,
    the same baseline and the same formatting rules, or a number lifted from one into the other
    would silently change meaning.
    """
    frame = frame.merge(baseline, on=["regime", "fold"], how="left")
    frame["delta_k5"] = frame.aupr - frame.get("baseline_k5_aupr", np.nan)
    for column in (
        "cal_ece",
        "cal_dead_auroc",
        "cal_n_called_median",
        "cal_n_true_median",
        "cal_count_spread",
        "cal_call_f1",
        "cal_power_ratio",
    ):
        if column not in frame:
            frame[column] = np.nan

    done = frame.groupby("stage").size() if "stage" in frame else pd.Series(dtype=int)
    lines = [
        "# The decision rule — sweeping λ on `L = L_infonce + λ · L_bce`",
        "",
        "Generated by `scripts/run_lambda_sweep.py`. The question is `T38`: the model produced a",
        "ranking and the deliverable needs a **call**, and no defensible threshold existed.",
        "",
        "**The diagnosis.** `multi_positive_infonce` is a per-row softmax, so it depends only on",
        "differences *within* one protein's row and is exactly invariant to adding a constant to",
        "that row. Nothing in it ever says where a protein's scores should sit relative to",
        "another's. That is why the null anchor settled inside the negative cloud rather than",
        "between the classes, and why a global Platt fit on the validation slice could not rescue",
        "it either — no post-hoc step recovers an offset the objective never constrained.",
        "",
        "**The change.** A masked, unweighted binary cross-entropy term over the same cells,",
        "through a two-scalar calibration head, weighted by `λ`. It is not shift-invariant, and",
        "`sigmoid` of its output is a probability comparable across proteins. **`λ = 0` is the",
        "pure-ranking objective exactly**, and it is re-run inside this experiment rather than",
        "quoted from `training.md` — GPU training is not reproducible across processes (~0.003",
        "AUPR on `S1/fold-2`), so a control from another invocation would carry that into every",
        "delta.",
        "",
        "**What the columns mean.** `calls` is the median number of 8-mers named per protein",
        "under the parameter-free rule — keep the top `round(Σ p)`, which has no threshold to",
        "choose — against `true`, the median number that actually bind. `spread` is the gap",
        "between the most and fewest called within a fold: `T38` measured a global cut calling 56",
        "for one protein and 408 for another with comparable true counts, and a spread that",
        "tracks the truth's is the outcome that fixes it. **`dead AUROC`** is the one the project",
        "wants — separating variants that bind *nothing* from variants that still bind, using",
        "the predicted interaction power `Σ p` and nothing else. `power ratio` is its paired",
        "form: a dead variant's predicted power over its own wild type's, where 0 is *the",
        "mutation abolished binding* and 1 is *the mutation did nothing*, which is what copying",
        "the wild type scores by construction.",
        "",
        f"**{len(frame)} runs so far.**"
        + (
            "  Stages completed: " + ", ".join(f"`{k}` ({v})" for k, v in done.items())
            if len(done)
            else ""
        ),
        "",
    ]

    curve_stages = ["1-curve", "2-seeds", "3-arm"]
    diagnostic = frame[frame.stage.isin(curve_stages)] if "stage" in frame else frame
    if len(diagnostic):
        lines += [
            "## 1. The λ curve, on the three diagnostic folds",
            "",
            "`S1/fold-2`, `S2/fold-3` and `P3/all` — the folds the step-budget curve was measured",
            "on, so this curve and that one stand on the same ground. Averaged over folds and",
            "seeds; the ± is the across-seed spread where more than one seed has run.",
            "",
        ]
        for arm, group in diagnostic.groupby("arm", sort=True):
            lines += [f"### arm `{arm}`", "", *_lambda_table(group), ""]

        lines += ["### Per fold, so an average cannot hide a disagreement", ""]
        for keys, group in diagnostic.groupby(["arm", "regime", "fold"], sort=True):
            arm, regime, fold = keys
            lines += [f"**`{arm}` · {regime}/`{fold}`**", "", *_lambda_table(group), ""]
        lines += _reading(diagnostic)

    # **Not selected by stage.** A run that appears in more than one stage of the plan is trained
    # once and carries the label of the first stage that wanted it — the `λ = 0` diagnostic folds
    # belong to `1-curve` and to `4-grid` both — so filtering the grid table by stage name would
    # silently drop three folds from it. The λ values and the base seed select it instead.
    base_seed = int(config["model"]["seed"])
    grid = (
        frame[frame.bce_weight.isin(grid_lambdas) & (frame.model_seed == base_seed)]
        if grid_lambdas
        else frame.iloc[0:0]
    )
    if len(grid):
        lines += [
            "",
            "## 2. The full 19-fold grid, at the λ values that got it",
            "",
            "Every regime, so the λ effect can be read where it matters — `S2`, unseen protein",
            "components — rather than only where the model does well. `k=5` is the",
            "nearest-neighbour baseline on the same held-out domains ([`nn_baseline.md`]"
            "(nn_baseline.md)).",
            "",
            "**`folds` says how complete each row is** — a λ whose grid stage has not run yet",
            "shows only the three diagnostic folds, and the count is how you can tell.",
            "",
            "| arm | regime | λ | folds | AUPR | `k=5` | delta | ECE | calls | true | dead AUROC |",
            "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
        ]
        for keys, group in grid.groupby(["arm", "regime", "bce_weight"], sort=True):
            arm, regime, lam = keys
            lines.append(
                f"| `{arm}` | {regime} | {lam:g} | {len(group)} | "
                f"**{_fmt(group.aupr.mean())}** | {_fmt(group.baseline_k5_aupr.mean())} | "
                f"{group.delta_k5.mean():+.4f} | {_fmt(group.cal_ece.mean())} | "
                f"{_fmt(group.cal_n_called_median.mean(), 0)} | "
                f"{_fmt(group.cal_n_true_median.mean(), 0)} | "
                f"**{_fmt(group.cal_dead_auroc.mean(), 3)}** |"
            )
        lines.append("")

    # `power_ratio` is absent entirely until a λ > 0 run has scored a dead variant whose wild
    # type stayed in training, so the column is checked for rather than assumed.
    if "power_ratio" in per_domain:
        dead = per_domain[np.isfinite(per_domain.power_ratio.to_numpy(dtype=float))]
    else:
        dead = per_domain.iloc[0:0]
    if len(dead):
        lines += [
            "",
            "## 3. The dead variants, one row each",
            "",
            "Variants whose binding measurably vanished, paired with their own wild type. One row",
            "per (fold, λ) occurrence, so `rows` exceeds `variants` wherever a variant is held out",
            "in more than one fold — the corpus holds 18 of them under `label_health`.",
            "`power ratio` near 0 is the model saying the mutation abolished binding; near 1 is",
            "the model saying it did nothing, which is what the nearest-neighbour baseline says",
            "by construction. `suppression` is the rank-based metric that predates this work and",
            "needs no threshold, kept beside it as a cross-check.",
            "",
            "| λ | rows | variants | mean power ratio | median | below 0.5 | mean suppression |",
            "|---:|---:|---:|---:|---:|---:|---:|",
            *[
                f"| {lam:g} | {len(g)} | {g.domain.nunique()} | {_fmt(g.power_ratio.mean(), 3)} | "
                f"{_fmt(g.power_ratio.median(), 3)} | "
                f"{_fmt(float((g.power_ratio < 0.5).mean()), 3)} | "
                f"{_fmt(np.nanmean(g.suppression), 3)} |"
                for lam, g in dead.groupby("bce_weight", sort=True)
            ],
            "",
        ]

    lines += [
        "",
        "## What this cannot say",
        "",
        "**The dead-variant numbers rest on 18 variants** (`label_health`'s `dead_variant`",
        "records, `no_evidence` excluded). At that `n` a detection AUROC is a direction and not a",
        "result, and the paired `power ratio` is the stronger read. `T30` is the open decision",
        "that would enlarge the set.",
        "",
        "**A calibrated probability is not a validated one.** ECE says the numbers are not lies;",
        "it does not say they are informative — a model predicting the base rate for every cell",
        "scores an excellent ECE and is useless. Read it beside `F1` and `calls`/`true`, never",
        "alone.",
        "",
        "**The calibration columns are blank at `λ = 0` by construction**, not by omission: the",
        "head takes no gradient there, so its output is `sigmoid` of an initialisation.",
        "",
        "## Configuration",
        "",
        "```yaml",
        f"width: {config['model']['width']}",
        f"dna: {config['model']['dna']}",
        f"protein: {config['model']['protein']}",
        f"training: {config['training']}",
        "```",
        "",
        _provenance_line(frame),
        "",
        "```bash",
        "python scripts/run_lambda_sweep.py --budget-hours 55   # resumable, see the docstring",
        "```",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
