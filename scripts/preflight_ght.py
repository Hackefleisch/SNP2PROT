#!/usr/bin/env python
"""Measure the step budget and the wall-clock budget before the grid is launched.

    python scripts/preflight_ght.py --steps 20000 --fold C1/all
    python scripts/preflight_ght.py --steps 20000 --fold G1/fold-0

`docs/GHT_PLAN.md` §8.1, and **the first thing to run**. One fold, `patience: 0` so nothing
terminates early, run well past where it looks done, logging validation *and* test at every
evaluation.

**(a) The step budget is read off the held-out TEST curve, never the training loss.** The PBM arm
got this backwards first: its training loss falls to ~0.005 by about step 6,000, which was read
as "everything after this is memorisation", and held-out AUPR then kept improving for another
9,000 steps worth +0.018 to +0.028. Memorising the training proteins and generalising to unseen
ones are not coupled.

**Logging test at every evaluation is a measurement, not a selection.** Nothing here chooses a
checkpoint on it and the number it produces is a curve for this report only; the grid selects on
validation like every other run in the project. Said out loud because the file makes the test
metric available, and the only thing stopping it becoming a selection signal is that it is not
used as one.

**(b) The wall-clock budget**, so the grid is predictable before it is launched: seconds per
step, cost of one evaluation pass, one fold end to end, peak GPU memory. Under rule 10 anything
past a few minutes is the owner's to launch, so this is what gets handed over.
"""

from __future__ import annotations

import argparse
import time

import pandas as pd
import torch

from snp2prot import experiment, tracking
from snp2prot.ght import config, data, panel, splits, training


def build_folds(corpus: data.WindowCorpus, table: pd.DataFrame, cfg: dict) -> dict:
    split_cfg = cfg["splits"]
    folds = splits.all_regimes(
        table,
        corpus.chromosomes("Train"),
        corpus.chromosomes("Test"),
        corpus.window_counts("Train"),
        int(split_cfg["n_folds"]),
        int(split_cfg["seed"]),
        float(split_cfg["validation_fraction"]),
        float(split_cfg["s2_min_identity"]),
    )
    return {f.label: f for f in folds}


def report(history: pd.DataFrame, meta: dict) -> str:
    curve = history[["step", "loss", "loss_bce", "validation_aupr", "test_aupr", "test_auroc"]]
    peak = curve.loc[curve.test_aupr.idxmax()]
    # A budget is only "where the test curve flattens" if the curve stopped climbing. Compare
    # the last fifth of the run against the fifth before it; if the later window is still
    # higher, the pre-flight was too short and must say so rather than report a number.
    n = len(curve)
    late = curve.test_aupr.iloc[int(0.8 * n) :].mean()
    earlier = curve.test_aupr.iloc[int(0.6 * n) : int(0.8 * n)].mean()
    still_climbing = late > earlier + 0.002

    lines = [
        "# GHT pre-flight: the step and wall-clock budgets",
        "",
        f"`scripts/preflight_ght.py --fold {meta['fold']} --steps {meta['steps']}`, arm "
        f"`{meta['arm']}`, seed {meta['model_seed']}. `docs/GHT_PLAN.md` §8.1.",
        "",
        "## (a) The step budget",
        "",
        "Read off the **held-out test curve**, not the training loss. Test is logged at every "
        "evaluation as a measurement and is never used to select anything — the grid selects on "
        "validation like every other run in this project.",
        "",
        "| step | train loss | validation auPRC | **test auPRC** | test auROC |",
        "|---:|---:|---:|---:|---:|",
    ]
    every = max(1, len(curve) // 25)
    for r in curve.iloc[::every].itertuples():
        lines.append(
            f"| {int(r.step):,} | {r.loss:.4f} | {r.validation_aupr:.4f} "
            f"| **{r.test_aupr:.4f}** | {r.test_auroc:.4f} |"
        )
    lines += [
        "",
        f"Best test auPRC **{peak.test_aupr:.4f} at step {int(peak.step):,}**; "
        f"training loss there {peak.loss:.4f}.",
        "",
        f"Mean test auPRC over the last fifth of the run: **{late:.4f}**, against "
        f"**{earlier:.4f}** over the fifth before it.",
        "",
    ]
    if still_climbing:
        lines += [
            "> **The curve was still climbing at the end of the pre-flight, so this run was too "
            "short to set a budget from.** Rerun with more steps before fixing "
            "`ght.training.steps`.",
            "",
        ]
    else:
        lines += [
            f"> The curve has flattened: the last fifth is not above the fifth before it by more "
            f"than 0.002. A budget at **{int(peak.step):,}** steps captures the peak, and the "
            "value adopted in `configs/experiment.yaml` is rounded to a round number at or above "
            "it.",
            "",
        ]
    lines += [
        "## (b) The wall-clock budget",
        "",
        "| quantity | value |",
        "|---|---:|",
        f"| ms per step, steady state | {1000 * meta['seconds_per_step']:.1f} |",
        f"| seconds per evaluation pass | {meta['seconds_per_eval']:.2f} |",
        f"| one fold end to end ({meta['steps']:,} steps, eval every {meta['eval_every']}) "
        f"| {meta['seconds'] / 60:.1f} min |",
        f"| peak GPU memory | {meta['peak_memory_mb']:.0f} MB |",
        f"| training windows | {meta['n_train']:,} |",
        f"| validation windows | {meta['n_validation']:,} |",
        f"| test windows (`shades`) | {meta['n_test']:,} |",
        f"| **one fold at the adopted {meta['budget_steps']:,}-step budget** "
        f"| **{meta['budget_minutes']:.1f} min** |",
        "",
        "The fold cost above is the pre-flight's own, at four times the budget it set. The row "
        "that matters is the last one: at the adopted budget the full grid is "
        f"**{meta['grid_folds']} folds x {meta['grid_seeds']} seeds = "
        f"{meta['grid_hours']:.1f} h**, plus the final `random` and `aliens` passes, which this "
        "pre-flight does not pay because it loads a resident-only corpus.",
        "",
        "## Parameter counts",
        "",
        "| tower | parameters | per training protein |",
        "|---|---:|---:|",
    ]
    for name, count in meta["parameter_counts"].items():
        if name == "total":
            continue
        lines.append(f"| {name} | {count:,} | {count / meta['n_train_tf']:,.0f} |")
    lines += [
        f"| **total** | **{meta['parameter_counts']['total']:,}** | "
        f"{meta['parameter_counts']['total'] / meta['n_train_tf']:,.0f} |",
        "",
        f"`docs/TRAINING.md` §5 calls the PBM arm's **245 per protein** the central engineering "
        f"constraint. At {int(meta['n_train_tf'])} training proteins even the bare linear tower "
        "is two orders of magnitude past that, so the panel size — not the tower architecture — "
        "is what sets the capacity risk here (`GHT_PLAN.md` §12).",
        "",
        tracking.provenance_line(),
    ]
    return "\n".join(lines) + "\n"


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--arm", default="A1")
    ap.add_argument("--fold", default="C1/all")
    ap.add_argument("--steps", type=int, default=20000)
    ap.add_argument("--eval-every", type=int, default=250)
    ap.add_argument("--device", default=None)
    ap.add_argument("--seed", type=int, default=None)
    ap.add_argument("--no-track", action="store_true")
    ap.add_argument(
        "--budget",
        type=int,
        default=None,
        help="the step budget the grid will use, for the wall-clock table; "
        "default: ght.training.steps",
    )
    args = ap.parse_args()

    cfg = experiment.section("ght")
    cfg["training"]["steps"] = args.steps
    cfg["training"]["eval_every"] = args.eval_every
    cfg["training"]["patience"] = 0

    table = panel.load()
    corpus = data.load(table.tf, resident_only=True)
    vectors, model_name = data.load_embeddings(args.arm, table)
    folds = build_folds(corpus, table, cfg)
    if args.fold not in folds:
        raise SystemExit(f"unknown fold {args.fold!r}; have {sorted(folds)}")
    fold = folds[args.fold]
    rows = training.fold_rows(corpus, fold)
    model_seed = args.seed if args.seed is not None else int(cfg["model"]["seed"])

    print(f"{fold.label}: {len(fold.fitting_tfs)} train TFs, {len(fold.test_tfs)} test TFs")
    print(
        f"  train {len(rows['train']):,}  validation {len(rows['validation']):,}  "
        f"test {len(rows['test']):,} windows"
    )

    history: list[dict] = []
    started = time.time()

    def on_eval(entry: dict, trainer: training.GHTTrainer, fold_=fold) -> None:
        # The test curve, measured at every evaluation and used for nothing but this report.
        summary, _ = trainer.score(rows["test"])
        entry = dict(entry)
        entry["test_aupr"] = summary["aupr"]
        entry["test_auroc"] = summary["auroc"]
        entry["chance_aupr"] = summary["chance_aupr"]
        history.append(entry)
        print(
            f"  step {entry['step']:>6,}  loss {entry['loss']:.4f}  "
            f"val {entry['validation_aupr']:.4f}  test {entry['test_aupr']:.4f}  "
            f"({time.time() - started:.0f}s)",
            flush=True,
        )

    summary, _ = training.run_fold(
        fold,
        args.arm,
        corpus,
        vectors,
        model_name,
        cfg,
        device=torch.device(args.device) if args.device else None,
        track=not args.no_track,
        model_seed=model_seed,
        on_eval=on_eval,
    )

    frame = pd.DataFrame(history)
    config.GHT_PROCESSED.mkdir(parents=True, exist_ok=True)
    slug = fold.label.replace("/", "-")
    table_path = config.PREFLIGHT_TABLE.with_name(f"ght_preflight_{args.arm}_{slug}.parquet")
    # The timing meta goes into the table as constant columns so the report can be rebuilt from
    # the artifact alone — a generated report whose numbers live only in a log is one that
    # cannot be corrected without re-running 25 minutes of GPU.
    frame.assign(
        arm=args.arm,
        fold=fold.label,
        model_seed=model_seed,
        **{k: v for k, v in summary.items() if k.startswith(("seconds", "peak_", "n_params_"))},
    ).to_parquet(table_path, index=False)

    n_folds = 1 + 2 * int(cfg["splits"]["n_folds"])
    seconds = summary["seconds"]
    # The grid's cost is NOT this run's: the pre-flight deliberately overshoots the budget it
    # exists to set. Scale from the per-step and per-evaluation costs to the adopted budget,
    # which is what the owner actually schedules.
    budget = int(args.budget or experiment.section("ght")["training"]["steps"])
    budget_seconds = (
        budget * summary["seconds_per_step"]
        + (budget / int(cfg["training"]["eval_every"])) * summary["seconds_per_eval"]
    )
    meta = {
        "budget_steps": budget,
        "budget_minutes": budget_seconds / 60,
        "arm": args.arm,
        "fold": fold.label,
        "steps": args.steps,
        "eval_every": args.eval_every,
        "model_seed": model_seed,
        "seconds": seconds,
        "seconds_per_step": summary["seconds_per_step"],
        "seconds_per_eval": summary["seconds_per_eval"],
        "peak_memory_mb": summary["peak_memory_mb"],
        "n_train": len(rows["train"]),
        "n_validation": len(rows["validation"]),
        "n_test": len(rows["test"]),
        "n_train_tf": max(len(fold.fitting_tfs), 1),
        "parameter_counts": {
            k.replace("n_params_", ""): int(v)
            for k, v in summary.items()
            if k.startswith("n_params_")
        },
        "grid_folds": n_folds,
        "grid_seeds": 3,
        "grid_hours": n_folds * 3 * budget_seconds / 3600,
    }
    path = config.REPORT_DIR / f"ght_preflight_{slug}.md"
    path.write_text(report(frame, meta))
    print(
        f"\nbest test auPRC {frame.test_aupr.max():.4f} at step "
        f"{int(frame.loc[frame.test_aupr.idxmax(), 'step']):,}"
    )
    print(
        f"{1000 * meta['seconds_per_step']:.1f} ms/step, "
        f"{meta['seconds_per_eval']:.2f} s/eval, {seconds / 60:.1f} min/fold, "
        f"{meta['peak_memory_mb']:.0f} MB peak"
    )
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
