#!/usr/bin/env python
"""Sweep `model.bce_weight` — the decision-rule experiment `T38` asked for.

    python scripts/run_lambda_sweep.py [--budget-hours 55] [--dry-run]

**What this answers.** `L = L_infonce + λ · L_bce` (`snp2prot.models.loss.hybrid_loss`). The
InfoNCE term is a per-row softmax and is therefore exactly invariant to a per-protein score
offset, so it cannot produce a threshold and nothing applied afterwards can recover one. Binary
cross-entropy is not shift-invariant. `λ` decides how much of each, and three things have to be
measured across it rather than argued about:

1. **does the ranking survive** — macro AUPR against the same folds and the same
   nearest-neighbour baseline every other report uses. `λ = 0` is the pure-ranking objective
   exactly, and it is re-run *inside* this experiment rather than read off `reports/training.md`:
   GPU training is not reproducible across processes (measured, ~0.003 AUPR on `S1/fold-2`), so a
   control from another invocation would carry that noise into every delta;
2. **does calibration arrive** — ECE, and the spread of per-protein call counts, which is where
   `T38` recorded the failure: one global cut called 56 8-mers for one protein and 408 for
   another whose true counts were comparable;
3. **can it see a dead variant** — `cal_dead_auroc`, the protein-level detection score. This is
   the capability the project wants and the one no rank-based metric can express.

**Priority-ordered, budget-guarded and resumable, because it runs unattended.** The job list is
enumerated in the order the questions matter (`STAGES` below), a run is only started if the
observed mean run time says it fits in the remaining budget, and every finished run is appended
to `data/processed/lambda_sweep_folds.parquet` before the next one starts. Re-running the script
skips whatever is already in that file, so a kill, a crash or a reboot costs one run. A partial
sweep is still readable: the report is rewritten after **every** run, and each stage is a
complete answer to one question.

⚠️ **This is a multi-day job and is the owner's to run** (`CLAUDE.md` rule 10). At the measured
9 minutes per run the default 55 h budget covers about 340 runs. `--dry-run` prints the plan and
the arithmetic without training anything.
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np
import pandas as pd
from run_grid import baseline_by_fold, write_calibration_report
from train import load_everything

from snp2prot import experiment, splits, training
from snp2prot.config import PROCESSED_DIR, REPORTS_DIR

DEFAULT_OUT = REPORTS_DIR / "calibration.md"
SUMMARY_OUT = PROCESSED_DIR / "lambda_sweep_folds.parquet"
PER_DOMAIN_OUT = PROCESSED_DIR / "lambda_sweep_domains.parquet"

#: The λ grid, logarithmic and placed by measurement rather than by taste. At initialisation the
#: InfoNCE term is 10.44 and the BCE term 0.0166 (measured 2026-08-28 on `A1`/`P3/all`), a ratio
#: of 631 — so `λ = 1` is nearly the control and the two terms are comparable around `λ ≈ 600`.
#: The grid brackets that by more than an order of magnitude each way, because the ratio does not
#: hold still: InfoNCE falls to ~0.005 by step 6,000 while the BCE term plateaus, so a λ that is
#: negligible at step 1 can dominate by step 10,000.
LAMBDAS = (0.0, 1.0, 3.0, 10.0, 30.0, 100.0, 300.0, 1000.0, 3000.0, 10000.0)

#: The three folds the step-budget curve was measured on (`configs/experiment.yaml`,
#: `training.steps`), reused here so the λ curve and that curve are on the same ground: one
#: regime where the model does well, the regime the project hangs on, and the variant transfer.
DIAGNOSTIC = ("S1/fold-2", "S2/fold-3", "P3/all")

#: The λ values that get the full 19-fold treatment. Chosen blind — the sweep runs unattended, so
#: it cannot wait for stage 1 to pick one — and deliberately spanning the region where the BCE
#: term goes from a tenth of the total to several times it.
FULL_GRID_LAMBDAS = (0.0, 100.0)
EXTRA_GRID_LAMBDAS = (30.0, 300.0)

#: How many runs may fail back to back before the sweep gives up. One failure is a flake worth
#: skipping over on an unattended job; five in a row is a fault, and continuing would spend the
#: remaining budget writing tracebacks.
MAX_CONSECUTIVE_FAILURES = 5


def _jobs(name: str, arms, folds, lambdas, seeds) -> list[dict]:
    """One stage's jobs. `folds=None` means "every fold", filled in later by `expand`, which is
    the only place that knows what the corpus produced."""
    folds = [None] if folds is None else folds
    return [
        {"stage": name, "arm": arm, "fold": fold, "bce_weight": lam, "model_seed": seed}
        for arm in arms
        for lam in lambdas
        for fold in folds
        for seed in seeds
    ]


def plan(config: dict) -> list[dict]:
    """Every run this sweep would like to do, most informative first.

    The order is the point. Each stage answers one question completely, so a budget that runs out
    half way through leaves a readable experiment rather than a scattering of cells — and the
    stages that would only *refine* an answer (more seeds, the second arm, the wider λ grid on
    every fold) come after the ones that establish it.
    """
    base = int(config["model"]["seed"])
    seeds3 = (base, base + 1, base + 2)
    s2 = [f"S2/fold-{i}" for i in range(int(config["splits"]["n_folds"]))]
    return [
        # 1. Does the λ curve exist at all? One arm, one seed, the three diagnostic folds.
        *_jobs("1-curve", ["A1"], DIAGNOSTIC, LAMBDAS, [base]),
        # 2. Is it robust to the initialisation? The same curve at two more seeds — `T37`'s
        #    question asked of this experiment, because a λ effect smaller than the seed spread
        #    is not an effect.
        *_jobs("2-seeds", ["A1"], DIAGNOSTIC, LAMBDAS, seeds3[1:]),
        # 3. Does it hold on the other embedding arm?
        *_jobs("3-arm", ["A4"], DIAGNOSTIC, LAMBDAS, [base]),
        # 4. The full regime table at the control and at a mid λ, so the deliverable report has
        #    every fold and not just three.
        *_jobs("4-grid", ["A1"], None, FULL_GRID_LAMBDAS, [base]),
        *_jobs("5-grid-A4", ["A4"], None, FULL_GRID_LAMBDAS, [base]),
        # 6. Two more λ across the whole grid: a curve of three points per fold rather than two.
        *_jobs("6-grid-wide", ["A1"], None, EXTRA_GRID_LAMBDAS, [base]),
        # 7. `S2` in depth — the regime that decides whether any of this transfers — with seeds.
        *_jobs("7-s2-seeds", ["A1"], s2, (*FULL_GRID_LAMBDAS, *EXTRA_GRID_LAMBDAS), seeds3[1:]),
        *_jobs("8-grid-wide-A4", ["A4"], None, EXTRA_GRID_LAMBDAS, [base]),
    ]


def expand(jobs: list[dict], folds_by_arm: dict[str, list[str]]) -> list[dict]:
    """Replace a `None` fold with every fold of that arm, keeping the stage order intact.

    **Duplicates are dropped here, and they are not hypothetical.** The stages overlap by design
    — `1-curve` runs `λ = 30` on `P3/all` and `6-grid-wide` wants every fold at `λ = 30`, which
    is the same run — and the resume check only compares against what is already on disk, so
    within one invocation the same configuration would be trained twice. On the 2026-08-28 sweep
    it was: **32 of 312 runs, about 4.3 h**, and every duplicated cell was double-weighted in any
    naive average over the summary table. The first stage that wants a configuration keeps it,
    which preserves the priority ordering.
    """
    out, seen = [], set()
    for job in jobs:
        wanted = (
            [job]
            if job["fold"] is not None
            else [dict(job, fold=label) for label in folds_by_arm[job["arm"]]]
        )
        for one in wanted:
            key = (one["arm"], one["fold"], one["bce_weight"], one["model_seed"])
            if key in seen:
                continue
            seen.add(key)
            out.append(one)
    return out


def load_done(path: Path) -> tuple[pd.DataFrame, set]:
    """Runs already on disk, and their keys. An unreadable file is a stop, never a silent reset:
    overwriting a half-finished multi-day sweep because a parquet was truncated is not a
    recovery."""
    if not path.exists():
        return pd.DataFrame(), set()
    frame = pd.read_parquet(path)
    # The key uses the `regime/fold` label a job carries; the frame keeps its own two columns
    # untouched, because every report downstream groups on them separately.
    keys = {
        (str(r.arm), f"{r.regime}/{r.fold}", float(r.bce_weight), int(r.model_seed))
        for r in frame.itertuples()
    }
    return frame, keys


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--budget-hours",
        type=float,
        default=55.0,
        help="stop before starting a run that would not finish inside this (default 55)",
    )
    ap.add_argument("--steps", type=int, default=None, help="override training.steps")
    ap.add_argument("--device", default=None)
    ap.add_argument("--no-track", action="store_true", help="do not write to mlruns/")
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--dry-run", action="store_true", help="print the plan and the arithmetic")
    ap.add_argument(
        "--estimate-seconds",
        type=float,
        default=560.0,
        help="assumed run cost before any run has been timed; measured 2026-08-28 at ~540 s",
    )
    args = ap.parse_args()

    config = experiment.load()
    if args.steps is not None:
        config["training"]["steps"] = args.steps
    device = training.device_for(args.device)
    budget = args.budget_hours * 3600.0

    # Fold labels are a property of the corpus, not of the arm, but they are read per arm so
    # that a mismatch between an arm's embeddings and the corpus surfaces here rather than
    # 40 hours in.
    cache: dict[str, tuple] = {}

    def inputs(arm: str):
        if arm not in cache:
            cache[arm] = load_everything(arm)
        return cache[arm]

    domains, _, _, dist, _ = inputs("A1")
    labels = [f.label for f in splits.all_regimes(domains, dist)]
    jobs = expand(plan(config), {"A1": labels, "A4": labels})

    done_frame, done = load_done(SUMMARY_OUT)
    todo = [j for j in jobs if (j["arm"], j["fold"], j["bce_weight"], j["model_seed"]) not in done]

    print(f"{len(jobs)} runs planned, {len(done)} already on disk, {len(todo)} to go")
    per_stage = pd.Series([j["stage"] for j in todo]).value_counts().sort_index()
    for stage, n in per_stage.items():
        print(f"  {stage:<16} {n:>4} runs  ~{n * args.estimate_seconds / 3600:.1f} h")
    print(
        f"budget {args.budget_hours:.1f} h; at {args.estimate_seconds:.0f} s per run that is "
        f"about {int(budget / args.estimate_seconds)} runs"
    )
    if args.dry_run:
        return

    started = time.time()
    summaries: list[dict] = [] if done_frame.empty else done_frame.to_dict("records")
    per_domain: list[pd.DataFrame] = []
    if PER_DOMAIN_OUT.exists():
        per_domain.append(pd.read_parquet(PER_DOMAIN_OUT))
    observed: list[float] = []
    skipped = 0
    failures = 0

    for i, job in enumerate(todo, 1):
        elapsed = time.time() - started
        # The estimate becomes the observed mean as soon as there is one, so a slower machine or
        # a larger fold shortens the plan rather than overrunning it.
        expected = float(np.mean(observed)) if observed else args.estimate_seconds
        if elapsed + expected > budget:
            skipped = len(todo) - i + 1
            print(
                f"\nSTOPPING: {elapsed / 3600:.1f} h used, next run needs ~{expected / 60:.0f} min "
                f"and the budget is {args.budget_hours:.1f} h. {skipped} runs not attempted."
            )
            break

        arm, label = job["arm"], job["fold"]
        d, matrix, vectors, distance, trainable = inputs(arm)
        fold = {f.label: f for f in splits.all_regimes(d, distance)}[label]
        config["model"]["bce_weight"] = job["bce_weight"]

        run_started = time.time()
        try:
            summary, rows = training.run_fold(
                fold,
                arm,
                d,
                matrix,
                vectors,
                distance,
                config,
                trainable,
                device=device,
                track=not args.no_track,
                model_seed=job["model_seed"],
            )
        except Exception as failure:  # noqa: BLE001
            # **One bad run must not end forty hours of work.** The sweep is unattended, so a
            # failure is recorded and skipped rather than raised — but a *systematic* failure
            # would otherwise burn the whole budget printing tracebacks, so a run of them stops
            # it. A skipped run is simply absent from the parquet, which is exactly what a
            # resume looks for, so re-running the script retries it.
            failures += 1
            print(
                f"[{i}/{len(todo)}] FAILED  {job['stage']} {arm} {label} "
                f"lam={job['bce_weight']:g} seed={job['model_seed']}: "
                f"{type(failure).__name__}: {failure}",
                flush=True,
            )
            if failures >= MAX_CONSECUTIVE_FAILURES:
                print(
                    f"\nSTOPPING: {failures} runs failed in a row. This is a fault, not a "
                    "flake — read the traceback above before restarting.",
                    flush=True,
                )
                break
            continue
        failures = 0
        observed.append(time.time() - run_started)
        summary["stage"] = job["stage"]
        rows["stage"] = job["stage"]
        summaries.append(summary)
        per_domain.append(rows)

        # Written before the next run starts, so a kill costs exactly one run.
        frame = pd.DataFrame(summaries)
        SUMMARY_OUT.parent.mkdir(parents=True, exist_ok=True)
        frame.to_parquet(SUMMARY_OUT, index=False)
        domains_frame = pd.concat(per_domain, ignore_index=True)
        domains_frame.to_parquet(PER_DOMAIN_OUT, index=False)
        write_calibration_report(
            args.out,
            frame,
            domains_frame,
            config,
            baseline_by_fold(),
            grid_lambdas=(*FULL_GRID_LAMBDAS, *EXTRA_GRID_LAMBDAS),
        )

        print(
            f"[{i}/{len(todo)}] {job['stage']:<14} {arm} {label:<18} "
            f"lam={job['bce_weight']:<8g} seed={job['model_seed']}  "
            f"AUPR {summary['aupr']:.4f}  "
            f"ECE {summary.get('cal_ece', float('nan')):.4f}  "
            f"dead {summary.get('cal_dead_auroc', float('nan')):.3f}  "
            f"calls {summary.get('cal_n_called_median', float('nan')):.0f}"
            f"/{summary.get('cal_n_true_median', float('nan')):.0f}  "
            f"({(time.time() - started) / 3600:.2f} h)",
            flush=True,
        )

    print(
        f"\n{len(observed)} runs in {(time.time() - started) / 3600:.2f} h"
        + (f", {failures} failed" if failures else "")
        + (f", mean {np.mean(observed) / 60:.1f} min each" if observed else "")
        + (f", {skipped} left for a later invocation" if skipped else "")
    )
    print(f"wrote {args.out}, {SUMMARY_OUT} and {PER_DOMAIN_OUT}")


if __name__ == "__main__":
    main()
