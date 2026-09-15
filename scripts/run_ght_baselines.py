#!/usr/bin/env python
"""The bar: the PWM baseline in Setting 1 and motif transfer by DBD identity in Setting 2.

    python scripts/run_ght_baselines.py            # ~10 min on the GPU

`docs/GHT_PLAN.md` §9. Two passes over `pwms.zip`:

**Pass 1, selection.** Every TF's own PWMs are scanned over that TF's **training**-chromosome
windows and the single best-scoring motif is kept, by auPRC. Best-hit and sum-occupancy are
selected independently. Nothing here ever looks at a test chromosome.

**Pass 2, scoring.** Every TF's test-chromosome windows are scanned against the union of the 33
selected motifs — one pass, then any combination of them is a column selection. Setting 1 reads
the TF's own column; Setting 2 reads the columns of its nearest training neighbours and never its
own.

ArChIPelago is cited, not rerun (`D11`): we follow their positives, negatives and splits, so
their published numbers apply.
"""

from __future__ import annotations

import argparse
import time

import numpy as np
import pandas as pd
import torch

from snp2prot import experiment, tracking
from snp2prot.evaluation import metrics
from snp2prot.ght import baselines, config, data, panel, splits, training


def score_block(truth: np.ndarray, scores: np.ndarray) -> dict:
    return {
        "n": int(len(truth)),
        "n_pos": int(truth.sum()),
        "positive_rate": float(truth.mean()),
        "aupr": metrics.average_precision(truth, scores),
        "auroc": metrics.auroc(truth, scores),
    }


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--device", default=None)
    ap.add_argument("--top-k", type=int, default=5)
    args = ap.parse_args()

    device = torch.device(args.device) if args.device else None
    cfg = experiment.section("ght")
    table = panel.load()
    corpus = data.load(table.tf)
    pwms = baselines.read_pwms(table.tf)
    print(f"{len(table)} TFs, {sum(len(v) for v in pwms.values()):,} PWMs")
    started = time.time()

    # --- pass 1: select one motif per TF on its own training chromosomes ---------------------
    train_chromosomes = tuple(corpus.chromosomes("Train"))
    test_chromosomes = tuple(corpus.chromosomes("Test"))
    selected: dict[str, dict] = {}
    for tf in table.tf:
        rows = corpus.rows((tf,), train_chromosomes, training.TRAIN_NEGATIVES)
        truth = (corpus.label[rows] == 1).astype(np.int64)
        best, total = baselines.scan(corpus.tokens[rows], pwms[tf], device=device)
        i_best, aupr_best = baselines.select(best, truth)
        i_total, aupr_total = baselines.select(total, truth)
        use_total = aupr_total > aupr_best
        chosen = pwms[tf][i_total if use_total else i_best]
        selected[tf] = {
            "tf": tf,
            "pwm": chosen.name,
            "source": chosen.source,
            "width": chosen.width,
            "aggregation": "sum_occupancy" if use_total else "best_hit",
            "train_aupr": max(aupr_best, aupr_total),
            "n_pwms": len(pwms[tf]),
            "matrix": chosen.matrix,
        }
        print(
            f"  {tf:<9} {len(pwms[tf]):>4} PWMs -> {chosen.name.split('@')[-1][:28]:<28} "
            f"{selected[tf]['aggregation']:<13} train auPRC {selected[tf]['train_aupr']:.3f} "
            f"({time.time() - started:.0f}s)",
            flush=True,
        )

    order = list(table.tf)
    bank = [
        baselines.PWM(selected[tf]["pwm"], tf, selected[tf]["source"], selected[tf]["matrix"])
        for tf in order
    ]
    aggregation = np.array([selected[tf]["aggregation"] for tf in order])

    # --- pass 2: score every TF's test windows against every selected motif ------------------
    identity = baselines.identity_matrix(table)
    folds = {
        f.label: f
        for f in splits.all_regimes(
            table,
            list(train_chromosomes),
            list(test_chromosomes),
            corpus.window_counts("Train"),
            int(cfg["splits"]["n_folds"]),
            int(cfg["splits"]["seed"]),
            float(cfg["splits"]["validation_fraction"]),
            float(cfg["splits"]["s2_min_identity"]),
        )
    }
    tf_regime_folds = {k: v for k, v in folds.items() if v.regime in splits.TF_REGIMES}

    rows_out = []
    for position, tf in enumerate(order):
        for negatives in training.EVALUATION_NEGATIVES:
            rows = corpus.rows((tf,), test_chromosomes, negatives)
            truth = (corpus.label[rows] == 1).astype(np.int64)
            best, total = baselines.scan(corpus.tokens[rows], bank, device=device)
            picked = np.where(aggregation == "sum_occupancy", 1, 0)
            scores = np.where(picked[None, :] == 1, total, best)
            standard = baselines.standardise(scores)

            # Setting 1: the TF's own motif, selected on its own training chromosomes.
            rows_out.append(
                {
                    "regime": splits.CHROMOSOME_REGIME,
                    "fold": "all",
                    "method": "pwm",
                    "tf": tf,
                    "negset": negatives,
                    "neighbour": tf,
                    "identity": 1.0,
                    **score_block(truth, scores[:, position]),
                }
            )

            # Setting 2: the nearest training neighbours' motifs, never the TF's own.
            for fold in tf_regime_folds.values():
                if tf not in fold.test_tfs:
                    continue
                train_positions = np.array([order.index(t) for t in fold.train_tfs])
                picks, values = baselines.neighbours(
                    identity, [position], train_positions, args.top_k
                )
                usable = picks[0] >= 0
                if not usable.any():
                    continue
                columns = train_positions[picks[0][usable]]
                weights = np.clip(values[0][usable], 0.0, None)
                for k in (1, args.top_k):
                    take = min(k, len(columns))
                    w = weights[:take]
                    w = w / w.sum() if w.sum() > 0 else np.full(take, 1 / take)
                    blended = standard[:, columns[:take]] @ w
                    rows_out.append(
                        {
                            "regime": fold.regime,
                            "fold": fold.name,
                            "method": f"nn{k}",
                            "tf": tf,
                            "negset": negatives,
                            "neighbour": order[columns[0]],
                            "identity": float(values[0][0]),
                            **score_block(truth, blended),
                        }
                    )
        print(f"  scored {tf} ({time.time() - started:.0f}s)", flush=True)

    frame = pd.DataFrame(rows_out)
    code = tracking.code_version()
    frame = frame.assign(code_commit=code["commit"], code_dirty=code["dirty"])
    config.GHT_PROCESSED.mkdir(parents=True, exist_ok=True)
    frame.to_parquet(config.BASELINE_TABLE, index=False)
    pd.DataFrame(
        [{k: v for k, v in s.items() if k != "matrix"} for s in selected.values()]
    ).to_parquet(config.GHT_PROCESSED / "ght_selected_pwms.parquet", index=False)

    print(f"\n{len(frame):,} baseline rows in {time.time() - started:.0f}s")
    summary = (
        frame[frame.negset == training.TRAIN_NEGATIVES]
        .groupby(["regime", "method"])
        .agg(aupr=("aupr", "mean"), auroc=("auroc", "mean"), n=("tf", "size"))
    )
    print(summary.to_string())
    print(f"wrote {config.BASELINE_TABLE}")


if __name__ == "__main__":
    main()
