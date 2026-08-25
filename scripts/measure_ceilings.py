#!/usr/bin/env python
"""What the protein representation allows, before any question of how well a model trains.

    python scripts/measure_ceilings.py [--arm A1 --arm A4] [--fold S2/fold-0 ...]

`T35`. The three measurements [`docs/ML_RESULTS.md`](../docs/ML_RESULTS.md) §4 rests on, which
between them decide **which link is at fault** when the two-tower model loses to nearest-neighbour
lookup: the width of the shared space, the tower's strength, or the protein representation itself.

Roughly three minutes. Reads only cached artifacts — the 8-mer matrix, the embeddings and the
distances — and writes `reports/representation_ceiling.md`.

This is also the harness for `T34`: a new pooling scheme is compared by rebuilding
`embeddings/<arm>.npz` and re-running this, which reports it against the same ceilings on the
same folds.
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np
import pandas as pd

from snp2prot import corpus, distances, embeddings, experiment, splits
from snp2prot.config import PROCESSED_DIR, REPORTS_DIR
from snp2prot.data.matrix import KmerMatrix
from snp2prot.evaluation import ceilings

DEFAULT_OUT = REPORTS_DIR / "representation_ceiling.md"
GRID = PROCESSED_DIR / "training_folds.parquet"
C1_SET = PROCESSED_DIR / "c1_variants.parquet"

#: The folds the probes run on. Not all 19 — the probes cost minutes each and the question is
#: about the hard regimes, where the model and the baseline disagree most.
DEFAULT_FOLDS = ("S2/fold-0", "S2/fold-2", "S2/fold-4", "P1/Homeodomain")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--arm", action="append", dest="arms", help="repeatable; default A1 and A4")
    ap.add_argument("--fold", action="append", dest="folds", help="repeatable")
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = ap.parse_args()
    arms = args.arms or ["A1", "A4"]
    wanted = tuple(args.folds) if args.folds else DEFAULT_FOLDS

    started = time.time()
    domains = corpus.domains()
    matrix = KmerMatrix.load()
    dist = distances.DomainDistances.load()
    folds = {f.label: f for f in splits.all_regimes(domains, dist)}
    missing = [f for f in wanted if f not in folds]
    if missing:
        raise SystemExit(f"unknown folds: {missing}\nhave: {', '.join(folds)}")

    cfg = experiment.load()
    trainable = corpus.trainable(domains).to_numpy()
    excluded = domains.dbd_seq.isin(_c1_domains()).to_numpy()

    print("rank ceiling (oracle) ...", flush=True)
    ranks = ceilings.rank_ceiling(matrix)
    for rank, score in ranks.items():
        print(f"  rank {rank:>4}: {score:.4f}")

    rows = []
    for arm in arms:
        vectors = embeddings.DomainEmbeddings.load(arm)
        if list(vectors.domains) != list(domains.dbd_seq):
            raise SystemExit(f"{arm} embeddings were built for a different domain set")
        for label in wanted:
            fold = folds[label]
            pool = fold.train[trainable[fold.train] & ~excluded[fold.train]]
            inner = splits.Fold(fold.regime, fold.name, fold.test, pool, fold.held_out)
            train, validation = splits.validation_split(
                inner, domains, dist,
                float(cfg["splits"]["validation"]["fraction"]),
                int(cfg["splits"]["seed"]),
                str(cfg["splits"]["validation"]["grouping"]),
            )
            linear = ceilings.ridge_probe(matrix, vectors.vectors, train, validation, fold.test)
            kernel = ceilings.kernel_probe(matrix, vectors.vectors, train, validation, fold.test)
            rows.append(
                {
                    "arm": arm, "regime": fold.regime, "fold": fold.name,
                    "n_train": len(train), "n_test": len(fold.test),
                    "linear_ceiling": linear["ceiling"], "linear_honest": linear.get("honest"),
                    "kernel_ceiling": kernel["ceiling"], "kernel_honest": kernel.get("honest"),
                }
            )
            print(
                f"  {arm} {label:<16} linear {linear['ceiling']:.4f}"
                f" ({linear.get('honest', float('nan')):.4f})"
                f"  rbf {kernel['ceiling']:.4f} ({kernel.get('honest', float('nan')):.4f})"
                f"  [{time.time() - started:.0f}s]",
                flush=True,
            )

    frame = pd.DataFrame(rows)
    write_report(args.out, ranks, frame, cfg)
    print(f"wrote {args.out} in {time.time() - started:.0f}s")


def _c1_domains() -> set[str]:
    if not C1_SET.exists():
        raise SystemExit(f"no C1 set at {C1_SET}\nrun scripts/run_nn_baseline.py first")
    return set(pd.read_parquet(C1_SET).domain)


def _model_and_baseline() -> pd.DataFrame:
    """The grid's own numbers, so the ceilings are read next to what was achieved."""
    if not GRID.exists():
        return pd.DataFrame()
    return pd.read_parquet(GRID)[["arm", "regime", "fold", "aupr", "baseline_aupr"]]


def _fmt(value, places: int = 4) -> str:
    if value is None or not np.isfinite(value):
        return "n/a"
    return f"{value:.{places}f}"


def write_report(path: Path, ranks: dict, frame: pd.DataFrame, cfg: dict) -> None:
    grid = _model_and_baseline()
    merged = frame.merge(grid, on=["arm", "regime", "fold"], how="left") if len(grid) else frame

    lines = [
        "# What the protein representation allows",
        "",
        "Generated by `scripts/measure_ceilings.py` (`T35`). These are the measurements",
        "[`ML_RESULTS.md`](../docs/ML_RESULTS.md) §4 rests on. They answer one question —",
        "**when the",
        "two-tower model loses to nearest-neighbour lookup, which link is at fault?** — by",
        "bounding",
        "each link separately.",
        "",
        "## 1. Is the shared space wide enough?",
        "",
        "Best macro AUPR from a rank-`r` approximation of the E-score matrix. An **oracle**: the",
        "factors are fitted to the matrix being scored, so this bounds what a bilinear model of",
        "that",
        "width could represent and says nothing about generalisation.",
        "",
        "| rank | macro AUPR |",
        "|---:|---:|",
    ]
    for rank, score in ranks.items():
        marker = "  ← the configured `model.width`" if rank == int(cfg["model"]["width"]) else ""
        lines.append(f"| {rank} | {score:.4f}{marker} |")

    lines += [
        "",
        "## 2. Is the protein tower too weak, and would non-linearity help?",
        "",
        "**Linear** is a ridge regression from the pooled embedding straight onto all 32,896",
        "8-mers",
        "— unconstrained output, no shared space, no DNA tower — so it strictly upper-bounds any",
        "model whose protein side is a linear map of the same vector, ours included. **RBF** is",
        "the",
        "same in kernel form, which isolates non-linearity from width and depth.",
        "",
        "Each is given twice. **ceiling** picks its regularisation on the *test* fold and is an",
        "oracle — an upper bound on the model class, not an achievable score. **honest** picks",
        "it on",
        "the fold's own validation slice, which is what a real run would get. A bound that is too",
        "tight cannot rule a class out, which is why the oracle is the one to compare against.",
        "",
        "| arm | fold | linear ceiling | (honest) | RBF ceiling | (honest) | two-tower "
        "| NN baseline |",
        "|---|---|---:|---:|---:|---:|---:|---:|",
    ]
    for row in merged.itertuples():
        model = _fmt(getattr(row, "aupr", None))
        base = _fmt(getattr(row, "baseline_aupr", None))
        lines.append(
            f"| `{row.arm}` | {row.regime}/`{row.fold}` | {_fmt(row.linear_ceiling)} | "
            f"{_fmt(row.linear_honest)} | {_fmt(row.kernel_ceiling)} | {_fmt(row.kernel_honest)} | "
            f"**{model}** | {base} |"
        )

    lines += [
        "",
        "## How to read it",
        "",
        "- **two-tower ≈ linear ceiling** — the tower is already extracting what a linear map of",
        "  this embedding can extract. Widening or deepening it is not the lever.",
        "- **RBF ≈ linear** — non-linearity in the protein map is not the missing ingredient",
        "  either.",
        "- **rank ceiling far above everything else** — the shared space is not the constraint.",
        "",
        "All three together leave the **protein representation itself**: what a mean-pooled",
        "protein-LM vector does and does not carry. That is the conclusion `ML_RESULTS.md` §4.4",
        "draws, and `T34` is the experiment that tests it.",
        "",
        "```bash",
        "python scripts/measure_ceilings.py",
        "```",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
