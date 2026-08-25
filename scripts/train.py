#!/usr/bin/env python
"""Train the two-tower model on one arm and one fold.

    python scripts/train.py --arm A1 --fold P3/all [--steps 500] [--no-track]

Phase 7, build step 1. This is the single-run entry point, for iterating on the design and for
checking a fold against a number with a known answer; `scripts/run_grid.py` runs the whole grid.

**`P3/all` is the fold to try first.** With no variant in training, the nearest-neighbour
baseline copies each held-out variant's own wild type and scores **0.928** — so it is the one
fold where a wrong answer is obvious immediately (`reports/nn_baseline.md`).

The design is `docs/TRAINING.md`: multi-positive InfoNCE over the complete 8-mer axis with a
null anchor, batches drawn on the protein axis, everything resident on the device.
"""

from __future__ import annotations

import argparse

import pandas as pd

from snp2prot import corpus, distances, embeddings, experiment, splits, training
from snp2prot.config import PROCESSED_DIR
from snp2prot.data.matrix import KmerMatrix

C1_SET = PROCESSED_DIR / "c1_variants.parquet"


def load_everything(arm: str):
    """The cached artifacts every run reads, and the two training-pool masks."""
    domains = corpus.domains()
    matrix = KmerMatrix.load()
    vectors = embeddings.DomainEmbeddings.load(arm)
    dist = distances.DomainDistances.load()

    order = list(domains.dbd_seq)
    for name, table in (("8-mer matrix", matrix), ("distances", dist), (arm, vectors)):
        if list(table.domains) != order:
            raise SystemExit(f"the cached {name} was built for a different domain set — rebuild it")

    trainable = corpus.trainable(domains).to_numpy()
    excluded = domains.dbd_seq.isin(_c1_domains()).to_numpy()
    return domains, matrix, vectors, dist, trainable, excluded


def _c1_domains() -> set[str]:
    """The C1 evaluation set, which is never trained on (`D6`)."""
    if not C1_SET.exists():
        raise SystemExit(f"no C1 evaluation set at {C1_SET}\nrun scripts/run_nn_baseline.py first")
    return set(pd.read_parquet(C1_SET).domain)


def add_common_arguments(ap: argparse.ArgumentParser) -> None:
    """Flags both entry points share. `--arm` is not among them: this script takes exactly one
    and `run_grid.py` takes a repeatable list, so each declares its own."""
    ap.add_argument("--steps", type=int, default=None, help="override training.steps")
    ap.add_argument("--device", default=None, help="default: cuda when available")
    ap.add_argument("--no-track", action="store_true", help="do not write to mlruns/")


def apply_overrides(config: dict, args) -> dict:
    if args.steps is not None:
        config["training"]["steps"] = args.steps
    return config


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    add_common_arguments(ap)
    ap.add_argument("--arm", default="A1", help="which protein embedding arm")
    ap.add_argument("--fold", default="P3/all", help="regime/fold, e.g. S2/fold-0")
    args = ap.parse_args()

    config = apply_overrides(experiment.load(), args)
    domains, matrix, vectors, dist, trainable, excluded = load_everything(args.arm)

    folds = {f.label: f for f in splits.all_regimes(domains, dist)}
    if args.fold not in folds:
        raise SystemExit(f"unknown fold {args.fold!r}; have {', '.join(folds)}")
    fold = folds[args.fold]

    print(f"{args.arm} ({vectors.model}) on {fold.label} — {fold.held_out}")
    summary, per_domain = training.run_fold(
        fold,
        args.arm,
        domains,
        matrix,
        vectors,
        dist,
        config,
        trainable,
        excluded,
        device=training.device_for(args.device),
        track=not args.no_track,
    )

    print(
        f"  trained {int(summary['steps_run'])} steps in {summary['seconds']:.0f}s, "
        f"best at {int(summary['best_step'])}"
    )
    print(f"  validation AUPR {summary['validation_aupr']:.4f}")
    print(f"  test AUPR       {summary['aupr']:.4f}  (median {summary['aupr_median']:.4f})")
    print(
        f"  test AUROC {summary['auroc']:.3f}  Spearman {summary['spearman']:.3f}  "
        f"P@50 {summary['precision_at_50']:.3f}"
    )
    print(
        f"  parameters: protein {int(summary['n_params_protein_tower']):,}  "
        f"dna {int(summary['n_params_dna_tower']):,}"
    )
    worst = per_domain.nsmallest(3, "aupr")[["gene", "family", "n_pos", "aupr"]]
    print(f"  worst-predicted held-out domains:\n{worst.to_string(index=False)}")


if __name__ == "__main__":
    main()
