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

from snp2prot import corpus, distances, embeddings, experiment, splits, training
from snp2prot.data.matrix import KmerMatrix


def load_everything(arm: str):
    """The cached artifacts every run reads, and the training-pool mask.

    One mask, not two. The C1 set used to be masked out of training as well; it is not any more,
    because the nearest-neighbour baseline never applied the same mask and the comparison against
    it was therefore between two different training pools — see `snp2prot.training.run_fold`.
    """
    domains = corpus.domains()
    matrix = KmerMatrix.load()
    vectors = embeddings.DomainEmbeddings.load(arm)
    dist = distances.DomainDistances.load()

    try:
        corpus.require_aligned(
            domains.dbd_seq, matrix=matrix.domains, distances=dist.domains, **{arm: vectors.domains}
        )
    except ValueError as mismatch:
        raise SystemExit(str(mismatch)) from mismatch

    trainable = corpus.trainable(domains).to_numpy()
    return domains, matrix, vectors, dist, trainable


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
    domains, matrix, vectors, dist, trainable = load_everything(args.arm)

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
        device=training.device_for(args.device),
        track=not args.no_track,
    )

    print(
        f"  trained {int(summary['steps_run'])} steps in {summary['seconds']:.0f}s, "
        f"best at {int(summary['best_step'])}"
    )
    print(f"  validation AUPR {summary['validation_aupr']:.4f}")
    chance, p95 = summary["chance_aupr"], summary["null_p95"]
    lift = summary["aupr"] / chance
    verdict = (
        "AT CHANCE"
        if summary["aupr"] <= p95
        else f"{lift:.0f}x chance"
        if lift >= 10
        else f"{lift:.1f}x chance"
    )
    print(f"  test AUPR       {summary['aupr']:.4f}  (median {summary['aupr_median']:.4f})")
    print(f"  chance {chance:.4f}, null 95th {p95:.4f}  ->  {verdict}")
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
