#!/usr/bin/env python
"""Train every GHT fold at every seed and write the per-fold and per-TF tables.

    python scripts/run_ght_grid.py --seeds 3                  # the full grid
    python scripts/run_ght_grid.py --regime C1 --seeds 1      # just Setting 1
    python scripts/run_ght_grid.py --resume                   # skip (arm, fold, seed) already done

`docs/GHT_PLAN.md` §8-§9. Eleven folds — `C1/all`, `G1/fold-0..4`, `G2/fold-0..4` — at
`--seeds` model seeds each. Wall-clock comes from `reports/ght_preflight_*.md`; under rule 10
the full grid is the owner's to launch, and this is the command.

**Resumable.** Each completed run appends to `data/processed/ght/ght_folds.parquet` and
`ght_tfs.parquet` immediately, so an interrupted grid restarts where it stopped rather than from
zero. Deleting those two files restarts it.
"""

from __future__ import annotations

import argparse
import time

import pandas as pd
import torch

from snp2prot import experiment
from snp2prot.ght import config, data, panel, splits, training


def existing() -> pd.DataFrame:
    if config.FOLD_TABLE.exists():
        return pd.read_parquet(config.FOLD_TABLE)
    return pd.DataFrame(columns=["arm", "regime", "fold", "model_seed"])


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--arm", action="append", default=None, help="default: A1")
    ap.add_argument("--regime", action="append", default=None, help="C1, G1, G2; default: all")
    ap.add_argument("--seeds", type=int, default=3)
    ap.add_argument(
        "--protein-mode",
        action="append",
        default=None,
        help="real | constant | shuffled (repeatable); default: real. "
        "`constant` is the protein-blind control and is what makes a Setting 1 number readable",
    )
    ap.add_argument(
        "--tower",
        default=None,
        choices=["pooled", "residue"],
        help="override ght.model.protein.tower; `residue` is the GHT_PLAN section 6 experiment",
    )
    ap.add_argument(
        "--pbm-protein",
        default=None,
        help="start the protein tower from a PBM checkpoint's tower "
        "(data/processed/checkpoints/*.pt) — rung 5 of the fallback ladder",
    )
    ap.add_argument(
        "--freeze-protein",
        action="store_true",
        help="with --pbm-protein: hold the PBM tower fixed so only the DNA tower adapts",
    )
    ap.add_argument("--steps", type=int, default=None, help="override ght.training.steps")
    ap.add_argument("--device", default=None)
    ap.add_argument("--resume", action="store_true")
    ap.add_argument("--no-track", action="store_true")
    args = ap.parse_args()

    arms = args.arm or ["A1"]
    modes = args.protein_mode or ["real"]
    cfg = experiment.section("ght")
    if args.steps:
        cfg["training"]["steps"] = args.steps
    if args.tower:
        cfg["model"]["protein"]["tower"] = args.tower
    tower = str(cfg["model"]["protein"].get("tower", "pooled"))
    if args.pbm_protein:
        cfg["model"]["protein"]["init_from"] = args.pbm_protein
        cfg["model"]["protein"]["freeze"] = bool(args.freeze_protein)
        # The transfer is recorded in `protein_mode`, not in a separate column: every row of the
        # results table already has to say what the protein tower was fed, and a run started from
        # the PBM arm's tower is a different answer to that same question.
        modes = [f"pbm_{'frozen' if args.freeze_protein else 'finetuned'}"]

    table = panel.load()
    corpus = data.load(table.tf)
    folds = splits.all_regimes(
        table,
        corpus.chromosomes("Train"),
        corpus.chromosomes("Test"),
        corpus.window_counts("Train"),
        int(cfg["splits"]["n_folds"]),
        int(cfg["splits"]["seed"]),
        float(cfg["splits"]["validation_fraction"]),
        float(cfg["splits"]["s2_min_identity"]),
    )
    if args.regime:
        folds = [f for f in folds if f.regime in set(args.regime)]

    columns = ["arm", "regime", "fold", "model_seed", "protein_mode", "protein_tower"]
    done = existing() if args.resume else pd.DataFrame(columns=columns)
    for column in ("protein_mode", "protein_tower"):
        if column not in done:
            done[column] = "real" if column == "protein_mode" else "pooled"
    seen = set(
        zip(
            done.arm,
            done.regime,
            done.fold,
            done.model_seed,
            done.protein_mode,
            done.protein_tower,
            strict=True,
        )
    )
    base_seed = int(cfg["model"]["seed"])
    seeds = [base_seed + i for i in range(args.seeds)]

    total = len(arms) * len(folds) * len(seeds) * len(modes)
    print(
        f"{total} runs: {len(arms)} arms x {len(folds)} folds x {len(seeds)} seeds "
        f"x {len(modes)} protein modes, {tower} tower"
    )
    print(f"steps {cfg['training']['steps']:,}; skipping {len(seen)} already recorded")
    started = time.time()
    n = 0

    for arm in arms:
        if tower == "residue":
            base, mask, model_name = data.load_residue_embeddings(arm, table)
        else:
            base, model_name = data.load_embeddings(arm, table)
            mask = None
        for mode in modes:
            # The ablation is applied to the embedding table, before any fold sees it, so every
            # fold of one mode is fed the same mis-assignment rather than a fresh one.
            # A `pbm_*` mode names where the TOWER came from, not what it is fed, so the
            # embeddings pass through unablated.
            vectors = data.ablate(
                base, "real" if mode.startswith("pbm_") else mode, seed=int(cfg["splits"]["seed"])
            )
            for fold in folds:
                for seed in seeds:
                    n += 1
                    key = (arm, fold.regime, fold.name, seed, mode, tower)
                    if key in seen:
                        print(f"[{n}/{total}] {arm} {fold.label} seed{seed} {mode} — done")
                        continue
                    print(f"[{n}/{total}] {arm} {fold.label} seed{seed} {mode}/{tower}", flush=True)
                    summary, per_tf = training.run_fold(
                        fold,
                        arm,
                        corpus,
                        vectors,
                        model_name,
                        cfg,
                        device=torch.device(args.device) if args.device else None,
                        track=not args.no_track,
                        model_seed=seed,
                        protein_mask=mask,
                        protein_mode=mode,
                    )
                    # A fresh Trainer per run allocates its own GPU mirror of the resident
                    # windows; 33 of them in one process is the kind of thing that OOMs an 8 GB
                    # card through fragmentation rather than through peak use.
                    if torch.cuda.is_available():
                        torch.cuda.empty_cache()
                    _append(config.FOLD_TABLE, pd.DataFrame([summary]))
                    if not per_tf.empty:
                        _append(config.TF_TABLE, per_tf)
                    print(
                        f"      auPRC {summary['aupr']:.4f} "
                        f"(chance {summary['chance_aupr']:.4f}), "
                        f"auROC {summary['auroc']:.4f}, "
                        f"random auROC {summary.get('auroc_random', float('nan')):.4f}, "
                        f"aliens auROC {summary.get('auroc_aliens', float('nan')):.4f}, "
                        f"{summary['seconds'] / 60:.1f} min  "
                        f"[{(time.time() - started) / 60:.0f} min elapsed]",
                        flush=True,
                    )

    print(f"\ndone in {(time.time() - started) / 60:.0f} min")
    print(f"wrote {config.FOLD_TABLE} and {config.TF_TABLE}")


def _append(path, frame: pd.DataFrame) -> None:
    """Write-through after every run, so an interrupted grid keeps what it earned."""
    path.parent.mkdir(parents=True, exist_ok=True)
    if path.exists():
        frame = pd.concat([pd.read_parquet(path), frame], ignore_index=True)
    frame.to_parquet(path, index=False)


if __name__ == "__main__":
    main()
