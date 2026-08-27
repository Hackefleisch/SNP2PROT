"""Train one fold and dump the shared protein-DNA space, for `ML_PLAN.md` 5.2's figure.

Reproduces `training.run_fold`'s path exactly — same training pool, same validation split, same
seed — but keeps the trained model, which `run_fold` does not return, so both towers can be
projected into the space they share.

    .venv/bin/python docs/talk/dump_shared_space.py     # ~3 min on the GPU

Writes `data/processed/talk_shared_space.npz` — 32 MB, so it lives with the other
cached modelling artifacts rather than in `docs/` — which `make_shared_space_figure.py` reads.
Prints the fold's test AUPR next to the recorded one, so a drift in the code shows up
immediately rather than silently changing the picture: run on 2026-08-26 against commit
9390e41 it gave **0.8727** where `reports/training.md` records **0.8754**.

**This is a figure script, not a result.** A UMAP of a contrastively-trained space shows
clusters almost by construction (`ML_PLAN.md` 5.2) — it illustrates, it does not evidence.
"""

from __future__ import annotations

import sys
import time

import numpy as np
import torch

sys.path.insert(0, "src")
sys.path.insert(0, "scripts")

from train import load_everything  # noqa: E402

from snp2prot import experiment, splits, training  # noqa: E402
from snp2prot.config import PROCESSED_DIR  # noqa: E402
from snp2prot.evaluation import metrics  # noqa: E402

ARM, FOLD = "A4", "P3/all"

domains, matrix, vectors, dist, trainable, excluded = load_everything(ARM)
config = experiment.load()
folds = {f.label: f for f in splits.all_regimes(domains, dist)}
fold = folds[FOLD]

split_cfg = config["splits"]
seed = int(split_cfg["seed"])
pool = fold.train[trainable[fold.train] & ~excluded[fold.train]]
inner = splits.Fold(fold.regime, fold.name, fold.test, pool, fold.held_out)
train_rows, validation_rows = splits.validation_split(
    inner,
    domains,
    dist,
    float(split_cfg["validation"]["fraction"]),
    seed,
    str(split_cfg["validation"]["grouping"]),
)
print(
    f"{ARM} on {fold.label}: train {len(train_rows)}, val {len(validation_rows)}, "
    f"test {len(fold.test)}",
    flush=True,
)

t0 = time.time()
trainer = training.Trainer(matrix, vectors, config, device=training.device_for(None), seed=seed)
result = trainer.train(train_rows, validation_rows, seed=seed)
print(
    f"  {result.steps_run} steps, best at {result.best_step}, {time.time() - t0:.0f}s", flush=True
)

predicted = trainer.predict(fold.test)
scored = [
    metrics.score_domain(
        matrix.label[row], matrix.escore[row], predicted[i], domain=str(matrix.domains[row])
    )
    for i, row in enumerate(fold.test)
]
summary = metrics.macro_average(scored)
print(f"  test AUPR {summary['aupr']:.4f}   (reports/training.md records 0.8754)", flush=True)

# project both towers into the shared space
trainer.model.eval()
with torch.no_grad():
    prot = trainer.model.protein(trainer.proteins).float().cpu().numpy()
    table = trainer.model.dna_table(trainer.tokens).float().cpu().numpy()

np.savez_compressed(
    PROCESSED_DIR / "talk_shared_space.npz",
    protein=prot,
    dna=table[:-1],
    null=table[-1],
    domains=matrix.domains,
    kmers=matrix.kmers,
    aupr=summary["aupr"],
)
print("  wrote _shared_space.npz", prot.shape, table.shape, flush=True)
