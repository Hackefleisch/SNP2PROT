#!/usr/bin/env python
"""Embed the GHT panel's domains with a protein language model and cache the vectors.

    python scripts/build_ght_embeddings.py --arm A1                 # ~10 s on the GPU
    python scripts/build_ght_embeddings.py --arm A1 --per-residue   # the un-pooled cache

Phase A of `docs/GHT_PLAN.md` §11. The same inference `scripts/build_embeddings.py` runs for the
PBM corpus — same code, via `snp2prot.embeddings.encode` — over a different set of 33 domains,
written under `data/processed/ght/embeddings/`.

Separate files rather than one shared table because the two corpora are different objects: the
PBM table is keyed by `corpus.domains()` order and a GHT panel domain need not appear in it at
all. The **model** is shared, which is the point; the cache is not.
"""

from __future__ import annotations

import argparse
import time

import numpy as np

from snp2prot import embeddings
from snp2prot.ght import config, panel


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--arm", default="A1", choices=sorted(embeddings.ARMS))
    ap.add_argument("--batch-tokens", type=int, default=8192)
    ap.add_argument("--device", default=None, help="default: cuda when available")
    ap.add_argument("--per-residue", action="store_true")
    args = ap.parse_args()

    table = panel.load()
    sequences = [str(s) for s in table.dbd_seq]
    print(f"{len(sequences)} domains, {min(map(len, sequences))}-{max(map(len, sequences))} aa")
    print(f"model {embeddings.ARMS[args.arm]}")
    started = time.time()

    result = embeddings.encode(
        sequences,
        arm=args.arm,
        device=args.device,
        batch_tokens=args.batch_tokens,
        per_residue=args.per_residue,
    )
    config.EMBEDDING_DIR.mkdir(parents=True, exist_ok=True)

    if args.per_residue:
        offsets = np.zeros(len(result) + 1, dtype=np.int64)
        np.cumsum([len(b) for b in result], out=offsets[1:])
        cache = embeddings.ResidueEmbeddings(
            arm=args.arm,
            model=embeddings.ARMS[args.arm],
            domains=np.asarray(sequences, dtype=np.str_),
            offsets=offsets,
            vectors=np.concatenate(result, axis=0),
        )
        path = cache.save(config.residue_embedding_table(args.arm))
        print(f"{cache.width}-d over {offsets[-1]:,} residues in {len(result)} domains")
    else:
        cache = embeddings.DomainEmbeddings(
            arm=args.arm,
            model=embeddings.ARMS[args.arm],
            domains=np.asarray(sequences, dtype=np.str_),
            vectors=result,
        )
        path = cache.save(config.embedding_table(args.arm))
        norms = np.linalg.norm(result, axis=1)
        print(f"{cache.width}-d, norms {norms.min():.2f}-{norms.max():.2f}")
    print(f"wrote {path} ({path.stat().st_size / 1e6:.0f} MB) in {time.time() - started:.0f}s")


if __name__ == "__main__":
    main()
