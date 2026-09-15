#!/usr/bin/env python
"""Embed every domain with a protein language model and cache the pooled vectors.

    python scripts/build_embeddings.py --arm A1 [--batch-tokens 8192] [--device cuda]
    python scripts/build_embeddings.py --arm A1 --per-residue      # the un-pooled cache (T36)

Phase 7, build step 1 (`docs/ML_PLAN.md` §9). A one-off inference over the corpus's 1,338
canonical domains, mean-pooled to one vector per domain and written to
`data/processed/embeddings/<arm>.npz`. The language model never runs again — the training loop
reads this file (`ML_PLAN.md` §3.1).

**The inference itself lives in `snp2prot.embeddings.encode`**, not here, since 2026-09-15: the
GHT arm embeds its own 33-domain panel with the same model and the same pooling, and a copy of
the loader in a second script is how two arms quietly stop being comparable.

**Weights are cached under `data/external/models/`**, not in the user's home directory, so that
a downloaded artifact the build depends on sits with the Pfam HMMs and the UniProt cache and can
carry a `PROVENANCE.md` row like everything else.

Order and identity: rows come out in `corpus.domains()` order and the domain list is written
into the file, so an embedding table built either side of a rebuild cannot be paired with a
mismatched matrix silently.

**`--per-residue` writes the un-pooled cache instead** — `<arm>_residues.npz`, every residue's
vector rather than their mean, about 600 MB. It is the same forward pass; the pooled run simply
averages the result and throws the rest away. It exists because that averaging was measured to
destroy the single-residue signal the project's central claim depends on
(`docs/ML_RESULTS.md` §9.3), and `scripts/check_unpooling.py` reads this file to find out whether
the signal survives before the mean is taken.
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np

from snp2prot import corpus, embeddings
from snp2prot.config import embedding_table, residue_embedding_table


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--arm", default="A1", choices=sorted(embeddings.ARMS))
    ap.add_argument("--batch-tokens", type=int, default=8192)
    ap.add_argument("--device", default=None, help="default: cuda when available")
    ap.add_argument("--out", type=Path, default=None)
    ap.add_argument(
        "--per-residue",
        action="store_true",
        help="write the un-pooled cache <arm>_residues.npz instead of the pooled table (T36)",
    )
    args = ap.parse_args()

    domains = corpus.domains()
    sequences = [str(s) for s in domains.dbd_seq]
    print(f"{len(sequences)} domains, {min(map(len, sequences))}-{max(map(len, sequences))} aa")
    print(f"model {embeddings.ARMS[args.arm]}")
    started = time.time()

    def progress(done: int, total: int) -> None:
        if done % 200 < 32:
            print(f"  {done}/{total}  ({time.time() - started:.0f}s)", flush=True)

    result = embeddings.encode(
        sequences,
        arm=args.arm,
        device=args.device,
        batch_tokens=args.batch_tokens,
        per_residue=args.per_residue,
        progress=progress,
    )

    if args.per_residue:
        offsets = np.zeros(len(result) + 1, dtype=np.int64)
        np.cumsum([len(b) for b in result], out=offsets[1:])
        table = embeddings.ResidueEmbeddings(
            arm=args.arm,
            model=embeddings.ARMS[args.arm],
            domains=np.asarray(sequences, dtype=np.str_),
            offsets=offsets,
            vectors=np.concatenate(result, axis=0),
        )
        path = table.save(args.out or residue_embedding_table(args.arm))
        print(f"{table.width}-d over {offsets[-1]:,} residues in {len(result)} domains")
    else:
        table = embeddings.DomainEmbeddings(
            arm=args.arm,
            model=embeddings.ARMS[args.arm],
            domains=np.asarray(sequences, dtype=np.str_),
            vectors=result,
        )
        path = table.save(args.out or embedding_table(args.arm))
        norms = np.linalg.norm(result, axis=1)
        print(f"{table.width}-d, norms {norms.min():.2f}-{norms.max():.2f}")
    print(f"wrote {path} ({path.stat().st_size / 1e6:.0f} MB) in {time.time() - started:.0f}s")


if __name__ == "__main__":
    main()
