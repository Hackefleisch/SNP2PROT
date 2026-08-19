#!/usr/bin/env python
"""Turn the merged table into the dense domain x 8-mer matrix modelling reads.

    python scripts/build_matrix.py [--out data/processed/kmer_matrix.npz]

Phase 7, step 0. One pass over `data/processed/training.parquet`, writing an `.npz` holding
the E-score matrix, the label matrix and the two index vectors. About 40 seconds and 250 MB
on disk; see `snp2prot.data.matrix` for why the shape change is worth a file.

Row order is `corpus.domains()`, which is `dbd_seq`-sorted, and is written into the artifact
so that a matrix and a distance matrix built either side of a rebuild cannot be paired
silently.
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np

from snp2prot import corpus
from snp2prot.config import KMER_MATRIX, MERGED_TABLE
from snp2prot.data import matrix


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--table", type=Path, default=MERGED_TABLE)
    ap.add_argument("--out", type=Path, default=KMER_MATRIX)
    args = ap.parse_args()

    started = time.time()
    domains = corpus.domains()
    print(f"{len(domains)} domains, {domains.dbd_family.nunique()} families")

    m = matrix.build(domains.dbd_seq.to_numpy(), args.table)
    n_domains, n_kmers = m.shape
    print(f"matrix {n_domains} x {n_kmers} = {n_domains * n_kmers:,} cells")

    counts = {int(v): int(c) for v, c in zip(*np.unique(m.label, return_counts=True), strict=True)}
    total = n_domains * n_kmers
    for value, name in ((1, "positive"), (0, "negative"), (-1, "no-call")):
        n = counts.get(value, 0)
        print(f"  {name:>9}: {n:>12,}  ({n / total:.3%})")
    print(f"  E-score range: {np.nanmin(m.escore):.3f} .. {np.nanmax(m.escore):.3f}")

    path = m.save(args.out)
    print(f"wrote {path} ({path.stat().st_size / 1e6:.0f} MB) in {time.time() - started:.0f}s")


if __name__ == "__main__":
    main()
