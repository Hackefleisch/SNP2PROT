#!/usr/bin/env python
"""Align every pair of the corpus's domains and cache the result.

    python scripts/build_distances.py [--processes 24]

Phase 7, step 0. 894,453 pairs at 0.17 ms is 2.5 minutes on one core and well under a minute
on this machine's 24; the output is a 4 MB `.npz` that three consumers read — the
nearest-neighbour baseline, the `S2` split grouping and the nearest-neighbour-identity
histogram (`docs/ML_PLAN.md` §8.1).

Reports the connected components the `S2` regime holds out, which is the measurement `T27`
turned on: single linkage is only safe here if it does not chain the corpus into a blob.
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np

from snp2prot import corpus, distances, experiment, thresholds
from snp2prot.config import DISTANCE_MATRIX


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--processes", type=int, default=None, help="default: every core")
    ap.add_argument("--out", type=Path, default=DISTANCE_MATRIX)
    args = ap.parse_args()

    min_overlap = float(thresholds.load()["cluster"]["min_overlap"])
    min_identity = float(experiment.section("splits")["s2_min_identity"])

    domains = corpus.domains()
    n = len(domains)
    print(f"{n} domains, {n * (n - 1) // 2:,} pairs")

    started = time.time()
    last = [0.0]

    def progress(done: int) -> None:
        now = time.time()
        if now - last[0] > 10:
            last[0] = now
            print(f"  {done}/{n - 1} rows, {now - started:.0f}s", flush=True)

    d = distances.build(domains.dbd_seq.to_numpy(), processes=args.processes, progress=progress)
    print(f"aligned in {time.time() - started:.0f}s")

    identity = d.identity()
    off = ~np.eye(n, dtype=bool)
    guard = d.comparable(min_overlap) & off
    print(f"  pairs passing the {min_overlap:.0%} overlap guard: {guard.sum() // 2:,}")
    print(f"  identity, guarded pairs: median {np.nanmedian(identity[guard]):.3f}")
    print(f"  identity, all pairs:     median {np.nanmedian(identity[off]):.3f}")

    families = domains.dbd_family.to_numpy()
    labels = distances.connected_components(d, min_identity, min_overlap, families)
    sizes = np.bincount(labels)
    print(
        f"  S2 components at >= {min_identity:.0%} identity, >= {min_overlap:.0%} overlap: "
        f"{len(sizes)}, largest {sizes.max()} domains"
    )
    merged = sum(1 for c in range(len(sizes)) if domains.wt_id[labels == c].nunique() > 1)
    print(f"  components merging more than one cluster: {merged}")
    # Single linkage at a loose floor chains, which is the safe direction for a split but has
    # to be visible: a component whose members are mostly NOT within the floor of each other
    # is a chain of near-relatives, not a group of near-identical domains.
    largest = np.flatnonzero(labels == int(np.argmax(sizes)))
    inside = identity[np.ix_(largest, largest)][~np.eye(len(largest), dtype=bool)]
    print(
        f"  largest component: {len(largest)} domains, median internal identity "
        f"{np.nanmedian(inside):.3f}, {np.nanmean(inside >= min_identity):.0%} of its pairs "
        f"at the floor"
    )

    path = d.save(args.out)
    print(f"wrote {path} ({path.stat().st_size / 1e6:.1f} MB)")


if __name__ == "__main__":
    main()
