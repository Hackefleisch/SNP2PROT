#!/usr/bin/env python
"""Press the Pfam-A HMM library into HMMER's binary format. Run once after downloading it.

    python scripts/press_pfam.py

`Pfam-A.hmm` is 2.2 GB of text. Unpressed, every call to `snp2prot.domains.scan` has to
parse all 30,134 profiles from scratch — 19.5 s per source against roughly 0.1 s of actual
scanning, so a 30-source build spent about ten minutes re-reading one unchanging file.

Pressing writes four binary siblings (`.h3f`, `.h3i`, `.h3m`, `.h3p`) holding profiles
already optimized for the search pipeline. `HMMFile` picks them up automatically — nothing
in the codebase has to be told they exist — and the load drops to ~0.5 s with each scan
after it costing ~0.6 s.

This changes no result. The domain calls are identical either way, which is checked by
`tests/test_domains.py::test_pressed_and_unpressed_libraries_agree`.

Costs about 2.5 GB of disk next to the library, all of it git-ignored and derived. Delete
the four `.h3*` files to undo this completely.
"""

from __future__ import annotations

import sys
import time

import pyhmmer
from pyhmmer.plan7 import HMMFile

from snp2prot.config import EXTERNAL_DIR

PFAM = EXTERNAL_DIR / "pfam" / "Pfam-A.hmm"
SUFFIXES = (".h3f", ".h3i", ".h3m", ".h3p")


def main() -> None:
    if not PFAM.exists():
        raise SystemExit(f"Pfam library not found: {PFAM}")

    existing = [PFAM.with_suffix(PFAM.suffix + s) for s in SUFFIXES]
    if all(p.exists() for p in existing):
        with HMMFile(PFAM) as f:
            if f.is_pressed():
                total = sum(p.stat().st_size for p in existing)
                print(f"already pressed ({total / 1e9:.1f} GB alongside {PFAM.name})")
                return

    print(f"pressing {PFAM} ({PFAM.stat().st_size / 1e9:.1f} GB) ...")
    t = time.perf_counter()
    with HMMFile(PFAM) as f:
        pyhmmer.hmmer.hmmpress(f, PFAM)
    written = sum(p.stat().st_size for p in existing if p.exists())
    print(f"done in {time.perf_counter() - t:.0f}s — wrote {written / 1e9:.1f} GB")

    with HMMFile(PFAM) as f:
        if not f.is_pressed():
            print("WARNING: the library still does not report as pressed", file=sys.stderr)
            raise SystemExit(1)
    print("scans will now load the library in ~0.5s instead of ~19.5s")


if __name__ == "__main__":
    main()
