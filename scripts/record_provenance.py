#!/usr/bin/env python
"""Append a row to PROVENANCE.md for a file under data/raw/.

Usage:
    python scripts/record_provenance.py data/raw/<source>/<file> \
        --url <where it came from> --accession <id> --desc "one line"

Computes size and sha256 itself; the descriptive fields are yours to supply, because a
wrong description is worse than a missing one.
"""

from __future__ import annotations

import argparse
import hashlib
from datetime import datetime, timezone
from pathlib import Path

from snp2prot.config import PROVENANCE_FILE, RAW_DIR


def sha256(path: Path, chunk: int = 1 << 20) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        while block := fh.read(chunk):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("path", type=Path, help="file under data/raw/")
    ap.add_argument("--url", default="", help="exact URL the file came from")
    ap.add_argument("--accession", default="", help="dataset accession, if any")
    ap.add_argument("--desc", default="", help="one line: what this file contains")
    args = ap.parse_args()

    path = args.path.resolve()
    if not path.is_file():
        raise SystemExit(f"not a file: {path}")
    try:
        rel = path.relative_to(RAW_DIR)
    except ValueError:
        raise SystemExit(f"{path} is not under {RAW_DIR}") from None

    source = rel.parts[0]
    row = " | ".join(
        [
            "",
            source,
            f"`{rel}`",
            args.url or "TODO",
            args.accession or "-",
            datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
            f"{path.stat().st_size:,}",
            f"`{sha256(path)}`",
            args.desc or "TODO",
            "",
        ]
    ).strip()

    print(row)
    print(f"\nAppend the line above to the 'Downloaded files' table in {PROVENANCE_FILE}")


if __name__ == "__main__":
    main()
