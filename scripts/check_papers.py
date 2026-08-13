#!/usr/bin/env python
"""Report which source-paper PDFs are present in docs/papers/ and which are missing.

The expected filenames are parsed out of the manifest table in docs/papers/README.md, so
the README stays the single source of truth and this script cannot drift from it.

    python scripts/check_papers.py
    python scripts/check_papers.py --phase 1     # only what a given phase needs
"""

from __future__ import annotations

import argparse
import re

from snp2prot.config import DOCS_DIR

PAPERS_DIR = DOCS_DIR / "papers"
MANIFEST = PAPERS_DIR / "README.md"

#: `| `name.pdf` | ! P1 | Paper cite | [doi](url) |`
ROW = re.compile(
    r"^\|\s*`(?P<file>[^`]+\.pdf)`\s*\|\s*(?P<prio>!?)\s*P(?P<phase>\d)\+?\s*\|"
    r"\s*(?P<cite>[^|]+?)\s*\|"
)


def expected() -> list[dict[str, str]]:
    if not MANIFEST.exists():
        raise SystemExit(f"manifest not found: {MANIFEST}")
    out = []
    for line in MANIFEST.read_text().splitlines():
        if m := ROW.match(line):
            out.append(m.groupdict())
    if not out:
        raise SystemExit(f"no manifest rows parsed from {MANIFEST} — has the table changed?")
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--phase", type=int, help="only papers needed by this phase or earlier")
    args = ap.parse_args()

    rows = expected()
    if args.phase is not None:
        rows = [r for r in rows if int(r["phase"]) <= args.phase]

    present = {p.name for p in PAPERS_DIR.glob("*.pdf")}
    missing = []

    for r in rows:
        mark = "ok     " if r["file"] in present else "MISSING"
        flag = "!" if r["prio"] else " "
        print(f"{mark} {flag} P{r['phase']}  {r['file']:<44} {r['cite']}")
        if r["file"] not in present:
            missing.append(r)

    unlisted = sorted(present - {r["file"] for r in expected()})
    for name in unlisted:
        print(f"UNLISTED  --  {name:<44} not in the manifest — rename it or add a row")

    blocking = [r for r in missing if r["prio"]]
    print(f"\n{len(rows) - len(missing)}/{len(rows)} present, {len(blocking)} blocking")
    if blocking:
        print("blocking the next phase: " + ", ".join(r["file"] for r in blocking))


if __name__ == "__main__":
    main()
