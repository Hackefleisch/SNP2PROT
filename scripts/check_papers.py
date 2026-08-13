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
INBOX_DIR = DOCS_DIR / "papers_inbox"
MANIFEST = PAPERS_DIR / "README.md"

#: `| `name.pdf` | ! P1 | Paper cite | [doi](url) |` — the DOI cell may be empty, and the
#: extension is not assumed to be .pdf (one supplement is a .docx).
ROW = re.compile(
    r"^\|\s*`(?P<file>[^`]+\.\w+)`\s*\|\s*(?P<prio>!?)\s*P(?P<phase>\d)\+?\s*\|"
    r"\s*(?P<cite>[^|]*?)\s*\|"
)


def is_ancillary(filename: str) -> bool:
    """Supplements and standalone figures are optional context; missing ones never block."""
    return "_supp" in filename or "_fig" in filename


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

    present = {p.name for p in PAPERS_DIR.iterdir() if p.name != "README.md"}
    missing = []

    for r in rows:
        ok = r["file"] in present
        mark = "ok     " if ok else "MISSING"
        flag = "!" if r["prio"] and not is_ancillary(r["file"]) else " "
        print(f"{mark} {flag} P{r['phase']}  {r['file']:<56} {r['cite']}")
        if not ok:
            missing.append(r)

    unlisted = sorted(present - {r["file"] for r in expected()})
    for name in unlisted:
        print(f"UNLISTED  --  {name:<56} not in the manifest — rename it or add a row")

    blocking = [r for r in missing if r["prio"] and not is_ancillary(r["file"])]
    n_anc = sum(1 for r in rows if is_ancillary(r["file"]))
    print(
        f"\n{len(rows) - len(missing)}/{len(rows)} present "
        f"({len(rows) - n_anc} main, {n_anc} supplementary), "
        f"{len(unlisted)} unlisted, {len(blocking)} blocking"
    )
    if blocking:
        print("blocking the next phase: " + ", ".join(r["file"] for r in blocking))

    waiting = sorted(f.name for f in INBOX_DIR.glob("*") if f.name != "README.md")
    if waiting:
        print(f"\n{len(waiting)} file(s) waiting in {INBOX_DIR.name}/ to be identified and filed:")
        for name in waiting:
            print(f"  {name}")

    print("\nknown gaps are documented at the end of docs/papers/README.md")


if __name__ == "__main__":
    main()
