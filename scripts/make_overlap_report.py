#!/usr/bin/env python
"""Regenerate `reports/overlap.md` — the cross-source label agreement check.

    python scripts/make_overlap_report.py

Corpus-wide rather than per-source, so it gets its own entry point alongside
`build_clusters.py` rather than living in `make_reports.py --source`.

Two passes, because the corpus is 45M rows and only ~3M of them are replicated: pass one
reads the `dbd_seq` column of every interim table to find the domains more than one source
measured, pass two reads back only those rows. Holding all 19 frames at once would cost tens
of gigabytes for a report about 3% of them.
"""

from __future__ import annotations

import glob
from collections import Counter
from pathlib import Path

import pandas as pd
import pyarrow.compute as pc
import pyarrow.dataset as ds

from snp2prot import reports
from snp2prot.config import INTERIM_DIR, PROJECT_ROOT, REPORTS_DIR

#: Columns the report needs. Reading `dbd_seq` for 45M rows is already the expensive part;
#: pulling the other 17 columns as well would double it for nothing.
COLUMNS = ["source_dataset", "dbd_seq", "dbd_family", "wt_id", "dna_seq", "label", "raw_score"]


def main() -> None:
    # Discovered from disk rather than from the parser REGISTRY, which also lists sources
    # screened out and never parsed (GD09, PP15, ...). The protein-side table lives under
    # the same tree and is not a source.
    paths = sorted(glob.glob(str(INTERIM_DIR / "*" / "*.parquet")))
    tables = {Path(p).parent.name: Path(p) for p in paths if Path(p).parent.name != "proteins"}
    if not tables:
        raise SystemExit("nothing parsed yet\nrun scripts/build_dataset.py --all")

    counts: Counter[str] = Counter()
    for source, path in tables.items():
        seqs = ds.dataset(path).to_table(columns=["dbd_seq"]).column("dbd_seq").unique()
        counts.update(seqs.to_pylist())
        print(f"  {source}: {len(seqs):,} domains")

    shared = [seq for seq, n in counts.items() if n > 1]
    print(f"{len(shared):,} domains measured by more than one source")
    if not shared:
        print("nothing to compare")

    frames: dict[str, pd.DataFrame] = {}
    for source, path in tables.items():
        table = ds.dataset(path).to_table(columns=COLUMNS, filter=pc.field("dbd_seq").isin(shared))
        if table.num_rows:
            frames[source] = table.to_pandas()

    REPORTS_DIR.mkdir(parents=True, exist_ok=True)
    out = REPORTS_DIR / "overlap.md"
    out.write_text(reports.overlap_report(frames))
    print(f"wrote {out.relative_to(PROJECT_ROOT)}")


if __name__ == "__main__":
    main()
