#!/usr/bin/env python
"""Merge every parsed source into the single training table.

    python scripts/build_merged.py [--out data/processed/training.parquet] [--dry-run]

Phase 5. The parsers produce one schema-conforming table per source and are forbidden to know
about each other (rule 7); everything that requires seeing two sources at once happens here,
and it is exactly one thing:

**One record per domain survives.** 47 domains are measured by more than one source and their
labels do not agree — the median pair on 46% of the 8-mers either called positive. Keeping
both would hand a model identical input with two different labels, inside a single cluster
where no split can separate them. The record with more positives wins, unless one candidate
belongs to a variant series, in which case the series wins, because a wild type has to be
measured by the same lab and array design as the variants it is the reference for
(`snp2prot.merge`, `docs/DECISIONS.md` 2026-08-18).

**What this step deliberately does not do:**

- *drop records with no positive evidence.* 34 records have no positive 8-mer and no control
  to interpret that against, and 20 more are variants that measurably lost binding. Which of
  those a model sees is a training decision, applied at featurization through
  `label_health.usable`, not baked into the dataset (`T21`).
- *touch a label, a threshold or a score.* The merged table is a concatenation minus the
  losing duplicates. Every cutoff was applied per experiment at parse time and is recorded in
  each row.
- *drop the no-call band.* `label == -1` rows stay: they are the 8-mers between the two
  cutoffs, and excluding them is likewise a modelling choice.

The result is written once and read many times, so it is a single Parquet file rather than a
directory: the split code reads it whole, and per-source provenance survives in
`source_dataset`.
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq

from snp2prot import label_health, merge, schema
from snp2prot.config import PROCESSED_DIR, REPORTS_DIR, source_tables

DEFAULT_OUT = PROCESSED_DIR / "training.parquet"
BATCH_ROWS = 2_000_000


def dropped_records(records: pd.DataFrame) -> set[tuple[str, str]]:
    """The `(dbd_seq, source_dataset)` pairs that lose their duplicate contest."""
    resolution = merge.resolve(records)
    losers = resolution[~resolution.keep]
    return set(zip(losers.dbd_seq, losers.source_dataset, strict=True))


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument("--dry-run", action="store_true", help="report the plan, write nothing")
    args = ap.parse_args()

    tables = source_tables()
    if not tables:
        raise SystemExit("nothing parsed yet\nrun scripts/build_dataset.py --all")
    records = label_health.load()
    drop = dropped_records(records)
    print(f"{len(records)} records across {len(tables)} sources; {len(drop)} dropped as duplicates")

    started = time.time()
    writer = None
    kept_rows = seen_rows = 0
    per_source: dict[str, int] = {}
    try:
        for source, path in sorted(tables.items()):
            for batch in pq.ParquetFile(path).iter_batches(batch_size=BATCH_ROWS):
                df = batch.to_pandas()
                seen_rows += len(df)
                if drop:
                    pairs = zip(df["dbd_seq"], df["source_dataset"], strict=True)
                    df = df[[(seq, src) not in drop for seq, src in pairs]]
                if df.empty:
                    continue
                df = schema.coerce(df)
                kept_rows += len(df)
                per_source[source] = per_source.get(source, 0) + len(df)
                if args.dry_run:
                    continue
                table = pa.Table.from_pandas(df, preserve_index=False)
                if writer is None:
                    args.out.parent.mkdir(parents=True, exist_ok=True)
                    writer = pq.ParquetWriter(args.out, table.schema, compression="zstd")
                writer.write_table(table)
            print(f"  {source:<16} {per_source.get(source, 0):>10,} rows", flush=True)
    finally:
        if writer is not None:
            writer.close()

    print(f"\n{seen_rows:,} rows in -> {kept_rows:,} rows out ({seen_rows - kept_rows:,} dropped)")
    print(f"took {time.time() - started:.0f}s")
    if args.dry_run:
        print("\n--dry-run: nothing written")
        return
    print(f"wrote {args.out} ({args.out.stat().st_size / 1e9:.2f} GB)")

    report = REPORTS_DIR / "merge.md"
    report.write_text(render(per_source, records, drop, kept_rows, args.out))
    print(f"wrote {report}")


def render(per_source, records, drop, kept_rows, out: Path) -> str:
    kept = records[
        [(s, d) not in drop for s, d in zip(records.dbd_seq, records.source_dataset, strict=True)]
    ]
    lines = [
        "# The merged training table",
        "",
        "Regenerated by `scripts/build_merged.py`. Do not hand-edit.",
        "",
        f"`{out.relative_to(out.parents[2])}` — **{kept_rows:,} rows**, one per "
        "(domain, 8-mer), from the sources below. Every cutoff was applied per experiment at "
        "parse time and travels in the row; nothing here re-thresholds anything.",
        "",
        "## What the merge step did",
        "",
        f"- **{len(drop)} of {len(records)} domain records dropped** as duplicate measurements "
        "of one domain: more positives wins, unless one candidate belongs to a variant series "
        "(`docs/DECISIONS.md`, 2026-08-18). Every resolution is listed in `overlap.md`.",
        "- Nothing else. No label, threshold or score was touched, no-call rows were kept, and "
        "records with no positive evidence were kept and left flagged in `label_health.md` for "
        "the training filter to decide.",
        "",
        "## Rows by source",
        "",
        "| source | domains kept | rows |",
        "|---|---:|---:|",
    ]
    for source, rows in sorted(per_source.items(), key=lambda kv: -kv[1]):
        n = int((kept.source_dataset == source).sum())
        lines.append(f"| `{source}` | {n:,} | {rows:,} |")
    lines += [
        f"| **total** | **{len(kept):,}** | **{kept_rows:,}** |",
        "",
        "## Label composition",
        "",
        "| label | rows | share |",
        "|---|---:|---:|",
    ]
    composition = (
        ("binding (1)", "n_pos"),
        ("non-binding (0)", "n_neg"),
        ("no call (-1)", "n_gray"),
    )
    for name, col in composition:
        n = int(kept[col].sum())
        lines.append(f"| {name} | {n:,} | {n / kept_rows:.3%} |")
    calls = int(kept.n_pos.sum() + kept.n_neg.sum())
    lines += [
        "",
        f"Negative-to-positive ratio **{int(kept.n_neg.sum() / kept.n_pos.sum())}:1** on calls "
        f"made ({int(kept.n_pos.sum()):,} of {calls:,}). That imbalance is the assay, not a "
        "defect: a universal PBM scores every protein against every 8-mer and a transcription "
        "factor binds a small minority of them.",
        "",
    ]
    return "\n".join(lines)


if __name__ == "__main__":
    main()
