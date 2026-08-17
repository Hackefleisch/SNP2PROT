#!/usr/bin/env python
"""Run one or more source parsers, validate, and checkpoint to data/interim/.

    python scripts/build_dataset.py --source BAR15A
    python scripts/build_dataset.py --all

Nothing is written unless the schema validator passes. The merge step (Phase 5) is a
separate entry point — this script never combines sources.
"""

from __future__ import annotations

import argparse

from snp2prot import schema
from snp2prot.config import interim_dir, interim_table
from snp2prot.parsers import REGISTRY


def build(source: str) -> None:
    if source not in REGISTRY:
        raise SystemExit(f"no parser registered for {source!r}; have {sorted(REGISTRY) or 'none'}")

    df = schema.coerce(REGISTRY[source]())
    report = schema.validate(df, source=source)
    print(report)
    report.raise_if_failed()

    interim_dir(source, create=True)
    out = interim_table(source)
    df.to_parquet(out, index=False)
    print(f"wrote {len(df):,} rows -> {out}")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    g = ap.add_mutually_exclusive_group(required=True)
    g.add_argument("--source", help="source_dataset key, e.g. BAR15A")
    g.add_argument("--all", action="store_true", help="every registered parser")
    args = ap.parse_args()

    sources = sorted(REGISTRY) if args.all else [args.source]
    failed: list[tuple[str, str]] = []
    for src in sources:
        # With --all, one source that yields nothing must not abort the rest: an accession
        # whose every construct is rejected by the domain policy is a normal outcome.
        try:
            build(src)
        except Exception as exc:  # noqa: BLE001 - reported, not swallowed
            failed.append((src, f"{type(exc).__name__}: {exc}"))
            print(f"SKIPPED {src}: {exc}")
            # A source that stops yielding must not leave last build's output behind. Two
            # accessions went to zero when canonicalisation landed, and their stale Parquet
            # was still being read as part of the corpus.
            stale = interim_table(src)
            if stale.exists():
                stale.unlink()
                print(f"  removed stale {stale}")
        if not args.all:
            raise SystemExit(0)

    print(f"\n{len(sources) - len(failed)}/{len(sources)} sources built")
    for src, why in failed:
        print(f"  no output: {src}  ({why})")


if __name__ == "__main__":
    main()
