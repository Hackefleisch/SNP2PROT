#!/usr/bin/env python
"""Assign cluster ids across the whole corpus, and drop redundant short-window entries.

    python scripts/build_clusters.py [--max-edits 5] [--dry-run]

Clustering is corpus-wide by nature, so it cannot live in a parser -- rule 7 keeps a parser to
one source. This runs after `build_dataset.py`, reads every interim table, and rewrites
`wt_id` from the cluster the domain belongs to.

Two things happen here, in this order.

**Subsumption.** Canonicalisation removes length variation wherever a construct covered the
window or a reference could be resolved, but an engineered chimera has neither and keeps a
short window. If some *other* source assayed the same domain and did produce a full window,
the short copy describes the same protein at a different length -- exactly the pollution the
canonical form exists to prevent. Domains that are zero edits apart under free terminal gaps
are the same domain, so all but the longest are dropped.

**Clustering.** CD-HIT greedy incremental over what remains (see `snp2prot.clusters`). Clusters
are for building splits only: every member stays its own training row, and `n_mut_from_wt` is
recomputed against the cluster representative so one coordinate frame covers the cluster.
"""

from __future__ import annotations

import argparse
import glob

import pandas as pd

from snp2prot import align, clusters
from snp2prot.config import INTERIM_DIR, REPORTS_DIR

#: Marks a wt_id as cluster-derived rather than lineage-derived, and makes the rewrite
#: idempotent -- see `load_domains`.
CLUSTER_PREFIX = "C:"


def load_domains() -> pd.DataFrame:
    cols = ["dbd_seq", "dbd_family", "wt_id", "source_dataset"]
    frames = [
        pd.read_parquet(f, columns=cols).drop_duplicates()
        for f in sorted(glob.glob(str(INTERIM_DIR / "*" / "*.parquet")))
        if "/proteins/" not in f
    ]
    df = pd.concat(frames, ignore_index=True).drop_duplicates("dbd_seq").reset_index(drop=True)
    # This script rewrites wt_id in place, so a re-run -- or a resumed one after an
    # interruption -- can meet its own output. Strip the marker so the lineage name underneath
    # is recovered rather than nested into "C:C:BAR15A:ARX".
    df["wt_id"] = df.wt_id.astype("string").str.removeprefix(CLUSTER_PREFIX)
    return df


def subsumed(domains: pd.DataFrame) -> dict[str, str]:
    """Short-window domain -> the longer domain that supersedes it.

    Zero edits apart with free terminal gaps means one is the other with different padding.
    """
    seqs = {r.dbd_seq: r.dbd_seq for r in domains.itertuples()}
    fams = {r.dbd_seq: r.dbd_family for r in domains.itertuples()}
    out: dict[str, str] = {}
    for c in clusters.cluster(seqs, max_edits=0, families=fams):
        for m in c.members:
            if m != c.representative and len(m) != len(c.representative):
                out[m] = c.representative
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--max-edits", type=int, default=5)
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    domains = load_domains()
    print(f"{len(domains)} distinct canonical domains")

    drop = subsumed(domains)
    print(f"redundant short-window domains superseded by a fuller copy: {len(drop)}")
    kept = domains[~domains.dbd_seq.isin(drop)].reset_index(drop=True)

    seqs = {r.dbd_seq: r.dbd_seq for r in kept.itertuples()}
    fams = {r.dbd_seq: r.dbd_family for r in kept.itertuples()}
    cs = clusters.cluster(seqs, max_edits=args.max_edits, families=fams)
    rep = clusters.assign(cs)
    multi = sum(1 for c in cs if c.size > 1)
    print(f"{len(cs)} clusters at <= {args.max_edits} edits; {multi} hold more than one domain")
    print(f"  largest cluster: {max(c.size for c in cs)} domains")

    # A readable id: the representative's own lineage wt_id, so a cluster stays recognisable.
    # But that name is NOT unique -- `ROG18A:FoxJ3` is the lineage name of six different
    # representatives, because one gene folder held six engineered chimeras that now sit in
    # six clusters. Left alone they collapse back into one id and the clustering is undone,
    # for precisely the constructs it matters most for. Collisions get a suffix.
    lineage = dict(zip(kept.dbd_seq, kept.wt_id, strict=True))
    used: dict[str, int] = {}
    rep_name: dict[str, str] = {}
    for c in sorted(cs, key=lambda c: c.representative):
        base = lineage[c.representative]
        used[base] = used.get(base, 0) + 1
        rep_name[c.representative] = base if used[base] == 1 else f"{base}#{used[base]}"
    cluster_id = {seq: f"{CLUSTER_PREFIX}{rep_name[r]}" for seq, r in rep.items()}
    if any(v > 1 for v in used.values()):
        n = sum(v - 1 for v in used.values() if v > 1)
        print(f"  {n} cluster ids disambiguated with a suffix (shared lineage name)")

    if args.dry_run:
        print("\n--dry-run: nothing written")
        return

    # Per DISTINCT domain, not per row: the corpus has ~1,350 domains and 39 million rows,
    # and aligning once per row would be four orders of magnitude of wasted work.
    profile = {}
    for seq, r in rep.items():
        if r == seq:
            profile[seq] = (0, "")
        else:
            e = align.edit_profile(r, seq)
            profile[seq] = (e.n_edits, e.positions_str)
    print(f"edit profiles computed for {len(profile)} distinct domains")

    rows_before = rows_after = 0
    for path in sorted(glob.glob(str(INTERIM_DIR / "*" / "*.parquet"))):
        if "/proteins/" in path:
            continue
        df = pd.read_parquet(path)
        rows_before += len(df)
        df = df[~df.dbd_seq.isin(drop)]
        if df.empty:
            print(f"  {path}: every domain superseded; left unchanged")
            continue
        df["wt_id"] = df.dbd_seq.map(cluster_id).astype("string")
        df["n_mut_from_wt"] = df.dbd_seq.map(lambda s: profile[s][0]).astype("Int32")
        df["mut_positions"] = df.dbd_seq.map(lambda s: profile[s][1]).astype("string")
        rows_after += len(df)
        df.to_parquet(path, index=False)
    print(f"\nrows {rows_before:,} -> {rows_after:,}")

    inv = REPORTS_DIR / "clusters.md"
    sizes = pd.Series([c.size for c in cs]).value_counts().sort_index()
    lines = [
        "# Cluster inventory — whole corpus",
        "",
        "Regenerated by `scripts/build_clusters.py`. Do not hand-edit.",
        "",
        f"- algorithm: CD-HIT greedy incremental, `<= {args.max_edits}` edits to the rep",
        f"- domains: {len(kept):,} (after dropping {len(drop)} superseded short-window copies)",
        f"- clusters: {len(cs):,}, of which {multi} hold more than one domain",
        "",
        "| domains per cluster | clusters |",
        "|---:|---:|",
        *[f"| {k} | {v} |" for k, v in sizes.items()],
    ]
    inv.write_text("\n".join(lines) + "\n")
    print(f"wrote {inv}")


if __name__ == "__main__":
    main()
