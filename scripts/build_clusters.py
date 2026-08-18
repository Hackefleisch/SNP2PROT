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
recomputed against the cluster's reference -- its medoid -- so one coordinate frame covers
the cluster and a wild type is not described as a mutant of one of its own variants.
"""

from __future__ import annotations

import argparse

import pandas as pd

from snp2prot import align, clusters, thresholds
from snp2prot.config import CLUSTER_TABLE, REPORTS_DIR, source_tables

#: Marks a wt_id as cluster-derived rather than lineage-derived, and makes the rewrite
#: idempotent -- see `load_domains`.
CLUSTER_PREFIX = "C:"


def load_domains() -> pd.DataFrame:
    cols = ["dbd_seq", "dbd_family", "wt_id", "source_dataset"]
    frames = [pd.read_parquet(f, columns=cols).drop_duplicates() for f in source_tables().values()]
    df = pd.concat(frames, ignore_index=True).drop_duplicates("dbd_seq").reset_index(drop=True)
    # This script rewrites wt_id in place, so a re-run -- or a resumed one after an
    # interruption -- can meet its own output. Strip the marker so the lineage name underneath
    # is recovered rather than nested into "C:C:BAR15A:ARX".
    df["wt_id"] = df.wt_id.astype("string").str.removeprefix(CLUSTER_PREFIX)
    return df


def subsumed(domains: pd.DataFrame) -> dict[str, str]:
    """Short-window domain -> the longer domain that supersedes it.

    Zero edits apart with free terminal gaps means one is the other with different padding —
    but only if the shorter one is *wholly* contained, which is what `min_overlap=1.0`
    demands. Without it, two unrelated domains whose six aligned residues happen to match
    would be called one domain and one of them silently dropped from the corpus.
    """
    seqs = {r.dbd_seq: r.dbd_seq for r in domains.itertuples()}
    fams = {r.dbd_seq: r.dbd_family for r in domains.itertuples()}
    out: dict[str, str] = {}
    for c in clusters.cluster(seqs, max_edits=0, families=fams, min_overlap=1.0):
        for m in c.members:
            if m != c.seed and len(m) != len(c.seed):
                out[m] = c.seed
    return out


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--max-edits", type=int, default=5)
    ap.add_argument(
        "--min-overlap",
        type=float,
        default=None,
        help="fraction of the shorter domain that must align; default from thresholds.yaml",
    )
    ap.add_argument("--dry-run", action="store_true")
    args = ap.parse_args()

    cfg = thresholds.load()["cluster"]
    min_overlap = float(cfg["min_overlap"]) if args.min_overlap is None else args.min_overlap

    domains = load_domains()
    print(f"{len(domains)} distinct canonical domains")

    drop = subsumed(domains)
    print(f"redundant short-window domains superseded by a fuller copy: {len(drop)}")
    kept = domains[~domains.dbd_seq.isin(drop)].reset_index(drop=True)

    seqs = {r.dbd_seq: r.dbd_seq for r in kept.itertuples()}
    fams = {r.dbd_seq: r.dbd_family for r in kept.itertuples()}
    cs = clusters.cluster(seqs, max_edits=args.max_edits, families=fams, min_overlap=min_overlap)
    rep = clusters.assign(cs)
    # The seed decided membership; the reference is the medoid and decides the coordinate
    # frame every member's mut_positions are expressed in (`docs/DECISIONS.md`, `D3`).
    reference = clusters.references(cs)
    reference_of = {c.seed: reference[c.seed] for c in cs}
    reseated = sum(1 for c in cs if reference_of[c.seed] != c.seed)
    multi = sum(1 for c in cs if c.size > 1)
    print(
        f"{len(cs)} clusters at <= {args.max_edits} edits and >= {min_overlap:.0%} overlap; "
        f"{multi} hold more than one domain"
    )
    print(f"  largest cluster: {max(c.size for c in cs)} domains")
    print(f"  {reseated} clusters take a reference other than their seed")

    # A readable id: the reference's own lineage wt_id, so a cluster stays recognisable.
    # But that name is NOT unique -- `ROG18A:FoxJ3` is the lineage name of six different
    # references, because one gene folder held six engineered chimeras that now sit in
    # six clusters. Left alone they collapse back into one id and the clustering is undone,
    # for precisely the constructs it matters most for. Collisions get a suffix.
    lineage = dict(zip(kept.dbd_seq, kept.wt_id, strict=True))
    used: dict[str, int] = {}
    rep_name: dict[str, str] = {}
    for c in sorted(cs, key=lambda c: reference_of[c.seed]):
        base = lineage[reference_of[c.seed]]
        used[base] = used.get(base, 0) + 1
        rep_name[reference_of[c.seed]] = base if used[base] == 1 else f"{base}#{used[base]}"
    cluster_id = {seq: f"{CLUSTER_PREFIX}{rep_name[reference_of[r]]}" for seq, r in rep.items()}
    if any(v > 1 for v in used.values()):
        n = sum(v - 1 for v in used.values() if v > 1)
        print(f"  {n} cluster ids disambiguated with a suffix (shared lineage name)")

    if args.dry_run:
        print("\n--dry-run: nothing written")
        return

    # Per DISTINCT domain, not per row: the corpus has ~1,350 domains and 39 million rows,
    # and aligning once per row would be four orders of magnitude of wasted work.
    profile = {}
    over = []
    for seq, r in rep.items():
        ref = reference_of[r]
        if ref == seq:
            profile[seq] = (0, "")
            continue
        e = align.edit_profile(ref, seq)
        profile[seq] = (e.n_edits, e.positions_str)
        # Members are within max_edits of the SEED by construction; of the reference only
        # within 2 * max_edits, by the triangle inequality. Checked, never assumed.
        if e.n_edits > args.max_edits:
            over.append((rep_name[ref], e.n_edits))
    print(f"edit profiles computed for {len(profile)} distinct domains")
    if over:
        print(f"  WARNING: {len(over)} domains sit more than {args.max_edits} edits from their")
        print("           cluster reference -- the medoid is farther from an outlier than the")
        print(
            "           seed was. Worst: "
            + ", ".join(f"{n} at {d}" for n, d in sorted(over, key=lambda x: -x[1])[:5])
        )

    rows_before = rows_after = 0
    # Accumulated while the rows are already in memory: the inventory needs a row count per
    # cluster, and a second pass over 45 million rows to get one would cost more than
    # everything else in this script put together.
    constructs = []
    for source, path in source_tables().items():
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

        counts = df.groupby(["dbd_seq", "protein_id"], dropna=False).size()
        per_construct = (
            df[["dbd_seq", "protein_id", "dbd_family", "wt_id", "n_mut_from_wt"]]
            .drop_duplicates(["dbd_seq", "protein_id"])
            .set_index(["dbd_seq", "protein_id"])
        )
        per_construct["n_rows"] = counts
        per_construct["source_dataset"] = source
        # Carried so the inventory can record what membership was decided against, which is
        # not the same domain as the coordinate frame (`docs/DECISIONS.md`, `D3`).
        per_construct["seed"] = per_construct.index.get_level_values("dbd_seq").map(rep)
        constructs.append(per_construct.reset_index())
    print(f"\nrows {rows_before:,} -> {rows_after:,}")

    inv_table = clusters.inventory(pd.concat(constructs, ignore_index=True))
    CLUSTER_TABLE.parent.mkdir(parents=True, exist_ok=True)
    inv_table.to_parquet(CLUSTER_TABLE, index=False)
    print(f"wrote {CLUSTER_TABLE} ({len(inv_table):,} clusters)")

    inv = REPORTS_DIR / "clusters.md"
    sizes = pd.Series([c.size for c in cs]).value_counts().sort_index()
    with_variants = inv_table[inv_table.n_variants > 0]
    by_family = (
        inv_table[inv_table.n_domains > 1].groupby("dbd_family").size().sort_values(ascending=False)
    )
    lines = [
        "# Cluster inventory — whole corpus",
        "",
        "Regenerated by `scripts/build_clusters.py`. Do not hand-edit.",
        "",
        f"- algorithm: CD-HIT greedy incremental, `<= {args.max_edits}` edits to the rep",
        f"- a member must also align over `>= {min_overlap:.0%}` of the shorter sequence, "
        "which is what stops free terminal gaps from calling two unrelated domains "
        "near-identical (`docs/DECISIONS.md`, `T25`)",
        f"- domains: {len(kept):,}",
        f"- clusters: {len(cs):,}, of which {multi} hold more than one domain",
        f"- variant domains: {int(inv_table.n_variants.sum()):,} across "
        f"{len(with_variants):,} clusters",
        "",
        # Deliberately not "after dropping N superseded copies": N describes the state of
        # `data/interim/` when the script ran (a fresh build has copies to drop, a re-run has
        # none), not the dataset. A committed report whose diff depends on how many times a
        # step was re-run is a report that cannot be trusted as evidence.
        "Short-window duplicates are removed before clustering, so these counts are the same"
        " whether the corpus was just built or the pass was repeated. Per-cluster figures are"
        " in `data/interim/clusters/clusters.parquet`; read them with `snp2prot.clusters.load`"
        " rather than by grouping the row tables.",
        "",
        "| domains per cluster | clusters |",
        "|---:|---:|",
        *[f"| {k} | {v} |" for k, v in sizes.items()],
        "",
        "## Multi-domain clusters by family",
        "",
        "| family | clusters |",
        "|---|---:|",
        *[f"| `{k}` | {v} |" for k, v in by_family.items()],
        "",
        "## The largest clusters",
        "",
        "| `wt_id` | family | domains | variants | max edits | sources |",
        "|---|---|---:|---:|---:|---:|",
        *[
            f"| `{r.wt_id}` | `{r.dbd_family}` | {r.n_domains} | {r.n_variants} |"
            f" {r.max_edits} | {r.n_sources} |"
            for r in inv_table.head(15).itertuples()
        ],
        "",
    ]
    inv.write_text("\n".join(lines))
    print(f"wrote {inv}")


if __name__ == "__main__":
    main()
