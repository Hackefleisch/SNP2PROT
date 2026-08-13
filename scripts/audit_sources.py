#!/usr/bin/env python
"""Cross-check every UniPROBE source against the invariants past bugs violated.

    python scripts/audit_sources.py            # all registered accessions
    python scripts/audit_sources.py --source GR09

Every defect found in this project so far was noticed because a distribution or a structural
assumption looked wrong, never because something raised. This script turns each of those
observations into a check that runs over all sources at once:

1. every archive file is claimed by exactly one construct
   -- HOXD13's `_8mers_11111111.txt`, silently dropped by a `_8mers.txt` glob
2. every construct on a page is emitted or rejected with a reason
   -- ROG18A's six chimeras collapsed into one protein
3. no two constructs share a detail-page key
   -- GR09's bare `Insert sequence` label, which gave 88 genes one empty key
4. the E-score column is identifiable in every file
   -- GR09/SCI09 reading median intensity; RAD13A has no E-score at all
5. scores look like E-scores: bounded, signed, centred near zero
   -- the 4:1 label ratio with raw_score reaching 776,106
6. every domain has exactly 32,896 rows -- a partially-read or duplicated experiment
7. the 8-mer set matches every other source -- a source on a different array design
8. no cluster is a wild outlier in size -- the 45-domain `GR09:` cluster

It reads archives, detail pages and the existing interim tables. It does **not** re-scan Pfam,
so it is cheap to run before a rebuild rather than after.
"""

from __future__ import annotations

import argparse
import zipfile
from dataclasses import dataclass, field

import pandas as pd

from snp2prot.config import interim_table, raw_dir
from snp2prot.parsers import _uniprobe, uniprobe_panels

EXPECTED_8MERS = 32896


@dataclass
class Finding:
    source: str
    check: str
    severity: str  # "ERROR" or "warn"
    detail: str


@dataclass
class SourceAudit:
    source: str
    findings: list[Finding] = field(default_factory=list)
    stats: dict = field(default_factory=dict)

    def error(self, check: str, detail: str) -> None:
        self.findings.append(Finding(self.source, check, "ERROR", detail))

    def warn(self, check: str, detail: str) -> None:
        self.findings.append(Finding(self.source, check, "warn", detail))


def audit_structure(acc: str, a: SourceAudit) -> None:
    """Checks 1-4: archive coverage, construct coverage, key collisions, readable E-scores."""
    panel = uniprobe_panels.PANELS.get(acc)
    archive = raw_dir(acc) / (panel.archive if panel else f"{acc}_contig8mers.zip")
    details_dir = raw_dir(acc) / "details"
    if not archive.exists() or not details_dir.exists():
        a.warn("inputs", f"missing archive or details ({archive.name})")
        return

    pages = _uniprobe.load_details(details_dir)
    experiments = uniprobe_panels.group_experiments(archive)
    all_members = {m for ms in experiments.values() for m in ms}
    a.stats["experiment_files"] = len(all_members)
    a.stats["detail_pages"] = len(pages)

    # 3: a page key collision silently fuses proteins.
    seen: dict[str, str] = {}
    for gene, page in pages.items():
        for key in page.inserts:
            if not key:
                a.error("page-keys", f"{gene}: empty construct key")
            elif key in seen and seen[key] != gene:
                a.error("page-keys", f"construct key {key!r} used by {seen[key]} and {gene}")
            seen[key] = gene
    a.stats["constructs"] = len(seen)

    # 1 + 2: coverage both ways, classified by cause. Only a name mismatch is a defect --
    # the other two are genuine gaps in what UniPROBE publishes, and the parser reports them.
    no_page, no_sequence, mismatch = [], [], []
    for gene, members in experiments.items():
        page = pages.get(gene)
        if page is None:
            no_page.append(gene)
        elif not page.inserts:
            no_sequence.append(gene)
        else:
            by_construct = uniprobe_panels.split_by_construct(members, list(page.inserts))
            if sum(len(v) for v in by_construct.values()) < len(members):
                mismatch.append(gene)
    if mismatch:
        a.error(
            "construct-names",
            f"{len(mismatch)} gene(s) whose experiments match no construct on their own "
            f"page, so the data is dropped silently: {mismatch[:3]}",
        )
    if no_page:
        a.warn(
            "no-sequence",
            f"{len(no_page)} archive gene(s) have no detail page — usually protein "
            f"complexes, which fail admission condition 1 regardless: {no_page[:3]}",
        )
    if no_sequence:
        a.warn(
            "no-sequence",
            f"{len(no_sequence)} gene(s) have a page carrying no sequence: {no_sequence[:3]}",
        )

    # 4: can an E-score actually be read from every file?
    unreadable_by_gene: dict[str, int] = {}
    total_by_gene: dict[str, int] = {}
    with zipfile.ZipFile(archive) as z:
        for gene, members in experiments.items():
            for m in members:
                total_by_gene[gene] = total_by_gene.get(gene, 0) + 1
                try:
                    _uniprobe.read_8mer_table(z, m)
                except _uniprobe.EscoreColumnError:
                    unreadable_by_gene[gene] = unreadable_by_gene.get(gene, 0) + 1
    # Unreadable files only cost data when a gene has no readable file left. SCI09 ships
    # 20-column combined files alongside per-replicate ones; rejecting the former loses nothing.
    starved = [g for g, n in unreadable_by_gene.items() if n == total_by_gene[g]]
    n_unreadable = sum(unreadable_by_gene.values())
    if starved:
        a.error(
            "escore-column",
            f"{len(starved)} gene(s) have no file with a readable E-score: {starved[:3]}",
        )
    elif n_unreadable:
        a.warn(
            "escore-column",
            f"{n_unreadable} file(s) unreadable but every gene keeps a readable replicate",
        )


def audit_table(acc: str, a: SourceAudit, reference_8mers: set[str] | None) -> set[str] | None:
    """Checks 5-8, from the parsed table."""
    path = interim_table(acc)
    if not path.exists():
        a.warn("interim", "not parsed")
        return reference_8mers

    df = pd.read_parquet(path, columns=["wt_id", "dbd_seq", "dna_seq", "label", "raw_score"])
    a.stats["rows"] = len(df)
    a.stats["domains"] = df["dbd_seq"].nunique()

    # 5: an E-score is bounded, signed and sits near zero for most 8-mers.
    s = df["raw_score"]
    if s.min() < -0.5 or s.max() > 0.5:
        a.error("score-range", f"raw_score outside [-0.5, 0.5]: [{s.min():.4g}, {s.max():.4g}]")
    if s.min() > -0.3:
        a.warn("score-shape", f"minimum only {s.min():.3f}; an E-score should reach near -0.5")
    if abs(s.median()) > 0.1:
        a.warn("score-shape", f"median {s.median():.3f}; expected near 0")
    pos_rate = float((df["label"] == 1).mean())
    a.stats["pos_rate"] = pos_rate
    if pos_rate > 0.05:
        a.error("positive-rate", f"{pos_rate:.1%} positive — implausible for a PBM at E>=0.45")

    # 6: the design is fully crossed, so every domain sees every 8-mer exactly once.
    per_domain = df.groupby("dbd_seq").size()
    odd = per_domain[per_domain != EXPECTED_8MERS]
    if len(odd):
        a.error(
            "rows-per-domain",
            f"{len(odd)} domain(s) not at {EXPECTED_8MERS:,} rows, e.g. {odd.head(2).to_dict()}",
        )

    # 7: one shared DNA axis across the whole corpus.
    kmers = set(df["dna_seq"].unique())
    if len(kmers) != EXPECTED_8MERS:
        a.error("8mer-set", f"{len(kmers):,} distinct 8-mers, expected {EXPECTED_8MERS:,}")
    if reference_8mers is not None and kmers != reference_8mers:
        a.error(
            "8mer-set",
            f"8-mer set differs from the first source by {len(kmers ^ reference_8mers)} k-mers",
        )

    # 8: one cluster far larger than the rest means a key collapsed.
    sizes = df.drop_duplicates("dbd_seq").groupby("wt_id").size()
    if len(sizes) > 2:
        biggest, median = int(sizes.max()), float(sizes.median())
        if biggest > max(12, 6 * median):
            a.error(
                "cluster-size",
                f"largest cluster has {biggest} domains vs median {median:g}: {sizes.idxmax()!r}",
            )
    return reference_8mers if reference_8mers is not None else kmers


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--source")
    args = ap.parse_args()

    sources = [args.source] if args.source else sorted(uniprobe_panels.PANELS) + ["BAR15A"]
    audits, reference = [], None
    for acc in sources:
        a = SourceAudit(acc)
        try:
            audit_structure(acc, a)
            reference = audit_table(acc, a, reference)
        except Exception as exc:  # noqa: BLE001 - an audit must never abort the sweep
            a.error("audit", f"{type(exc).__name__}: {exc}")
        audits.append(a)
        st = a.stats
        flag = "".join(sorted({f.severity[0].upper() for f in a.findings}))
        print(
            f"  {acc:8s} files={st.get('experiment_files', '-'):>4} "
            f"constructs={st.get('constructs', '-'):>4} "
            f"domains={st.get('domains', '-'):>4} "
            f"rows={st.get('rows', 0):>10,} {flag}"
        )

    print("\n" + "=" * 78)
    errors = [f for a in audits for f in a.findings if f.severity == "ERROR"]
    warns = [f for a in audits for f in a.findings if f.severity == "warn"]
    for f in errors:
        print(f"ERROR  {f.source:8s} [{f.check}] {f.detail}")
    for f in warns:
        print(f"warn   {f.source:8s} [{f.check}] {f.detail}")
    print(f"\n{len(errors)} error(s), {len(warns)} warning(s) across {len(sources)} source(s)")
    raise SystemExit(1 if errors else 0)


if __name__ == "__main__":
    main()
