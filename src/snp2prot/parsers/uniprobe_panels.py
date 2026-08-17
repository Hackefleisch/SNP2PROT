"""UniPROBE family panels: one protein per experiment, no engineered variants.

These accessions contribute breadth across protein sequence space rather than depth within a
cluster — the complement to BAR15A, which is all depth and little breadth. Every protein is
its own reference, so each becomes a `wt_id` cluster of size one with `n_mut_from_wt = 0`.
That is a real limitation, not a bookkeeping artefact, and the cluster inventory reports it.

| accession | citation                        | contributes            |
|-----------|---------------------------------|------------------------|
| `Cell08`  | Berger et al., Cell 2008        | ~168 mouse homeodomains|
| `EMBO10`  | Wei et al., EMBO J 2010         | ETS family             |
| `PNAS13`  | Nakagawa et al., PNAS 2013      | forkhead family        |
| `SCI09`   | Badis et al., Science 2009      | 104 mouse TFs, broad   |
| `GR09`    | Zhu et al., Genome Res 2009     | yeast TFs              |
| `MAR17A`  | Mariani et al., Cell Syst 2017  | bZIP and forkhead      |
| `SHO18A`  | Shokri et al., Cell Reports 2019| mixed                  |
| `ROG18A`  | Rogers et al., Mol Cell 2019    | forkhead, replicated   |

Every UniPROBE accession that publishes a `<ACC>_contig8mers.zip` is registered here, so the
corpus spans one uniform DNA space: all 32,896 non-redundant 8-mers, identical across every
source. Each construct is screened by the domain policy, so a panel contributes only those
whose signal is attributable to a single continuous domain — yields are typically well under
half the gene count, and some accessions contribute nothing at all.

`dbd_seq` is the detail page's Pfam-trimmed **DNA binding domain** field where present. That
differs from BAR15A, which must use the clone insert because its variant sequences exist only
in that form. `dbd_source` records which was used; the two are not interchangeable and are
reconciled at merge time (open item #11).
"""

from __future__ import annotations

import zipfile
from dataclasses import dataclass
from pathlib import Path

import pandas as pd

# `references` is aliased: the name is already taken inside parse_panel by the
# gene -> reference-sequence map used for clustering variants.
from snp2prot import align, canonical, domains, schema, thresholds
from snp2prot import references as reference_db
from snp2prot.config import raw_dir
from snp2prot.parsers import _uniprobe

DBD_SOURCE = "pfam_hmmer_padded"


@dataclass(frozen=True)
class Panel:
    accession: str
    citation: str
    archive: str

    @property
    def source_file(self) -> str:
        return f"{self.accession}/{self.archive}"


PANELS: dict[str, Panel] = {
    p.accession: p
    # Not registered: MIZ19A and MMB08 advertise a contig8mers.zip the server does not
    # have; GB11, LAI20A and NAR10 publish detail pages with no sequence at all, so the
    # assayed construct is unknown and nothing can be attributed to a domain. See
    # PROVENANCE.md for all five. Path10 publishes only the enriched end of each table
    # (341-1,391 rows cut at E >= 0.25, not the full 32,896), so it carries no negatives
    # and its top hits would be mislabelled as non-binding. GD09 is excluded twice over:
    # its tables are one 8-mer short of the full set, and its only protein has no Pfam
    # hit above the gathering threshold anywhere in its 294-residue construct.
    for p in (
        Panel("Cell08", "Berger et al., Cell 2008", "Cell08_contig8mers.zip"),
        Panel("EMBO10", "Wei et al., EMBO J 2010", "EMBO10_contig8mers.zip"),
        Panel("PNAS13", "Nakagawa et al., PNAS 2013", "PNAS13_contig8mers.zip"),
        Panel("SCI09", "Badis et al., Science 2009", "SCI09_contig8mers.zip"),
        Panel("GR09", "Zhu et al., Genome Res 2009", "GR09_contig8mers.zip"),
        Panel("MAR17A", "Mariani et al., Cell Systems 2017", "MAR17A_contig8mers.zip"),
        Panel("SHO18A", "Shokri et al., Cell Reports 2019", "SHO18A_contig8mers.zip"),
        Panel("ROG18A", "Rogers et al., Mol Cell 2019", "ROG18A_contig8mers.zip"),
        Panel("CB11", "Helfer et al., Curr Biol 2011", "CB11_contig8mers.zip"),
        Panel("CR09", "Scharer et al., Cancer Res 2009", "CR09_contig8mers.zip"),
        Panel("Cell09", "Grove et al., Cell 2009", "Cell09_contig8mers.zip"),
        Panel("DEV12", "Busser et al., Development 2012", "DEV12_contig8mers.zip"),
        Panel("GD09", "Lesch et al., Genes Dev 2009", "GD09_contig8mers.zip"),
        Panel("GD12", "Peterson et al., Genes Dev 2012", "GD12_contig8mers.zip"),
        Panel("GD13", "Soruco et al., Genes Dev 2013", "GD13_contig8mers.zip"),
        Panel("KUR17A", "Li et al., Nature 2017", "KUR17A_contig8mers.zip"),
        Panel("LIN14B", "Lindemose et al, Nucleic Acids Res. 2014", "LIN14B_contig8mers.zip"),
        Panel("LIU18A", "Liu et al., Cell 2018", "LIU18A_contig8mers.zip"),
        Panel("LIU18B", "Liu et al., eLife 2018", "LIU18B_contig8mers.zip"),
        Panel("MBE14", "Cheatle Jarvela et al., Mol Biol Evol 2014", "MBE14_contig8mers.zip"),
        Panel("NAR11", "De Masi et al., NAR 2011", "NAR11_contig8mers.zip"),
        Panel("NBT06", "Berger et al., Nat Biotech 2006", "NBT06_contig8mers.zip"),
        Panel("PNAS08", "De Silva et al., PNAS 2008", "PNAS08_contig8mers.zip"),
        Panel("PNAS12", "Busser et al., PNAS 2012", "PNAS12_contig8mers.zip"),
        Panel("PO10", "Del Bianco et al., PLoS ONE 2010", "PO10_contig8mers.zip"),
        Panel("PP15", "Lehti-Shiu et al., PP 2015", "PP15_contig8mers.zip"),
        Panel("Path10", "Campbell et al., PLoS Pathog 2010", "Path10_contig8mers.zip"),
        Panel("RAD13A", "Radke et al., PNAS 2013", "RAD13A_contig8mers.zip"),
        Panel("STI21B", "Stielow et al., Science Advances 2021", "STI21B_contig8mers.zip"),
    )
}


@dataclass
class SkippedProtein:
    gene: str
    reason: str


def group_experiments(archive: Path) -> dict[str, list[str]]:
    """Gene -> its experiment members. The gene is the first path element in every panel."""
    with zipfile.ZipFile(archive) as z:
        members = [n for n in z.namelist() if _uniprobe.CONTIG_8MER_RE.search(n)]
    out: dict[str, list[str]] = {}
    for m in members:
        if "/" in m:
            gene = m.split("/")[0]
        else:
            # PNAS08's archive is flat -- the file sits at the top level with no gene folder,
            # so the gene has to come from the filename. Its detail pages are named for the
            # stem exactly ("Cgd2_3490" for "Cgd2_3490_contig8mers.txt").
            gene = _uniprobe.CONTIG_8MER_RE.sub("", m)
        out.setdefault(gene, []).append(m)
    return {k: sorted(v) for k, v in sorted(out.items())}


def split_by_construct(members: list[str], constructs: list[str]) -> dict[str, list[str]]:
    """Split one gene's experiments across the distinct constructs on its detail page.

    A gene folder does not always hold one protein. ROG18A's `FoxJ3_N3/` holds six chimeras
    (`FoxJ3_N3_6aa`, `FoxJ3_N3_loop`, ...), NAR11's `HLH-1/` holds a point-mutant series
    (`HLH-1_L13R`, `_L13T`, `_L13V`), LIN14B's `ANAC092/` holds a domain construct and a
    full-length one. Each is a different protein, and treating them as replicates of one
    sequence would average distinct proteins together and attach the wrong sequence to the
    measurements.

    The construct is encoded as a path element, so members are matched against the detail
    page's own construct names, longest first (`HLH-1_L13R` must win over `HLH-1`). Members
    matching nothing fall back to the single-construct case.
    """
    if len(constructs) <= 1:
        return {constructs[0]: members} if constructs else {}

    out: dict[str, list[str]] = {}
    ordered = sorted(constructs, key=len, reverse=True)
    for m in members:
        parts = set(m.split("/"))
        hit = next((c for c in ordered if c in parts), None)
        if hit is not None:
            out.setdefault(hit, []).append(m)
    return {k: sorted(v) for k, v in out.items()}


def parse_panel(
    accession: str,
    genes: list[str] | None = None,
    threshold_path: str | Path | None = None,
) -> tuple[pd.DataFrame, list[SkippedProtein]]:
    """Parse one panel. Returns the frame and the proteins that had to be skipped."""
    panel = PANELS[accession]
    cfg = thresholds.for_assay("pbm", threshold_path)
    pos_cut, neg_cut = float(cfg["positive"]), float(cfg["negative"])
    if not cfg.get("per_experiment", False):
        raise ValueError("UniPROBE panels binarize per experiment; thresholds.yaml disables it")

    root = raw_dir(accession)
    archive = root / panel.archive
    details = _uniprobe.load_details(root / "details")
    experiments = group_experiments(archive)

    frames: list[pd.DataFrame] = []
    skipped: list[SkippedProtein] = []

    # Resolve each gene's DBD, then group genes that share one. Within a panel, paralogues
    # routinely have byte-identical domains — Cell08 has 11 such pairs (Evx1/Evx2, Lhx2/Lhx9,
    # Pitx2/Pitx3, ...). Their measurements are independent experiments on the same protein
    # sequence, so they are reconciled exactly like replicates. Keeping only the first would
    # discard real data and, worse, hide whether the two agree.
    domain_cfg = thresholds.load(threshold_path)["domain"]
    pad = canonical.padding(domain_cfg)
    # Every construct on every detail page, not one per gene: a gene folder can hold several
    # distinct engineered proteins.
    constructs = {
        name: seq
        for g, page in details.items()
        if not genes or g in genes
        for name, seq in page.inserts.items()
    }
    hits = domains.scan(constructs)

    groups: dict[str, dict] = {}
    unmatched: list[str] = []
    for gene, members in experiments.items():
        if genes and gene not in genes:
            continue
        page = details.get(gene)
        if page is None:
            # Not a download gap: UniPROBE publishes no detail page for these at all.
            # Most are protein complexes -- Myc_Max, Kay_Jra (Fos/Jun), Da_Twi, the C.
            # elegans HLH-2 heterodimers, the CSL/NOTCH/MAML ternary complexes. Two or three
            # different chains form one binding unit, so no single dbd_seq is responsible
            # for the measurement and they fail the first admission condition anyway.
            skipped.append(
                SkippedProtein(gene, "UniPROBE publishes no sequence (often a protein complex)")
            )
            continue
        if not page.inserts:
            skipped.append(SkippedProtein(gene, "detail page has no clone insert sequence"))
            continue

        by_construct = split_by_construct(members, list(page.inserts))
        placed = sum(len(v) for v in by_construct.values())
        if placed < len(members):
            unmatched.append(f"{gene}: {len(members) - placed} experiment(s) matched no construct")

        for construct, cmembers in by_construct.items():
            # The construct that was on the array is what gets annotated -- never the
            # protein's own domain annotation, which describes the full-length protein.
            insert = page.inserts[construct]
            call = domains.call_domain(insert, hits[construct], domain_cfg)
            if not call.ok:
                skipped.append(
                    SkippedProtein(construct, f"rejected: {call.rejection} ({call.family})")
                )
                continue
            # Stored sequence is canonical, so the same domain is one string however much
            # flank this particular lab happened to clone.
            h = call.hits[0]
            ref = (
                None
                if canonical.covers(insert, h.start, h.end, pad)
                else reference_db.resolve(page.protein_id)
            )
            cc = canonical.canonicalise(insert, h.start, h.end, reference=ref, pad=pad)
            if not cc.ok:
                skipped.append(SkippedProtein(construct, f"not canonicalisable: {cc.rejection}"))
                continue
            g = groups.setdefault(
                cc.sequence,
                {
                    "genes": [],
                    "members": [],
                    "call": call,
                    "page": page,
                    "gene": gene,
                    "inserts": {},
                },
            )
            g["genes"].append(construct)
            g["members"].extend(cmembers)
            g["inserts"][construct] = insert

    for u in unmatched:
        skipped.append(SkippedProtein("<unmatched>", u))

    # Constructs sharing a canonical domain but differing OUTSIDE it carry a signal the model
    # can never see -- `Cell09:HLH-25` and `HLH-27` are distinct genes with identical 76 aa
    # domains and four substitutions just beyond. Bin the group.
    for dbd in [d for d, g in groups.items() if len(g["inserts"]) > 1]:
        if canonical.conflicting_constructs(groups[dbd]["inserts"]):
            skipped.append(
                SkippedProtein(
                    "/".join(sorted(groups[dbd]["inserts"])),
                    "constructs share a canonical domain but differ outside it; discarded",
                )
            )
            del groups[dbd]

    # Cluster engineered variants with their reference. A gene folder holding several
    # constructs is a variant series: NAR11's HLH-1 carries L13R/L13T/L13V, ROG18A's FoxN3
    # carries chimeras. The construct named after the gene is the reference; the rest are
    # variants of it. Lengths need not match: distance is alignment-based, with terminal
    # gaps free so that padding clipped by a short construct is not counted as an indel.
    references = {g["gene"]: seq for seq, g in groups.items() if g["gene"] in g["genes"]}

    def find_reference(gene: str, seq: str) -> tuple[str, str] | None:
        """The reference this construct varies from, if there is an unambiguous one."""
        ref = references.get(gene)
        if ref is not None:
            return gene, ref
        # Sibling folders: ROG18A puts FoxN3 and its chimeras FoxN3_J3_* side by side, so
        # the parent is the longest reference whose name prefixes this one. Equal length is
        # no longer required -- distance comes from an alignment, and the padded envelopes of
        # a protein and its own chimera routinely differ by however much padding fitted.
        cands = [(g, r) for g, r in references.items() if g != gene and gene.startswith(g + "_")]
        return max(cands, key=lambda c: len(c[0])) if cands else None

    with zipfile.ZipFile(archive) as z:
        for dbd_seq, g in groups.items():
            group_genes = sorted(g["genes"])
            found = find_reference(g["gene"], dbd_seq)
            if found is not None:
                ref_gene, ref_seq = found
                profile = align.edit_profile(ref_seq, dbd_seq)
                wt_id = f"{accession}:{ref_gene}"
                n_mut = profile.n_edits
                mut_positions = profile.positions_str
            else:
                # No same-length reference: this construct is its own reference.
                wt_id = f"{accession}:{'/'.join(group_genes)}"
                n_mut, mut_positions = 0, ""
            if len(group_genes) > 1:
                skipped.append(
                    SkippedProtein(
                        "/".join(group_genes),
                        f"{len(group_genes)} genes share one DBD sequence; their "
                        f"{len(g['members'])} experiments were reconciled as replicates",
                    )
                )
            page = g["page"]
            dna_seq, label, raw_score = _uniprobe.reconcile_replicates(
                z, sorted(g["members"]), pos_cut, neg_cut
            )
            frames.append(
                _uniprobe.build_frame(
                    dna_seq=dna_seq,
                    label=label,
                    raw_score=raw_score,
                    dbd_seq=dbd_seq,
                    dbd_family=g["call"].family,
                    dbd_source=DBD_SOURCE,
                    wt_id=wt_id,
                    n_mut_from_wt=n_mut,
                    mut_positions=mut_positions,
                    protein_id=page.protein_id,
                    species=page.species,
                    pos_cut=pos_cut,
                    neg_cut=neg_cut,
                    source_dataset=accession,
                    source_file=panel.source_file,
                )
            )

    if not frames:
        raise ValueError(f"{accession}: no experiments parsed")
    return schema.coerce(pd.concat(frames, ignore_index=True)), skipped


def make_parser(accession: str):
    """A zero-argument `parse()` for one accession, as the parser registry expects."""

    def parse() -> pd.DataFrame:
        return parse_panel(accession)[0]

    parse.__name__ = f"parse_{accession.lower()}"
    parse.__doc__ = f"Parse the {accession} panel ({PANELS[accession].citation})."
    return parse
