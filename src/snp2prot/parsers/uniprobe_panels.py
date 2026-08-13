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

The last five were added to rebuild family breadth after the domain policy left the corpus
79.5% homeodomain. Every one is screened by the same policy, so a panel contributes only the
constructs whose signal is attributable to a single continuous domain.

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

from snp2prot import domains, schema, thresholds
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
    for p in (
        Panel("Cell08", "Berger et al., Cell 2008", "Cell08_contig8mers.zip"),
        Panel("EMBO10", "Wei et al., EMBO J 2010", "EMBO10_contig8mers.zip"),
        Panel("PNAS13", "Nakagawa et al., PNAS 2013", "PNAS13_contig8mers.zip"),
        Panel("SCI09", "Badis et al., Science 2009", "SCI09_contig8mers.zip"),
        Panel("GR09", "Zhu et al., Genome Res 2009", "GR09_contig8mers.zip"),
        Panel("MAR17A", "Mariani et al., Cell Systems 2017", "MAR17A_contig8mers.zip"),
        Panel("SHO18A", "Shokri et al., Cell Reports 2019", "SHO18A_contig8mers.zip"),
        Panel("ROG18A", "Rogers et al., Mol Cell 2019", "ROG18A_contig8mers.zip"),
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
        out.setdefault(m.split("/")[0], []).append(m)
    return {k: sorted(v) for k, v in sorted(out.items())}


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
    constructs = {
        g: next(iter(p.inserts.values()))
        for g, p in details.items()
        if p.inserts and (not genes or g in genes)
    }
    hits = domains.scan(constructs)

    groups: dict[str, dict] = {}
    for gene, members in experiments.items():
        if genes and gene not in genes:
            continue
        page = details.get(gene)
        if page is None:
            skipped.append(SkippedProtein(gene, "no detail page downloaded"))
            continue
        if not page.inserts:
            skipped.append(SkippedProtein(gene, "detail page has no clone insert sequence"))
            continue
        # The construct that was on the array is what gets annotated -- never the protein's
        # own domain annotation, which describes the full-length protein and would say
        # "Homeobox, POU" for a construct that in fact carries both domains.
        construct = next(iter(page.inserts.values()))
        call = domains.call_domain(construct, hits[gene], domain_cfg)
        if not call.ok:
            skipped.append(SkippedProtein(gene, f"rejected: {call.rejection} ({call.family})"))
            continue
        g = groups.setdefault(
            call.sequence, {"genes": [], "members": [], "call": call, "page": page}
        )
        g["genes"].append(gene)
        g["members"].extend(members)

    with zipfile.ZipFile(archive) as z:
        for dbd_seq, g in groups.items():
            group_genes = sorted(g["genes"])
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
                    # No engineered variants: each distinct domain is its own reference.
                    wt_id=f"{accession}:{'/'.join(group_genes)}",
                    n_mut_from_wt=0,
                    mut_positions="",
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
