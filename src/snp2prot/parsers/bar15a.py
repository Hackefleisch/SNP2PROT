"""Barrera et al. 2016 (UniPROBE accession BAR15A) — universal PBM 8-mer E-scores.

41 human TF reference alleles plus their point-mutant variants, every allele scored over the
same 32,896 non-redundant 8-mers. This is the densest matched WT/mutant resource in the
project: within a gene, alleles differ by a single residue over an identical DNA space.

Two sources, both under `data/raw/BAR15A/`:

- `BAR15A_contig8mers.zip` — one `*_8mers.txt` per experiment, laid out as
  `GENE/GENE_ALLELE/GENE_ALLELE_REP/GENE_ALLELE_REP_8mers.txt`, columns
  `8-mer, 8-mer(revcomp), E-score, Median, Z-score`.
- `details/<GENE>.html` — the UniPROBE detail page, which carries the clone insert sequence
  for every allele plus the Pfam domain, Swiss-Prot accession and species.

Three things here are not obvious and are load-bearing:

1. **Eight of HOXD13's files are named `*_8mers_11111111.txt`, not `*_8mers.txt`.** A glob on
   the latter silently drops three alleles entirely (I297V, N298S, Q325K) and one replicate
   each from five more. `CONTIG_8MER_RE` matches both spellings.
2. **`dbd_seq` is the deposited clone insert, not a reconstructed Pfam domain.** Only 20 of
   the 41 detail pages carry a separate "DNA binding domain" field, so the insert is the one
   representation available for every allele — and it is what was actually on the array.
3. **Mutation bookkeeping is derived from the sequences, never from the allele name.** Four
   alleles disagree with their own names (see `ANOMALIES`); trusting the name would put the
   wrong residue at the wrong position in exactly the way the brief warns about.
"""

from __future__ import annotations

import re
import zipfile
from pathlib import Path

import pandas as pd

from snp2prot import schema, thresholds
from snp2prot.config import raw_dir
from snp2prot.parsers import _uniprobe

SOURCE = "BAR15A"
DBD_SOURCE = "uniprobe_clone_insert"

ARCHIVE = "BAR15A_contig8mers.zip"
DETAILS = "details"

#: Shared with the other UniPROBE accessions; see `_uniprobe.CONTIG_8MER_RE`.
CONTIG_8MER_RE = _uniprobe.CONTIG_8MER_RE
ALLELE_NAME_RE = re.compile(r"^([A-Z])(\d+)([A-Z])$")

#: Alleles whose deposited insert contradicts their own name. Characterised in Phase 1 by
#: deriving each gene's insert->protein offset from its self-consistent variants (all 41 genes
#: give a single unambiguous offset) and re-testing these against it. Reported, not silently
#: repaired — see reports/validation_BAR15A.md and reports/OPEN_ITEMS.md.
ANOMALIES: dict[str, str] = {
    "PITX2_T114P": (
        "deposited insert is byte-identical to REF; the T->P substitution at insert "
        "position 43 is simply absent, so no variant protein can be constructed"
    ),
    "POU3F4_A237G": (
        "carries two substitutions, A237G (intended) and A312V, the latter being POU3F4's "
        "other deposited allele; the clone appears to be a double mutant"
    ),
    "PAX6_R26G": (
        "off by one: name says R26G, insert changes Q27->G. Position 26 is R in both REF "
        "and variant"
    ),
    "PROP1_R112Q": (
        "right position and right wild-type residue, wrong substitution: name says R112Q, "
        "insert gives R112M"
    ),
}

#: Only this one cannot be represented at all — with no difference from REF it would collide
#: with the reference allele's own rows on `pair_id`. The other three are kept, with their
#: mutation bookkeeping taken from the sequence rather than the name.
UNUSABLE = {"PITX2_T114P"}


def load_metadata(details_dir: Path | None = None) -> dict[str, _uniprobe.DetailPage]:
    """Parse every `details/<GENE>.html` page. Every gene must expose a REF insert."""
    pages = _uniprobe.load_details(details_dir or (raw_dir(SOURCE) / DETAILS))
    for gene, page in pages.items():
        if "REF" not in page.inserts:
            raise ValueError(f"{gene}: no REF insert sequence found")
    return pages


def find_experiments(archive: Path | None = None) -> dict[str, list[str]]:
    """Map `GENE_ALLELE` -> list of member paths inside the archive, one per replicate."""
    archive = archive or (raw_dir(SOURCE) / ARCHIVE)
    with zipfile.ZipFile(archive) as z:
        members = [n for n in z.namelist() if CONTIG_8MER_RE.search(n)]
    out: dict[str, list[str]] = {}
    for m in members:
        out.setdefault(m.split("/")[1], []).append(m)
    return {k: sorted(v) for k, v in sorted(out.items())}


def _mutation_bookkeeping(ref: str, variant: str) -> tuple[int, str]:
    """`(n_mut_from_wt, mut_positions)` from the sequences themselves, 1-based in the DBD."""
    if len(ref) != len(variant):
        raise ValueError(f"insert length differs from REF: {len(variant)} vs {len(ref)}")
    positions = [i for i, (r, v) in enumerate(zip(ref, variant, strict=True), 1) if r != v]
    return len(positions), ",".join(str(p) for p in positions)


def parse(genes: list[str] | None = None, threshold_path: str | Path | None = None) -> pd.DataFrame:
    """Parse BAR15A into the unified schema.

    `genes` restricts the work to a subset, for development only — a partial table must not
    be written to `data/interim/`.
    """
    cfg = thresholds.for_assay("pbm", threshold_path)
    pos_cut, neg_cut = float(cfg["positive"]), float(cfg["negative"])
    if not cfg.get("per_experiment", False):
        raise ValueError("BAR15A binarizes per experiment; thresholds.yaml disables it")

    meta = load_metadata()
    experiments = find_experiments()
    archive = raw_dir(SOURCE) / ARCHIVE
    source_file = f"{SOURCE}/{ARCHIVE}"

    frames: list[pd.DataFrame] = []
    with zipfile.ZipFile(archive) as z:
        for allele_full, members in experiments.items():
            gene, _, allele = allele_full.partition("_")
            if genes and gene not in genes:
                continue
            if allele_full in UNUSABLE:
                continue
            page = meta[gene]
            insert = page.inserts[allele]
            n_mut, mut_positions = _mutation_bookkeeping(page.inserts["REF"], insert)

            dna_seq, label, raw_score = _uniprobe.reconcile_replicates(z, members, pos_cut, neg_cut)
            frames.append(
                _uniprobe.build_frame(
                    dna_seq=dna_seq,
                    label=label,
                    raw_score=raw_score,
                    dbd_seq=insert,
                    dbd_family=page.family,
                    dbd_source=DBD_SOURCE,
                    wt_id=f"{SOURCE}:{gene}",
                    n_mut_from_wt=n_mut,
                    mut_positions=mut_positions,
                    protein_id=page.protein_id,
                    species=page.species,
                    pos_cut=pos_cut,
                    neg_cut=neg_cut,
                    source_dataset=SOURCE,
                    source_file=source_file,
                )
            )

    if not frames:
        raise ValueError("no experiments parsed")
    return schema.coerce(pd.concat(frames, ignore_index=True))
