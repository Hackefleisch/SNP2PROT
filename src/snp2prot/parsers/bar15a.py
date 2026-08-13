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
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from snp2prot import schema, thresholds
from snp2prot.config import raw_dir

SOURCE = "BAR15A"
ASSAY = "PBM"
SCORE_TYPE = "pbm_escore"
DBD_SOURCE = "uniprobe_clone_insert"
DNA_CONTEXT = "core_only"  # a PBM 8-mer score aggregates over many flanking contexts

ARCHIVE = "BAR15A_contig8mers.zip"
DETAILS = "details"

#: Matches both spellings of the contiguous-8-mer file. See note 1 in the module docstring.
CONTIG_8MER_RE = re.compile(r"_8mers(_11111111)?\.txt$")

#: `<dt>LABEL</dt> <dd> <kbd>SEQUENCE</kbd>` — the detail pages put `<dd>` and `<kbd>` on
#: separate lines, so the whitespace between them is not optional.
SEQ_BLOCK_RE = re.compile(r"<dt>([^<]*?)</dt>\s*<dd>\s*<kbd>(.*?)</kbd>", re.S)
FIELD_RE = re.compile(r"<dt>\s*(Domain|Swiss-Prot|Species)\s*</dt>\s*<dd>(.*?)</dd>", re.S)
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


@dataclass(frozen=True)
class GeneMeta:
    """Everything the detail page tells us about one gene."""

    gene: str
    family: str
    protein_id: str
    species: str
    inserts: dict[str, str]  # allele -> clone insert sequence

    @property
    def ref(self) -> str:
        return self.inserts["REF"]


def _clean_sequence(html_fragment: str) -> str:
    """Strip the numbered, space-broken sequence formatting down to residues."""
    text = re.sub(r"<[^>]+>", " ", html_fragment).replace("&nbsp;", " ")
    return "".join(re.findall(r"[A-Z]", text))


def _clean_field(html_fragment: str) -> str:
    return re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", html_fragment)).strip()


def load_metadata(details_dir: Path | None = None) -> dict[str, GeneMeta]:
    """Parse every `details/<GENE>.html` page into a `GeneMeta`."""
    details_dir = details_dir or (raw_dir(SOURCE) / DETAILS)
    pages = sorted(details_dir.glob("*.html"))
    if not pages:
        raise FileNotFoundError(f"no detail pages under {details_dir}")

    out: dict[str, GeneMeta] = {}
    for page in pages:
        html = page.read_text(encoding="utf8", errors="ignore")
        fields = {m.group(1): _clean_field(m.group(2)) for m in FIELD_RE.finditer(html)}
        inserts = {}
        for m in SEQ_BLOCK_RE.finditer(html):
            label = m.group(1).strip()
            if label.endswith("Insert Sequence"):
                allele_full = label[: -len(" Insert Sequence")].strip()
                _, _, allele = allele_full.partition("_")
                inserts[allele] = _clean_sequence(m.group(2))
        if "REF" not in inserts:
            raise ValueError(f"{page.name}: no REF insert sequence found")
        out[page.stem] = GeneMeta(
            gene=page.stem,
            family=fields.get("Domain", ""),
            protein_id=fields.get("Swiss-Prot", ""),
            species=fields.get("Species", ""),
            inserts=inserts,
        )
    return out


def find_experiments(archive: Path | None = None) -> dict[str, list[str]]:
    """Map `GENE_ALLELE` -> list of member paths inside the archive, one per replicate."""
    archive = archive or (raw_dir(SOURCE) / ARCHIVE)
    with zipfile.ZipFile(archive) as z:
        members = [n for n in z.namelist() if CONTIG_8MER_RE.search(n)]
    out: dict[str, list[str]] = {}
    for m in members:
        out.setdefault(m.split("/")[1], []).append(m)
    return {k: sorted(v) for k, v in sorted(out.items())}


def _read_experiment(z: zipfile.ZipFile, member: str) -> pd.DataFrame:
    """One experiment's 8-mer table, as `dna_seq` + `escore`."""
    with z.open(member) as fh:
        df = pd.read_csv(fh, sep="\t", usecols=[0, 2], names=["dna_seq", "escore"], skiprows=1)
    df["dna_seq"] = df["dna_seq"].str.strip().str.upper()
    return df


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
            gm = meta[gene]
            insert = gm.inserts[allele]
            n_mut, mut_positions = _mutation_bookkeeping(gm.ref, insert)

            # Binarize each replicate independently, then reconcile. Pooling E-scores across
            # replicates would cross array designs, whose scales differ (brief section 4).
            labels, scores, keys = [], [], None
            for member in members:
                exp = _read_experiment(z, member)
                exp = exp.sort_values("dna_seq", kind="stable").reset_index(drop=True)
                if keys is None:
                    keys = exp["dna_seq"]
                elif not keys.equals(exp["dna_seq"]):
                    raise ValueError(f"{member}: 8-mer set differs from the first replicate")
                e = exp["escore"].to_numpy(dtype=float)
                lab = np.full(len(e), schema.LABEL_GRAY, dtype=np.int8)
                lab[e >= pos_cut] = schema.LABEL_BIND
                lab[e <= neg_cut] = schema.LABEL_NONBIND
                labels.append(lab)
                scores.append(e)

            stacked = np.vstack(labels)
            agreed = (stacked == stacked[0]).all(axis=0)
            label = np.where(agreed, stacked[0], schema.LABEL_GRAY).astype(np.int8)
            raw_score = np.vstack(scores).mean(axis=0)

            frame = pd.DataFrame(
                {
                    "dna_seq": keys.to_numpy(),
                    "label": label,
                    "raw_score": raw_score,
                }
            )
            frame["dbd_seq"] = insert
            frame["dbd_family"] = gm.family
            frame["dbd_source"] = DBD_SOURCE
            frame["wt_id"] = f"{SOURCE}:{gene}"
            frame["n_mut_from_wt"] = n_mut
            frame["mut_positions"] = mut_positions
            frame["protein_id"] = gm.protein_id
            frame["species"] = gm.species
            frame["dna_context"] = DNA_CONTEXT
            frame["score_type"] = SCORE_TYPE
            frame["threshold_pos"] = pos_cut
            frame["threshold_neg"] = neg_cut
            frame["assay"] = ASSAY
            frame["stringency"] = ""
            frame["neg_provenance"] = np.where(
                label == schema.LABEL_NONBIND, "assayed_unbound", None
            )
            frame["source_dataset"] = SOURCE
            frame["source_file"] = source_file
            frames.append(frame)

    if not frames:
        raise ValueError("no experiments parsed")
    return schema.coerce(pd.concat(frames, ignore_index=True))
