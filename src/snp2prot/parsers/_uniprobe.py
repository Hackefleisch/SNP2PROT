"""Machinery shared by every UniPROBE accession.

UniPROBE is one database but not one format. Across the four accessions parsed so far the
archive layout, the 8-mer filename and even the presence of a header row all differ:

| accession | member path                                   | header row |
|-----------|-----------------------------------------------|------------|
| `BAR15A`  | `GENE/GENE_ALLELE/GENE_ALLELE_REP/..._8mers.txt` (and `_8mers_11111111.txt`) | yes |
| `Cell08`  | `GENE/EXPID/GENE_EXPID_contig8mers.txt`       | yes        |
| `EMBO10`  | `GENE/GENE_8mers.txt`                         | **no**     |
| `PNAS13`  | `GENE/GENE_8mers_11111111.txt`                | yes        |

So the header is **detected**, never assumed — `EMBO10`'s first line is data, and skipping it
would silently drop the `AAAAAAAA` 8-mer from every experiment in that panel while leaving
row counts one short in a way nothing else would catch.

Detail pages vary too: the UniProt accession is labelled `Swiss-Prot` on some pages and
`Uniprot` on others, and clone sequences appear either as `<ALLELE> Insert Sequence`
(BAR15A, which needs them because variants exist only in that form) or as
`Clone <ID> insert sequence` (the single-protein panels).
"""

from __future__ import annotations

import re
import zipfile
from dataclasses import dataclass, field
from pathlib import Path

import numpy as np
import pandas as pd

from snp2prot import schema

ASSAY = "PBM"
SCORE_TYPE = "pbm_escore"
#: A PBM 8-mer score aggregates over many flanking contexts, so no specific flank was observed.
DNA_CONTEXT = "core_only"

#: Every contiguous-8-mer spelling seen so far. `_8mers_11111111.txt` is the gapped-k-mer
#: naming for the all-positions-contiguous pattern; `_contig8mers.txt` is Cell08's spelling.
CONTIG_8MER_RE = re.compile(r"(_8mers(_11111111)?|_contig8mers)\.txt$")


#: `<dt>LABEL</dt> <dd> <kbd>SEQUENCE</kbd>` — `<dd>` and `<kbd>` sit on separate lines.
SEQ_BLOCK_RE = re.compile(r"<dt>([^<]*?)</dt>\s*<dd>\s*<kbd>(.*?)</kbd>", re.S)
FIELD_RE = re.compile(r"<dt>\s*(Domain|Swiss-Prot|Uniprot|Species)\s*</dt>\s*<dd>(.*?)</dd>", re.S)


def strip_html(fragment: str) -> str:
    return re.sub(r"\s+", " ", re.sub(r"<[^>]+>", " ", fragment)).strip()


def clean_sequence(fragment: str) -> str:
    """Strip UniPROBE's numbered, space-broken sequence formatting down to residues."""
    text = re.sub(r"<[^>]+>", " ", fragment).replace("&nbsp;", " ")
    return "".join(re.findall(r"[A-Z]", text))


@dataclass
class DetailPage:
    """What one UniPROBE detail page tells us about one protein."""

    gene: str
    family: str = ""
    protein_id: str = ""
    species: str = ""
    dbd: str = ""
    #: allele label -> clone insert sequence. `REF` for a reference allele.
    inserts: dict[str, str] = field(default_factory=dict)


def parse_detail_page(html: str, gene: str) -> DetailPage:
    """Parse one detail page, tolerating both label vocabularies."""
    fields = {m.group(1): strip_html(m.group(2)) for m in FIELD_RE.finditer(html)}
    page = DetailPage(
        gene=gene,
        family=fields.get("Domain", ""),
        # Pages use one label or the other, never both.
        protein_id=fields.get("Swiss-Prot") or fields.get("Uniprot", ""),
        species=fields.get("Species", ""),
    )
    for m in SEQ_BLOCK_RE.finditer(html):
        label, seq = m.group(1).strip(), clean_sequence(m.group(2))
        low = label.lower()
        if low.endswith("dna binding domain"):
            page.dbd = seq
        elif low.endswith("insert sequence"):
            if low.startswith("clone "):
                page.inserts.setdefault("REF", seq)
            else:
                allele_full = label[: -len(" Insert Sequence")].strip()
                _, _, allele = allele_full.partition("_")
                page.inserts[allele or "REF"] = seq
    return page


def load_details(details_dir: Path) -> dict[str, DetailPage]:
    pages = sorted(details_dir.glob("*.html"))
    if not pages:
        raise FileNotFoundError(f"no detail pages under {details_dir}")
    return {
        p.stem: parse_detail_page(p.read_text(encoding="utf8", errors="ignore"), p.stem)
        for p in pages
    }


def read_8mer_table(z: zipfile.ZipFile, member: str) -> pd.DataFrame:
    """One experiment's 8-mer table as `dna_seq` + `escore`.

    Neither the header nor the column count is assumed. EMBO10's files have no header row at
    all; Cell08's have a seven-column descriptive one ("enrichment score", "pvalue", "Qvalue")
    where the others have five terse ones. Columns 0 and 2 are the 8-mer and its E-score in
    every layout seen, so they are selected positionally and renamed.
    """
    with z.open(member) as fh:
        first = fh.readline().decode("utf8", errors="ignore")
    has_header = first.lower().lstrip().startswith("8-mer")

    with z.open(member) as fh:
        df = pd.read_csv(fh, sep="\t", header=0 if has_header else None, usecols=[0, 2])
    df.columns = ["dna_seq", "escore"]
    df["dna_seq"] = df["dna_seq"].astype(str).str.strip().str.upper()
    df["escore"] = pd.to_numeric(df["escore"], errors="coerce")
    return df.dropna(subset=["escore"])


def binarize(escores: np.ndarray, pos_cut: float, neg_cut: float) -> np.ndarray:
    """Per-experiment binarization. Never call this on scores pooled across experiments."""
    lab = np.full(len(escores), schema.LABEL_GRAY, dtype=np.int8)
    lab[escores >= pos_cut] = schema.LABEL_BIND
    lab[escores <= neg_cut] = schema.LABEL_NONBIND
    return lab


def reconcile_replicates(
    z: zipfile.ZipFile, members: list[str], pos_cut: float, neg_cut: float
) -> tuple[pd.Series, np.ndarray, np.ndarray]:
    """Binarize each replicate independently, then combine.

    Replicates can sit on different array designs whose intensity scales differ by an order of
    magnitude, so E-scores are never averaged *before* thresholding. Replicates that disagree
    on a given 8-mer fall to the gray band rather than being resolved by majority or by mean.

    Returns `(dna_seq, label, mean_escore)`.
    """
    labels, scores, keys = [], [], None
    for member in members:
        exp = read_8mer_table(z, member).sort_values("dna_seq", kind="stable")
        exp = exp.reset_index(drop=True)
        if keys is None:
            keys = exp["dna_seq"]
        elif not keys.equals(exp["dna_seq"]):
            raise ValueError(f"{member}: 8-mer set differs from the first replicate")
        e = exp["escore"].to_numpy(dtype=float)
        labels.append(binarize(e, pos_cut, neg_cut))
        scores.append(e)

    stacked = np.vstack(labels)
    agreed = (stacked == stacked[0]).all(axis=0)
    label = np.where(agreed, stacked[0], schema.LABEL_GRAY).astype(np.int8)
    return keys, label, np.vstack(scores).mean(axis=0)


def build_frame(
    *,
    dna_seq: pd.Series,
    label: np.ndarray,
    raw_score: np.ndarray,
    dbd_seq: str,
    dbd_family: str,
    dbd_source: str,
    wt_id: str,
    n_mut_from_wt: int,
    mut_positions: str,
    protein_id: str,
    species: str,
    pos_cut: float,
    neg_cut: float,
    source_dataset: str,
    source_file: str,
) -> pd.DataFrame:
    """Assemble one experiment's rows. Scalars broadcast, so each string is stored once."""
    frame = pd.DataFrame({"dna_seq": dna_seq.to_numpy(), "label": label, "raw_score": raw_score})
    frame["dbd_seq"] = dbd_seq
    frame["dbd_family"] = dbd_family
    frame["dbd_source"] = dbd_source
    frame["wt_id"] = wt_id
    frame["n_mut_from_wt"] = n_mut_from_wt
    frame["mut_positions"] = mut_positions
    frame["protein_id"] = protein_id
    frame["species"] = species
    frame["dna_context"] = DNA_CONTEXT
    frame["score_type"] = SCORE_TYPE
    frame["threshold_pos"] = pos_cut
    frame["threshold_neg"] = neg_cut
    frame["assay"] = ASSAY
    frame["stringency"] = ""
    frame["neg_provenance"] = np.where(label == schema.LABEL_NONBIND, "assayed_unbound", None)
    frame["source_dataset"] = source_dataset
    frame["source_file"] = source_file
    return frame
