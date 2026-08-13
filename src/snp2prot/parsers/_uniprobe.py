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

#: A universal PBM E-score is a rank statistic bounded to this interval by construction.
#: The bound is what lets the E-score column be identified when a file has no header.
ESCORE_MIN, ESCORE_MAX = -0.5, 0.5

ASSAY = "PBM"
SCORE_TYPE = "pbm_escore"
#: A PBM 8-mer score aggregates over many flanking contexts, so no specific flank was observed.
DNA_CONTEXT = "core_only"

#: Every contiguous-8-mer spelling seen so far. `_8mers_11111111.txt` is the gapped-k-mer
#: naming for the all-positions-contiguous pattern; `_contig8mers.txt` is Cell08's spelling.
CONTIG_8MER_RE = re.compile(r"(_8mers(_11111111)?|_contig8mers)\.txt$")


#: `<dt>LABEL</dt> <dd> <kbd>SEQUENCE</kbd>` — `<dd>` and `<kbd>` sit on separate lines.
ESCORE_HEADER_RE = re.compile(r"e[-_ ]?score|enrichment", re.I)

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
    #: construct name -> clone insert sequence, keyed by the FULL name the page gives
    #: ("ARX_L343Q", "FoxJ3_N3_6aa", "HLH-1_L13R"). The full name matters: UniPROBE archives
    #: encode the construct as a path element, so keying by anything shorter makes distinct
    #: engineered proteins impossible to tell apart and silently merges them.
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
            # The label names the construct ("ARX_L343Q Insert Sequence"), but not always:
            # it can name the plasmid ("Clone pTH3418 insert sequence") or nothing at all
            # (GR09 uses a bare "Insert sequence"). Both fall back to the gene, which is the
            # construct in those panels. An empty key would fuse every such gene into one
            # cluster -- GR09 alone would collapse 88 proteins.
            name = label[: -len(" Insert Sequence")].strip()
            if low.startswith("clone ") or not name:
                page.inserts.setdefault(gene, seq)
            else:
                page.inserts[name] = seq
    return page


def load_details(details_dir: Path) -> dict[str, DetailPage]:
    pages = sorted(details_dir.glob("*.html"))
    if not pages:
        raise FileNotFoundError(f"no detail pages under {details_dir}")
    return {
        p.stem: parse_detail_page(p.read_text(encoding="utf8", errors="ignore"), p.stem)
        for p in pages
    }


class EscoreColumnError(ValueError):
    """Raised when a file's E-score column cannot be identified unambiguously."""


def read_8mer_table(z: zipfile.ZipFile, member: str) -> pd.DataFrame:
    """One experiment's 8-mer table as `dna_seq` + `escore`.

    **The E-score column is identified by its value range, not by position.** Position is not
    safe across UniPROBE: the column order and count both vary, and at least one accession
    has no E-score column at all.

    | accession | layout |
    |---|---|
    | `BAR15A`, `EMBO10`, `PNAS13` | 5 cols, E-score at index 2 |
    | `Cell08` | 7 cols with a descriptive header, E-score ("enrichment score") at index 2 |
    | `GR09`, `SCI09` | headerless, median intensity at index 2, **E-score at index 3** |
    | `RAD13A` | 4 cols, `8-mer / 8-mer / Median / Z-score` — **no E-score at all** |

    A universal PBM E-score is a rank statistic bounded to [-0.5, 0.5] that always takes
    negative values across 32,896 8-mers, since most are unbound. Intensities run to
    hundreds of thousands, z-scores past 1, and p/q-values are non-negative -- so those
    two properties together single the E-score out. A header naming it wins outright.

    Zero candidates means the file has no E-score; more than one means it concatenates
    several experiments side by side
    (SCI09 ships 20-column combined files alongside its 9-column per-replicate ones). Both
    are errors here rather than a silent guess.
    """
    with z.open(member) as fh:
        first = fh.readline().decode("utf8", errors="ignore")
    has_header = first.lower().lstrip().startswith("8-mer")

    with z.open(member) as fh:
        df = pd.read_csv(fh, sep="\t", header=0 if has_header else None, dtype=str)
    if df.shape[1] < 3:
        raise EscoreColumnError(f"{member}: only {df.shape[1]} columns")

    numeric = {}
    for i in range(1, df.shape[1]):
        vals = pd.to_numeric(df.iloc[:, i], errors="coerce")
        if vals.notna().sum() >= len(df) * 0.9:
            numeric[i] = vals

    # A named column wins outright: Cell08 calls it "enrichment score".
    if has_header:
        named = [i for i in numeric if ESCORE_HEADER_RE.search(str(df.columns[i]))]
        if len(named) == 1:
            return _finish(df, numeric[named[0]])

    # Otherwise: bounded to [-0.5, 0.5] AND actually taking negative values. The bound alone
    # is not enough -- p-values and q-values sit in [0, 0.5] too -- but an E-score over 32,896
    # 8-mers always runs negative, since most 8-mers are not bound, while a probability cannot.
    candidates = [
        (i, v)
        for i, v in numeric.items()
        if v.min() >= ESCORE_MIN and v.max() <= ESCORE_MAX and v.min() < 0
    ]

    if not candidates:
        raise EscoreColumnError(
            f"{member}: no column lies within [{ESCORE_MIN}, {ESCORE_MAX}] — this file "
            f"carries no E-score (columns: {list(df.columns)[:6]})"
        )
    if len(candidates) > 1:
        raise EscoreColumnError(
            f"{member}: {len(candidates)} columns look like E-scores (indices "
            f"{[i for i, _ in candidates]}); the file probably concatenates experiments"
        )

    return _finish(df, candidates[0][1])


def _finish(df: pd.DataFrame, escore: pd.Series) -> pd.DataFrame:
    out = pd.DataFrame(
        {"dna_seq": df.iloc[:, 0].astype(str).str.strip().str.upper(), "escore": escore}
    )
    return out.dropna(subset=["escore"])


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
    skipped: list[str] = []
    for member in members:
        try:
            exp = read_8mer_table(z, member)
        except EscoreColumnError as exc:
            # A file we cannot read an E-score from is dropped with its reason, rather than
            # silently contributing whatever column happened to sit in that position.
            skipped.append(str(exc))
            continue
        exp = exp.sort_values("dna_seq", kind="stable").reset_index(drop=True)
        if keys is None:
            keys = exp["dna_seq"]
        elif not keys.equals(exp["dna_seq"]):
            raise ValueError(f"{member}: 8-mer set differs from the first replicate")
        e = exp["escore"].to_numpy(dtype=float)
        labels.append(binarize(e, pos_cut, neg_cut))
        scores.append(e)

    if not labels:
        raise EscoreColumnError("; ".join(skipped) or "no readable replicate")

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
