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

from snp2prot.parsers._pbm import (
    ASSAY,
    DNA_CONTEXT,
    ESCORE_MAX,
    ESCORE_MIN,
    N_NONREDUNDANT_8MERS,
    SCORE_TYPE,
    EscoreColumnError,
    binarize,
    build_frame,
    check_complete,
    combine_replicates,
    escore_column,
)

#: Re-exported so callers and tests keep one import site for this source's machinery. The
#: definitions live in `_pbm` because they describe the assay, not UniPROBE's file layout.
__all__ = [
    "ASSAY",
    "CONTIG_8MER_RE",
    "DNA_CONTEXT",
    "ESCORE_MAX",
    "ESCORE_MIN",
    "N_NONREDUNDANT_8MERS",
    "SCORE_TYPE",
    "DetailPage",
    "EscoreColumnError",
    "binarize",
    "build_frame",
    "clean_sequence",
    "escore_column",
    "load_details",
    "parse_detail_page",
    "read_8mer_table",
    "reconcile_replicates",
    "strip_html",
]

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


def read_8mer_table(z: zipfile.ZipFile, member: str, require_complete: bool = True) -> pd.DataFrame:
    """One experiment's 8-mer table as `dna_seq` + `escore`.

    **The header is detected, never assumed.** EMBO10's first line is data, and skipping it
    would silently drop the `AAAAAAAA` 8-mer from every experiment in that panel while leaving
    row counts one short in a way nothing else would catch.

    The E-score column is then identified by `_pbm.escore_column`, because position is not
    safe across UniPROBE: the column order and count both vary, and at least one accession has
    no E-score column at all.

    | accession | layout |
    |---|---|
    | `BAR15A`, `EMBO10`, `PNAS13` | 5 cols, E-score at index 2 |
    | `Cell08` | 7 cols with a descriptive header, E-score ("enrichment score") at index 2 |
    | `GR09`, `SCI09` | headerless, median intensity at index 2, **E-score at index 3** |
    | `RAD13A` | 4 cols, `8-mer / 8-mer / Median / Z-score` — **no E-score at all** |
    """
    with z.open(member) as fh:
        first = fh.readline().decode("utf8", errors="ignore")
    has_header = first.lower().lstrip().startswith("8-mer")

    with z.open(member) as fh:
        df = pd.read_csv(fh, sep="\t", header=0 if has_header else None, dtype=str)
    if df.shape[1] < 3:
        raise EscoreColumnError(f"{member}: only {df.shape[1]} columns")

    escore = escore_column(df, has_header, member)
    if require_complete:
        check_complete(len(df), member)
    return _finish(df, escore)


def _finish(df: pd.DataFrame, escore: pd.Series) -> pd.DataFrame:
    out = pd.DataFrame(
        {"dna_seq": df.iloc[:, 0].astype(str).str.strip().str.upper(), "escore": escore}
    )
    return out.dropna(subset=["escore"])


def reconcile_replicates(
    z: zipfile.ZipFile,
    members: list[str],
    pos_cut: float,
    neg_cut: float,
    require_complete: bool = True,
) -> tuple[pd.Series, np.ndarray, np.ndarray]:
    """Read each replicate, then combine them (see `_pbm.combine_replicates`).

    A file we cannot read an E-score from is dropped with its reason rather than silently
    contributing whatever column happened to sit in that position.
    """
    tables, skipped = [], []
    for member in members:
        try:
            tables.append((member, read_8mer_table(z, member, require_complete)))
        except EscoreColumnError as exc:
            skipped.append(str(exc))
    if not tables:
        raise EscoreColumnError("; ".join(skipped) or "no readable replicate")
    return combine_replicates(tables, pos_cut, neg_cut)
