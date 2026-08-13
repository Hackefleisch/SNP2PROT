"""The protein-side companion table.

The main table stores one thing about the protein: `dbd_seq`, the padded DNA-binding domain.
That is the unit the labels are attributable to (see `docs/DOMAIN_POLICY.md`), but it is not
the only representation worth having — comparing embeddings computed on the bare domain, the
padded domain, the assayed construct and the full-length protein needs all four available
from one build.

So this table carries them, keyed by `dbd_seq`, and the main table stays exactly as the brief
specifies. Everything downstream is a join.

| column | meaning |
|---|---|
| `dbd_seq` | join key; the padded domain as stored in the main table |
| `dbd_bare` | the same domain without the padding |
| `pfam_family`, `pfam_start`, `pfam_end` | unpadded Pfam envelope, 1-based in the construct |
| `dbd_start`, `dbd_end` | padded envelope, 1-based in the construct |
| `construct_seq` | what was physically assayed |
| `full_seq` | canonical UniProt sequence, where the accession resolved |
| `dbd_start_protein`, `dbd_end_protein` | padded envelope mapped onto `full_seq` |
| `protein_id`, `gene`, `species`, `source_dataset` | provenance |
"""

from __future__ import annotations

from dataclasses import asdict, dataclass

import pandas as pd

COLUMNS = (
    "dbd_seq",
    "dbd_bare",
    "pfam_family",
    "pfam_start",
    "pfam_end",
    "dbd_start",
    "dbd_end",
    "construct_seq",
    "full_seq",
    "dbd_start_protein",
    "dbd_end_protein",
    "protein_id",
    "gene",
    "species",
    "source_dataset",
)


@dataclass
class ProteinRecord:
    dbd_seq: str
    dbd_bare: str
    pfam_family: str
    pfam_start: int
    pfam_end: int
    dbd_start: int
    dbd_end: int
    construct_seq: str
    full_seq: str | None
    dbd_start_protein: int | None
    dbd_end_protein: int | None
    protein_id: str
    gene: str
    species: str
    source_dataset: str


def locate_in_protein(construct: str, full: str, start: int, end: int) -> tuple[int, int] | None:
    """Map a construct-relative span onto the full-length protein.

    Exact substring match only. Clone constructs routinely carry vector-derived residues at
    their termini that are absent from the canonical sequence, so a partial match is
    attempted from the domain outwards; if even that fails the mapping is left null rather
    than guessed at.
    """
    if not full:
        return None
    offset = full.find(construct)
    if offset >= 0:
        return offset + start, offset + end
    # Fall back to locating the stored domain itself, which is interior to the construct.
    domain = construct[start - 1 : end]
    offset = full.find(domain)
    if offset >= 0:
        return offset + 1, offset + len(domain)
    return None


def to_frame(records: list[ProteinRecord]) -> pd.DataFrame:
    df = pd.DataFrame([asdict(r) for r in records])
    return df[list(COLUMNS)].drop_duplicates("dbd_seq").reset_index(drop=True)
