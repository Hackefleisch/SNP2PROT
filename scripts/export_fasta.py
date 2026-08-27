"""Export the corpus's sequences as FASTA for structure prediction — both axes.

`--what protein` writes the protein axis, `--what dna` the DNA axis, and the default writes
both. Neither ever alters a sequence: what goes into the file is what the dataset stores.

## The protein axis

One record per domain of the merged table — 1,338 of them — carrying `dbd_seq` **verbatim**.
That is the padded Pfam envelope the dataset stores (`docs/DOMAIN_POLICY.md`); nothing here
extends, trims or re-pads it, so a predicted structure's residue *i* is `dbd_seq[i - 1]` and
the join back to the dataset is the sequence itself.

Ordering and the `priority=` tag answer the collaborators' question of what to fold first.
The 281 domains that share a CD-HIT cluster with at least one other domain are the dense part
of the protein axis — a wild type and the single-residue variants around it — and they are
the records where a structure has to resolve differences of one residue. They come first,
grouped by cluster with the cluster reference ahead of its variants, and are tagged
`priority=high`. Everything else follows, tagged `priority=standard`.

Names are `<tier><rank>_<gene>_<sha8>`: sortable, unique, and readable in a directory listing,
with `<sha8>` the first 8 hex of the SHA-256 of `dbd_seq` — stable across rebuilds because it
depends on nothing but the sequence. The sidecar TSV is the full mapping back.

## The DNA axis

One record per 8-mer of the vocabulary every protein is scored against — 32,896 of them,
in the column order of `kmer_matrix.npz`. The dataset stores each reverse-complement pair
**once**, under whichever of the two strands sorts first, so a record stands for both strands
and its `rc=` field names the one that is not written. The 256 palindromes are their own
reverse complement and are flagged. Nothing is padded onto an 8-mer here — padding is a
modelling decision and the validator rejects it in the stored table (`CLAUDE.md`, rule 6).

    .venv/bin/python scripts/export_fasta.py        # ~15 s, both axes
"""

from __future__ import annotations

import argparse
import hashlib
import re
from pathlib import Path

import pandas as pd

from snp2prot import align, clusters, config, corpus
from snp2prot.data import matrix

#: FASTA line width for the sequence body.
WRAP = 60

#: Longest gene fragment kept in a record name; a few identifiers are 30+ characters and the
#: name only has to be readable, not complete — the TSV carries the untruncated value.
MAX_GENE_CHARS = 24

OUT_DIR = config.PROCESSED_DIR / "structures"

DNA_MANIFEST_COLUMNS = (
    "name",
    "matrix_index",
    "dna_seq",
    "revcomp",
    "palindrome",
    "n_domains_pos",
    "n_domains_neg",
    "n_domains_nocall",
    "max_escore",
)

MANIFEST_COLUMNS = (
    "name",
    "priority",
    "sha256",
    "gene",
    "protein_id",
    "dbd_family",
    "source_dataset",
    "wt_id",
    "role",
    "cluster_size",
    "n_edits_from_reference",
    "edit_positions",
    "verdict",
    "n_pos",
    "length",
    "dbd_seq",
)


def _safe(value: str) -> str:
    """Reduce a free-text identifier to something safe in a filename and a FASTA name."""
    cleaned = re.sub(r"[^A-Za-z0-9]+", "-", str(value)).strip("-")
    return cleaned or "unnamed"


def build(domains: pd.DataFrame, inventory: pd.DataFrame) -> pd.DataFrame:
    """The manifest: one row per domain, ordered as the FASTA is written."""
    reference = dict(zip(inventory.wt_id, inventory.reference, strict=True))

    rows = []
    for row in domains.itertuples(index=False):
        ref = reference.get(row.wt_id)
        # Edits are counted against the cluster reference, in the reference's frame, so every
        # member of a cluster shares one coordinate system (`CLAUDE.md`, conventions).
        profile = align.edit_profile(ref, row.dbd_seq) if ref else None
        rows.append(
            {
                "priority": "high" if row.cluster_size > 1 else "standard",
                "sha256": hashlib.sha256(row.dbd_seq.encode()).hexdigest(),
                "gene": row.gene,
                "protein_id": row.protein_id,
                "dbd_family": row.dbd_family,
                "source_dataset": row.source_dataset,
                "wt_id": row.wt_id,
                "role": "variant" if row.is_variant else "reference",
                "cluster_size": int(row.cluster_size),
                "n_edits_from_reference": None if profile is None else profile.n_edits,
                "edit_positions": "" if profile is None else ",".join(map(str, profile.positions)),
                "verdict": row.verdict,
                "n_pos": int(row.n_pos),
                "length": len(row.dbd_seq),
                "dbd_seq": row.dbd_seq,
            }
        )

    out = pd.DataFrame(rows)
    # High priority first, then whole clusters together with the reference ahead of its
    # variants; the rest by family so a folding run can be split along family lines.
    out["_tier"] = (out.priority != "high").astype(int)
    out["_role"] = (out.role != "reference").astype(int)
    out = out.sort_values(
        ["_tier", "wt_id", "_role", "gene", "dbd_family"],
        kind="stable",
    ).reset_index(drop=True)

    tiers = {"high": "HP", "standard": "SD"}
    names, rank = [], {"high": 0, "standard": 0}
    for row in out.itertuples(index=False):
        rank[row.priority] += 1
        gene = _safe(row.gene)[:MAX_GENE_CHARS]
        names.append(f"{tiers[row.priority]}{rank[row.priority]:04d}_{gene}_{row.sha256[:8]}")
    out.insert(0, "name", names)
    return out.loc[:, list(MANIFEST_COLUMNS)]


COMPLEMENT = str.maketrans("ACGT", "TGCA")


def revcomp(seq: str) -> str:
    return seq.translate(COMPLEMENT)[::-1]


def build_dna(km: matrix.KmerMatrix) -> pd.DataFrame:
    """The DNA manifest: one row per 8-mer, in the matrix's column order."""
    kmers = [str(k) for k in km.kmers]
    label = km.label
    escore = km.escore
    out = pd.DataFrame(
        {
            "matrix_index": range(len(kmers)),
            "dna_seq": kmers,
            "revcomp": [revcomp(k) for k in kmers],
            "n_domains_pos": (label == 1).sum(axis=0),
            "n_domains_neg": (label == 0).sum(axis=0),
            "n_domains_nocall": (label == -1).sum(axis=0),
            "max_escore": escore.max(axis=0).round(5),
        }
    )
    out.insert(0, "name", [f"K{i + 1:05d}_{k}" for i, k in enumerate(kmers)])
    out.insert(4, "palindrome", out.dna_seq == out.revcomp)
    return out.loc[:, list(DNA_MANIFEST_COLUMNS)]


def dna_header(row: pd.Series) -> str:
    """The DNA description line. `rc=` is the strand the dataset does not store separately."""
    return (
        f">{row.name} idx={row.matrix_index} rc={row.revcomp} "
        f"palindrome={'yes' if row.palindrome else 'no'} len={len(row.dna_seq)} "
        f"bound_by={row.n_domains_pos} max_escore={row.max_escore:.5f}"
    )


def header(row: pd.Series) -> str:
    """The FASTA description line: the name, then space-free `key=value` metadata."""
    fields = [
        f"priority={row.priority}",
        f"gene={_safe(row.gene)}",
        f"family={_safe(row.dbd_family)}",
        f"source={row.source_dataset}",
        f"cluster={row.wt_id}",
        f"role={row.role}",
        f"cluster_size={row.cluster_size}",
        f"len={row.length}",
        f"sha256={row.sha256[:16]}",
    ]
    if row.role == "variant":
        fields.insert(6, f"edits={row.n_edits_from_reference}")
        fields.insert(7, f"edit_positions={row.edit_positions}")
    return f">{row.name} " + " ".join(fields)


def write_fasta(manifest: pd.DataFrame, path: Path, header_fn, seq_column: str) -> None:
    lines = []
    for row in manifest.itertuples(index=False):
        lines.append(header_fn(row))
        seq = getattr(row, seq_column)
        lines.extend(seq[i : i + WRAP] for i in range(0, len(seq), WRAP))
    path.write_text("\n".join(lines) + "\n")


def export_protein(out_dir: Path) -> None:
    manifest = build(corpus.domains(), clusters.load())
    fasta = out_dir / "dbd_domains.fasta"
    tsv = out_dir / "dbd_domains_manifest.tsv"
    write_fasta(manifest, fasta, header, "dbd_seq")
    manifest.to_csv(tsv, sep="\t", index=False)

    high = manifest[manifest.priority == "high"]
    print(f"{len(manifest)} protein sequences -> {fasta}")
    print(f"  high priority : {len(high)} in {high.wt_id.nunique()} clusters")
    print(f"  standard      : {len(manifest) - len(high)}")
    print(f"  lengths       : {manifest.length.min()}-{manifest.length.max()} aa")
    print(f"manifest -> {tsv}")


def export_dna(out_dir: Path) -> None:
    manifest = build_dna(matrix.KmerMatrix.load())
    fasta = out_dir / "dna_8mers.fasta"
    tsv = out_dir / "dna_8mers_manifest.tsv"
    write_fasta(manifest, fasta, dna_header, "dna_seq")
    manifest.to_csv(tsv, sep="\t", index=False)

    print(f"{len(manifest)} 8-mers -> {fasta}")
    print(f"  palindromes   : {int(manifest.palindrome.sum())}")
    print(f"  bound by >= 1 : {int((manifest.n_domains_pos > 0).sum())}")
    print(f"manifest -> {tsv}")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out-dir", type=Path, default=OUT_DIR)
    parser.add_argument("--what", choices=("protein", "dna", "all"), default="all")
    args = parser.parse_args()

    args.out_dir.mkdir(parents=True, exist_ok=True)
    if args.what in ("protein", "all"):
        export_protein(args.out_dir)
    if args.what in ("dna", "all"):
        export_dna(args.out_dir)


if __name__ == "__main__":
    main()
