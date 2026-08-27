"""The merged table as a dense domain x 8-mer matrix.

The merged table is a **complete** matrix — 1,338 domains x 32,896 8-mers, every cell
present, 44,014,848 = 1,338 x 32,896 exactly (`TODO.md` `N1`). It is *stored* as 44 million
rows in a 974 MB Parquet because that is the schema every parser writes and every validator
checks. Nothing downstream of the merge wants it in that shape: a split holds out domains, a
metric ranks one domain's 32,896 8-mers, and the nearest-neighbour baseline copies one row of
it. All three are one array indexing away once the table is a matrix, and a group-by over 44
million rows otherwise.

Two arrays, both resident:

| array | shape | dtype | size |
|---|---|---|---:|
| `escore` | 1,338 x 32,896 | `float32` | 168 MiB |
| `label` | 1,338 x 32,896 | `int8` | 42 MiB |

`float32` is not lossy here in any way that matters: E-scores live on a fixed [-0.5, 0.5]
scale with about four significant digits, and the cutoff comparison that produced `label` was
made at parse time in double precision and is stored, not recomputed.

**The label keeps its three values.** `1` binds, `0` does not, `-1` is the no-call band
between the two cutoffs — 1.7% of cells. `-1` is *absent evidence, not a negative*, so it is
carried in the matrix and excluded by a mask wherever a metric or a loss is computed
(`ML_PLAN.md` §3.1: two masks, not filters).

**Row order is `corpus.domains()` order and column order is the sorted 8-mer vocabulary.**
Neither is stored anywhere else, so both are written into the artifact — a matrix built before
a rebuild that changed the domain set would otherwise be silently misaligned with a distance
matrix built after it.

**The order is not checked here.** `load` cannot compare itself against the corpus without
importing `snp2prot.corpus`, and a data container should not depend on the tables that describe
it. The check is `corpus.require_aligned`, called wherever two of these artifacts are indexed
together: `snp2prot.training.Trainer` and `run_fold`, and each script after it loads them.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pyarrow.parquet as pq

from snp2prot.config import KMER_MATRIX, MERGED_TABLE

#: Cells with no measurement. Every cell of the merged table is present, so this is only
#: reachable if `build` is given a table that is not the complete matrix — which it checks.
MISSING = np.float32(np.nan)


@dataclass(frozen=True)
class KmerMatrix:
    """E-scores and labels for every (domain, 8-mer) pair, plus the two index vectors."""

    domains: np.ndarray  # (n_domains,) dbd_seq, in corpus.domains() order
    kmers: np.ndarray  # (n_kmers,) the 8-mer vocabulary, sorted
    escore: np.ndarray  # (n_domains, n_kmers) float32
    label: np.ndarray  # (n_domains, n_kmers) int8

    @property
    def shape(self) -> tuple[int, int]:
        return self.escore.shape

    def index(self) -> dict[str, int]:
        """`dbd_seq` -> row, for callers holding domain metadata rather than positions."""
        return {seq: i for i, seq in enumerate(self.domains)}

    def rows_for(self, domains) -> np.ndarray:
        """Row positions for an iterable of `dbd_seq`, in the order given."""
        lookup = self.index()
        missing = [d for d in domains if d not in lookup]
        if missing:
            raise KeyError(f"{len(missing)} domains are not in the matrix, e.g. {missing[0][:30]}")
        return np.array([lookup[d] for d in domains], dtype=np.int64)

    def scored(self) -> np.ndarray:
        """Boolean mask of cells a metric may use: everything outside the no-call band."""
        return self.label != -1

    def save(self, path: str | Path | None = None) -> Path:
        p = Path(path) if path else KMER_MATRIX
        p.parent.mkdir(parents=True, exist_ok=True)
        np.savez(
            p,
            domains=self.domains,
            kmers=self.kmers,
            escore=self.escore,
            label=self.label,
        )
        return p

    @classmethod
    def load(cls, path: str | Path | None = None) -> KmerMatrix:
        p = Path(path) if path else KMER_MATRIX
        if not p.exists():
            raise FileNotFoundError(f"no 8-mer matrix at {p}\nrun scripts/build_matrix.py")
        with np.load(p, allow_pickle=False) as z:
            return cls(z["domains"], z["kmers"], z["escore"], z["label"])


def build(domains: np.ndarray, table: str | Path | None = None) -> KmerMatrix:
    """Read the merged table row-group by row-group and scatter it into the matrix.

    `domains` fixes the row order; the 8-mer vocabulary is read from the table rather than
    generated, because which of a reverse-complement pair represents it is a fact about how
    the arrays were designed and not one to re-derive here.

    **The vocabulary comes from row group 0 alone**, which holds about a million rows against a
    32,896-wide vocabulary, so it covers it many times over. It is safe rather than lucky: an
    8-mer outside the derived vocabulary makes `_codes` raise, and a vocabulary that came out
    short leaves cells unwritten, which the `-2` sweep below catches. Both failures are loud, so
    this does not pay for a second pass over 44 million rows to be sure.

    Streaming by row group keeps peak memory at one group (about 1 million rows) rather than
    the 44 million the file holds, and the scatter is done through each group's own
    dictionary — 32 domains and 32,896 8-mers hashed per group instead of 2 million strings.
    """
    path = Path(table) if table else MERGED_TABLE
    if not path.exists():
        raise FileNotFoundError(f"no merged table at {path}\nrun scripts/build_merged.py")

    reader = pq.ParquetFile(path)
    first = reader.read_row_group(0, columns=["dna_seq"])["dna_seq"].unique()
    kmers = np.sort(np.asarray(first.to_pylist(), dtype=np.str_))
    domain_row = {seq: i for i, seq in enumerate(domains)}
    kmer_col = {kmer: i for i, kmer in enumerate(kmers)}

    n_domains, n_kmers = len(domains), len(kmers)
    escore = np.full((n_domains, n_kmers), MISSING, dtype=np.float32)
    label = np.full((n_domains, n_kmers), -2, dtype=np.int8)  # -2 = never written

    for group in range(reader.num_row_groups):
        batch = reader.read_row_group(group, columns=["dbd_seq", "dna_seq", "raw_score", "label"])
        rows = _codes(batch["dbd_seq"], domain_row, "domain")
        cols = _codes(batch["dna_seq"], kmer_col, "8-mer")
        flat = rows * n_kmers + cols
        escore.reshape(-1)[flat] = batch["raw_score"].to_numpy(zero_copy_only=False)
        label.reshape(-1)[flat] = batch["label"].to_numpy(zero_copy_only=False)

    unwritten = int((label == -2).sum())
    if unwritten:
        raise ValueError(
            f"{unwritten} of {n_domains * n_kmers} cells were not in the table — the merged "
            "table is not the complete matrix it is asserted to be (TODO.md N1)"
        )
    # Fixed-width unicode, not object: an object array in an `.npz` needs `allow_pickle`
    # on load, and a cached artifact should never require that to be read.
    return KmerMatrix(np.asarray(domains, dtype=np.str_), kmers, escore, label)


def _codes(column, lookup: dict[str, int], what: str) -> np.ndarray:
    """Map a chunked string column onto positions, hashing each chunk's dictionary once."""
    encoded = column.dictionary_encode()
    out = np.empty(len(column), dtype=np.int64)
    at = 0
    for chunk in encoded.chunks:
        values = chunk.dictionary.to_pylist()
        try:
            mapping = np.array([lookup[v] for v in values], dtype=np.int64)
        except KeyError as exc:
            raise KeyError(f"{what} in the merged table but not in the index: {exc}") from exc
        indices = chunk.indices.to_numpy(zero_copy_only=False)
        out[at : at + len(chunk)] = mapping[indices]
        at += len(chunk)
    return out
