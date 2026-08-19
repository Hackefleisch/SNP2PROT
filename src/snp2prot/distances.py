"""All-vs-all distance between the corpus's domains, built once and read by three consumers.

`ML_PLAN.md` §8.1: this matrix is **shared infrastructure, not baseline overhead**. Three
things need exactly these numbers and no others —

1. the **nearest-neighbour baseline**, which copies the profile of the most identical
   training domain (`snp2prot.baselines.nn_lookup`);
2. the **`S2` split grouping**, whose connected components are this matrix thresholded
   (`snp2prot.splits`, `T27`);
3. the **nearest-neighbour-identity histogram** that explains the degradation curve
   (`ML_PLAN.md` §5.2).

894,453 pairs at 0.17 ms each is 2.5 minutes on one core, so it is computed in parallel and
cached rather than recomputed per experiment.

## What is stored, and why not just identity

`n_edits` and `n_aligned`, both integers, from `snp2prot.align.edit_profile` — the same
alignment the clustering uses: BLOSUM62, affine gaps, **free terminal gaps**. Identity is
`1 - n_edits / n_aligned` and overlap is `n_aligned / min(len_i, len_j)`, both derived on
demand. Storing the two counts rather than the two ratios keeps the artifact exact and lets a
caller change its mind about the guard without a rebuild.

**The overlap guard is not optional.** With terminal gaps free, the aligner may discard an
arbitrary prefix of one domain and suffix of the other at no cost and then report a handful of
edits over the sliver that survives: two unrelated Myb domains align on six residues and come
out 3 edits apart where an honest alignment says 59. That produced 31 false cluster
memberships before it was caught (`T25`, `docs/DECISIONS.md`). Identity computed over such a
window is not a small error, it is a different quantity — and under a family holdout, where
every candidate neighbour is cross-family by construction, it is the *typical* case rather
than the exception. Every consumer here applies `cluster.min_overlap`.

## Symmetry

`n_edits` and `n_aligned` are symmetric — the alignment score is, and both counts are read off
the alignment rather than off either sequence's frame. Only `mut_positions` is asymmetric, and
that is not stored here. The upper triangle is computed and mirrored.
"""

from __future__ import annotations

from dataclasses import dataclass
from multiprocessing import Pool
from pathlib import Path

import numpy as np

from snp2prot import align
from snp2prot.config import DISTANCE_MATRIX

#: Domains with no aligned residue at all have no identity to report, and no caller should
#: silently treat one as "0% identical" — that is a different statement from "not comparable".
UNCOMPARABLE = np.float32(np.nan)

_SEQUENCES: list[str] = []


@dataclass(frozen=True)
class DomainDistances:
    """The pairwise alignment counts, plus the domain order they are indexed by."""

    domains: np.ndarray  # (n,) dbd_seq, in corpus.domains() order
    n_edits: np.ndarray  # (n, n) int16, symmetric, 0 on the diagonal
    n_aligned: np.ndarray  # (n, n) int32, symmetric, len(domain) on the diagonal

    @property
    def n_domains(self) -> int:
        return len(self.domains)

    def identity(self) -> np.ndarray:
        """`1 - n_edits / n_aligned`: percent identity over the aligned region."""
        with np.errstate(invalid="ignore", divide="ignore"):
            out = 1.0 - self.n_edits / np.where(self.n_aligned > 0, self.n_aligned, np.nan)
        return out.astype(np.float32)

    def overlap(self) -> np.ndarray:
        """`n_aligned / min(len_i, len_j)`: how much of the shorter domain actually aligned."""
        lengths = np.array([len(s) for s in self.domains], dtype=np.float32)
        shorter = np.minimum(lengths[:, None], lengths[None, :])
        return (self.n_aligned / shorter).astype(np.float32)

    def comparable(self, min_overlap: float) -> np.ndarray:
        """Boolean mask of pairs whose alignment is trustworthy enough to rank on."""
        return self.overlap() >= min_overlap

    def save(self, path: str | Path | None = None) -> Path:
        p = Path(path) if path else DISTANCE_MATRIX
        p.parent.mkdir(parents=True, exist_ok=True)
        np.savez(p, domains=self.domains, n_edits=self.n_edits, n_aligned=self.n_aligned)
        return p

    @classmethod
    def load(cls, path: str | Path | None = None) -> DomainDistances:
        p = Path(path) if path else DISTANCE_MATRIX
        if not p.exists():
            raise FileNotFoundError(f"no distance matrix at {p}\nrun scripts/build_distances.py")
        with np.load(p, allow_pickle=False) as z:
            return cls(z["domains"], z["n_edits"], z["n_aligned"])


def _init(sequences: list[str]) -> None:
    global _SEQUENCES
    _SEQUENCES = sequences


def _row(i: int) -> tuple[int, np.ndarray, np.ndarray]:
    """One row of the upper triangle: domain `i` against every domain after it."""
    reference = _SEQUENCES[i]
    rest = _SEQUENCES[i + 1 :]
    edits = np.zeros(len(rest), dtype=np.int16)
    aligned = np.zeros(len(rest), dtype=np.int32)
    for k, other in enumerate(rest):
        profile = align.edit_profile(reference, other)
        edits[k] = profile.n_edits
        aligned[k] = round(profile.overlap * min(len(reference), len(other)))
    return i, edits, aligned


def build(domains, processes: int | None = None, progress=None) -> DomainDistances:
    """Align every pair of `domains` (an iterable of `dbd_seq`) and return the matrices.

    `progress` is called with the number of rows finished, for a script that wants to say
    something during the two minutes of work.
    """
    sequences = [str(s) for s in domains]
    n = len(sequences)
    edits = np.zeros((n, n), dtype=np.int16)
    aligned = np.zeros((n, n), dtype=np.int32)
    for i, seq in enumerate(sequences):
        aligned[i, i] = len(seq)

    # Row i costs n - i alignments, so the longest rows are handed out first; otherwise the
    # last worker starts row 0 while the others have nothing left to do.
    with Pool(processes=processes, initializer=_init, initargs=(sequences,)) as pool:
        for done, (i, row_edits, row_aligned) in enumerate(
            pool.imap_unordered(_row, range(n - 1), chunksize=4), start=1
        ):
            edits[i, i + 1 :] = row_edits
            edits[i + 1 :, i] = row_edits
            aligned[i, i + 1 :] = row_aligned
            aligned[i + 1 :, i] = row_aligned
            if progress is not None:
                progress(done)
    # Fixed-width unicode, not object: an object array in an `.npz` needs `allow_pickle`
    # on load, and a cached artifact should never require that to be read.
    return DomainDistances(np.asarray(sequences, dtype=np.str_), edits, aligned)


def connected_components(
    distances: DomainDistances,
    min_identity: float,
    min_overlap: float,
    families: np.ndarray | None = None,
) -> np.ndarray:
    """Component id per domain, joining any two domains at or above `min_identity`.

    The `S2` grouping (`ML_PLAN.md` §5, `T27`, `D5`). Single linkage is deliberate and is the
    opposite of the choice `snp2prot.clusters` makes for *clustering*: for **splitting**,
    over-grouping removes more from training than strictly necessary but never leaks, whereas
    in a *stored* cluster it would be a false claim about relatedness. That asymmetry is why
    chaining is acceptable here and was rejected there.

    **The floor is identity, not an edit count, and it is set in `configs/experiment.yaml`.**
    Grouping at `<= 5` edits — about 94% identity on a 77 aa domain — was measured on
    2026-08-19 to remove the near-twins (24% of held-out domains had a >= 90% identical
    neighbour under `S1`, 7% under that grouping) and to change the nearest-neighbour
    baseline's AUPR by 0.018, because copying a 70-90% identical neighbour already scores
    0.927. The floor moved to **0.5** on the owner's decision (`docs/DECISIONS.md`): the split
    is meant to be hard, and 0.5 is well below the 60-70% band Weirauch et al. 2014 set for
    transferring a motif by identity.

    It chains, and that is understood rather than overlooked: at 0.5 the largest component is
    273 domains whose *median* internal identity is 0.377, so it is a chain of near-relatives
    and not a clique. Since chaining only ever removes more from training, the cost is
    statistical power, not validity.

    `families` restricts edges to within one Pfam family, which is how the clustering was
    built and what makes `P1`'s family holdout consistent with this grouping.
    """
    n = distances.n_domains
    identity = distances.identity()
    close = (identity >= min_identity) & np.isfinite(identity) & distances.comparable(min_overlap)
    if families is not None:
        close &= np.asarray(families)[:, None] == np.asarray(families)[None, :]

    parent = np.arange(n)

    def find(x: int) -> int:
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x

    for i, j in zip(*np.nonzero(np.triu(close, k=1)), strict=True):
        a, b = find(int(i)), find(int(j))
        if a != b:
            parent[max(a, b)] = min(a, b)

    roots = np.array([find(i) for i in range(n)])
    # Relabel to 0..k-1 in order of first appearance, so the ids do not depend on n.
    _, labels = np.unique(roots, return_inverse=True)
    return labels.astype(np.int32)
