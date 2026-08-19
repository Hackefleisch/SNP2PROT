"""Pooled protein-language-model embeddings: one fixed-width vector per domain.

**A precomputed lookup, not part of the model** (`docs/ML_PLAN.md` §3.1). 1,338 domains through
a PLM is a one-off inference cached to disk; the language model never runs inside the training
loop, which is what makes the experiment grid an afternoon rather than a week.

## One vector per domain, mean-pooled

Decided 2026-08-19 (§3.1), and the reason is **cross-arm comparability** rather than
convenience. A structural embedder produces one vector per structure, so if the sequence arms
were per-residue and the structure arms per-structure, `A1`/`A4` and `A2`/`A3` would feed
differently shaped towers — and the `A1` → `A2` comparison that carries claim **C2** would
confound modality with architecture. One vector per domain everywhere keeps the four arms
interchangeable, and it is also what TransBind does, which makes §8.2's baseline exactly
controlled rather than nearly so.

**Pooling is over the residues of the domain only** — the model's start-of-sequence and
end-of-sequence tokens are dropped. They carry sequence-level summary information that would be
identical in kind for every domain and would dilute the per-residue signal the pooling is
supposed to average.

## The cost lands on C1, so it is measured before it is trusted

Mean pooling dilutes a single-residue change: a variant differs from its reference by 1-5
residues out of a median 77. The naive bound is one part in 77, but that assumes only the
mutated residue's representation changes, and a PLM is **contextual** — a substitution shifts
its neighbours' representations too, which can only help. The true attenuation is empirical,
which is what `scripts/check_pooling.py` measures before any training run.

## The arms

| arm | model | width |
|---|---|---:|
| `A1` | ESM-2 650M (`esm2_t33_650M_UR50D`) — the most standard PLM available | 1280 |
| `A4` | ESM-DBP — ESM-2 domain-adapted to DNA-binding proteins (Zeng et al. 2024) | 1280 |

Equal width is not a coincidence and is not required either: every tower ends in a learned
projection to a common `D` (§4.2). It does mean `A1` and `A4` enter that projection on equal
terms, so the `A1` → `A4` delta isolates the pretraining corpus and nothing else.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np

from snp2prot.config import embedding_table

#: Arm -> the checkpoint it is built from. `A2`/`A3` are structural and are not built here.
ARMS = {
    "A1": "esm2_t33_650M_UR50D",
    "A4": "esm_dbp",
}


@dataclass(frozen=True)
class DomainEmbeddings:
    """One pooled vector per domain, indexed by the same order `corpus.domains()` returns."""

    arm: str
    model: str
    domains: np.ndarray  # (n,) dbd_seq
    vectors: np.ndarray  # (n, width) float32

    @property
    def width(self) -> int:
        return self.vectors.shape[1]

    def rows_for(self, domains) -> np.ndarray:
        lookup = {seq: i for i, seq in enumerate(self.domains)}
        missing = [d for d in domains if d not in lookup]
        if missing:
            raise KeyError(f"{len(missing)} domains have no embedding, e.g. {missing[0][:30]}")
        return np.array([lookup[d] for d in domains], dtype=np.int64)

    def save(self, path: str | Path | None = None) -> Path:
        p = Path(path) if path else embedding_table(self.arm)
        p.parent.mkdir(parents=True, exist_ok=True)
        np.savez(
            p,
            arm=np.asarray(self.arm),
            model=np.asarray(self.model),
            domains=self.domains,
            vectors=self.vectors,
        )
        return p

    @classmethod
    def load(cls, arm: str, path: str | Path | None = None) -> DomainEmbeddings:
        p = Path(path) if path else embedding_table(arm)
        if not p.exists():
            raise FileNotFoundError(
                f"no {arm} embeddings at {p}\nrun scripts/build_embeddings.py --arm {arm}"
            )
        with np.load(p, allow_pickle=False) as z:
            return cls(str(z["arm"]), str(z["model"]), z["domains"], z["vectors"])


def cosine_distance(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Row-wise cosine distance between two equally shaped stacks of vectors.

    Cosine rather than Euclidean because PLM embedding norms vary with sequence length, and a
    variant is by construction almost the same length as its reference — a norm difference
    would report as a distance without meaning one.
    """
    a = a.astype(np.float64)
    b = b.astype(np.float64)
    denominator = np.linalg.norm(a, axis=-1) * np.linalg.norm(b, axis=-1)
    with np.errstate(invalid="ignore", divide="ignore"):
        return 1.0 - (a * b).sum(axis=-1) / np.where(denominator > 0, denominator, np.nan)


def pairwise_cosine_distance(vectors: np.ndarray) -> np.ndarray:
    """All-vs-all cosine distance for a stack of vectors."""
    normalised = vectors.astype(np.float64)
    norms = np.linalg.norm(normalised, axis=1, keepdims=True)
    normalised = normalised / np.where(norms > 0, norms, np.nan)
    return 1.0 - normalised @ normalised.T
