"""Nearest-neighbour lookup: copy the most identical training domain's E-score profile.

**The bar, and not a strawman.** Transferring a motif to an uncharacterised TF by DBD
sequence identity is how CIS-BP does inference, and Weirauch et al. 2014 established the
per-family identity thresholds for exactly that — a paper that is one of this dataset's own 19
sources. It is the standard method in the field and the first question a biologist asks of any
model here: *isn't this just copying the most similar protein's motif?* Any model that does not
beat it under `S2` has learned nothing transferable (`CLAUDE.md`; `docs/ML_PLAN.md` §8.1).

**It copies the neighbour's binary calls, not its E-scores.** It used to copy the continuous
profile, and that was wrong twice over.

*It is not a like-for-like comparison.* Measured 2026-08-26, same neighbour and same overlap
guard, the continuous form scored **0.786 against 0.472** on `S1/fold-0` and **0.391 against
0.145** on `S2/fold-0`. Between 29% and 76% of the old baseline's AUPR came from the ordering
*within* the copied profile — information the model is never given, since it trains on labels.
Scoring a binary-trained model against a continuous-profile lookup compares inputs, not methods.

*And the ordering is not a quantity the assay reports.* A universal-PBM E-score is a
rank-enrichment statistic on a fixed [-0.5, 0.5] scale, read by the field at a cutoff and stored
here at 0.45 / 0.35 (`configs/thresholds.yaml`). It measures statistical enrichment against
background, **not graded affinity**, so ranking 8-mers by it reads a precision into the number
that is not there. That is why the dataset stores a binary label, and why the whole PBM corpus —
training, baseline and evaluation alike — is binary.

**Selection is percent identity over the aligned domain, under the overlap guard.** Identity
is `1 - n_edits / n_aligned` from `snp2prot.distances`, and a candidate is only eligible if the
alignment covered at least `cluster.min_overlap` of the shorter domain. The guard is mandatory
rather than tidy: with free terminal gaps the aligner can park almost all of both sequences in
gaps that cost nothing and report a handful of edits over the sliver that survives, which is
how 31 false cluster memberships were formed before it was caught (`T25`). Under `P1` every
candidate is cross-family by construction, so that degenerate case is the *typical* one there,
not an edge case. Measured on this corpus every domain has at least 128 guarded candidates,
so the guard never empties the pool — `MEAN_PROFILE` below is defined for completeness and
does not fire.

**The training pool is the caller's to choose.** `fit_predict` copies from exactly the rows it
is given, so a caller that has dropped the `no_evidence` records (`label_health.usable`) gets a
baseline that never copies an unverifiable silent profile.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from snp2prot.data.matrix import KmerMatrix
from snp2prot.distances import DomainDistances

#: What a held-out domain is predicted when no training domain clears the overlap guard: the
#: mean training profile, i.e. the corpus's average 8-mer preference. The lookup has no
#: neighbour to copy and says so; falling back to the best *unguarded* candidate instead would
#: copy the profile of whichever unrelated domain the aligner mangled most favourably.
MEAN_PROFILE = "mean"


def _calls(matrix: KmerMatrix) -> np.ndarray:
    """The neighbour's *positive calls* as 1.0 and everything else as 0.0.

    Its own no-call band included: "not called a binder" is the prediction being transferred,
    and the held-out domain's gray band is masked by the metric rather than by this.
    """
    return (matrix.label == 1).astype(np.float32)


@dataclass(frozen=True)
class Prediction:
    """Predicted profiles for the held-out domains, and which neighbour produced each."""

    #: (n_test, n_kmers) float32 — the copied (or averaged) E-score profile.
    profile: np.ndarray
    #: One row per held-out domain: the neighbour chosen and how close it was.
    neighbours: pd.DataFrame

    @property
    def n_without_neighbour(self) -> int:
        return int((self.neighbours.neighbour == MEAN_PROFILE).sum())


def choose(
    distances: DomainDistances,
    test_rows: np.ndarray,
    train_rows: np.ndarray,
    min_overlap: float,
    k: int = 1,
) -> tuple[np.ndarray, np.ndarray]:
    """The `k` most identical guarded training rows for each test row, best first.

    Returns positions into `train_rows` and their identities; a test row with fewer than `k`
    eligible candidates is padded with `-1` and `nan`, which `fit_predict` reads as "use what
    there is".

    Ties in identity break towards the lower training row, which is the alphabetically earlier
    `dbd_seq` — arbitrary, but fixed, so a rerun reproduces the same neighbour.

    **A stable sort, not `argpartition`.** `argpartition` is faster and does not promise an order
    among equal values, so the rule above was simply untrue: measured across the 19 folds, 9.4%
    of held-out domains have a tied best identity and in **51% of those the pick was not the
    lowest training row**. It moved the baseline by up to 0.008 AUPR and would move again on a
    different numpy version or a reordered training pool. A stable sort over a 1,338-wide row
    costs microseconds, which is not a price worth paying an unstated rule for.
    """
    block = distances.identity()[np.ix_(test_rows, train_rows)]
    guarded = distances.comparable(min_overlap)[np.ix_(test_rows, train_rows)]
    eligible = np.where(guarded & np.isfinite(block), block, -np.inf)

    k = min(k, eligible.shape[1])
    top = np.argsort(-eligible, axis=1, kind="stable")[:, :k]
    values = np.take_along_axis(eligible, top, axis=1)

    empty = ~np.isfinite(values)
    top = np.where(empty, -1, top)
    return top, np.where(empty, np.nan, values)


def fit_predict(
    matrix: KmerMatrix,
    distances: DomainDistances,
    test_rows: np.ndarray,
    train_rows: np.ndarray,
    min_overlap: float,
    k: int = 1,
) -> Prediction:
    """Predict every held-out domain's 8-mer profile by copying its nearest training domain.

    `k = 1` is the primary form — the pure lookup table, and the thing to beat. Above 1 the
    profile is the identity-weighted mean of the top `k` neighbours, which is a slightly
    stronger bar and one extra line (`docs/ML_PLAN.md` §8.1).

    """
    test_rows = np.asarray(test_rows, dtype=np.int64)
    train_rows = np.asarray(train_rows, dtype=np.int64)
    if not len(train_rows):
        raise ValueError("empty training pool: there is nothing to copy from")

    picks, identities = choose(distances, test_rows, train_rows, min_overlap, k)
    escore = _calls(matrix)
    fallback = escore[train_rows].mean(axis=0)
    overlap = distances.overlap()

    profile = np.empty((len(test_rows), escore.shape[1]), dtype=np.float32)
    records = []
    for i, row in enumerate(test_rows):
        usable = picks[i] >= 0
        if not usable.any():
            profile[i] = fallback
            records.append((matrix.domains[row], MEAN_PROFILE, np.nan, -1, np.nan))
            continue
        chosen = train_rows[picks[i][usable]]
        weights = identities[i][usable].astype(np.float64)
        # Identity is bounded below by 0 in practice but the formula does not guarantee it;
        # a negative weight would be nonsense, so the floor is explicit.
        weights = np.clip(weights, 0.0, None)
        total = weights.sum()
        weights = weights / total if total > 0 else np.full(len(chosen), 1 / len(chosen))
        profile[i] = np.average(escore[chosen], axis=0, weights=weights).astype(np.float32)
        best = int(chosen[0])
        records.append(
            (
                matrix.domains[row],
                matrix.domains[best],
                float(identities[i][0]),
                int(distances.n_edits[row, best]),
                float(overlap[row, best]),
            )
        )

    neighbours = pd.DataFrame(
        records, columns=["domain", "neighbour", "identity", "n_edits", "overlap"]
    )
    return Prediction(profile, neighbours)
