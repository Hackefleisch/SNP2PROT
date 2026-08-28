"""What the *shared space* allows, independent of how well the model is trained.

## Why there is only one probe here, and what it can honestly say

This module used to hold three: a rank ceiling, a ridge probe and an RBF kernel probe, which
between them were read as a decomposition — *when the two-tower model loses to nearest-neighbour
lookup, which link is at fault?* Two of the three were removed on 2026-08-26, because the
decomposition did not follow from what they measure.

**Every probe of this kind is a *lower* bound on its model class, never an upper one.** Each
fits one estimator and reports that estimator's AUPR; the class can always contain a better
member. `ridge_probe` made the strongest form of the mistake — it minimised *squared error on
E-scores* and was then scored by *AUPR*, two objectives that demonstrably disagree (on
`S2/fold-0` the MSE-optimal ridge scored 0.4324 where the AUPR-optimal one scored 0.4736) — and
was documented as "strictly upper-bounds any model whose protein side is a linear map". It was
violated on **5 of the 8 rows of its own report**, by up to 68%: the two-tower model *is* a
linear protein map, fitted with a ranking loss instead of least squares, and on `P1` that is
worth 0.0588 against the ridge's 0.0286.

The asymmetry decides what survives:

> A lower bound licenses a conclusion only when the number comes out **high**. It can rule a
> component *out* as the constraint — "a solution this good exists, so this is not what stops
> us". It can never rule one *in*: a low number means only that this particular fit was poor.

`rank_ceiling` survives on exactly that basis, and the other two did not. Their numbers came out
*low*, which is where a lower bound says nothing.

## What `rank_ceiling` measures, and why its height is the point

The two-tower's score matrix is `(B x D) @ (D x K)`, so it is rank `D` at most, whatever the
towers do. `rank_ceiling` truncates the SVD of the binary label matrix itself and asks whether a
rank `D` matrix can order these 8-mers at all.

It is **not** the AUPR-optimal rank-`D` matrix — the truncated SVD is optimal in Frobenius norm,
which is again not the metric — so it understates what rank `D` can do. That is the safe
direction: it is a lower bound, and it comes out **high** (0.9167 at `D = 256` against the best
number any model here has reached, 0.875 on the easiest fold and 0.15-0.59 under `S2`). A
representable solution far better than anything trained therefore exists at the configured width,
and the shared space is not the binding constraint. That conclusion needs only the lower bound,
so it holds.

## What replaced the other two

The question they were meant to answer — *would a stronger or a non-linear protein tower help* —
is answered by **achievable numbers, not bounds**: train with `model.protein.hidden` set and read
the AUPR off the same folds (`TODO.md` `T36`). An ablation says "changing this helps by X", which
is both sound and actionable, where a probe could only ever say "some estimator I chose got Y".
"""

from __future__ import annotations

import numpy as np

from snp2prot.data.matrix import KmerMatrix
from snp2prot.evaluation import metrics


def macro_aupr(matrix: KmerMatrix, rows: np.ndarray, predicted: np.ndarray) -> float:
    """Mean per-protein AUPR, scored exactly as every other number in the project is."""
    scores = []
    for i, row in enumerate(np.asarray(rows)):
        scored = matrix.label[row] != -1
        binary = (matrix.label[row][scored] == 1).astype(np.int64)
        if binary.sum():
            scores.append(metrics.average_precision(binary, predicted[i][scored]))
    return float(np.mean(scores)) if scores else float("nan")


def rank_ceiling(matrix: KmerMatrix, ranks=(16, 64, 256, 512)) -> dict[int, float]:
    """Macro AUPR of a rank-`r` approximation of the **binary label** matrix.

    An **oracle and a lower bound at once**: the factors are fitted to the very matrix being
    scored, so it says nothing about generalisation; and it is the Frobenius-optimal rank-`r`
    approximation rather than the AUPR-optimal one, so the true maximum over rank `r` is higher
    still. Both work in the same direction — the number is a floor on what rank `r` can express,
    and it is the *height* of that floor that licenses the reading (see the module docstring).

    **The label matrix, not the E-score matrix** (changed 2026-08-28). Factorising the E-scores
    asked whether a rank-`r` matrix can reproduce a *continuous* quantity the project does not
    otherwise use and does not predict — a PBM E-score is a rank-enrichment statistic the field
    reads at a cutoff, so its ordering is not a target (`D7`). Factorising the calls asks the
    question that matches the task. It also raises the answer, from 0.9167 to **0.9379** at
    `r = 256`, so the reading below is unchanged and slightly stronger.

    The no-call band is filled with each 8-mer's mean call rather than treated as a negative,
    which is the same "mask, do not filter" rule the metrics follow.

    Computed through the domain-side Gram matrix, which is 1,338 x 1,338 rather than the full
    SVD of a 1,338 x 32,896.
    """
    calls = (matrix.label == 1).astype(np.float64)
    gray = matrix.label == -1
    if gray.any():
        column_mean = np.where(gray, np.nan, calls)
        calls = np.where(gray, np.nanmean(column_mean, axis=0, keepdims=True), calls)
    scores = calls
    mean = scores.mean(axis=0, keepdims=True)
    centred = scores - mean
    eigenvalues, vectors = np.linalg.eigh(centred @ centred.T)
    vectors = vectors[:, np.argsort(-eigenvalues)]

    rows = np.arange(len(matrix.domains))
    out = {}
    for rank in ranks:
        basis = vectors[:, :rank]
        out[int(rank)] = macro_aupr(matrix, rows, basis @ (basis.T @ centred) + mean)
    return out
