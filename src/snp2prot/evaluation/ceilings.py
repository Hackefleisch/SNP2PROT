"""What the *representation* allows, independent of how well the model is trained.

`ML_RESULTS.md` §4 rests on three measurements, and this module is where they live so that they
regenerate like everything else in the project rather than surviving as numbers someone once
computed in a shell.

The question they answer together is: **when the two-tower model loses to nearest-neighbour
lookup, which link is at fault?** Three links can be measured separately.

| link | measured by | what a low number would mean |
|---|---|---|
| the shared space is too narrow | `rank_ceiling` | no bilinear model this wide can fit it |
| the protein tower is too weak | `ridge_probe` | a better tower would help |
| the tower needs to be non-linear | `kernel_probe` | non-linearity is the missing piece |

## Two numbers per probe, and the difference matters

Each probe is reported twice:

- **ceiling** — the regularisation strength chosen by looking at the *test* fold. This is an
  oracle and is not an achievable score; it is an **upper bound on the model class**, which is
  exactly what is wanted when the purpose is to rule a class out.
- **honest** — the same strength chosen on the fold's own validation slice, which is what a real
  run would get.

Quoting only the ceiling would overstate the alternatives; quoting only the honest number would
understate them, and a bound that is too tight cannot rule anything out. Both are cheap.

## The rank ceiling is an oracle too, and more so

`rank_ceiling` truncates the SVD of the E-score matrix itself, so it sees every held-out label.
It answers one narrow question — *could a rank-`D` bilinear score matrix represent this data at
all* — and nothing else. It is a ceiling on the architecture, never a target.
"""

from __future__ import annotations

import numpy as np

from snp2prot.data.matrix import KmerMatrix
from snp2prot.evaluation import metrics

#: Regularisation strengths swept by the linear probe. Wide, because the right scale depends on
#: the embedding's norm and differs between arms.
RIDGE_LAMBDAS = (1e0, 1e1, 1e2, 1e3, 1e4)
#: The kernel probe solves in sample space, where the scale is completely different.
KERNEL_LAMBDAS = (1e-4, 1e-3, 1e-2, 1e-1)
#: Bandwidths as multiples of the median squared pairwise distance — the usual heuristic.
KERNEL_BANDWIDTHS = (0.5, 1.0, 2.0)


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
    """Best macro AUPR from a rank-`r` approximation of the E-score matrix.

    An **oracle**: the factors are fitted to the very matrix being scored, so this bounds what
    any bilinear model of that width could represent and says nothing about generalisation.
    Computed through the domain-side Gram matrix, which is 1,338 x 1,338 rather than the full
    SVD of a 1,338 x 32,896.
    """
    scores = matrix.escore.astype(np.float64)
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


def _solve_and_score(matrix, design, targets, train, evaluate, lambdas, identity_size):
    """Fit a ridge for every lambda and score each on `evaluate`. Returns {lambda: aupr}."""
    gram = design[train].T @ design[train]
    cross = design[train].T @ targets[train]
    eye = np.eye(identity_size)
    out = {}
    for lam in lambdas:
        weights = np.linalg.solve(gram + lam * eye, cross)
        out[lam] = macro_aupr(matrix, evaluate, design[evaluate] @ weights)
    return out


def ridge_probe(
    matrix: KmerMatrix,
    embeddings: np.ndarray,
    train: np.ndarray,
    validation: np.ndarray,
    test: np.ndarray,
    lambdas=RIDGE_LAMBDAS,
) -> dict[str, float]:
    """The best **linear-in-embedding** predictor of the whole E-score profile.

    Unconstrained in output — it maps the protein vector straight onto all 32,896 8-mers with no
    shared space and no DNA tower — so it strictly upper-bounds any model whose protein side is a
    linear map of the same vector, the two-tower model included.
    """
    targets = matrix.escore.astype(np.float64)
    design = np.hstack([embeddings.astype(np.float64), np.ones((len(embeddings), 1))])

    on_test = _solve_and_score(matrix, design, targets, train, test, lambdas, design.shape[1])
    best_by_test = max(on_test, key=on_test.get)
    out = {"ceiling": on_test[best_by_test], "ceiling_lambda": best_by_test}

    if len(validation):
        on_validation = _solve_and_score(
            matrix, design, targets, train, validation, lambdas, design.shape[1]
        )
        chosen = max(on_validation, key=on_validation.get)
        out |= {"honest": on_test[chosen], "honest_lambda": chosen}
    return out


def kernel_probe(
    matrix: KmerMatrix,
    embeddings: np.ndarray,
    train: np.ndarray,
    validation: np.ndarray,
    test: np.ndarray,
    lambdas=KERNEL_LAMBDAS,
    bandwidths=KERNEL_BANDWIDTHS,
) -> dict[str, float]:
    """The same, but non-linear: RBF kernel ridge, solved in sample space.

    The comparison against `ridge_probe` is the whole point — it isolates whether *non-linearity*
    in the protein map is the missing ingredient, separately from the tower's width or depth.
    """
    targets = matrix.escore.astype(np.float64)
    unit = embeddings.astype(np.float64)
    unit = unit / np.linalg.norm(unit, axis=1, keepdims=True)
    squared = (np.sum(unit**2, 1)[:, None] + np.sum(unit**2, 1)[None, :] - 2 * unit @ unit.T).clip(
        min=0
    )
    median = float(np.median(squared[np.triu_indices(len(unit), 1)]))

    on_test: dict[tuple, float] = {}
    on_validation: dict[tuple, float] = {}
    for bandwidth in bandwidths:
        gram = np.exp(-squared[np.ix_(train, train)] / (bandwidth * median))
        eye = np.eye(len(train))
        for lam in lambdas:
            weights = np.linalg.solve(gram + lam * eye, targets[train])
            k_test = np.exp(-squared[np.ix_(test, train)] / (bandwidth * median))
            on_test[(bandwidth, lam)] = macro_aupr(matrix, test, k_test @ weights)
            if len(validation):
                k_val = np.exp(-squared[np.ix_(validation, train)] / (bandwidth * median))
                on_validation[(bandwidth, lam)] = macro_aupr(matrix, validation, k_val @ weights)

    best = max(on_test, key=on_test.get)
    out = {
        "ceiling": on_test[best],
        "ceiling_bandwidth": best[0],
        "ceiling_lambda": best[1],
    }
    if on_validation:
        chosen = max(on_validation, key=on_validation.get)
        out |= {
            "honest": on_test[chosen],
            "honest_bandwidth": chosen[0],
            "honest_lambda": chosen[1],
        }
    return out
