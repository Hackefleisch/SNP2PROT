"""The two evaluation settings of `GHT_PLAN.md` §9, and the validation carve for each.

**This arm has two holdout axes where the PBM arm has one**, and which of them a regime varies
is the whole experiment:

| regime | held out | what it tests | who else can play |
|---|---|---|---|
| `C1` | the **test chromosomes**, all TFs in training | the DNA encoder | PWM, ArChIPelago |
| `G1` | **random TFs**, 5-fold | an unseen protein | 1NN, 5NN |
| `G2` | **components at >= 0.5 identity** | an unseen neighbourhood | 1NN, 5NN |

`C1` is Setting 1 and the goal there is **parity, not victory**: a PWM fitted on that TF's own
training peaks is a strong, specific model and beating it is not the claim. `G1`/`G2` are
Setting 2, where a per-TF method cannot run at all because it was never fitted on the held-out
protein — that is the capability the talk exists to demonstrate.

## Both axes are held out in Setting 2, not just the protein

A `G1` fold tests held-out TFs on the benchmark's **test** chromosomes, never on the training
chromosomes the other TFs were trained on. Holding out only the protein would leave a held-out
TF's evaluation window at a locus some training TF was trained on — a co-bound promoter, a
repeat family — and a model that memorised "this locus is bound" would score for it. Costing
nothing and removing the question is the right trade, and it is why `Fold.digest` hashes the
chromosome set alongside the TF set (`GHT_PLAN.md` §8.2).

## Validation mirrors the test task where it can

`C1` carves whole **chromosomes** out of the training side; the TF regimes carve whole **TFs**.
Neither ever touches a test chromosome. The one regime where the mirror is imperfect is `C1` at
n = 33 proteins — validation there says "held-out loci", which is exactly what the test says, so
in fact `C1` is the clean one and the TF regimes inherit the PBM arm's usual mismatch.
"""

from __future__ import annotations

import hashlib
from dataclasses import dataclass

import numpy as np
import pandas as pd

from snp2prot import distances as distances_module
from snp2prot import thresholds
from snp2prot.splits import _balanced_groups

CHROMOSOME_REGIME = "C1"
TF_REGIMES = ("G1", "G2")


@dataclass(frozen=True)
class GHTFold:
    """One train/test division, as TF names and chromosome sets rather than row positions.

    Rows are resolved against a corpus at use time, so a fold survives a rebuild that changes
    how many windows a TF has — which row `i` is would not.
    """

    regime: str
    name: str
    test_tfs: tuple[str, ...]
    train_tfs: tuple[str, ...]
    test_chromosomes: tuple[str, ...]
    train_chromosomes: tuple[str, ...]
    #: Chromosomes carved out of `train_chromosomes` for validation; empty for a TF regime.
    validation_chromosomes: tuple[str, ...] = ()
    #: TFs carved out of `train_tfs` for validation; empty for `C1`.
    validation_tfs: tuple[str, ...] = ()
    held_out: str = ""

    @property
    def label(self) -> str:
        return f"{self.regime}/{self.name}"

    @property
    def fitting_tfs(self) -> tuple[str, ...]:
        """Training TFs minus the validation carve — what the optimiser actually sees."""
        carved = set(self.validation_tfs)
        return tuple(t for t in self.train_tfs if t not in carved)

    @property
    def fitting_chromosomes(self) -> tuple[str, ...]:
        carved = set(self.validation_chromosomes)
        return tuple(c for c in self.train_chromosomes if c not in carved)

    def digest(self) -> str:
        """A hash of **both** held-out axes.

        A TF list alone does not pin this split: the same 7 TFs evaluated on the training
        chromosomes would be a different and much easier experiment. `tracking.start_run` refuses
        a run without a digest for exactly this class of ambiguity.
        """
        payload = "|".join(
            [
                self.regime,
                self.name,
                ",".join(sorted(self.test_tfs)),
                ",".join(sorted(self.test_chromosomes)),
            ]
        )
        return hashlib.sha256(payload.encode()).hexdigest()[:12]


def _validation_chromosomes(
    chromosomes: list[str], weights: dict[str, int], fraction: float, seed: int
) -> tuple[str, ...]:
    """Whole chromosomes covering about `fraction` of the training windows, smallest-first.

    Smallest-first because the training chromosomes run from chr1 (248 Mb) to chr22, and taking
    one large one overshoots badly: a greedy fill from the small end lands within a few percent
    of the target without ever splitting a chromosome, which is the thing that must not happen.
    """
    rng = np.random.default_rng(seed)
    order = rng.permutation(len(chromosomes))
    ordered = [chromosomes[i] for i in order]
    ordered.sort(key=lambda c: weights.get(c, 0))
    target = fraction * sum(weights.get(c, 0) for c in chromosomes)
    taken: list[str] = []
    total = 0
    for chrom in ordered:
        if total >= target:
            break
        taken.append(chrom)
        total += weights.get(chrom, 0)
    return tuple(sorted(taken))


def chromosome_fold(
    tfs: list[str],
    train_chromosomes: list[str],
    test_chromosomes: list[str],
    weights: dict[str, int],
    validation_fraction: float,
    seed: int,
) -> GHTFold:
    """Setting 1: every TF trains, the benchmark's test chromosomes are held out."""
    validation = _validation_chromosomes(train_chromosomes, weights, validation_fraction, seed)
    return GHTFold(
        regime=CHROMOSOME_REGIME,
        name="all",
        test_tfs=tuple(sorted(tfs)),
        train_tfs=tuple(sorted(tfs)),
        test_chromosomes=tuple(sorted(test_chromosomes)),
        train_chromosomes=tuple(sorted(train_chromosomes)),
        validation_chromosomes=validation,
        held_out=f"{len(test_chromosomes)} chromosomes, all {len(tfs)} TFs in training",
    )


def _tf_folds(
    regime: str,
    tfs: list[str],
    groups: np.ndarray,
    train_chromosomes: list[str],
    test_chromosomes: list[str],
    n_folds: int,
    seed: int,
    validation_fraction: float,
    held_out: str,
) -> list[GHTFold]:
    """`n_folds` TF holdouts over `groups`, every member of a group landing in one fold."""
    rng = np.random.default_rng(seed)
    ids, inverse = np.unique(groups, return_inverse=True)
    sizes = np.bincount(inverse, minlength=len(ids))
    assignment = _balanced_groups(sizes, min(n_folds, len(ids)), rng)
    fold_of = assignment[inverse]

    names = np.asarray(tfs)
    folds = []
    for f in range(int(fold_of.max()) + 1):
        test = sorted(names[fold_of == f].tolist())
        train = sorted(names[fold_of != f].tolist())
        # The validation carve is whole TFs out of the training side, drawn with a seed that
        # depends on the fold so two folds do not carve the same proteins.
        carve = max(1, int(round(validation_fraction * len(train))))
        picked = np.random.default_rng([seed, f]).permutation(len(train))[:carve]
        folds.append(
            GHTFold(
                regime=regime,
                name=f"fold-{f}",
                test_tfs=tuple(test),
                train_tfs=tuple(train),
                test_chromosomes=tuple(sorted(test_chromosomes)),
                train_chromosomes=tuple(sorted(train_chromosomes)),
                validation_tfs=tuple(sorted(train[i] for i in picked)),
                held_out=f"{len(test)} TFs: {', '.join(test)}",
            )
        )
    return folds


def identity_components(panel: pd.DataFrame, min_identity: float) -> np.ndarray:
    """Connected components of the panel's domains at `min_identity`, within a family.

    The `S2` grouping of the PBM arm applied to 33 proteins instead of 1,338
    (`snp2prot.distances.connected_components`). At this size most components are singletons —
    the panel is 15 families over 33 proteins — so `G2` is close to `G1` with the near-duplicate
    pairs forced apart, which is exactly what it is for.
    """
    sequences = list(panel.dbd_seq)
    distances = distances_module.build(sequences, processes=1)
    min_overlap = float(thresholds.load()["cluster"]["min_overlap"])
    return distances_module.connected_components(
        distances, min_identity, min_overlap, families=panel.dbd_family.to_numpy()
    )


def all_regimes(
    panel: pd.DataFrame,
    train_chromosomes: list[str],
    test_chromosomes: list[str],
    weights: dict[str, int],
    n_folds: int,
    seed: int,
    validation_fraction: float,
    s2_min_identity: float,
) -> list[GHTFold]:
    """Every fold of every regime, in report order."""
    tfs = list(panel.tf)
    folds = [
        chromosome_fold(
            tfs, train_chromosomes, test_chromosomes, weights, validation_fraction, seed
        )
    ]
    folds += _tf_folds(
        "G1",
        tfs,
        np.arange(len(tfs)),
        train_chromosomes,
        test_chromosomes,
        n_folds,
        seed,
        validation_fraction,
        "random TFs",
    )
    components = identity_components(panel, s2_min_identity)
    folds += _tf_folds(
        "G2",
        tfs,
        components,
        train_chromosomes,
        test_chromosomes,
        n_folds,
        seed,
        validation_fraction,
        f"identity components at >= {s2_min_identity}",
    )
    return folds
