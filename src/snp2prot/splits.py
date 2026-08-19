"""Split protocols. Five regimes, always reported separately.

A single pooled number hides exactly the failure this project is designed to detect (brief
§5b): a model that memorises per-protein profiles scores well wherever a near-relative of the
held-out domain stayed in training, and the only way to see that is to vary how far the
holdout reaches and report each distance separately.

Every regime holds out **domains**, never (domain, 8-mer) pairs. Splitting pairs would put
the same protein on both sides of the split and measure almost nothing; and the DNA-side
holdout was dropped on 2026-08-19, because the talk is about the protein axis
(`docs/ML_PLAN.md` §5).

| regime | held out | tests |
|---|---|---|
| **S1** | random domains, 5-fold | the ceiling; the optimistic case |
| **S2** | connected components at >= 0.5 identity, 5-fold | genuinely unseen proteins |
| **P1** | the whole `Homeodomain` family | can an unseen *fold* be predicted |
| **P2** | the non-homeodomain variant-bearing clusters | does mutation sensitivity transfer |
| **P3** | a fraction of the 173 variants | how much mutation supervision is needed |

**S2 groups by connected component, not by cluster** (`T27`, settled 2026-08-19,
`docs/DECISIONS.md` §2). Two domains of the same family are joined by an edge when they reach
`splits.s2_min_identity` over an alignment covering `cluster.min_overlap` of the shorter one,
and the whole component is held out together. CD-HIT's greedy assignment compares each
sequence only to cluster *seeds*, so two near-identical domains can land in different clusters
and a cluster-level holdout would leave one of them in training — which is the leak the
baseline exists to detect. Single linkage is right here and wrong in `build_clusters.py`, and
the asymmetry is the point: over-grouping removes more from training than strictly necessary
but never leaks, whereas in a *stored* cluster it would be a false claim about relatedness.

**The floor is 0.5 identity, and it is deliberately far looser than a cluster** (`D5`, decided
2026-08-19). Grouping at `<= 5` edits — the cluster rule, about 94% identity on a 77 aa domain
— was measured and cost the nearest-neighbour baseline 0.018 AUPR against a random split,
because copying a 70-90% identical neighbour already scores 0.927. At 0.5 the split is hard by
construction: nothing within reach of the identity band Weirauch et al. 2014 set for
transferring a motif stays in training. It chains, which is the safe direction here and only
ever costs statistical power.

**The grouping is evaluation-time only.** Nothing here writes to `data/`; `wt_id`,
`mut_positions` and the cluster inventory are untouched.

**The 1,057 singleton clusters stay in training in every regime.** They are the backbone that
keeps the protein axis at full width; only the designated group is ever removed
(`docs/ML_PLAN.md` §6.1).

## Fold membership is recorded, not reconstructed

A metric is meaningless without knowing which domains were held out, and "regime S2, fold 3,
seed 20260819" does not pin that down — the assignment also depends on the corpus, which
grows. Every `Fold` therefore carries a `digest` of its held-out domains, which is what goes
into a run record (`docs/ML_PLAN.md` §9.2: log the split, not just the hyperparameters).
"""

from __future__ import annotations

import hashlib
from collections.abc import Iterator, Sequence
from dataclasses import dataclass

import numpy as np
import pandas as pd

from snp2prot import distances as distances_module
from snp2prot import experiment, thresholds
from snp2prot.distances import DomainDistances


@dataclass(frozen=True)
class Fold:
    """One train/test division, as positions into the domain table it was built from."""

    regime: str
    #: Distinguishes folds within a regime: `"fold-2"`, `"all"`, `"half:seed-1"`.
    name: str
    test: np.ndarray
    train: np.ndarray
    #: What the regime removed, for the report: `"Homeodomain"`, `"86 of 173 variants"`.
    held_out: str = ""

    @property
    def label(self) -> str:
        return f"{self.regime}/{self.name}"

    def digest(self, domains: pd.DataFrame) -> str:
        """A short hash of the held-out domain sequences: the run's record of this split."""
        held = sorted(domains.dbd_seq.to_numpy()[self.test])
        return hashlib.sha256("\n".join(held).encode()).hexdigest()[:12]


def _fold(regime: str, name: str, test: np.ndarray, n: int, held_out: str = "") -> Fold:
    test = np.sort(np.asarray(test, dtype=np.int64))
    train = np.setdiff1d(np.arange(n, dtype=np.int64), test, assume_unique=False)
    return Fold(regime, name, test, train, held_out)


def _balanced_groups(sizes: Sequence[int], n_folds: int, rng: np.random.Generator) -> np.ndarray:
    """Assign groups to folds, largest first into the emptiest fold.

    Random assignment of *groups* leaves folds of very different sizes when the groups are —
    here one component of 8 domains against 1,057 of one. Largest-first into the emptiest fold
    is the standard greedy balance; the shuffle only breaks ties between equal-sized groups, so
    the seed still matters and the result is still random.
    """
    order = rng.permutation(len(sizes))
    order = order[np.argsort(-np.asarray(sizes)[order], kind="stable")]
    load = np.zeros(n_folds, dtype=np.int64)
    assignment = np.empty(len(sizes), dtype=np.int64)
    for g in order:
        f = int(np.argmin(load))
        assignment[g] = f
        load[f] += sizes[g]
    return assignment


def grouped_folds(groups: np.ndarray, n_folds: int, seed: int, regime: str) -> list[Fold]:
    """5-fold CV over groups of domains: every domain of a group lands in one fold."""
    rng = np.random.default_rng(seed)
    ids, inverse = np.unique(groups, return_inverse=True)
    sizes = np.bincount(inverse, minlength=len(ids))
    fold_of_group = _balanced_groups(sizes, n_folds, rng)
    fold_of_domain = fold_of_group[inverse]
    n = len(groups)
    return [
        _fold(
            regime,
            f"fold-{f}",
            np.flatnonzero(fold_of_domain == f),
            n,
            held_out=f"{int((fold_of_domain == f).sum())} domains",
        )
        for f in range(n_folds)
    ]


def s1_random(domains: pd.DataFrame, n_folds: int, seed: int) -> list[Fold]:
    """S1 — 5-fold CV holding out random *domains*.

    Optimistic by construction, and that is its job: variants of one cluster can land either
    side of the split, so it measures the ceiling rather than generalisation.
    """
    return grouped_folds(np.arange(len(domains)), n_folds, seed, "S1")


def s2_components(
    domains: pd.DataFrame,
    distances: DomainDistances,
    n_folds: int,
    seed: int,
    min_identity: float | None = None,
    min_overlap: float | None = None,
) -> list[Fold]:
    """S2 — 5-fold CV holding out whole connected components of the domain graph."""
    if min_identity is None:
        min_identity = float(experiment.section("splits")["s2_min_identity"])
    if min_overlap is None:
        min_overlap = float(thresholds.load()["cluster"]["min_overlap"])
    components = distances_module.connected_components(
        distances, min_identity, min_overlap, domains.dbd_family.to_numpy()
    )
    return grouped_folds(components, n_folds, seed, "S2")


def p1_family(domains: pd.DataFrame, family: str) -> list[Fold]:
    """P1 — hold out one whole Pfam family, `Homeodomain` by default.

    Unaffected by the cluster-versus-component question: clustering runs `same_family_only`,
    so every cluster and every component lies inside one family and a family holdout removes
    both members of a near-identical pair together.
    """
    test = np.flatnonzero(domains.dbd_family.to_numpy() == family)
    if not len(test):
        raise ValueError(f"no domains of family {family!r}")
    return [_fold("P1", family, test, len(domains), held_out=f"all {len(test)} {family} domains")]


def variant_clusters(domains: pd.DataFrame, exclude_family: str | None = None) -> np.ndarray:
    """`wt_id`s of the clusters that hold at least one variant, optionally excluding a family."""
    variant_bearing = domains.loc[domains.is_variant, ["wt_id", "dbd_family"]].drop_duplicates()
    if exclude_family is not None:
        variant_bearing = variant_bearing[variant_bearing.dbd_family != exclude_family]
    return variant_bearing.wt_id.to_numpy()


def p2_variant_clusters(domains: pd.DataFrame, exclude_family: str) -> list[Fold]:
    """P2 — hold out the variant-bearing clusters outside one family, references included.

    Scoped to the variant-bearing clusters and not to their whole families: the 1,057
    singletons stay, so the model *has* seen these folds and just never a mutated member of
    one. That makes P2 a mutation-transfer test rather than a second family holdout, and the
    deliberately gentler mirror of P1 (`docs/ML_PLAN.md` §6.1).
    """
    ids = variant_clusters(domains, exclude_family)
    test = np.flatnonzero(domains.wt_id.isin(ids).to_numpy())
    n_variants = int(domains.is_variant.to_numpy()[test].sum())
    return [
        _fold(
            "P2",
            f"non-{exclude_family}",
            test,
            len(domains),
            held_out=f"{len(ids)} clusters, {len(test)} domains, {n_variants} of them variants",
        )
    ]


def p3_variants(
    domains: pd.DataFrame, fraction: float, seed: int, name: str, repeat: int = 0
) -> Fold:
    """P3 — hold out a fraction of the variants, leaving every wild type in training.

    At `fraction = 1.0` this is the zero-shot extreme, and it is what makes the
    nearest-neighbour baseline the exact null hypothesis for claim C1: with no variant in
    training, a held-out variant's nearest training domain is its own wild type one edit away,
    so the baseline predicts precisely "the mutation has no effect" (`docs/ML_PLAN.md` §6.1).
    """
    variants = np.flatnonzero(domains.is_variant.to_numpy())
    n = max(1, int(round(fraction * len(variants))))
    if n >= len(variants):
        test = variants
    else:
        # Seeded per repeat, so the three draws of `half` differ and the report can show the
        # spread rather than a sample of one.
        test = np.random.default_rng(seed + repeat).choice(variants, size=n, replace=False)
    label = name if fraction >= 1.0 else f"{name}:draw-{repeat}"
    held_out = f"{len(test)} of {len(variants)} variants"
    return _fold("P3", label, test, len(domains), held_out=held_out)


def all_regimes(
    domains: pd.DataFrame, distances: DomainDistances, config: dict | None = None
) -> Iterator[Fold]:
    """Every fold of every regime, in the order the report presents them."""
    cfg = experiment.section("splits") if config is None else config
    n_folds, seed = int(cfg["n_folds"]), int(cfg["seed"])

    yield from s1_random(domains, n_folds, seed)
    yield from s2_components(domains, distances, n_folds, seed, float(cfg["s2_min_identity"]))
    yield from p1_family(domains, cfg["p1_family"])
    yield from p2_variant_clusters(domains, cfg["p1_family"])
    for name, fraction in cfg["p3"]["fractions"].items():
        for repeat in range(int(cfg["p3"]["repeats"][name])):
            yield p3_variants(domains, float(fraction), seed, name, repeat)
