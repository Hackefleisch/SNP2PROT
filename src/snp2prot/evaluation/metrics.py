"""Retrieval metrics, computed per protein and macro-averaged.

**Positives first, everywhere** (`docs/ML_PLAN.md` §5.3, decided 2026-08-19). The dataset is
466:1 negative to positive, so any metric that rewards correct negatives measures the class
balance rather than the model: calling every pair non-binding scores 99.8% accuracy and is
worthless. AUPR is the primary number, precision@k and recall-at-fixed-precision are the
readable secondaries, and AUROC is reported only because some readers expect it.

**Nothing here scores against the raw E-score.** A Spearman correlation against it was reported
until 2026-08-26 and was removed: a universal-PBM E-score is a rank-enrichment statistic read by
the field at a cutoff, not a graded affinity, so correlating with its ordering reads a precision
into the assay that is not there. The label is the measurement. `spearman` survives in this
module as a utility — `scripts/check_pooling.py` correlates embedding displacement against edit
count with it — but it is not a model metric.

**Per protein, then macro-averaged.** For each held-out domain, rank all 32,896 8-mers, score
that ranking, and average across domains — every protein counts once whatever its positive
count. Never pooled: the question is *"does this recover this protein's motif"*, which is a
retrieval problem per protein, and pooling hides exactly the per-protein failures the `S2`
regime exists to expose. The per-domain values are returned alongside the mean, because the
distribution says far more than its average when positive counts run from 1 to 236.

Two things this module refuses to paper over:

**A domain with no positives has no AUPR.** 54 records have no positive 8-mer (`T21`), 20 of
them variants that measurably lost binding. Their AUPR is not zero and not one — it does not
exist, so it is `nan` here and the macro-average skips it and reports how many it skipped.

**But "undefined" and "failed" are different, and only the first may be skipped.** That
distinction was lost once: `recall_at_precision` returned `nan` for a ranking that never reached
the target precision, which is a failure and not an absence, and the macro-average dropped those
domains — reporting the nearest-neighbour baseline's `P1` recall as 0.130 when over all
positive-bearing domains it is 0.008. A failure now scores 0 and stays in the mean. Every skip
this module still makes is counted per metric in `macro_average`, so a figure can never again be
an average over an unstated subset.

**The no-call band is excluded, not counted as negative.** `label == -1` is absent evidence:
the 8-mer fell between the two cutoffs, or replicates disagreed. Ranking metrics see only
`label in (0, 1)` cells.

**Ties are not broken, they are collapsed.** Sorting equal scores by array position and then
scoring the result is how a metric silently rewards an ordering the model never expressed: a
prediction that scores every 8-mer identically comes out at AUPR 1.0 if the positives happen to
sit first. So a run of tied scores is treated as **one threshold** — precision and recall are
read at the end of the run (the definition scikit-learn's `average_precision_score` uses), and
precision@k takes the expected value over the tie group straddling the cut. A constant
prediction then scores the positive rate, which is what it deserves. On this corpus the largest
tie group in a real E-score profile is 2 of 32,896 8-mers, so this changes no number measured
here; it stops a coarser prediction from being flattered later.
"""

from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np

#: The default `k` values for precision@k, overridden from `configs/experiment.yaml`.
PRECISION_AT = (10, 50, 100)


def average_precision(labels: np.ndarray, scores: np.ndarray) -> float:
    """Area under the precision-recall curve, as the mean precision at each positive.

    The same quantity as scikit-learn's `average_precision_score`, implemented here rather
    than adding a dependency for twelve lines.
    """
    if labels.sum() == 0:
        return float("nan")
    precision, recall = _pr_curve(labels, scores)
    return float((np.diff(np.r_[0.0, recall]) * precision).sum())


def auroc(labels: np.ndarray, scores: np.ndarray) -> float:
    """Rank-based AUROC, ties averaged. Never a headline number at 466:1 (§5.3)."""
    n_pos = int(labels.sum())
    n_neg = len(labels) - n_pos
    if n_pos == 0 or n_neg == 0:
        return float("nan")
    ranks = _average_ranks(scores)
    return float((ranks[labels == 1].sum() - n_pos * (n_pos + 1) / 2) / (n_pos * n_neg))


def precision_at_k(labels: np.ndarray, scores: np.ndarray, k: int) -> float:
    """Of the top `k` 8-mers this ranking calls, the fraction that bind."""
    if k > len(scores):
        return float("nan")
    cut = np.partition(scores, len(scores) - k)[len(scores) - k]
    above = scores > cut
    tied = scores == cut
    # The tie group straddling the cut contributes its expected share, not whichever of its
    # members the sort happened to place first.
    share = (k - int(above.sum())) / int(tied.sum())
    return float((labels[above].sum() + labels[tied].sum() * share) / k)


def recall_at_precision(labels: np.ndarray, scores: np.ndarray, target: float) -> float:
    """Recall at the deepest cut whose precision still reaches `target`.

    **A ranking that never reaches the target scores 0, not `nan`.** It used to return `nan`,
    on the reading that "no operating point exists" is a different statement from "recall 0" —
    and `macro_average` then dropped those domains from the mean, so the figure was an average
    over whichever domains happened to succeed. Measured on the nearest-neighbour baseline
    (2026-08-26), that inflated the published `P1` number from **0.008 to 0.130**, an average
    over 24 of 412 positive-bearing domains; `S2/fold-4` went 0.172 to 0.623.

    Zero is also the right value on its own terms. At depth 1 the precision is 1.0 whenever the
    top-ranked 8-mer is a positive, so a ranking that never reaches 0.5 has put a negative first
    *and* failed to recover at any depth. The recall attainable at precision >= target is then
    exactly none — a failure with a value, not a quantity that does not exist.

    `nan` remains for the one genuinely undefined case: a domain with no positive 8-mer, where
    there is no recall to measure at any precision. `aupr` and `auroc` are `nan` for the same
    kind of reason and keep that behaviour.
    """
    if labels.sum() == 0:
        return float("nan")
    precision, recall = _pr_curve(labels, scores)
    reached = np.flatnonzero(precision >= target)
    return float(recall[reached[-1]]) if len(reached) else 0.0


def spearman(a: np.ndarray, b: np.ndarray) -> float:
    """Spearman rank correlation, ties averaged.

    **Not a model metric.** Nothing in the evaluation path scores a prediction against the raw
    E-score any more (see the module docstring). This stays because `scripts/check_pooling.py`
    needs a rank correlation for a different question — how far a pooled embedding moves against
    how many residues changed — and importing scipy for twelve lines is not worth it.
    """
    ra, rb = _average_ranks(a), _average_ranks(b)
    ra = ra - ra.mean()
    rb = rb - rb.mean()
    denominator = np.sqrt((ra**2).sum() * (rb**2).sum())
    return float((ra * rb).sum() / denominator) if denominator > 0 else float("nan")


def _pr_curve(labels: np.ndarray, scores: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Precision and recall at every distinct score, deepest cut last.

    One point per distinct score rather than one per element: a run of tied scores is a single
    threshold, and reading precision inside it would be reading an order the prediction never
    expressed.
    """
    order = np.argsort(-scores, kind="stable")
    hits = np.cumsum(labels[order])
    depth = np.arange(1, len(order) + 1)
    sorted_scores = scores[order]
    last_of_run = np.r_[sorted_scores[1:] != sorted_scores[:-1], True]
    hits, depth = hits[last_of_run], depth[last_of_run]
    return hits / depth, hits / labels.sum()


def _average_ranks(x: np.ndarray) -> np.ndarray:
    """Ascending ranks starting at 1, tied values sharing their mean rank."""
    order = np.argsort(x, kind="stable")
    ranks = np.empty(len(x), dtype=np.float64)
    ranks[order] = np.arange(1, len(x) + 1)
    values = x[order]
    start = 0
    for stop in np.flatnonzero(np.r_[values[1:] != values[:-1], True]) + 1:
        if stop - start > 1:
            ranks[order[start:stop]] = (start + stop + 1) / 2
        start = stop
    return ranks


def chance_aupr(labels: np.ndarray) -> float:
    """The macro AUPR a random ranking scores on these domains: the mean positive rate.

    **Not a formality — it is what makes an AUPR comparable between folds.** A per-protein AUPR
    is anchored to the protein's own positive rate, and the regimes hold out different domain
    mixes: across the 19 folds this ranges from 0.0015 to 0.0034, a factor of 2.3. Two folds
    reporting 0.05 are not reporting the same thing, and `P1`'s 0.0035 against a 0.0030 null is
    a different statement from `S2/fold-4`'s 0.15 against 0.0015.

    `labels` is `(n_domains, n_kmers)` in `{1, 0, -1}`. Domains with no positive are skipped,
    exactly as `macro_average` skips them.
    """
    scored = labels != -1
    n = scored.sum(axis=1)
    positive = ((labels == 1) & scored).sum(axis=1)
    keep = (positive > 0) & (n > 0)
    return float(np.mean(positive[keep] / n[keep])) if keep.any() else float("nan")


def random_baseline(labels: np.ndarray, repeats: int = 200, seed: int = 0) -> dict[str, float]:
    """The distribution of macro AUPR under a random ranking: `{random_mean, random_p95}`.

    **Not to be confused with the model's null anchor** (`snp2prot.models.two_tower`), which is a
    learned vector inside the shared space. This is a property of the held-out labels alone: what
    a predictor with no information scores on them.

    `chance_aupr` gives the null's centre but not its width, and a ratio near 1 is an eyeball
    rather than a test. This samples the whole macro statistic so a result can be called
    indistinguishable from random or not.

    **Sampled by drawing the positives' positions, not by shuffling a prediction.** The two are
    the same null — a shuffled ranking *is* a uniformly random ordering — but a shuffle costs an
    argsort over 32,896 elements per domain per repeat, while the positions of `n_pos` positives
    in a random permutation can be drawn directly and average precision read off them in
    `O(n_pos)`. Measured equal to the brute-force shuffle to three decimals, and about a thousand
    times faster: one fold at 200 repeats is 1.2 s rather than 20 minutes.

    Note the null mean sits slightly *above* the positive rate — average precision is biased
    upward at small `n_pos` — which is why the band is sampled rather than assumed.
    """
    rng = np.random.default_rng(seed)
    scored = labels != -1
    n = scored.sum(axis=1)
    positive = ((labels == 1) & scored).sum(axis=1)
    keep = (positive > 0) & (n > 0)
    n, positive = n[keep], positive[keep]
    if not len(n):
        return {"random_mean": float("nan"), "random_p95": float("nan")}

    samples = np.empty(repeats, dtype=np.float64)
    for r in range(repeats):
        samples[r] = np.mean(
            [
                np.mean(np.arange(1, k + 1) / (np.sort(rng.choice(total, k, replace=False)) + 1))
                for total, k in zip(n, positive, strict=True)
            ]
        )
    return {"random_mean": float(samples.mean()), "random_p95": float(np.percentile(samples, 95))}


def suppression(
    variant_scores: np.ndarray, reference_scores: np.ndarray, reference_labels: np.ndarray
) -> float:
    """Of the sites the wild type binds, the fraction the model scores *lower* for the variant.

    **The metric for a domain that binds nothing.** 20 records in a test fold have no positive
    8-mer — every one a `dead_variant`, a variant whose binding measurably vanished (`T21`) —
    and AUPR, AUROC and recall-at-precision are all undefined for them. They are also the
    sharpest evidence for claim C1, so leaving them unscored throws away the best case in the
    corpus.

    Counting "predicted positives" would need a decision threshold, and the model has none: the
    null anchor was intended as one but sits inside the negative cloud, with a median of 15,159
    of 32,460 8-mers scoring above it (measured 2026-08-26). So this compares the variant with
    its own wild type down the **protein axis** instead, which needs no threshold and no
    continuous assay value — only the model's two predictions against each other:

    - **1.0** — every site the wild type binds falls in the variant: the model saw that the
      mutation abolished binding.
    - **0.5** — the sites moved at random.
    - **0.0** — the variant is predicted exactly like its wild type.

    **The comparison is between ranks, not scores.** Comparing raw scores would hand 1.0 to any
    model that simply scores the variant lower everywhere, which is a global offset and not
    sensitivity to a mutation. Within-profile rank is invariant to that, so an untargeted
    downward shift scores 0.5 as it should.

    The nearest-neighbour baseline scores **exactly 0** whenever the wild type is in its training
    pool, because it then predicts the variant by copying that wild type — the null hypothesis
    for C1 in the literal sense rather than by interpretation. Where the wild type is held out
    too (`P1`, `S2`) the baseline copies something else and the number stops being interpretable,
    which is why the callers only score a variant whose reference the model was allowed to see.

    `nan` when the reference binds nothing either — the 2 of 20 dead records that are their own
    cluster's reference and simply do not bind.
    """
    positives = reference_labels == 1
    if not positives.any():
        return float("nan")
    variant_rank = _average_ranks(np.asarray(variant_scores, dtype=np.float64))
    reference_rank = _average_ranks(np.asarray(reference_scores, dtype=np.float64))
    return float((variant_rank[positives] < reference_rank[positives]).mean())


@dataclass
class DomainScores:
    """One held-out domain's ranking, scored."""

    domain: str
    n_pos: int
    aupr: float
    auroc: float
    recall_at_precision: float
    precision_at: dict[int, float] = field(default_factory=dict)


def score_domain(
    labels: np.ndarray,
    truth_escore: np.ndarray,
    predicted: np.ndarray,
    domain: str = "",
    precision_at: tuple[int, ...] = PRECISION_AT,
    precision_target: float = 0.5,
) -> DomainScores:
    """Score one predicted ranking over all 32,896 8-mers against one domain's measurements.

    `labels` carries the three-valued label, `truth_escore` the continuous score behind it.
    Ranking metrics see only the cells outside the no-call band.
    """
    scored = labels != -1
    binary = (labels[scored] == 1).astype(np.int64)
    ranked = predicted[scored]
    return DomainScores(
        domain=domain,
        n_pos=int(binary.sum()),
        aupr=average_precision(binary, ranked),
        auroc=auroc(binary, ranked),
        recall_at_precision=recall_at_precision(binary, ranked, precision_target),
        precision_at={k: precision_at_k(binary, ranked, k) for k in precision_at},
    )


def macro_average(scores: list[DomainScores]) -> dict[str, float]:
    """Mean over domains of each metric, skipping the ones where it is undefined.

    `n_scored` is how many domains the AUPR average is over and `n_undefined` how many were
    skipped for having no positive 8-mer — reported rather than absorbed, because a regime
    that holds out the 20 dead variants would otherwise show a mean over a silently smaller
    set (`ML_PLAN.md` §5.3).
    """
    if not scores:
        return {}
    aupr = np.array([s.aupr for s in scores])
    columns: dict[str, list] = {
        "aupr": list(aupr),
        "auroc": [s.auroc for s in scores],
        "recall_at_precision": [s.recall_at_precision for s in scores],
        **{
            f"precision_at_{k}": [s.precision_at[k] for s in scores] for k in scores[0].precision_at
        },
    }
    out: dict[str, float] = {
        "n_domains": float(len(scores)),
        "n_scored": float(np.isfinite(aupr).sum()),
        "n_undefined": float((~np.isfinite(aupr)).sum()),
        "aupr_median": float(np.nanmedian(aupr)) if np.isfinite(aupr).any() else float("nan"),
    }
    for name, values in columns.items():
        array = np.asarray(values, dtype=np.float64)
        out[name] = _nanmean(array)
        # How many domains each mean is actually over. `n_scored` covers AUPR and used to be
        # taken to cover everything, which is how a figure over a sixth of the domains was read
        # as a figure over all of them.
        out[f"n_scored_{name}"] = float(np.isfinite(array).sum())
    return out


def _nanmean(values) -> float:
    array = np.asarray(values, dtype=np.float64)
    return float(np.nanmean(array)) if np.isfinite(array).any() else float("nan")
