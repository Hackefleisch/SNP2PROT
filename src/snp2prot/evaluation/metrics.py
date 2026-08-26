"""Retrieval metrics, computed per protein and macro-averaged.

**Positives first, everywhere** (`docs/ML_PLAN.md` §5.3, decided 2026-08-19). The dataset is
466:1 negative to positive, so any metric that rewards correct negatives measures the class
balance rather than the model: calling every pair non-binding scores 99.8% accuracy and is
worthless. AUPR is the primary number, precision@k and recall-at-fixed-precision are the
readable secondaries, Spearman against the raw E-score is the one that does not throw away
the no-call band's information, and AUROC is reported only because some readers expect it.

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
Their Spearman against the raw E-score is perfectly well defined, and is the number to read
for them.

**The no-call band is excluded, not counted as negative.** `label == -1` is absent evidence:
the 8-mer fell between the two cutoffs, or replicates disagreed. Ranking metrics see only
`label in (0, 1)` cells; Spearman uses every cell, because it reads the continuous score
underneath the label rather than the label.

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

    `nan` when the ranking never reaches that precision at any depth — which is a different
    statement from "recall 0" and is reported as such.
    """
    if labels.sum() == 0:
        return float("nan")
    precision, recall = _pr_curve(labels, scores)
    reached = np.flatnonzero(precision >= target)
    return float(recall[reached[-1]]) if len(reached) else float("nan")


def spearman(a: np.ndarray, b: np.ndarray) -> float:
    """Spearman rank correlation, ties averaged."""
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


def null_aupr(labels: np.ndarray, repeats: int = 200, seed: int = 0) -> dict[str, float]:
    """The distribution of macro AUPR under a random ranking: `{null_mean, null_p95}`.

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
        return {"null_mean": float("nan"), "null_p95": float("nan")}

    samples = np.empty(repeats, dtype=np.float64)
    for r in range(repeats):
        samples[r] = np.mean(
            [
                np.mean(np.arange(1, k + 1) / (np.sort(rng.choice(total, k, replace=False)) + 1))
                for total, k in zip(n, positive, strict=True)
            ]
        )
    return {"null_mean": float(samples.mean()), "null_p95": float(np.percentile(samples, 95))}


@dataclass
class DomainScores:
    """One held-out domain's ranking, scored."""

    domain: str
    n_pos: int
    aupr: float
    auroc: float
    spearman: float
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
    Ranking metrics see only the cells outside the no-call band; Spearman sees all of them.
    """
    scored = labels != -1
    binary = (labels[scored] == 1).astype(np.int64)
    ranked = predicted[scored]
    return DomainScores(
        domain=domain,
        n_pos=int(binary.sum()),
        aupr=average_precision(binary, ranked),
        auroc=auroc(binary, ranked),
        spearman=spearman(truth_escore, predicted),
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
    out: dict[str, float] = {
        "n_domains": float(len(scores)),
        "n_scored": float(np.isfinite(aupr).sum()),
        "n_undefined": float((~np.isfinite(aupr)).sum()),
        "aupr": _nanmean(aupr),
        "aupr_median": float(np.nanmedian(aupr)) if np.isfinite(aupr).any() else float("nan"),
        "auroc": _nanmean([s.auroc for s in scores]),
        "spearman": _nanmean([s.spearman for s in scores]),
        "recall_at_precision": _nanmean([s.recall_at_precision for s in scores]),
    }
    for k in scores[0].precision_at:
        out[f"precision_at_{k}"] = _nanmean([s.precision_at[k] for s in scores])
    return out


def _nanmean(values) -> float:
    array = np.asarray(values, dtype=np.float64)
    return float(np.nanmean(array)) if np.isfinite(array).any() else float("nan")
