#!/usr/bin/env python
"""Score the nearest-neighbour lookup baseline in every split regime.

    python scripts/run_nn_baseline.py [--top-k] [--out reports/nn_baseline.md]

**Phase 7, step 0 — this runs before any model is trained** (`docs/ML_PLAN.md` §8.1, §9). It
does two jobs at once:

1. it sets the bar. Every later figure carries this number as a horizontal line, because a
   degradation curve is not diagnostic on its own: a family-holdout AUPR of 0.22 means the
   model is a lookup table if the baseline scores 0.21, and is the talk if the baseline scores
   0.05. Same number, opposite conclusions.
2. it validates the splits before any training run is spent on them. If the lookup scores near
   its `S1` ceiling under `S2`, the split groups leak — which is exactly what `T27` was about,
   and this is the cheapest possible check that connected-component grouping fixed it.

Reads three cached artifacts and nothing else: the 8-mer matrix, the distance matrix and the
cluster inventory. About 4 minutes, almost all of it in the per-domain metrics.

The training pool excludes the `no_evidence` records in every regime, matching what a model
would be trained on (`label_health.usable`, `T21`): the baseline must not be allowed to copy a
profile a model would never have seen.
"""

from __future__ import annotations

import argparse
import time
from pathlib import Path

import numpy as np
import pandas as pd

from snp2prot import corpus, distances, experiment, splits, thresholds, tracking
from snp2prot.baselines import nn_lookup
from snp2prot.config import PROCESSED_DIR, REPORTS_DIR
from snp2prot.data.matrix import KmerMatrix
from snp2prot.evaluation import c1, metrics

DEFAULT_OUT = REPORTS_DIR / "nn_baseline.md"
#: Every held-out domain's own numbers, for the figures §5.2 asks for: the report can
#: only carry summaries, and "which proteins fail" is the question worth plotting.
PER_DOMAIN_OUT = PROCESSED_DIR / "nn_baseline_domains.parquet"
#: The C1 evaluation set (`D6`): the held-out variants the wild-type copy fails on. Fixed
#: before any model runs, so every arm is scored on the same variants.
C1_OUT = PROCESSED_DIR / "c1_variants.parquet"

#: Nearest-neighbour identity bands. The boundary that matters is 0.7: Weirauch et al.
#: 2014 set per-family identity thresholds for motif transfer in the 60-70% range, so a
#: held-out domain above it is one the incumbent method would simply have looked up.
IDENTITY_BANDS = (0.0, 0.3, 0.5, 0.7, 0.9, 1.01)


def score_fold(
    fold: splits.Fold,
    domains: pd.DataFrame,
    matrix: KmerMatrix,
    dist: distances.DomainDistances,
    trainable: np.ndarray,
    min_overlap: float,
    k: int,
    precision_at: tuple[int, ...],
    precision_target: float,
) -> tuple[dict, pd.DataFrame]:
    """One fold: choose neighbours, copy profiles, score every held-out domain."""
    pool = fold.train[trainable[fold.train]]
    prediction = nn_lookup.fit_predict(matrix, dist, fold.test, pool, min_overlap, k)

    scores = [
        metrics.score_domain(
            matrix.label[row],
            matrix.escore[row],
            prediction.profile[i],
            domain=str(matrix.domains[row]),
            precision_at=precision_at,
            precision_target=precision_target,
        )
        for i, row in enumerate(fold.test)
    ]
    summary = metrics.macro_average(scores)
    summary |= {
        "regime": fold.regime,
        "fold": fold.name,
        "digest": fold.digest(domains),
        "n_train": float(len(pool)),
        "identity": float(prediction.neighbours.identity.median()),
        "n_without_neighbour": float(prediction.n_without_neighbour),
    }

    per_domain = pd.DataFrame(
        {
            "regime": fold.regime,
            "fold": fold.name,
            "domain": [s.domain for s in scores],
            "family": domains.dbd_family.to_numpy()[fold.test],
            "gene": domains.gene.to_numpy()[fold.test],
            "is_variant": domains.is_variant.to_numpy()[fold.test],
            "verdict": domains.verdict.to_numpy()[fold.test],
            "n_pos": [s.n_pos for s in scores],
            "aupr": [s.aupr for s in scores],
            "auroc": [s.auroc for s in scores],
            "spearman": [s.spearman for s in scores],
            "neighbour_identity": prediction.neighbours.identity.to_numpy(),
            "neighbour_edits": prediction.neighbours.n_edits.to_numpy(),
        }
    )
    return summary, per_domain


def positive_rate(matrix: KmerMatrix, rows: np.ndarray) -> float:
    """The AUPR a random ranking would score: the macro-average of each domain's own rate."""
    rates = []
    for row in rows:
        scored = matrix.label[row] != -1
        n = int(scored.sum())
        if n:
            rates.append(float((matrix.label[row] == 1).sum() / n))
    return float(np.mean(rates)) if rates else float("nan")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    ap.add_argument(
        "--top-k",
        action="store_true",
        help="also run the identity-weighted top-k form, the secondary bar (ML_PLAN 8.1)",
    )
    args = ap.parse_args()

    cfg = experiment.load()
    min_overlap = float(thresholds.load()["cluster"]["min_overlap"])
    precision_at = tuple(int(k) for k in cfg["metrics"]["precision_at"])
    precision_target = float(cfg["metrics"]["recall_at_precision"])
    top_k = int(cfg["baseline"]["nn_lookup"]["top_k"])

    domains = corpus.domains()
    matrix = KmerMatrix.load()
    dist = distances.DomainDistances.load()
    if list(matrix.domains) != list(domains.dbd_seq) or list(dist.domains) != list(domains.dbd_seq):
        raise SystemExit(
            "the cached matrices were built for a different domain set\n"
            "rerun scripts/build_matrix.py and scripts/build_distances.py"
        )
    trainable = corpus.trainable(domains).to_numpy()
    print(f"{len(domains)} domains, {trainable.sum()} of them trainable")

    ks = [1, top_k] if args.top_k else [1]
    started = time.time()
    summaries: list[dict] = []
    per_domain: list[pd.DataFrame] = []
    folds = list(splits.all_regimes(domains, dist))
    for fold in folds:
        for k in ks:
            summary, rows = score_fold(
                fold,
                domains,
                matrix,
                dist,
                trainable,
                min_overlap,
                k,
                precision_at,
                precision_target,
            )
            summary["k"] = float(k)
            rows["k"] = k
            summaries.append(summary)
            per_domain.append(rows)
        print(
            f"  {fold.label:<20} k=1 AUPR {summaries[-len(ks)]['aupr']:.4f}  "
            f"({time.time() - started:.0f}s)",
            flush=True,
        )

    frame = pd.DataFrame(summaries)
    domains_frame = pd.concat(per_domain, ignore_index=True)
    PER_DOMAIN_OUT.parent.mkdir(parents=True, exist_ok=True)
    domains_frame.to_parquet(PER_DOMAIN_OUT, index=False)
    c1.select(domains_frame).to_parquet(C1_OUT, index=False)
    write_report(args.out, frame, domains_frame, folds, matrix, cfg, min_overlap, ks)
    print(f"wrote {PER_DOMAIN_OUT} and {C1_OUT}")
    print(f"wrote {args.out} in {time.time() - started:.0f}s")


def _fmt(value: float, places: int = 4) -> str:
    return "n/a" if not np.isfinite(value) else f"{value:.{places}f}"


def write_report(
    path: Path,
    frame: pd.DataFrame,
    per_domain: pd.DataFrame,
    folds: list[splits.Fold],
    matrix: KmerMatrix,
    cfg: dict,
    min_overlap: float,
    ks: list[int],
) -> None:
    """The phase-boundary deliverable: one table per regime plus the diagnostics."""
    primary = frame[frame.k == 1]
    lines = [
        "# Nearest-neighbour lookup — the baseline every later number is read against",
        "",
        "Generated by `scripts/run_nn_baseline.py`. **This is step 0 of the modelling plan**",
        "(`docs/ML_PLAN.md` §8.1, §9): it runs before any model is trained, because a",
        "degradation curve cannot be interpreted without it and because it is the cheapest",
        "available check that the split regimes hold out what they claim to.",
        "",
        "For each held-out domain the baseline copies the **E-score profile of the most",
        "identical training domain**, ranks all 32,896 8-mers by it, and is scored per protein",
        "and macro-averaged. Identity is `1 - n_edits / n_aligned` over a BLOSUM62 alignment",
        f"with free terminal gaps, and a candidate must align over at least {min_overlap:.0%} of",
        "the shorter domain to be eligible.",
        "",
        "## The regimes",
        "",
        "| regime | fold | test | train | AUPR | median | AUROC | Spearman | P@10 | P@50 "
        "| R@P0.5 | NN identity |",
        "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in primary.itertuples():
        lines.append(
            f"| {row.regime} | `{row.fold}` | {int(row.n_domains)} | {int(row.n_train)} | "
            f"**{_fmt(row.aupr)}** | {_fmt(row.aupr_median)} | {_fmt(row.auroc, 3)} | "
            f"{_fmt(row.spearman, 3)} | {_fmt(row.precision_at_10, 3)} | "
            f"{_fmt(row.precision_at_50, 3)} | {_fmt(row.recall_at_precision, 3)} | "
            f"{_fmt(row.identity, 3)} |"
        )

    lines += [
        "",
        "`AUPR` is the macro-average over held-out domains and `median` the median of the same",
        "per-domain values; they differ because the distribution is skewed, and the plan asks",
        "for the distribution rather than the mean alone (§5.3). `NN identity` is the median",
        "percent identity between a held-out domain and the neighbour that was copied — the",
        "quantity that explains the differences between regimes.",
        "",
        "## The random baseline for each regime",
        "",
        "A per-protein AUPR is not comparable across regimes on its own: a random ranking",
        "scores a domain's own positive rate, and regimes hold out different domain mixes.",
        "",
        "| regime | fold | random AUPR | NN AUPR | ratio |",
        "|---|---|---:|---:|---:|",
    ]
    by_label = {f.label: f for f in folds}
    for row in primary.itertuples():
        fold = by_label[f"{row.regime}/{row.fold}"]
        chance = positive_rate(matrix, fold.test)
        ratio = row.aupr / chance if chance > 0 else float("nan")
        lines.append(
            f"| {row.regime} | `{row.fold}` | {_fmt(chance)} | {_fmt(row.aupr)} | "
            f"{_fmt(ratio, 1)}x |"
        )

    lines += [
        "",
        "## Does each regime hold out what it claims to?",
        "",
        "The reason to run this before training anything (§8.1). A regime is only as good as",
        "the distance it puts between a held-out domain and the closest thing left in",
        "training, and that distance is measurable directly. Weirauch et al. 2014 — one of",
        "this dataset's own sources — set per-family identity thresholds for motif transfer",
        "in the 60-70% range, so a held-out domain whose nearest training neighbour is above",
        "that is one the incumbent method would simply have looked up.",
        "",
        "| regime | median NN identity | >= 0.9 | >= 0.7 | < 0.5 |",
        "|---|---:|---:|---:|---:|",
    ]
    for regime, group in per_domain[per_domain.k == 1].groupby("regime", sort=False):
        identity = group.neighbour_identity.to_numpy()
        n = len(identity)
        lines.append(
            f"| {regime} | {_fmt(float(np.nanmedian(identity)), 3)} | "
            f"{np.nansum(identity >= 0.9) / n:.0%} | {np.nansum(identity >= 0.7) / n:.0%} | "
            f"{np.nansum(identity < 0.5) / n:.0%} |"
        )

    lines += [
        "",
        "And the same domains' AUPR, banded by how close that neighbour was. This is the",
        "quantitative companion §5.2 asks for: it *explains* a degradation curve rather than",
        "gesturing at it, and it is what a model's own numbers have to be read against band",
        "by band.",
        "",
        "| NN identity | n | AUPR | n with no positive |",
        "|---|---:|---:|---:|",
    ]
    banded = per_domain[per_domain.k == 1]
    for low, high in zip(IDENTITY_BANDS[:-1], IDENTITY_BANDS[1:], strict=True):
        inside = banded[(banded.neighbour_identity >= low) & (banded.neighbour_identity < high)]
        if not len(inside):
            continue
        aupr = inside.aupr.to_numpy()
        lines.append(
            f"| {low:.1f} – {min(high, 1.0):.1f} | {len(inside)} | "
            f"{_fmt(float(np.nanmean(aupr)))} | {int((~np.isfinite(aupr)).sum())} |"
        )

    lines += [
        "",
        "## Per-domain distribution, by regime",
        "",
        "| regime | n scored | n with no positive | AUPR p10 | p25 | median | p75 | p90 |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for regime, group in per_domain[per_domain.k == 1].groupby("regime", sort=False):
        aupr = group.aupr.to_numpy()
        finite = aupr[np.isfinite(aupr)]
        q = np.percentile(finite, [10, 25, 50, 75, 90]) if len(finite) else [np.nan] * 5
        lines.append(
            f"| {regime} | {len(finite)} | {len(aupr) - len(finite)} | "
            + " | ".join(_fmt(v, 3) for v in q)
            + " |"
        )

    lines += [
        "",
        "## Held-out variants, by family",
        "",
        "The variant regimes carry claim C1, and 84 of the 173 variants are homeodomain — so a",
        "pooled variant number is a homeodomain number wearing a general claim's clothes",
        "(`docs/ML_PLAN.md` §6.1). Under `P3/all` no variant is in training, so the nearest",
        "training neighbour of every held-out variant **is its own wild type**: this row is the",
        "hypothesis *the mutation has no effect*, stated as a number.",
        "",
        "| family | n | AUPR | median NN identity | median edits to neighbour |",
        "|---|---:|---:|---:|---:|",
    ]
    p3 = per_domain[(per_domain.k == 1) & (per_domain.fold == "all")]
    for family, group in p3.groupby("family", sort=False):
        if len(group) < 3:
            continue
        lines.append(
            f"| {family} | {len(group)} | {_fmt(np.nanmean(group.aupr))} | "
            f"{_fmt(group.neighbour_identity.median(), 3)} | "
            f"{int(group.neighbour_edits.median())} |"
        )
    other = p3.groupby("family").filter(lambda g: len(g) < 3)
    if len(other):
        lines.append(
            f"| *{other.family.nunique()} families with < 3 variants* | {len(other)} | "
            f"{_fmt(np.nanmean(other.aupr))} | {_fmt(other.neighbour_identity.median(), 3)} | "
            f"{int(other.neighbour_edits.median())} |"
        )

    dead = p3[p3.verdict == "dead_variant"]
    if len(dead):
        lines += [
            "",
            "### The dead variants",
            "",
            f"{len(dead)} of the held-out variants measurably **lost** binding and have no",
            "positive 8-mer at all, so their AUPR does not exist and they are excluded from",
            "every average above. They are the sharpest C1 evidence in the corpus: the baseline",
            "copies their wild type, which is maximally wrong for them, and a model that does",
            "the same has failed in precisely the way this project exists to detect (`T30`).",
            "Spearman against their own E-scores is the number that is defined for them.",
            "",
            "- median Spearman of the copied wild-type profile: "
            f"**{_fmt(dead.spearman.median(), 3)}**",
            "- median identity to the copied neighbour: "
            f"{_fmt(dead.neighbour_identity.median(), 3)}",
            f"- families: {', '.join(sorted(dead.family.unique()))}",
        ]

    selected = c1.select(per_domain)
    variants = per_domain[(per_domain.k == 1) & (per_domain.fold == "all")]
    cut = float(cfg["c1_set"]["max_baseline_aupr"])
    lines += [
        "",
        "## The C1 evaluation set",
        "",
        "**The number above that C1 has to be read against is `P3/all`, and it is close to",
        "saturated**: mean AUPR "
        f"{_fmt(float(primary[primary.fold == 'all'].aupr.iloc[0]))}, median "
        f"{_fmt(float(primary[primary.fold == 'all'].aupr_median.iloc[0]), 3)}. For most of the",
        "173 variants the wild-type profile simply is the right",
        "answer, because most single substitutions do not measurably change which 8-mers a",
        "domain binds. So claim C1 does not live in that mean — it lives in the variants where",
        "the wild type is wrong, and pooled with the rest they are invisible (`D6`, decided",
        "2026-08-19).",
        "",
        f"The set is every held-out variant whose wild-type copy scores below **{cut}**, plus",
        "every variant whose AUPR does not exist at all because it has no positive 8-mer — a",
        "`dead_variant`, which measurably lost binding and has a control in its own series.",
        "",
        f"**{len(selected)} of {len(variants)} variants**: "
        f"{int((selected.c1_reason == c1.LOST_BINDING).sum())} that lost binding entirely and "
        f"{int((selected.c1_reason == c1.POORLY_PREDICTED).sum())} whose profile changed enough",
        "that their own wild type does not predict it. Written to",
        "`data/processed/c1_variants.parquet`, fixed before any model runs so every arm is",
        "scored on the same variants.",
        "",
        "| family | variants | in the set | of those, lost binding | share |",
        "|---|---:|---:|---:|---:|",
    ]
    for row in c1.summarise(selected, variants).itertuples():
        if not row.n_in_set:
            continue
        lines.append(
            f"| {row.family} | {row.n_variants} | {row.n_in_set} | {row.n_dead} | "
            f"{row.fraction:.0%} |"
        )
    lines += [
        "",
        "**Read this as a caveat, not a result.** The set is *defined* by the baseline scoring",
        'badly on it, so the baseline scores badly on it by construction and "the model beats',
        "the baseline here\" is close to vacuous. What makes it mean something is the model's",
        "**absolute** number on the set, reported next to its number on the complement: a model",
        "that had merely learned to distrust wild types everywhere would gain here and lose",
        "there, which is a shifted prior and not C1.",
        "",
        f"Note also that {selected.family.nunique()} families carry the whole set, so any C1 "
        "claim from it is a claim about those folds and has to be worded that way.",
        "",
    ]

    if len(ks) > 1:
        lines += [
            "",
            "## The top-k form",
            "",
            "Identity-weighted mean of the top neighbours instead of a single copy — the",
            "secondary, slightly stronger bar.",
            "",
            f"| regime | fold | AUPR k=1 | AUPR k={ks[1]} | delta |",
            "|---|---|---:|---:|---:|",
        ]
        merged = primary.merge(
            frame[frame.k == ks[1]], on=["regime", "fold"], suffixes=("_1", "_k")
        )
        for row in merged.itertuples():
            lines.append(
                f"| {row.regime} | `{row.fold}` | {_fmt(row.aupr_1)} | {_fmt(row.aupr_k)} | "
                f"{row.aupr_k - row.aupr_1:+.4f} |"
            )

    lines += [
        "",
        "## Reproducing this",
        "",
        "```bash",
        "python scripts/build_matrix.py        # domain x 8-mer arrays, ~5 s",
        "python scripts/build_distances.py     # all-vs-all domain alignment, ~15 s",
        "python scripts/run_nn_baseline.py     # this report",
        "```",
        "",
        f"Seed `{cfg['splits']['seed']}`, {cfg['splits']['n_folds']} folds, overlap guard "
        f"{min_overlap:.0%}. Fold membership is hashed rather than described, because a regime",
        "name and a seed do not pin down which domains were held out once the corpus changes",
        "(`docs/ML_PLAN.md` §9.2).",
        "",
        "| regime | fold | held out | digest |",
        "|---|---|---|---|",
    ]
    for row in primary.itertuples():
        fold = by_label[f"{row.regime}/{row.fold}"]
        lines.append(f"| {row.regime} | `{row.fold}` | {fold.held_out} | `{row.digest}` |")

    lines += ["", tracking.provenance_line()]

    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
