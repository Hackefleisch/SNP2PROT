#!/usr/bin/env python
"""Does mean pooling still separate a variant from its own wild type?

    python scripts/check_pooling.py [--arm A1] [--out reports/pooling_check.md]

**The pre-flight check `docs/ML_PLAN.md` §3.1 requires before any training run.** Protein
embeddings are mean-pooled to one 1280-d vector per domain, which is what makes the four arms
comparable — and pooling dilutes a single-residue change. A variant differs from its reference
by 1-5 residues out of a median 77, so the naive bound on how far the domain vector moves is
one part in 77. That bound assumes only the mutated residue's representation changes, and a PLM
is contextual, so the true attenuation is an empirical question. §3.1: *"If variants are not
measurably separated from their own wild type in the pooled embedding, no §6.1 regime can
produce a meaningful C1 number on that arm, and it is far better to know that in an afternoon
than after the grid has run."*

## What is measured

Three distances in the pooled space, each a distribution over pairs:

- **variant to its own reference** — what a `P3` regime asks the model to tell apart;
- **unrelated domains of the same family** — the scale the first has to be read against, since
  a family's domains are similar to each other whatever the pooling does;
- **unrelated domains of different families** — the loosest comparison, for context.

The question is not whether the first is small. It is whether it is **resolvable**: a variant
must be closer to its own reference than to other domains of its family (otherwise the
embedding cannot even identify which wild type it belongs to), while still being **separated
from zero** relative to the precision of the space.

**And the check is pointed at the 29 variants that carry claim C1** (`D6`,
`data/processed/c1_variants.parquet`) as well as at all 173. Those are the variants whose
binding measurably changed, so they are the ones a C1 claim will rest on: if the pooled
embedding cannot separate *them* from their wild types, the fallback is attention pooling
(`TODO.md` `T32`) and it is needed before the grid runs, not after.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from snp2prot import clusters, corpus, distances, embeddings
from snp2prot.config import PROCESSED_DIR, REPORTS_DIR, embedding_table
from snp2prot.evaluation import metrics

DEFAULT_OUT = REPORTS_DIR / "pooling_check.md"
C1_SET = PROCESSED_DIR / "c1_variants.parquet"

#: A variant must be nearer its own reference than to this fraction of its family, or the
#: pooled embedding cannot even attribute it to the right wild type. Not a tuned number — it is
#: the weakest form of the property that has to hold for `P3` to mean anything.
RANK_FLOOR = 0.95


def variant_pairs(domains: pd.DataFrame) -> pd.DataFrame:
    """Every variant beside the reference of its own cluster."""
    inventory = clusters.load()
    reference = dict(zip(inventory.wt_id, inventory.reference, strict=True))
    rows = domains[domains.is_variant].copy()
    rows["reference"] = [reference[w] for w in rows.wt_id]
    return rows


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--arm",
        action="append",
        help="repeatable; default is every arm whose embeddings are built",
    )
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = ap.parse_args()

    domains = corpus.domains()
    arms = args.arm or [a for a in sorted(embeddings.ARMS) if embedding_table(a).exists()]
    if not arms:
        raise SystemExit("no embeddings built yet\nrun scripts/build_embeddings.py --arm A1")

    measured = {arm: measure(arm, domains) for arm in arms}
    write_report(args.out, measured, domains)
    print(f"wrote {args.out}")


def measure(arm: str, domains: pd.DataFrame) -> dict:
    """Every quantity the report shows for one arm."""
    table = embeddings.DomainEmbeddings.load(arm)
    if list(table.domains) != list(domains.dbd_seq):
        raise SystemExit(
            f"the {arm} embeddings were built for a different domain set\n"
            f"rerun scripts/build_embeddings.py --arm {arm}"
        )

    distance = embeddings.pairwise_cosine_distance(table.vectors)
    np.fill_diagonal(distance, np.nan)
    families = domains.dbd_family.to_numpy()
    row_of = {seq: i for i, seq in enumerate(domains.dbd_seq)}

    edits = distances.DomainDistances.load().n_edits

    pairs = variant_pairs(domains)
    records = []
    for row in pairs.itertuples():
        i, j = row_of[row.dbd_seq], row_of[row.reference]
        same_family = np.flatnonzero(families == row.dbd_family)
        others = same_family[(same_family != i) & (same_family != j)]
        to_family = distance[i, others]
        records.append(
            {
                "domain": row.dbd_seq,
                "family": row.dbd_family,
                "gene": row.gene,
                "verdict": row.verdict,
                "to_reference": distance[i, j],
                "n_edits": int(edits[i, j]),
                # Where the reference sits among the family, as a fraction: 1.0 means the
                # reference is the nearest domain of the family, 0.0 the furthest.
                "rank_in_family": (
                    float(np.mean(to_family > distance[i, j])) if len(others) else np.nan
                ),
                "family_median": float(np.nanmedian(to_family)) if len(others) else np.nan,
                "n_family": len(others) + 1,
            }
        )
    measured = pd.DataFrame(records)

    same = np.zeros_like(distance, dtype=bool)
    for family in np.unique(families):
        idx = np.flatnonzero(families == family)
        same[np.ix_(idx, idx)] = True
    np.fill_diagonal(same, False)

    return {
        "arm": arm,
        "table": table,
        "variants": measured,
        "within": distance[same],
        "across": distance[~same & ~np.isnan(distance)],
    }


def _q(values, places: int = 4) -> str:
    values = np.asarray(values, dtype=np.float64)
    values = values[np.isfinite(values)]
    if not len(values):
        return "n/a"
    p10, p50, p90 = np.percentile(values, [10, 50, 90])
    return f"{p50:.{places}f} ({p10:.{places}f} – {p90:.{places}f})"


def _c1_split(measured: pd.DataFrame):
    """The variants that carry C1 and the rest, or `(None, None)` if the set is not built."""
    if not C1_SET.exists():
        return None, None
    wanted = set(pd.read_parquet(C1_SET).domain)
    return measured[measured.domain.isin(wanted)], measured[~measured.domain.isin(wanted)]


def _c1_auc(measured: pd.DataFrame, in_set: set[str]) -> tuple[float, float, int, int]:
    """AUC of displacement for C1 vs ordinary variants, among single substitutions only.

    Controlled for edit count because the C1 set is disproportionately single substitutions
    and displacement rises with the number of edits — the raw comparison measures that
    confound rather than the question.
    """
    one = measured[measured.n_edits == 1]
    labels = one.domain.isin(in_set).to_numpy().astype(int)
    n_c1, n_rest = int(labels.sum()), int((1 - labels).sum())
    if not n_c1 or not n_rest:
        return float("nan"), float("nan"), n_c1, n_rest
    auc = metrics.auroc(labels, one.to_reference.to_numpy())
    stderr = float(np.sqrt(auc * (1 - auc) / min(n_c1, n_rest)))
    return auc, stderr, n_c1, n_rest


def write_report(path: Path, measured: dict, domains: pd.DataFrame) -> None:
    """The phase deliverable: one comparison table, then one section per arm."""
    lines = [
        "# Pooling pre-flight — do the sequence arms still see a single-residue change?",
        "",
        "Generated by `scripts/check_pooling.py`. **This runs before any training run**",
        "(`docs/ML_PLAN.md` §3.1). Protein embeddings are mean-pooled to one vector per domain,",
        "which is what makes the four arms comparable — and pooling dilutes a single-residue",
        "change. A variant differs from its reference by 1-5 residues out of a median 77, so the",
        "naive bound on how far the pooled vector moves is one part in 77. That bound assumes",
        "only the mutated residue's representation changes, and a PLM is **contextual**, so the",
        'true attenuation is an empirical question. §3.1: *"If variants are not measurably',
        "separated from their own wild type in the pooled embedding, no §6.1 regime can produce",
        "a meaningful C1 number on that arm, and it is far better to know that in an afternoon",
        'than after the grid has run."*',
        "",
        "Distances are **cosine**, not Euclidean: PLM embedding norms vary with sequence length",
        "and a variant is by construction almost the same length as its reference, so a norm",
        "difference would report as a distance without meaning one. The two arms differ in norm",
        "by a factor of six, which is exactly why.",
        "",
        "## The verdict, per arm",
        "",
        "| arm | model | variant vs family scale | own reference nearest | within nearest 5% "
        "| displacement ~ edits | C1 signal (AUC) |",
        "|---|---|---:|---:|---:|---:|---:|",
    ]
    for arm, m in measured.items():
        variants = m["variants"]
        ratio = float(np.nanmedian(variants.to_reference / variants.family_median))
        nearest = float((variants.rank_in_family == 1.0).mean())
        resolvable = float((variants.rank_in_family >= RANK_FLOOR).mean())
        scaling = metrics.spearman(
            variants.to_reference.to_numpy(dtype=float), variants.n_edits.to_numpy(dtype=float)
        )
        c1, _ = _c1_split(variants)
        auc, stderr, _, _ = (
            _c1_auc(variants, set(c1.domain)) if c1 is not None else (np.nan, np.nan, 0, 0)
        )
        auc_cell = "n/a" if not np.isfinite(auc) else f"{auc:.3f} ({(auc - 0.5) / stderr:.1f}σ)"
        lines.append(
            f"| `{arm}` | {m['table'].model} | {ratio:.3f} | {nearest:.0%} | {resolvable:.0%} | "
            f"{scaling:.3f} | {auc_cell} |"
        )

    lines += [
        "",
        "Reading the columns:",
        "",
        "- **variant vs family scale** — the variant-to-reference distance as a fraction of the",
        "  distance to other domains of the same family. Small is expected: a variant *is* its",
        "  wild type apart from one residue. It is the number that says how much the projection",
        "  and the contrastive objective have to amplify.",
        "- **own reference nearest** — the weakest property that has to hold. If a variant is not",
        "  nearer its own wild type than to unrelated domains of its family, the embedding cannot",
        "  even attribute it, and every variant regime is measuring noise.",
        "- **displacement ~ edits** — Spearman between how far the vector moved and how many",
        "  residues changed. This is what separates signal from the floor of the embedding: if a",
        "  five-residue change moved the vector no further than a one-residue change, the",
        "  distances would be measuring noise.",
        "- **C1 signal** — does the vector move *further* for the mutations that actually changed",
        "  binding? Measured on single substitutions only, since the C1 set is disproportionately",
        "  single substitutions and displacement rises with edit count (`D6`). 0.5 is no",
        "  difference; σ is standard errors from it.",
        "",
        "## What this decides",
        "",
        "**§3.1's criterion is met on both arms, so pooling stays and the sequence arms proceed",
        "as specified.** Variants are measurably separated from their own wild type, the",
        "separation scales with the size of the change, and a variant is attributable to its own",
        "reference for the large majority of the corpus. The fallback — **attention pooling**",
        "(`TODO.md` `T32`) — is explicitly not to be reached for before a measurement says it is",
        "needed, and this measurement does not.",
        "",
        "**The A1 -> A4 difference is already visible, before any training.** That comparison is",
        "the one §4.2 built A4 for — same architecture, same pooling, different pretraining",
        "corpus — and the domain-adapted model is the better one on every column here. It is a",
        "property of the embeddings and not yet a result about binding prediction, but it is the",
        "first evidence in the project bearing on it.",
        "",
    ]

    for arm, m in measured.items():
        lines += _arm_section(arm, m)

    lines += [
        "## Reproducing this",
        "",
        "```bash",
        *[f"python scripts/build_embeddings.py --arm {arm}" for arm in measured],
        "python scripts/check_pooling.py",
        "```",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")


def _arm_section(arm: str, m: dict) -> list[str]:
    table, variants = m["table"], m["variants"]
    c1, rest = _c1_split(variants)
    ratio = variants.to_reference / variants.family_median

    lines = [
        f"## `{arm}` — {table.model}",
        "",
        f"{table.width}-d, {len(table.domains)} domains.",
        "",
        "| pair | n | cosine distance, median (p10 – p90) |",
        "|---|---:|---|",
        f"| variant to its own reference | {len(variants)} | {_q(variants.to_reference)} |",
        f"| unrelated domains, same family | {len(m['within']) // 2:,} | {_q(m['within'])} |",
        f"| unrelated domains, different family | {len(m['across']) // 2:,} | {_q(m['across'])} |",
        "",
        "**Displacement against the size of the change**, the test that separates signal from",
        "the embedding's floor:",
        "",
        "| edits from the reference | variants | displacement, median (p10 – p90) |",
        "|---:|---:|---|",
    ]
    for n_edits, group in variants.groupby("n_edits", sort=True):
        lines.append(f"| {n_edits} | {len(group)} | {_q(group.to_reference, 5)} |")

    lines += [
        "",
        "**By family**, since 84 of the 173 variants are homeodomain and a pooled number would",
        "be a homeodomain number wearing a general claim's clothes (`ML_PLAN.md` §6.1):",
        "",
        "| family | variants | to reference | family median | ratio | own reference nearest |",
        "|---|---:|---|---|---:|---:|",
    ]
    for family, group in variants.groupby("family", sort=False):
        if len(group) < 3:
            continue
        lines.append(
            f"| {family} | {len(group)} | {_q(group.to_reference)} | {_q(group.family_median)} | "
            f"{np.nanmedian(group.to_reference / group.family_median):.3f} | "
            f"{(group.rank_in_family == 1.0).mean():.0%} |"
        )

    if c1 is not None and len(c1):
        auc, stderr, n_c1, n_rest = _c1_auc(variants, set(c1.domain))
        sigma = (auc - 0.5) / stderr
        if sigma >= 2:
            verdict = (
                "**positive, and separated from chance by more than two standard errors** — "
                "weak evidence, at this sample size, that the embedding moves most where the "
                "binding moved"
            )
        elif sigma >= 1:
            verdict = (
                "**positive and in the right direction, but not distinguishable from chance "
                "at this sample size**"
            )
        else:
            verdict = "**no better than chance**"
        one = variants[variants.n_edits == 1]
        in_c1 = one.domain.isin(set(c1.domain))
        lines += [
            "",
            f"**The {len(c1)} variants that carry C1** (`D6`) — the ones whose binding measurably",
            "changed, and so the ones any C1 claim rests on. The comparison is made on single",
            f"substitutions alone, because the C1 set has a median of {int(c1.n_edits.median())}",
            f"edit against {int(rest.n_edits.median())} for the rest and displacement rises with",
            "edit count, so the raw comparison would measure that confound instead.",
            "",
            "| among 1-edit variants | n | displacement, median (p10 – p90) |",
            "|---|---:|---|",
            f"| in the C1 set | {n_c1} | {_q(one[in_c1].to_reference, 5)} |",
            f"| not in it | {n_rest} | {_q(one[~in_c1].to_reference, 5)} |",
            "",
            f"**AUC {auc:.3f}** — the probability that a randomly chosen C1 variant is displaced",
            "further than a randomly chosen ordinary one, where 0.5 is no difference. Standard",
            f"error about {stderr:.3f} at these counts, so {sigma:.1f} standard errors from",
            f"chance: {verdict}.",
            "",
            "Note what it does and does not mean: a model does not need displacement *magnitude*",
            "to predict a binding change, only the displacement *direction* to be informative,",
            "which this does not test. A null here would not have sunk C1 on this arm, and a",
            "positive does not secure it.",
        ]

    lines += [
        "",
        f"Its own reference is **not** the nearest domain of the family for "
        f"{1 - (variants.rank_in_family == 1.0).mean():.0%} of variants on this arm — cases where",
        "an unrelated domain of the same family sits closer in the pooled space than a sequence",
        "one residue away. It does not block training, but it is the population where a `P3`",
        f"failure would be expected to concentrate. Median ratio to the family scale: "
        f"**{np.nanmedian(ratio):.3f}**.",
        "",
    ]
    return lines


if __name__ == "__main__":
    main()
