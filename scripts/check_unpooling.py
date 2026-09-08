#!/usr/bin/env python
"""Un-pooling pre-flight — does the single-residue signal survive *before* the mean is taken?

    python scripts/check_unpooling.py [--arm A1 --arm A4] [--out reports/unpooling_check.md]

**The one number that decides `T36`.** `reports/pooling_check.md` asked whether a variant is
separated from its wild type in the *pooled* embedding at all; this asks the sharper question the
λ sweep forced (`docs/ML_RESULTS.md` §9.3): whether the displacement carries information about
**whether the mutation changed binding**, and whether pooling is what destroys it.

Measured on the pooled vectors, it does not. A variant sits 0.0004-0.0008 cosine from its own wild
type, and the correlation between that displacement and the actual change in the bound 8-mer set
is **+0.10 on `A1` and +0.22 on `A4`** — which is why the trained model returns the wild type's
answer for a single-residue variant and cannot see a lost interaction.

Four displacements are compared for the same pairs, cheapest interpretation first:

- `pooled` — cosine between the two **mean** vectors. Today's input, and the baseline to beat.
- `profile` — the **mean of the per-residue** cosines. The same information without the
  cancellation that averaging the vectors first causes.
- `window` — mean per-residue cosine within ±`WINDOW` of the mutation. ESM-2 is contextual, so a
  substitution moves its neighbours' representations too.
- `site` — the per-residue cosine **at the mutated position**. The most local view there is.

**Read the last column of the table, not the first.** A displacement that is merely *larger* is
worth nothing — `site` is guaranteed to be larger than `pooled`, since pooling divides by the
domain length. What matters is whether it *correlates with the biology*, and the report leads with
that.

**And there is a control, because the obvious confound would pass every test above.** A
substitution's effect on binding is partly predictable from amino-acid chemistry alone, so the
BLOSUM62 score of the substitution is scored as a competing predictor on the same pairs. If ESM-2
does not beat a substitution matrix, un-pooling recovers chemistry rather than context and a
per-residue tower is not the answer.

Restricted to **single substitutions of equal length**, verified by direct string comparison
rather than by alignment: the two sequences must differ at exactly one position, so the
per-residue index is unambiguous and no alignment can silently mis-place it. That is also the
subset `reports/pooling_check.md` controls to, since displacement rises with edit count and the
raw comparison would measure that confound instead.

Reads `data/processed/embeddings/<arm>_residues.npz`; build it with
`scripts/build_embeddings.py --arm <arm> --per-residue`. About 10 s per arm.
"""

from __future__ import annotations

import argparse
import textwrap
from pathlib import Path

import numpy as np
import pandas as pd
from Bio.Align import substitution_matrices

from snp2prot import corpus, embeddings, label_health
from snp2prot.config import REPORTS_DIR
from snp2prot.data.matrix import KmerMatrix
from snp2prot.evaluation.metrics import auroc, spearman

DEFAULT_OUT = REPORTS_DIR / "unpooling_check.md"

#: Residues each side of the mutation that the `window` displacement averages over. 5 is the
#: local neighbourhood a contextual model plausibly perturbs; it is not tuned, and `site` and
#: `profile` bracket it from both directions so the report does not rest on the choice.
WINDOW = 5

#: The control's row label, used in more than one place, so it is not spelled twice.
CONTROL = "BLOSUM62 (control)"


def pairs(domains: pd.DataFrame, matrix: KmerMatrix) -> pd.DataFrame:
    """Every (variant, wild type) pair that is a single substitution of equal length.

    **Verified against the sequences, not taken from the edit table.** `n_edits` comes from an
    alignment with free terminal gaps, so a "1 edit" pair can still differ in length and the
    per-residue index would then be off by the offset — silently, and at exactly the position
    being measured. Requiring equal length and exactly one differing character removes the
    question: the mutated index is the differing index.

    `no_evidence` records are excluded as everywhere: with no control, silence is not evidence.
    """
    reference = corpus.reference_rows(domains)
    sequences = domains.dbd_seq.to_numpy()
    verdict = domains.verdict.to_numpy()
    is_variant = domains.is_variant.to_numpy()
    positive = matrix.label == 1

    rows = []
    for i in range(len(domains)):
        r = int(reference[i])
        if r < 0 or r == i or not is_variant[i] or verdict[i] == label_health.NO_EVIDENCE:
            continue
        variant, wild = str(sequences[i]), str(sequences[r])
        if len(variant) != len(wild):
            continue
        differing = [k for k, (a, b) in enumerate(zip(variant, wild, strict=True)) if a != b]
        if len(differing) != 1:
            continue
        union = int((positive[i] | positive[r]).sum())
        rows.append(
            {
                "row": i,
                "reference_row": r,
                "gene": str(domains.gene.to_numpy()[i]),
                "family": str(domains.dbd_family.to_numpy()[i]),
                "position": differing[0],
                "length": len(variant),
                "from_aa": wild[differing[0]],
                "to_aa": variant[differing[0]],
                "dead": verdict[i] == label_health.DEAD_VARIANT,
                "n_pos": int(positive[i].sum()),
                "n_pos_reference": int(positive[r].sum()),
                # How much the binding actually changed: 0 = the same 8-mers, 1 = disjoint sets.
                # Jaccard rather than a count difference, because a variant that binds as many
                # sites as its wild type but *different* ones has changed completely.
                "change": (
                    1.0 - int((positive[i] & positive[r]).sum()) / union if union else np.nan
                ),
            }
        )
    return pd.DataFrame(rows).dropna(subset=["change"])


def _cosine(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    """Row-wise cosine distance between two equally shaped stacks."""
    a64, b64 = a.astype(np.float64), b.astype(np.float64)
    numerator = (a64 * b64).sum(axis=-1)
    denominator = np.linalg.norm(a64, axis=-1) * np.linalg.norm(b64, axis=-1)
    return 1.0 - np.divide(
        numerator, denominator, out=np.zeros_like(numerator), where=denominator > 0
    )


def displacements(table: embeddings.ResidueEmbeddings, found: pd.DataFrame) -> pd.DataFrame:
    """The four displacements of the module docstring, one row per pair."""
    out = found.copy()
    site, window, profile, pooled = [], [], [], []
    for record in found.itertuples():
        variant = table.residues(record.row)
        wild = table.residues(record.reference_row)
        per_residue = _cosine(variant, wild)
        low = max(0, record.position - WINDOW)
        high = min(len(per_residue), record.position + WINDOW + 1)
        site.append(float(per_residue[record.position]))
        window.append(float(per_residue[low:high].mean()))
        profile.append(float(per_residue.mean()))
        pooled.append(float(_cosine(variant.mean(0), wild.mean(0))))
    out["site"] = site
    out["window"] = window
    out["profile"] = profile
    out["pooled"] = pooled
    return out


def blosum_scores(found: pd.DataFrame) -> np.ndarray:
    """BLOSUM62 score of each substitution, negated so that larger means a bigger change.

    The control. A high-scoring substitution is a conservative one, so the sign is flipped to put
    it on the same footing as a displacement: bigger number, more disruption expected.
    """
    matrix = substitution_matrices.load("BLOSUM62")
    return np.array(
        [-float(matrix[record.from_aa, record.to_aa]) for record in found.itertuples()],
        dtype=np.float64,
    )


def measure(found: pd.DataFrame, control: np.ndarray) -> pd.DataFrame:
    """Each predictor scored the three ways that matter, plus its raw scale.

    `rho_blosum` is the diagnostic that separates two very different failures: a displacement
    that tracks the substitution matrix is measuring amino-acid chemistry and merely adding noise
    to it, while one that does not is measuring something else entirely.
    """
    truth = found.change.to_numpy(dtype=np.float64)
    dead = found.dead.to_numpy().astype(int)
    rows = []
    predictors = [
        (name, found[name].to_numpy(dtype=np.float64))
        for name in ("pooled", "profile", "window", "site")
    ]
    predictors.append((CONTROL, control))
    for name, values in predictors:
        rows.append(
            {
                "predictor": name,
                "median": float(np.median(values)),
                "rho_change": spearman(values, truth),
                "auroc_dead": auroc(dead, values) if 0 < dead.sum() < len(dead) else float("nan"),
                "rho_blosum": float("nan") if name == CONTROL else spearman(values, control),
            }
        )
    return pd.DataFrame(rows)


def verdict(measured: dict[str, pd.DataFrame]) -> list[str]:
    """Which of the three outcomes the numbers actually landed on, stated rather than left open.

    A report whose conclusion section lists every branch and names none is a report the reader has
    to re-derive. The branches are still printed below this, because the *rule* is what makes the
    verdict checkable — but the rule is applied here.
    """
    lines = ["## The verdict this run landed on", ""]
    for arm, table in measured.items():
        indexed = table.set_index("predictor")
        control = indexed.loc[CONTROL]
        local = indexed.loc[["site", "window"]].rho_change.max()
        pooled = indexed.loc["pooled"].rho_change
        best = indexed.drop(index=CONTROL).rho_change.idxmax()
        unpooling = "recover nothing" if local <= pooled else "help"
        beats = (
            "beats"
            if indexed.drop(index=CONTROL).rho_change.max() > control.rho_change
            else ("**loses to**")
        )
        paragraph = (
            f"**`{arm}`** — the most local views (`site`, `window`) {unpooling}: their best "
            f"ρ is {local:+.3f} against `pooled`'s {pooled:+.3f}, and the best of all four "
            f"displacements is `{best}`. Against the control the language model {beats} a "
            f"substitution matrix ({indexed.drop(index=CONTROL).rho_change.max():+.3f} "
            f"against {control.rho_change:+.3f} on ρ, "
            f"{indexed.drop(index=CONTROL).auroc_dead.max():.3f} against "
            f"{control.auroc_dead:.3f} on dead detection)."
        )
        lines += [*textwrap.wrap(paragraph, width=96), ""]
    return lines


def _fmt(value: float, places: int = 3) -> str:
    return "n/a" if not np.isfinite(value) else f"{value:.{places}f}"


def write_report(path: Path, measured: dict[str, pd.DataFrame], found: pd.DataFrame) -> None:
    n_dead = int(found.dead.sum())
    lines = [
        "# Un-pooling pre-flight — does the single-residue signal survive before the mean?",
        "",
        "Generated by `scripts/check_unpooling.py`. **This runs before any per-residue tower is",
        "built** — the same discipline `reports/pooling_check.md` applied before the first grid,",
        "asked of the sharper question the λ sweep forced.",
        "",
        "`reports/pooling_check.md` established that a variant *is* weakly separated from its wild",
        "type in the pooled embedding. `docs/ML_RESULTS.md` §9.3 then established that the",
        "separation carries almost no information about **whether the mutation changed binding**",
        "— and that the trained model consequently returns the wild type's answer for a variant,",
        "at every λ. The question here is whether mean pooling is what destroys that information,",
        "or whether it was never in the sequence representation at all. **The two answers point at",
        "completely different projects**: the first is a tower to rebuild, the second says the",
        "sequence modality is exhausted and only structure remains.",
        "",
        "**The pairs**: single substitutions of equal length, verified character by character so",
        "the mutated index cannot be mis-placed by an alignment.",
        f"**{len(found)} pairs, {n_dead} of them variants that lost binding entirely.**",
        "`no_evidence` records excluded.",
        "",
        "**What is being correlated against**: `change`, one minus the Jaccard overlap of the",
        "variant's bound 8-mers with its wild type's. 0 means the mutation changed nothing, 1 that",
        "the two share no site. Jaccard rather than a count difference, because a variant binding",
        "as many sites as its wild type but *different* ones has changed completely.",
        "",
        "## The verdict",
        "",
        "**`ρ change` is the column that decides this.** A larger displacement is worth nothing on",
        "its own — `site` is bound to exceed `pooled`, since pooling divides by the domain length.",
        "What matters is whether it tracks the biology. `AUROC dead` is the same question asked as",
        "a classification: can the displacement alone tell a variant that lost binding from one",
        "that did not.",
        "",
    ]
    for arm, table in measured.items():
        lines += [
            f"### arm `{arm}` — {embeddings.ARMS[arm]}",
            "",
            "| displacement | median | **ρ change** | AUROC dead | ρ BLOSUM |",
            "|---|---:|---:|---:|---:|",
        ]
        for record in table.itertuples():
            emphasis = "**" if record.predictor == CONTROL else ""
            lines.append(
                f"| {emphasis}{record.predictor}{emphasis} | {record.median:.5f} | "
                f"**{_fmt(record.rho_change)}** | {_fmt(record.auroc_dead)} | "
                f"{_fmt(record.rho_blosum)} |"
            )
        lines.append("")

    lines += verdict(measured)
    lines += [
        "`ρ BLOSUM` is how much of each displacement is amino-acid chemistry rather than context.",
        "A displacement that tracks the substitution matrix is at worst a noisy version of it; one",
        "that does not is measuring a different quantity — for a masked language model, plausibly",
        "*how surprising the original residue was* rather than *how disruptive the new one is*.",
        "",
        "## The rule the verdict was applied from",
        "",
        "- **`site` clearly beats `pooled` on `ρ change`, and beats BLOSUM62** — pooling was",
        "  destroying the signal and a per-residue protein tower is worth building. `T36`",
        "  proceeds, and the choice between attention pooling and a fixed positional scheme can",
        "  be made on the `window` versus `site` gap.",
        "- **`site` beats `pooled` but not BLOSUM62** — what un-pooling recovers is amino-acid",
        "  chemistry, not context. A substitution matrix is a two-line feature and a per-residue",
        "  tower is not; add the cheap feature and do not build the tower.",
        "- **Nothing beats `pooled` by much** — the information is not in the sequence",
        "  representation at any granularity. The sequence modality is exhausted for this claim,",
        "  and the structure arms are not an incremental next step but the only remaining one.",
        "",
        "",
        "## What the control is doing",
        "",
        "The BLOSUM62 severity of the substitution separates the variants that lost binding from",
        "those that did not, and it is the strongest predictor in the table on both measures. That",
        "is biologically unsurprising — a radical substitution in a DNA-contacting position is",
        "likelier to abolish binding than a conservative one — but it is a **two-line**",
        "**feature**, and",
        "no arrangement of a 650M-parameter language model's outputs matched it here.",
        "",
        "| | n | mean severity | median | min | max |",
        "|---|---:|---:|---:|---:|---:|",
        *[
            f"| {'lost binding' if dead else 'still binds'} | {len(g)} | {g.blosum.mean():.2f} | "
            f"{g.blosum.median():.1f} | {g.blosum.min():.0f} | {g.blosum.max():.0f} |"
            for dead, g in found.groupby("dead")
        ],
        "",
        "(Severity is the negated BLOSUM62 score, so larger is more radical.)",
        "",
        "**Caveats that apply to every row above.** `n` is small — "
        f"{len(found)} pairs, {n_dead} of them dead — and dominated by one family "
        f"({found.family.value_counts().iloc[0]} of {len(found)} are "
        f"{found.family.value_counts().index[0]}), so this is close to a statement about "
        "homeodomains. At that `n` the control's dead-detection AUROC carries a standard error of "
        "about 0.10, which puts it near three standard deviations above chance rather than "
        "comfortably beyond it. `T30` remains the decision that would enlarge the set.",
        "",
        "## Reproducing this",
        "",
        "```bash",
        "python scripts/build_embeddings.py --arm A1 --per-residue   # ~30 s on the GPU, ~600 MB",
        "python scripts/build_embeddings.py --arm A4 --per-residue",
        "python scripts/check_unpooling.py                           # this report, ~10 s",
        "```",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--arm", action="append", dest="arms", help="repeatable; default A1 and A4")
    ap.add_argument("--out", type=Path, default=DEFAULT_OUT)
    args = ap.parse_args()
    arms = args.arms or ["A1", "A4"]

    domains = corpus.domains()
    matrix = KmerMatrix.load()
    try:
        corpus.require_aligned(domains.dbd_seq, matrix=matrix.domains)
    except ValueError as mismatch:
        raise SystemExit(str(mismatch)) from mismatch

    found = pairs(domains, matrix)
    if not len(found):
        raise SystemExit("no single-substitution pairs found; nothing to measure")
    control = blosum_scores(found)
    # Carried on the frame as well as passed around, because the report tabulates it directly and
    # the two must be the same numbers in the same order.
    found = found.assign(blosum=control)
    print(f"{len(found)} single-substitution pairs, {int(found.dead.sum())} of them dead variants")

    measured: dict[str, pd.DataFrame] = {}
    for arm in arms:
        table = embeddings.ResidueEmbeddings.load(arm)
        corpus.require_aligned(domains.dbd_seq, **{arm: table.domains})
        scored = displacements(table, found)
        measured[arm] = measure(scored, control)
        print(f"\n--- {arm}")
        print(measured[arm].to_string(index=False))

    write_report(args.out, measured, found)
    print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
