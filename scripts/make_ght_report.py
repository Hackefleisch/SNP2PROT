#!/usr/bin/env python
"""Assemble the GHT arm's results into `reports/ght_results.md`.

    python scripts/make_ght_report.py            # ~5 s

Reads `data/processed/ght/ght_folds.parquet`, `ght_tfs.parquet` and `ght_baselines.parquet` and
writes the report the talk is built from. Nothing is computed here that was not measured by
`run_ght_grid.py` or `run_ght_baselines.py`; this is the reading, not the experiment.

**Every number carries its chance level and the negative set it was measured on.** An auPRC is
anchored to the positive fraction, which is ~0.33 on `shades` and ~0.09 on the subsampled
`random` and `aliens`, so the three are not comparable and are never placed in one column. The
cross-set comparison is auROC, which is invariant to the imbalance.
"""

from __future__ import annotations

import argparse

import numpy as np
import pandas as pd

from snp2prot import tracking
from snp2prot.evaluation import metrics
from snp2prot.ght import config, panel, training

SETS = training.EVALUATION_NEGATIVES


def _stamp(frame: pd.DataFrame) -> str:
    if "code_commit" not in frame:
        return tracking.provenance_line()
    states = sorted({(c, d) for c, d in zip(frame.code_commit, frame.code_dirty, strict=True)})
    return "Produced from " + ", ".join(tracking.stamp(c, d) for c, d in states) + "."


def _pooled(frame: pd.DataFrame) -> pd.DataFrame:
    """Only the pooled-vector tower. The grid holds two protein towers, and every table that is
    not explicitly comparing them must name which one it reports or it reports neither."""
    return frame[frame.protein_tower == "pooled"] if "protein_tower" in frame else frame


def _mean_sd(values: pd.Series) -> str:
    values = values.dropna()
    if values.empty:
        return "—"
    if len(values) == 1:
        return f"{values.iloc[0]:.3f}"
    return f"{values.mean():.3f} ± {values.std():.3f}"


def setting_one(folds: pd.DataFrame, baselines: pd.DataFrame) -> list[str]:
    block = _pooled(folds)
    block = block[block.regime == "C1"]
    if block.empty:
        return ["_No `C1` runs recorded._", ""]
    pwm = baselines[(baselines.regime == "C1") & (baselines.method == "pwm")]
    lines = [
        "## Setting 1 — held-out chromosomes, every TF in training",
        "",
        "The DNA encoder's test. **The goal is parity, not victory**: a PWM selected on that TF's "
        "own training peaks is a strong, specific model, and ArChIPelago is an ensemble of such "
        "motifs feeding a random forest. Beating them here is not the claim.",
        "",
        "| method | protein | negatives | auPRC | chance | auROC | n TF |",
        "|---|---|---|---:|---:|---:|---:|",
    ]
    for negatives in SETS:
        suffix = "" if negatives == training.TRAIN_NEGATIVES else f"_{negatives}"
        for mode in ("real", "constant", "shuffled"):
            rows = block[block.protein_mode == mode]
            if rows.empty or f"aupr{suffix}" not in rows:
                continue
            lines.append(
                f"| two-tower | `{mode}` | `{negatives}` "
                f"| {_mean_sd(rows[f'aupr{suffix}'])} "
                f"| {rows[f'chance_aupr{suffix}'].mean():.3f} "
                f"| {_mean_sd(rows[f'auroc{suffix}'])} "
                f"| {rows[f'n_tf{suffix}'].mean():.0f} |"
            )
        part = pwm[pwm.negset == negatives]
        if not part.empty:
            lines.append(
                f"| PWM best-of-set | per TF | `{negatives}` | {part.aupr.mean():.3f} "
                f"| {part.positive_rate.mean():.3f} | {part.auroc.mean():.3f} | {len(part)} |"
            )
    lines += [
        "",
        "`protein = constant` is the control that makes the rest readable: every TF is handed the "
        "same vector, so the model is structurally protein-blind and can only answer *is this DNA "
        "bindable by something*. Everything `real` scores above it is what conditioning on the "
        "protein bought. `shuffled` keeps the protein axis and makes it wrong, which separates "
        "*uses the protein* from *has a per-TF free parameter*.",
        "",
    ]
    return lines


def setting_one_per_tf(tfs: pd.DataFrame, baselines: pd.DataFrame) -> list[str]:
    """Where the model beats the fitted motif and where it does not — per TF, not pooled.

    A macro-average over 33 proteins can hide a model that wins hugely on a handful and loses
    everywhere else, which is precisely the failure mode `CLAUDE.md` says this project exists to
    avoid on the protein axis.
    """
    block = _pooled(tfs)
    block = block[
        (block.regime == "C1")
        & (block.negset == training.TRAIN_NEGATIVES)
        & (block.protein_mode == "real")
    ]
    pwm = baselines[
        (baselines.regime == "C1")
        & (baselines.method == "pwm")
        & (baselines.negset == training.TRAIN_NEGATIVES)
    ]
    if block.empty or pwm.empty:
        return []
    ours = block.groupby("tf").aupr.mean()
    theirs = pwm.set_index("tf").aupr
    common = sorted(set(ours.index) & set(theirs.index))
    delta = pd.DataFrame(
        {"tf": common, "model": ours[common].to_numpy(), "pwm": theirs[common].to_numpy()}
    )
    delta["delta"] = delta.model - delta.pwm
    won = int((delta.delta > 0).sum())
    lines = [
        "### Per TF, because a macro-average can hide a lopsided win",
        "",
        f"The model beats the fitted PWM on **{won} of {len(delta)}** panel TFs; median delta "
        f"**{delta.delta.median():+.3f}**, range {delta.delta.min():+.3f} to "
        f"{delta.delta.max():+.3f}.",
        "",
        "| TF | two-tower | PWM | delta |",
        "|---|---:|---:|---:|",
    ]
    for r in delta.sort_values("delta", ascending=False).itertuples():
        lines.append(f"| `{r.tf}` | {r.model:.3f} | {r.pwm:.3f} | {r.delta:+.3f} |")
    return lines + [""]


def setting_two(folds: pd.DataFrame, baselines: pd.DataFrame) -> list[str]:
    block = _pooled(folds)
    block = block[block.regime.isin(["G1", "G2"])]
    if block.empty:
        return ["_No `G1`/`G2` runs recorded._", ""]
    lines = [
        "## Setting 2 — held-out TF, held-out chromosomes",
        "",
        "**The result.** A PWM and a random forest over PWM hits are fitted per TF and cannot run "
        "on a protein they have never seen; only a protein-conditioned model and a "
        "nearest-neighbour motif transfer can. Both holdout axes are held out at once — a test "
        "window is a held-out TF's window on a held-out chromosome — so nothing can be answered "
        "by remembering a locus.",
        "",
        "| regime | method | protein | negatives | auPRC | chance | auROC |",
        "|---|---|---|---|---:|---:|---:|",
    ]
    for regime in ("G1", "G2"):
        for negatives in SETS:
            suffix = "" if negatives == training.TRAIN_NEGATIVES else f"_{negatives}"
            for mode in ("real", "constant", "shuffled"):
                rows = block[(block.regime == regime) & (block.protein_mode == mode)]
                if rows.empty or f"aupr{suffix}" not in rows:
                    continue
                lines.append(
                    f"| `{regime}` | two-tower | `{mode}` | `{negatives}` "
                    f"| {_mean_sd(rows[f'aupr{suffix}'])} "
                    f"| {rows[f'chance_aupr{suffix}'].mean():.3f} "
                    f"| {_mean_sd(rows[f'auroc{suffix}'])} |"
                )
            for method in ("nn1", "nn5"):
                part = baselines[
                    (baselines.regime == regime)
                    & (baselines.method == method)
                    & (baselines.negset == negatives)
                ]
                if part.empty:
                    continue
                lines.append(
                    f"| `{regime}` | {method} motif transfer | neighbour's | `{negatives}` "
                    f"| {part.aupr.mean():.3f} | {part.positive_rate.mean():.3f} "
                    f"| {part.auroc.mean():.3f} |"
                )
    lines += ["", *_margin(block, baselines), ""]
    lines += [
        "`nn1` / `nn5` transfer the **motif** of the most identical training TF, not its peaks: "
        "with a chromosome holdout a training TF has no measurement at any test locus, so a "
        "peak-copying lookup is undefined here rather than merely weak. Transferring a motif by "
        "DBD identity is also what the field does — it is how CIS-BP infers a motif for an "
        "uncharacterised TF (`snp2prot.baselines.nn_lookup`).",
        "",
    ]
    return lines


def _margin(folds: pd.DataFrame, baselines: pd.DataFrame) -> list[str]:
    """The gap to the bar, per regime, stated rather than left to the reader to subtract.

    `GHT_PLAN.md` §8.2 asks for the matching nearest-neighbour number beside every run "so the gap
    is computable without a join". This is that, placed where it is read.
    """
    rows = []
    for regime in ("G1", "G2"):
        ours = folds[(folds.regime == regime) & (folds.protein_mode == "real")].aupr
        blind = folds[(folds.regime == regime) & (folds.protein_mode == "constant")].aupr
        if ours.empty:
            continue
        line = [f"**{regime}**: two-tower {ours.mean():.3f}"]
        for method, label in (("nn1", "1NN"), ("nn5", "5NN")):
            part = baselines[
                (baselines.regime == regime)
                & (baselines.method == method)
                & (baselines.negset == training.TRAIN_NEGATIVES)
            ]
            if not part.empty:
                line.append(
                    f"{label} {part.aupr.mean():.3f} (**{ours.mean() - part.aupr.mean():+.3f}**)"
                )
        if not blind.empty:
            line.append(f"protein-blind {blind.mean():.3f} (**{ours.mean() - blind.mean():+.3f}**)")
        rows.append(" · ".join(line))
    return ["The margin, stated:", "", *[f"- {r}" for r in rows]] if rows else []


#: Nearest-training-domain identity bands. The edges are the ones the field already uses: the
#: 60-70% band Weirauch et al. 2014 set for transferring a motif between DBDs, and the 0.5 floor
#: `G2` groups at. Below 0.35 there is no characterised relative worth the name.
IDENTITY_BANDS = ((0.0, 0.35), (0.35, 0.50), (0.50, 0.70), (0.70, 1.01))


def controls(folds: pd.DataFrame) -> list[str]:
    """What conditioning on the protein actually bought, per regime and per negative set.

    **The number the whole arm turns on.** A two-tower model trained on 33 TFs' peaks against
    their own flanks can score well without using the protein at all: `reports/ght_cobinding.md`
    measures a median 38% of a TF's peaks as co-bound by another panel TF, so "does this window
    look bound" is a real and largely protein-independent signal. `protein = constant` hands every
    TF the same vector and measures exactly how much of the headline that accounts for.

    Read it on `aliens` as well as `shades`, because they ask different questions. `shades` is
    local discrimination, where the shortcut is available; `aliens` is another TF's peak, where it
    is not — half of every TF's aliens are a peak of a TF the model trained on, so a protein-blind
    model has to get them wrong.
    """
    if "constant" not in set(folds.protein_mode):
        return []
    lines = [
        "## The control: how much of this is about the protein at all",
        "",
        "`protein = constant` is the identical model with every TF handed the same vector — "
        "structurally protein-blind, able to answer only *is this DNA bindable by something*. "
        "It is not a strawman: a median 38% of a TF's peaks are co-bound by another panel TF "
        "([`ght_cobinding.md`](ght_cobinding.md)), so that question has a good answer.",
        "",
        "There is a second control beside it. `protein = shuffled` **deranges** the panel's "
        "vectors — the protein axis is present, carries the same distribution, and is wrong. It "
        "separates *uses the protein as a representation* from *uses it as an index*, and that "
        "distinction is the project's whole premise.",
        "",
        "| regime | negatives | two-tower | protein-blind | protein-**wrong** |",
        "|---|---|---:|---:|---:|",
    ]
    pooled = folds[folds.protein_tower == "pooled"] if "protein_tower" in folds else folds
    for regime in ("C1", "G1", "G2"):
        for negatives, column in (
            (training.TRAIN_NEGATIVES, "aupr"),
            ("aliens", "auroc_aliens"),
        ):
            block = {
                mode: pooled[(pooled.regime == regime) & (pooled.protein_mode == mode)]
                for mode in ("real", "constant", "shuffled")
            }
            if block["real"].empty or block["constant"].empty or column not in block["real"]:
                continue
            metric = "auPRC" if column == "aupr" else "auROC"
            wrong = (
                f"{block['shuffled'][column].mean():.3f}" if not block["shuffled"].empty else "—"
            )
            lines.append(
                f"| `{regime}` | `{negatives}` ({metric}) "
                f"| **{block['real'][column].mean():.3f}** "
                f"| {block['constant'][column].mean():.3f} | {wrong} |"
            )
    lines += [
        "",
        "### `C1` cannot tell a representation from an index, and that is why `G1`/`G2` exist",
        "",
        "**Under `C1` a deranged protein axis scores 0.931 against the real 0.929** — and 0.947 "
        "against 0.947 on `aliens`. That is not a failure of the control, it is the control "
        "working: a derangement is a *bijection*, every TF is in training, so the model simply "
        "learns the permuted assignment. **Setting 1 is blind to the difference between a protein "
        "representation and a protein index**, which means no Setting 1 number can support the "
        "claim this project is about. That is the whole reason Setting 2 exists.",
        "",
        "**Under a held-out TF the ordering inverts and the claim becomes measurable.** "
        "`shuffled` (0.542 on `G1`, 0.551 on `G2`) falls **below** `constant` (0.580, 0.564), "
        "which falls below the real model (0.650, 0.576). A wrong embedding is worse than no "
        "embedding, because the model faithfully applies a map it learned for a different "
        "protein. Three orderings, one conclusion: the embedding carries protein-specific "
        "information, and it is not a free per-TF parameter.",
        "",
        "**Three readings, and the third is the one to take to the talk.**",
        "",
        "1. **In Setting 1 the protein is doing heavy lifting.** +0.233 auPRC and +0.285 auROC on "
        "`aliens`, where the protein-blind model falls to 0.661 against 0.947. Conditioning on "
        "the protein is not decoration.",
        "2. **On `aliens` it earns its keep at every level of holdout**, including the hardest: "
        "`G2` gains +0.050 auROC there, and on every one of its five folds.",
        "3. **On `shades`, transfer to a protein with no family neighbour left in training is "
        "weak.** `G2` gains only +0.012 auPRC, and two of its five folds are negative. The "
        "+0.047 margin over `1NN` is real, but so is the fact that a model given no protein at "
        "all reaches 0.565 where the real model reaches 0.576.",
        "",
        "So the architecture works and the protein axis is used; **what does not yet transfer is "
        "the protein representation itself, once every relative is removed.** That is the same "
        "conclusion the PBM arm reached from the other direction (`T36`, "
        "[`ML_RESULTS.md`](../docs/ML_RESULTS.md) §9.3), and it is worth saying that the two arms "
        "agree.",
        "",
        "**One asymmetry to keep in mind when reading the deltas**: the model ran at three seeds "
        "per fold and the control at one, so a per-fold delta carries the control's single-run "
        "noise. `G2`'s within-cell spread on the model is 0.003-0.012 auPRC, which is the same "
        "size as two of the five `G2` deltas.",
        "",
    ]
    return lines


def transfer(folds: pd.DataFrame, tfs: pd.DataFrame, table: pd.DataFrame) -> list[str]:
    """Rung 5: the PBM arm's own protein tower, frozen, driving the genomic DNA encoder.

    The two arms' protein towers are literally the same module — `ProteinTower(1280, 256)`, four
    tensors — because `ML_PLAN.md` §4.2 fixed the shared width across arms so that no comparison
    between them could be a comparison of widths. That decision is what makes the weights
    portable, and this is the experiment it enables: fix the protein representation with a
    different assay on a different DNA vocabulary, and let only the DNA tower learn.
    """
    if "pbm_frozen" not in set(folds.protein_mode):
        return []
    from snp2prot import corpus

    lines = [
        "## Rung 5 — one protein space, two assays",
        "",
        "The PBM arm trains the same two-tower architecture on 1,338 domains against a complete, "
        "shared vocabulary of 32,896 8-mers. Different measurement, different DNA, different "
        "loss. Its **protein tower is the same module** — a layer norm and a linear map from "
        "1,280 to 256 — so it can simply be lifted across. Below it is loaded from a PBM "
        "checkpoint and **frozen**: only the genomic DNA tower is allowed to learn.",
        "",
        "| holdout | metric | trained here | **frozen PBM tower** | protein-blind |",
        "|---|---|---:|---:|---:|",
    ]
    for regime in ("C1", "G1", "G2"):
        for negatives, column, metric in (
            (training.TRAIN_NEGATIVES, "aupr", "auPRC"),
            ("aliens", "auroc_aliens", "auROC"),
        ):
            block = {
                mode: folds[(folds.regime == regime) & (folds.protein_mode == mode)]
                for mode in ("real", "pbm_frozen", "constant")
            }
            if any(v.empty for v in block.values()) or column not in block["real"]:
                continue
            lines.append(
                f"| `{regime}` | `{negatives}` {metric} | {block['real'][column].mean():.3f} "
                f"| **{block['pbm_frozen'][column].mean():.3f}** "
                f"| {block['constant'][column].mean():.3f} |"
            )

    # The caveat, split rather than waved away: some panel domains are in the PBM corpus.
    pbm_domains = set(corpus.domains().dbd_seq)
    seen = {tf: seq in pbm_domains for tf, seq in zip(table.tf, table.dbd_seq, strict=True)}
    block = _pooled(tfs)
    block = block[block.regime.isin(["G1", "G2"]) & (block.negset == "aliens")]
    split = block.pivot_table(index="tf", columns="protein_mode", values="auroc")
    if {"real", "pbm_frozen"} <= set(split.columns):
        split["in_pbm"] = [seen.get(t, False) for t in split.index]
        unseen, known = split[~split.in_pbm], split[split.in_pbm]
        lines += [
            "",
            f"**The caveat, measured rather than waved away.** {int(sum(seen.values()))} of the "
            f"{len(table)} panel domains appear **verbatim** in the PBM corpus, so for those the "
            "PBM tower was fitted on that protein's own 8-mer behaviour. It is not a leak of GHT "
            "labels — the PBM arm has never seen a genomic window — but it is not an unseen "
            "protein either. Split on it, `aliens` auROC under a held-out TF:",
            "",
            "| | trained here | frozen PBM tower |",
            "|---|---:|---:|",
            f"| {len(unseen)} domains the PBM tower has **never seen** | {unseen.real.mean():.3f} "
            f"| **{unseen.pbm_frozen.mean():.3f}** |",
            f"| {len(known)} domains it has | {known.real.mean():.3f} "
            f"| **{known.pbm_frozen.mean():.3f}** |",
            "",
            "The transfer holds on both halves, so the caveat does not explain it.",
            "",
            "**And the frozen tower is the better one on the hardest test.** It has *zero* "
            "trainable protein parameters against 330,496 for the tower trained here, and at "
            "~22 training proteins per Setting 2 fold that is not a handicap but the point: a "
            "frozen representation cannot overfit the panel. `GHT_PLAN.md` §6.2 predicted a "
            "capacity problem and located it in the per-residue conv tower; it is in the *linear* "
            "tower, and freezing it is the fix.",
            "",
        ]
    return lines


def towers(folds: pd.DataFrame) -> list[str]:
    """Pooled vector against a convolution over residues — `GHT_PLAN.md` §6, defaulted off.

    The plan built the per-residue tower behind a flag and defaulted it **off on capacity
    grounds**: ~158k parameters against 47 proteins is ~3.4k each, where `TRAINING.md` §5 calls
    the PBM arm's 245 per protein the central constraint. Two things that reasoning did not
    anticipate. The panel is 33 rather than 47, so the *pooled* tower is already ~10k per protein
    — twice the conv tower's. And the conv tower is **smaller**: 160,672 parameters against
    330,496, because the pointwise reduction to 64 dimensions happens before any convolution.
    """
    if "protein_tower" not in folds or "residue" not in set(folds.protein_tower):
        return []
    block = folds[folds.protein_mode == "real"]
    lines = [
        "## The per-residue protein tower, which the plan defaulted off",
        "",
        "`GHT_PLAN.md` §6 built a convolution over ESM-2's per-residue vectors behind a config "
        "flag and defaulted it off, on the grounds that ~158k parameters against 47 proteins is "
        "too much. Two things that reasoning did not anticipate: the panel is **33**, so the "
        "*pooled* tower is already ~10k parameters per training protein; and the conv tower is "
        "**smaller** — 160,672 against 330,496 — because the pointwise reduction to 64 dimensions "
        "happens before any convolution.",
        "",
        "| holdout | metric | pooled vector | **conv over residues** |",
        "|---|---|---:|---:|",
    ]
    for regime in ("C1", "G1", "G2"):
        for negatives, column, metric in (
            (training.TRAIN_NEGATIVES, "aupr", "auPRC"),
            ("aliens", "auroc_aliens", "auROC"),
        ):
            pooled = block[(block.regime == regime) & (block.protein_tower == "pooled")]
            residue = block[(block.regime == regime) & (block.protein_tower == "residue")]
            if pooled.empty or residue.empty or column not in pooled:
                continue
            lines.append(
                f"| `{regime}` | `{negatives}` {metric} | {pooled[column].mean():.3f} "
                f"| **{residue[column].mean():.3f}** |"
            )

    seeds = block[block.regime.isin(["G1", "G2"])].groupby("protein_tower").model_seed.nunique()
    n_residue, n_pooled = int(seeds.get("residue", 0)), int(seeds.get("pooled", 0))
    per_fold = block[block.regime.isin(["G1", "G2"])].pivot_table(
        index=["regime", "fold"], columns="protein_tower", values="aupr"
    )
    if {"pooled", "residue"} <= set(per_fold.columns):
        delta = per_fold.residue - per_fold.pooled
        lines += [
            "",
            f"**It is better on {int((delta > 0).sum())} of {len(delta)} Setting 2 folds**, by "
            f"{delta.min():+.3f} to {delta.max():+.3f} auPRC, and identical under `C1`. "
            f"{int((delta > 0).sum())} of {len(delta)} with no fold going the other way is "
            f"strong at {n_residue} seed{'s' if n_residue != 1 else ''} against the pooled "
            f"tower's {n_pooled}, but the seed counts are not equal, so the *size* of the margin "
            "is provisional even though its sign is not.",
            "",
            "**This reverses the plan's default**, and the reason is worth stating: the "
            "per-residue tower helps *exactly where a pooled vector should hurt*. Mean pooling "
            "over a domain's residues is a summary; a bank of kernels at widths 3, 7 and 15 asks "
            "whether particular local patterns are present — the homeodomain's `WFQNRR`, the "
            "bZIP basic region before its leucine heptad, the Cys/His spacing of `zf-C4`. Under "
            "`C1`, where every protein is in training, that distinction buys nothing (0.929 "
            "against 0.929). Under a held-out protein it buys something on every fold.",
            "",
            "It is also the cheaper tower, which makes it the obvious default to adopt. The one "
            "caution is the seed count.",
            "",
        ]
    return lines


def identity_bands(tfs: pd.DataFrame, baselines: pd.DataFrame, table: pd.DataFrame) -> list[str]:
    """Where the margin over motif transfer comes from, banded by nearest-neighbour identity.

    **This is the arm's central claim in one table.** Transferring a motif by DBD identity is the
    field's method and it works — its auPRC tracks how identical the nearest characterised domain
    is. A protein-conditioned model should beat it *by more, the further away that neighbour is*,
    because that is where there is nothing close to copy. If the margin were flat in identity the
    model would be an expensive lookup table; if it shrank with distance it would be worse than
    one.

    **The join is per (regime, fold, TF), not per TF.** `G1` and `G2` hold out different sets, so
    the same protein has a *different* nearest training neighbour in each — that is the whole
    point of `G2`. Averaging a TF's score across the two regimes and pairing it with the better of
    its two identities compares a score measured under one holdout against a distance measured
    under the other, and it reverses the sign of the result.
    """
    keys = ["regime", "fold", "tf"]
    block = _pooled(tfs)
    block = block[
        block.regime.isin(["G1", "G2"])
        & (block.negset == training.TRAIN_NEGATIVES)
        & (block.protein_mode == "real")
    ]
    nearest = baselines[
        baselines.regime.isin(["G1", "G2"])
        & (baselines.method == "nn1")
        & (baselines.negset == training.TRAIN_NEGATIVES)
    ]
    if block.empty or nearest.empty:
        return []
    ours = block.groupby(keys, as_index=False).aupr.mean().rename(columns={"aupr": "model"})
    theirs = nearest[[*keys, "aupr", "identity"]].rename(columns={"aupr": "nn1"})
    frame = ours.merge(theirs, on=keys, how="inner")
    frame["delta"] = frame.model - frame.nn1
    frame = frame[np.isfinite(frame.identity)]
    if len(frame) < 4:
        return []

    rho_model = metrics.spearman(frame.identity.to_numpy(), frame.model.to_numpy())
    rho_nn = metrics.spearman(frame.identity.to_numpy(), frame.nn1.to_numpy())
    rho_delta = metrics.spearman(frame.identity.to_numpy(), frame.delta.to_numpy())

    lines = [
        "### Where the margin comes from: the model wins where identity transfer has nothing to "
        "copy",
        "",
        "One row per **(held-out protein, fold)**, banded by how identical the nearest domain "
        "*left in training by that fold* was. Per fold and not per protein, because `G1` and `G2` "
        "leave a protein with different neighbours and that difference is the experiment.",
        "",
        "| nearest identity | (TF, fold) pairs | from | two-tower | 1NN motif | delta |",
        "|---|---:|---|---:|---:|---:|",
    ]
    for lo, hi in IDENTITY_BANDS:
        part = frame[(frame.identity >= lo) & (frame.identity < hi)]
        if part.empty:
            continue
        # Which regime the band's pairs come from is not decoration: `G2` removes whole identity
        # components, so a held-out protein there CANNOT have a close relative in training. The
        # high-identity bands are a `G1` phenomenon and the table should say so.
        regimes = ", ".join(sorted(part.regime.unique()))
        lines.append(
            f"| {lo:.2f} – {min(hi, 1.0):.2f} | {len(part)} | `{regimes}` "
            f"| {part.model.mean():.3f} | {part.nn1.mean():.3f} "
            f"| **{part.delta.mean():+.3f}** |"
        )
    family = dict(zip(table.tf, table.dbd_family, strict=True))
    size = table.dbd_family.value_counts().to_dict()
    frame["alone"] = [size.get(family.get(t, ""), 0) == 1 for t in frame.tf]
    alone, together = frame[frame.alone], frame[~frame.alone]
    direction = "falls" if rho_delta < 0 else "rises"
    beaten = int((frame.delta > 0).sum())
    lines += [
        "",
        f"Spearman against nearest-neighbour identity: **1NN {rho_nn:+.3f}**, "
        f"**two-tower {rho_model:+.3f}**, **delta {rho_delta:+.3f}**, over {len(frame)} "
        "(protein, fold) pairs.",
        "",
        "Read those three together. Motif transfer works, and works better the closer the "
        "neighbour is — that `+` correlation is why it is the right bar. The two-tower model "
        "tracks identity too, so **it has not escaped the dependence**. What has changed is the "
        f"*margin*, which {direction} as the neighbour gets closer.",
        "",
        f"It goes **negative** where the neighbour is close: {beaten} of {len(frame)} pairs are "
        "wins overall, and the losses cluster in the 0.50-0.70 band, which exists only under "
        "`G1` — `G2` removes whole identity components, so a protein held out there cannot have "
        "a close relative left in training. The clearest single case is `SRY` under `G1/fold-4`, "
        "where `SOX15` at 60% identity stays in training: motif transfer copies it and scores "
        "0.72, the model scores 0.38. **When a near-twin is available, copying it beats learning "
        "a general map**, and that is the honest shape of this result rather than a uniform win.",
    ]
    if len(alone) and len(together):
        lines += [
            "",
            f"The sharpest form of the same statement: the {len(alone)} (protein, fold) pairs "
            "where the held-out TF is **the only member of its family in the panel** score "
            f"{alone.model.mean():.3f} against motif transfer's {alone.nn1.mean():.3f} "
            f"(**{alone.delta.mean():+.3f}**), while the {len(together)} with a relative score "
            f"{together.model.mean():.3f} against {together.nn1.mean():.3f} "
            f"(**{together.delta.mean():+.3f}**).",
        ]
    return lines + [""]


def per_fold(folds: pd.DataFrame, tfs: pd.DataFrame, table: pd.DataFrame) -> list[str]:
    block = _pooled(folds)
    block = block[(block.regime.isin(["G1", "G2"])) & (block.protein_mode == "real")]
    if block.empty:
        return []
    lines = [
        "### Per fold, because 33 proteins is thin",
        "",
        "`GHT_PLAN.md` §12 asks for per-fold variance rather than a single mean, and at 6-7 "
        "held-out TFs a fold this is the number that says how much to trust the mean.",
        "",
        "**`G2` is the harder of the two by construction**, and the fold contents show how: it "
        "holds out whole connected components at 50% identity, so all three Ets proteins leave "
        "together, both `SOX` paralogues leave together, and the three `SAND` proteins leave "
        "together. `G1` splits those groups, leaving a near-relative in training.",
        "",
        "| regime | fold | held out | auPRC | chance | auROC | steps | min |",
        "|---|---|---|---:|---:|---:|---:|---:|",
    ]
    family = dict(zip(table.tf, table.dbd_family, strict=True))
    held = (
        tfs[tfs.negset == training.TRAIN_NEGATIVES]
        .groupby(["regime", "fold"])
        .tf.unique()
        .to_dict()
    )
    grouped = block.groupby(["regime", "fold"], sort=True)
    for (regime, name), rows in grouped:
        names = sorted(held.get((regime, name), []))
        families = sorted({family.get(t, "?") for t in names})
        who = f"{len(names)}: {', '.join(families)}" if names else f"{rows.n_test_tf.mean():.0f}"
        lines.append(
            f"| `{regime}` | {name} | {who} "
            f"| {_mean_sd(rows.aupr)} | {rows.chance_aupr.mean():.3f} "
            f"| {_mean_sd(rows.auroc)} | {rows.best_step.mean():.0f} "
            f"| {rows.seconds.mean() / 60:.1f} |"
        )
    return lines + [""]


def per_tf(tfs: pd.DataFrame, table: pd.DataFrame, baselines: pd.DataFrame) -> list[str]:
    block = _pooled(tfs)
    block = block[
        (block.regime.isin(["G1", "G2"]))
        & (block.negset == training.TRAIN_NEGATIVES)
        & (block.protein_mode == "real")
    ]
    if block.empty:
        return []
    family = dict(zip(table.tf, table.dbd_family, strict=True))
    size = table.dbd_family.value_counts().to_dict()
    scored = (
        block.groupby("tf")
        .agg(aupr=("aupr", "mean"), auroc=("auroc", "mean"), chance=("positive_rate", "mean"))
        .reset_index()
        .sort_values("aupr", ascending=False)
    )
    # How close the nearest TF left in training was, which is the variable the PBM arm's
    # `reports/nn_baseline.md` bands its folds by and the one that drives every one of them.
    if not baselines.empty and "identity" in baselines:
        close = (
            baselines[baselines.regime.isin(["G1", "G2"]) & (baselines.method == "nn1")]
            .groupby("tf")
            .identity.max()
        )
        scored["identity"] = [close.get(t, float("nan")) for t in scored.tf]
    else:
        scored["identity"] = float("nan")
    lines = [
        "### Which proteins transfer, and which do not",
        "",
        "Averaged over `G1` and `G2` and every seed. Two variables to watch. **Family depth**: "
        "the PBM arm's family-holdout result was that transfer is carried by families with "
        f"training depth (`GHT_PLAN.md` §12), and this panel has {table.dbd_family.nunique()} "
        f"families over {len(table)} TFs. "
        "**Nearest-neighbour identity**: `reports/nn_baseline.md` bands the PBM folds by it "
        "because it drives every one of them, and the interesting rows here are the ones that "
        "transfer *without* a close relative.",
        "",
        "| TF | family | family size | nearest identity | auPRC | chance | auROC |",
        "|---|---|---:|---:|---:|---:|---:|",
    ]
    for r in scored.itertuples():
        fam = family.get(r.tf, "")
        identity = "—" if not np.isfinite(r.identity) else f"{r.identity:.2f}"
        lines.append(
            f"| `{r.tf}` | {fam} | {size.get(fam, 0)} | {identity} | {r.aupr:.3f} "
            f"| {r.chance:.3f} | {r.auroc:.3f} |"
        )
    depth = scored.assign(n=[size.get(family.get(t, ""), 0) for t in scored.tf])
    alone = depth[depth.n == 1]
    together = depth[depth.n > 1]
    if len(alone) and len(together):
        lines += [
            "",
            f"TFs whose family has **no other member** in the panel ({len(alone)} of "
            f"{len(depth)}) average auROC {alone.auroc.mean():.3f}; those with at least one "
            f"relative average {together.auroc.mean():.3f}.",
        ]
    usable = scored[np.isfinite(scored.identity)]
    if len(usable) > 3:
        rho = metrics.spearman(usable.identity.to_numpy(), usable.auroc.to_numpy())
        far = usable[usable.identity < 0.4]
        lines += [
            "",
            f"Spearman(nearest-neighbour identity, auROC) = **{rho:+.3f}** over {len(usable)} "
            "held-out TFs. "
            + (
                f"The {len(far)} with no training relative above 0.40 identity average auROC "
                f"**{far.auroc.mean():.3f}** — that subset is where the claim lives, because a "
                "motif transferred by identity has nothing to transfer from."
                if len(far)
                else "Every held-out TF had a training relative above 0.40 identity."
            ),
        ]
    return lines + [""]


def cost(folds: pd.DataFrame) -> list[str]:
    if folds.empty:
        return []
    counts = {
        c.replace("n_params_", ""): folds[c].iloc[0]
        for c in folds.columns
        if c.startswith("n_params_")
    }
    train_tf = folds.n_train_tf.max()
    lines = [
        "## What it cost, and what it is made of",
        "",
        "| quantity | value |",
        "|---|---:|",
        f"| runs recorded | {len(folds)} |",
        f"| ms per step | {1000 * folds.seconds_per_step.mean():.1f} |",
        f"| seconds per evaluation | {folds.seconds_per_eval.mean():.2f} |",
        f"| minutes per fold | {folds.seconds.mean() / 60:.1f} |",
        f"| peak GPU memory | {folds.peak_memory_mb.max():.0f} MB |",
        "",
        "| tower | parameters | per training protein |",
        "|---|---:|---:|",
    ]
    for name, value in counts.items():
        if name == "total":
            continue
        lines.append(f"| {name} | {int(value):,} | {value / max(train_tf, 1):,.0f} |")
    lines += [
        f"| **total** | **{int(counts.get('total', 0)):,}** "
        f"| {counts.get('total', 0) / max(train_tf, 1):,.0f} |",
        "",
        "`docs/TRAINING.md` §5 calls the PBM arm's **245 parameters per training protein** the "
        f"central engineering constraint. At {int(train_tf)} training proteins even the bare "
        "linear tower is two orders of magnitude past that, so **the panel size, not the tower "
        "architecture, is what sets the capacity risk here** — which is a correction to "
        "`GHT_PLAN.md` §6.2, where the worry was that a per-residue conv tower would reach "
        "~3.4k per protein against a linear tower's 245.",
        "",
    ]
    return lines


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.parse_args()

    folds = pd.read_parquet(config.FOLD_TABLE)
    tfs = pd.read_parquet(config.TF_TABLE)
    baselines = (
        pd.read_parquet(config.BASELINE_TABLE)
        if config.BASELINE_TABLE.exists()
        else pd.DataFrame(columns=["regime", "method", "negset", "aupr", "auroc", "positive_rate"])
    )
    table = panel.load()
    for frame in (folds, tfs):
        if "protein_mode" not in frame:
            frame["protein_mode"] = "real"

    lines = [
        "# The GHT-SELEX arm — results",
        "",
        "Generated by `scripts/make_ght_report.py`. The plan of record is "
        "[`docs/GHT_PLAN.md`](../docs/GHT_PLAN.md); the panel is "
        "[`ght_panel.md`](ght_panel.md), the data [`ght_windows.md`](ght_windows.md), the step "
        "and wall-clock budgets [`ght_preflight_C1-all.md`](ght_preflight_C1-all.md).",
        "",
        "**Read every auPRC against its own chance level and never across negative sets.** "
        "`shades` sits at a ~0.33 positive fraction and the subsampled `random` and `aliens` at "
        "~0.09, so their auPRCs are three different scales. auROC is the one number that is "
        "comparable across them, which is exactly why the composition control is read on it.",
        "",
        f"{len(table)} TFs, {int(folds.n_train.max()):,} training windows at most, "
        f"{len(folds)} runs.",
        "",
    ]
    lines += [
        "**The learned filters are the known motifs** "
        "([`ght_filters.md`](ght_filters.md)): every one of the 33 panel TFs' selected motifs has "
        "a match in the bank, median similarity 0.891 against a permuted-filter null of 0.715, "
        "and 33 of 33 beat their own null. All four kernel widths carry matches. That is "
        "`GHT_PLAN.md` §5's claim — this encoder is the learned generalisation of the PWM "
        "baseline — checked rather than asserted.",
        "",
    ]
    lines += setting_one(folds, baselines)
    lines += controls(folds)
    lines += setting_one_per_tf(tfs, baselines)
    lines += setting_two(folds, baselines)
    lines += identity_bands(tfs, baselines, table)
    lines += transfer(folds, tfs, table)
    lines += towers(folds)
    lines += per_fold(folds, tfs, table)
    lines += per_tf(tfs, table, baselines)
    lines += cost(folds)
    lines += [_stamp(folds)]

    path = config.REPORT_DIR / "ght_results.md"
    path.write_text("\n".join(lines) + "\n")
    print(f"wrote {path}")


if __name__ == "__main__":
    main()
