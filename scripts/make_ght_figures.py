#!/usr/bin/env python
"""The talk's panels, from the tables `run_ght_grid.py` and `run_ght_baselines.py` wrote.

    python scripts/make_ght_figures.py            # ~5 s -> results/figures/ght_*.png

Four panels, one per thing the arm has to say:

| file | what it shows |
|---|---|
| `ght_preflight.png` | the step budget, read off the **held-out test** curve (`§8.1`) |
| `ght_settings.png` | Setting 1 parity and Setting 2 transfer, against the protein-blind control |
| `ght_negatives.png` | auROC by negative set — the composition and co-binding controls |
| `ght_per_tf.png` | which held-out proteins transfer, ordered, with family depth |

**Colour carries identity, never rank**, and the assignment is fixed: our model is always slot 1,
the protein-blind control always slot 2, a per-TF fitted baseline always slot 3, motif transfer
slots 4 and 5. A panel that drops a series does not repaint the survivors.

Every bar is direct-labelled. Three of the five slots sit below 3:1 contrast on a light surface,
so the value has to be legible without the colour — which is what a results figure wants anyway.
"""

from __future__ import annotations

import argparse

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from snp2prot import experiment  # noqa: E402
from snp2prot.evaluation import metrics  # noqa: E402
from snp2prot.ght import baselines, config, panel, training  # noqa: E402

#: The reference categorical palette, light mode, in its documented order — slots 1-5 of
#: `dataviz/references/palette.md`. Taken unchanged rather than re-picked: that order is the
#: configuration the skill's validator passes on the adjacent pairlist (worst CVD Delta E 9.1),
#: and re-ordering or substituting a hue would invalidate the check without re-running it.
SERIES = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4"]
SURFACE = "#fcfcfb"
INK = "#0b0b0b"
MUTED = "#52514e"
GRID = "#d8d7d2"

#: Method -> its fixed colour slot. Identity, not rank.
COLOUR = {
    "two-tower": SERIES[0],
    "protein-blind": SERIES[1],
    "PWM (per TF)": SERIES[2],
    "1NN motif": SERIES[3],
    "5NN motif": SERIES[4],
    "shuffled protein": SERIES[1],
}


#: The three protein configurations, in fixed order so a panel that lacks one does not shift the
#: others' colours.
MODEL_SERIES = (
    ("two-tower", "real"),
    ("protein-blind", "constant"),
    ("shuffled protein", "shuffled"),
)


def _style(ax) -> None:
    ax.set_facecolor(SURFACE)
    ax.spines[["top", "right"]].set_visible(False)
    ax.spines[["left", "bottom"]].set_color(GRID)
    ax.tick_params(colors=MUTED, labelsize=9)
    ax.yaxis.grid(True, color=GRID, linewidth=0.8)
    ax.set_axisbelow(True)


def _bars(ax, labels, values, errors=None, chance=None, ylabel="") -> None:
    x = np.arange(len(labels))
    colours = [COLOUR.get(label, SERIES[0]) for label in labels]
    # An error bar of exactly 0 draws a stray tick on a single-seed bar, which reads as a
    # measurement rather than as its absence.
    spread = errors if errors is not None and np.nansum(errors) > 0 else None
    ax.bar(x, values, width=0.62, color=colours, yerr=spread, ecolor=MUTED, capsize=3)
    for i, value in enumerate(values):
        if np.isfinite(value):
            ax.text(i, value + 0.012, f"{value:.3f}", ha="center", fontsize=9, color=INK)
    if chance is not None and np.isfinite(chance):
        ax.axhline(chance, color=MUTED, linewidth=1.2, linestyle=(0, (4, 3)))
        ax.text(
            -0.45,
            chance + 0.012,
            f"chance {chance:.3f}",
            ha="left",
            fontsize=8,
            color=MUTED,
        )
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=18, ha="right", fontsize=9, color=INK)
    ax.set_ylabel(ylabel, fontsize=9, color=MUTED)
    _style(ax)


def preflight_figure(out) -> None:
    files = sorted(config.GHT_PROCESSED.glob("ght_preflight_*.parquet"))
    if not files:
        return
    fig, axes = plt.subplots(
        1, len(files), figsize=(max(7.6, 5.4 * len(files)), 3.9), facecolor=SURFACE
    )
    axes = np.atleast_1d(axes)
    for ax, path in zip(axes, files, strict=True):
        curve = pd.read_parquet(path)
        ax.plot(curve.step, curve.test_aupr, color=SERIES[0], linewidth=2, label="held-out test")
        ax.plot(
            curve.step,
            curve.validation_aupr,
            color=SERIES[2],
            linewidth=2,
            linestyle=(0, (4, 3)),
            label="validation",
        )
        peak = curve.loc[curve.test_aupr.idxmax()]
        ax.plot([peak.step], [peak.test_aupr], "o", color=SERIES[0], markersize=7)
        ax.annotate(
            f"peak {peak.test_aupr:.3f}\nstep {int(peak.step):,}",
            (peak.step, peak.test_aupr),
            textcoords="offset points",
            xytext=(-10, -30),
            ha="right",
            fontsize=8,
            color=INK,
        )
        budget = int(experiment.section("ght")["training"]["steps"])
        ax.axvline(budget, color=MUTED, linewidth=1.2, linestyle=(0, (4, 3)))
        ax.text(
            budget * 1.03,
            curve.test_aupr.min() + 0.02,
            f"adopted budget\n{budget:,} steps",
            fontsize=8,
            color=MUTED,
        )
        # No second y-axis for the training loss: it lives on a different scale, and a dual-axis
        # chart is the one form this never draws. The loss is in the report's table instead.
        ax.set_title(str(curve.fold.iloc[0]), fontsize=10, color=INK, loc="left")
        ax.set_xlabel("step", fontsize=9, color=MUTED)
        ax.set_ylabel("auPRC (shades)", fontsize=9, color=MUTED)
        ax.legend(frameon=False, fontsize=8, labelcolor=INK)
        _style(ax)
    fig.suptitle(
        "The budget comes from the held-out curve, not the training loss",
        fontsize=11,
        color=INK,
        x=0.01,
        ha="left",
    )
    fig.tight_layout()
    fig.savefig(out / "ght_preflight.png", dpi=200, facecolor=SURFACE)
    plt.close(fig)


def _model(folds, regimes, mode, column, tower: str = "pooled"):
    """One configuration's mean and spread. **Pinned to one tower**, because the grid holds two
    and a bar that silently averages a pooled-vector run with a per-residue one is not a
    measurement of either."""
    rows = folds[folds.regime.isin(regimes) & (folds.protein_mode == mode)]
    if "protein_tower" in rows:
        rows = rows[rows.protein_tower == tower]
    return (float(rows[column].mean()), float(rows[column].std())) if len(rows) else (np.nan, 0.0)


def settings_figure(folds, baselines, out) -> None:
    negatives = training.TRAIN_NEGATIVES
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.2), facecolor=SURFACE)

    labels, values, errors = [], [], []
    for name, mode in MODEL_SERIES:
        mean, sd = _model(folds, ["C1"], mode, "aupr")
        if np.isfinite(mean):
            labels.append(name)
            values.append(mean)
            errors.append(sd)
    pwm = baselines[(baselines.regime == "C1") & (baselines.negset == negatives)]
    if not pwm.empty:
        labels.append("PWM (per TF)")
        values.append(float(pwm.aupr.mean()))
        errors.append(0.0)
    chance = folds[folds.regime == "C1"].chance_aupr.mean()
    _bars(axes[0], labels, values, errors, chance, "auPRC")
    axes[0].set_title(
        "Setting 1 — held-out chromosomes\nparity is the goal", fontsize=10, color=INK, loc="left"
    )

    labels, values, errors = [], [], []
    for name, mode in MODEL_SERIES:
        mean, sd = _model(folds, ["G1", "G2"], mode, "aupr")
        if np.isfinite(mean):
            labels.append(name)
            values.append(mean)
            errors.append(sd)
    for name, method in (("1NN motif", "nn1"), ("5NN motif", "nn5")):
        part = baselines[
            baselines.regime.isin(["G1", "G2"])
            & (baselines.method == method)
            & (baselines.negset == negatives)
        ]
        if not part.empty:
            labels.append(name)
            values.append(float(part.aupr.mean()))
            errors.append(0.0)
    chance = folds[folds.regime.isin(["G1", "G2"])].chance_aupr.mean()
    _bars(axes[1], labels, values, errors, chance, "auPRC")
    axes[1].set_title(
        "Setting 2 — held-out TF and chromosomes\nno per-TF method can run here",
        fontsize=10,
        color=INK,
        loc="left",
    )
    for ax in axes:
        ax.set_ylim(0, 1.05)
    fig.tight_layout()
    fig.savefig(out / "ght_settings.png", dpi=200, facecolor=SURFACE)
    plt.close(fig)


def negatives_figure(folds, out) -> None:
    """auROC by negative set. **auROC and not auPRC**, because the sets have different balances
    (~0.33 for `shades`, ~0.09 for the subsampled `random` and `aliens`) and auROC is the only
    number invariant to that. Reading auPRC across the three would be reading the imbalance."""
    sets = training.EVALUATION_NEGATIVES
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.0), facecolor=SURFACE)
    for ax, (title, regimes) in zip(
        axes,
        [("Setting 1 (C1)", ["C1"]), ("Setting 2 (G1 + G2)", ["G1", "G2"])],
        strict=True,
    ):
        width = 0.36
        x = np.arange(len(sets))
        for offset, (name, mode) in enumerate(
            (("two-tower", "real"), ("protein-blind", "constant"))
        ):
            values, errors = [], []
            for negatives in sets:
                suffix = "" if negatives == training.TRAIN_NEGATIVES else f"_{negatives}"
                mean, sd = _model(folds, regimes, mode, f"auroc{suffix}")
                values.append(mean)
                errors.append(sd)
            if not np.isfinite(values).any():
                continue
            position = x + (offset - 0.5) * (width + 0.02)
            ax.bar(
                position,
                values,
                width=width,
                color=COLOUR[name],
                label=name,
                yerr=errors,
                ecolor=MUTED,
                capsize=3,
            )
            for xi, value in zip(position, values, strict=True):
                if np.isfinite(value):
                    ax.text(xi, value + 0.012, f"{value:.3f}", ha="center", fontsize=8, color=INK)
        ax.axhline(0.5, color=MUTED, linewidth=1.2, linestyle=(0, (4, 3)))
        ax.text(-0.45, 0.515, "chance 0.500", ha="left", fontsize=8, color=MUTED)
        ax.set_xticks(x)
        ax.set_xticklabels(
            ["shades\n(local flanks)", "random\n(GC-matched)", "aliens\n(other TFs' peaks)"],
            fontsize=9,
            color=INK,
        )
        ax.set_ylim(0.4, 1.05)
        ax.set_ylabel("auROC", fontsize=9, color=MUTED)
        ax.set_title(title, fontsize=10, color=INK, loc="left")
        ax.legend(frameon=False, fontsize=8, labelcolor=INK, loc="lower right", ncols=2)
        _style(ax)
    fig.suptitle(
        "The negative set is the control: composition on random, the protein on aliens",
        fontsize=11,
        color=INK,
        x=0.01,
        ha="left",
    )
    fig.tight_layout()
    fig.savefig(out / "ght_negatives.png", dpi=200, facecolor=SURFACE)
    plt.close(fig)


def per_tf_figure(tfs, table, out) -> None:
    block = tfs[
        tfs.regime.isin(["G1", "G2"])
        & (tfs.negset == training.TRAIN_NEGATIVES)
        & (tfs.protein_mode == "real")
    ]
    if "protein_tower" in block:
        block = block[block.protein_tower == "pooled"]
    if block.empty:
        return
    size = table.dbd_family.value_counts().to_dict()
    family = dict(zip(table.tf, table.dbd_family, strict=True))
    scored = (
        block.groupby("tf")
        .agg(auroc=("auroc", "mean"), spread=("auroc", "std"))
        .reset_index()
        .sort_values("auroc")
    )
    scored["alone"] = [size.get(family.get(t, ""), 0) == 1 for t in scored.tf]

    fig, ax = plt.subplots(figsize=(7.2, 9.0), facecolor=SURFACE)
    y = np.arange(len(scored))
    colours = [SERIES[1] if a else SERIES[0] for a in scored.alone]
    ax.barh(y, scored.auroc, height=0.66, color=colours, xerr=scored.spread.fillna(0), ecolor=MUTED)
    for yi, value in zip(y, scored.auroc, strict=True):
        ax.text(value + 0.008, yi, f"{value:.2f}", va="center", fontsize=8, color=INK)
    ax.axvline(0.5, color=MUTED, linewidth=1.2, linestyle=(0, (4, 3)))
    ax.set_yticks(y)
    ax.set_yticklabels([f"{t}  ({family.get(t, '')})" for t in scored.tf], fontsize=8, color=INK)
    ax.set_xlim(0.35, 1.05)
    ax.set_xlabel("auROC on held-out chromosomes, held out as a protein", fontsize=9, color=MUTED)
    ax.set_title(
        "Which proteins transfer\nblue: a relative is in the panel   orange: the only member of "
        "its family",
        fontsize=10,
        color=INK,
        loc="left",
    )
    _style(ax)
    ax.xaxis.grid(True, color=GRID, linewidth=0.8)
    ax.yaxis.grid(False)
    fig.tight_layout()
    fig.savefig(out / "ght_per_tf.png", dpi=200, facecolor=SURFACE)
    plt.close(fig)


#: Sequence-logo base colours. Not the categorical palette: these four are a **domain
#: convention** readers of this field already have in their heads, and re-encoding them with the
#: chart palette would be correct by the rules and wrong for the reader. They are used nowhere
#: else on any panel, so they cannot be confused with a series.
BASE_COLOUR = {"A": "#1baf7a", "C": "#2a78d6", "G": "#eda100", "T": "#e34948"}


def _logo(ax, matrix: np.ndarray, title: str, scale_to: np.ndarray | None = None) -> None:
    """A sequence logo for a `(width, 4)` log-odds matrix, letters scaled by information.

    The matrix is turned into probabilities with a softmax over the four bases of each column —
    which is what a log-odds column means — and each column's total height is its information
    content, `2 - H` bits, the standard logo convention.

    **A learned filter's absolute scale is not identifiable**, so it is not drawn as if it were.
    The projection that follows the convolution can absorb any positive factor, and weight decay
    pushes the weights small: rendered raw, every learned filter is a flat line at 0 bits, which
    would say something about the optimiser's regularisation rather than about the motif. Passing
    `scale_to` rescales the filter to the standard deviation of the matrix it is being compared
    with. The Pearson correlation the pairing was made on is invariant to exactly this factor, so
    the rescaling changes no number — only whether the picture is legible.
    """
    if scale_to is not None:
        spread = float(np.std(matrix))
        if spread > 0:
            matrix = matrix * float(np.std(scale_to)) / spread
    from matplotlib.patches import PathPatch
    from matplotlib.textpath import TextPath
    from matplotlib.transforms import Affine2D

    centred = matrix - matrix.max(axis=1, keepdims=True)
    probability = np.exp(centred)
    probability /= probability.sum(axis=1, keepdims=True)
    with np.errstate(divide="ignore", invalid="ignore"):
        entropy = -np.nansum(probability * np.log2(probability), axis=1)
    information = np.clip(2.0 - entropy, 0, 2)

    for position in range(len(matrix)):
        order = np.argsort(probability[position])
        bottom = 0.0
        for index in order:
            height = float(probability[position, index] * information[position])
            if height < 0.01:
                continue
            letter = "ACGT"[index]
            font = {"family": "DejaVu Sans", "weight": "bold"}
            path = TextPath((0, 0), letter, size=1, prop=font)
            extent = path.get_extents()
            transform = (
                Affine2D()
                .translate(-extent.x0, -extent.y0)
                .scale(1 / extent.width * 0.92, 1 / extent.height * height)
                .translate(position + 0.04, bottom)
            )
            ax.add_patch(
                PathPatch(transform.transform_path(path), facecolor=BASE_COLOUR[letter], lw=0)
            )
            bottom += height
    ax.set_xlim(0, len(matrix))
    ax.set_ylim(0, 2.05)
    ax.set_yticks([0, 1, 2])
    ax.set_xticks([])
    ax.set_title(title, fontsize=9, color=INK, loc="left", pad=4)
    _style(ax)
    ax.yaxis.grid(False)


def filters_figure(out, n_pairs: int = 4) -> None:
    """The best-matching learned filter beside the motif it matches, for the top few TFs.

    The argument of `GHT_PLAN.md` §5 made visible: a PWM best-hit is one fixed filter plus a
    global maximum, and this encoder is a bank of learned ones. If that is true in content as
    well as in form, the filters should be readable as motifs.
    """
    matches = config.GHT_PROCESSED / "ght_filter_matches.parquet"
    found = sorted(config.CHECKPOINT_DIR.glob("ght_A1_C1_all_seed*.pt"))
    if not matches.exists() or not found:
        return
    frame = pd.read_parquet(matches)
    best = frame.loc[frame.groupby("tf").similarity.idxmax()].nlargest(n_pairs, "similarity")
    banks = baselines.learned_filters(found[0])
    pwms = baselines.read_pwms(list(best.tf))

    fig, axes = plt.subplots(2, n_pairs, figsize=(3.4 * n_pairs, 4.4), facecolor=SURFACE)
    for column, row in enumerate(best.itertuples()):
        learned = banks[row.filter_width][row.channel]
        known = next(p.matrix for p in pwms[row.tf] if p.name == row.pwm)
        # The encoder scans both strands and takes the maximum, so a filter's orientation is
        # arbitrary. Draw it in whichever one matches, and say so in the title.
        _, flipped, _ = baselines.best_alignment(learned, known)
        if flipped:
            learned = np.flip(learned, axis=(0, 1))
        strand = " · reverse strand" if flipped else ""
        _logo(
            axes[0, column],
            learned,
            f"learned filter w{row.filter_width}#{row.channel}{strand}",
            scale_to=known,
        )
        _logo(axes[1, column], known, f"{row.tf} · selected motif · r = {row.similarity:.2f}")
    axes[0, 0].set_ylabel("bits", fontsize=9, color=MUTED)
    axes[1, 0].set_ylabel("bits", fontsize=9, color=MUTED)
    fig.suptitle(
        "A convolution filter over one-hot DNA is a position weight matrix "
        "(rescaled to the motif it matches, and drawn on the strand that matches)",
        fontsize=11,
        color=INK,
        x=0.01,
        ha="left",
    )
    fig.tight_layout()
    fig.savefig(out / "ght_filters.png", dpi=200, facecolor=SURFACE)
    plt.close(fig)


def margin_figure(tfs, baselines, out) -> None:
    """The margin over motif transfer against how close the nearest training relative was.

    The arm's central claim as a picture: both methods degrade as the nearest characterised
    domain gets more distant, and the gap between them widens. One point per (held-out protein,
    fold) — `G1` and `G2` leave a protein with different neighbours, and that difference is the
    experiment, so they are not collapsed.
    """
    keys = ["regime", "fold", "tf"]
    block = tfs[
        tfs.regime.isin(["G1", "G2"])
        & (tfs.negset == training.TRAIN_NEGATIVES)
        & (tfs.protein_mode == "real")
    ]
    if "protein_tower" in block:
        block = block[block.protein_tower == "pooled"]
    nearest = baselines[
        baselines.regime.isin(["G1", "G2"])
        & (baselines.method == "nn1")
        & (baselines.negset == training.TRAIN_NEGATIVES)
    ]
    if block.empty or nearest.empty:
        return
    ours = block.groupby(keys, as_index=False).aupr.mean().rename(columns={"aupr": "model"})
    frame = ours.merge(
        nearest[[*keys, "aupr", "identity"]].rename(columns={"aupr": "nn1"}), on=keys
    )
    frame = frame[np.isfinite(frame.identity)]
    if len(frame) < 4:
        return

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.3), facecolor=SURFACE)
    ax = axes[0]
    for name, colour, column in (
        ("two-tower", COLOUR["two-tower"], "model"),
        ("1NN motif transfer", COLOUR["1NN motif"], "nn1"),
    ):
        ax.scatter(
            frame.identity,
            frame[column],
            s=34,
            color=colour,
            label=name,
            edgecolor=SURFACE,
            linewidth=0.8,
            zorder=3,
        )
    ax.set_xlabel("identity of the nearest domain left in training", fontsize=9, color=MUTED)
    ax.set_ylabel("auPRC on the held-out protein", fontsize=9, color=MUTED)
    ax.set_title("Both track identity", fontsize=10, color=INK, loc="left")
    ax.set_xlim(0, 1)
    ax.legend(frameon=False, fontsize=8, labelcolor=INK, loc="upper left")
    _style(ax)

    ax = axes[1]
    frame = frame.assign(delta=frame.model - frame.nn1)
    ax.axhline(0, color=MUTED, linewidth=1.2, linestyle=(0, (4, 3)))
    ax.scatter(
        frame.identity,
        frame.delta,
        s=34,
        color=COLOUR["two-tower"],
        edgecolor=SURFACE,
        linewidth=0.8,
        zorder=3,
    )
    # Band means as steps rather than a fitted line: the bands are what the report reports, and a
    # regression line through 27 points would imply a model of the relationship there is no
    # reason to assume.
    for lo, hi in ((0.0, 0.35), (0.35, 0.50), (0.50, 0.70), (0.70, 1.01)):
        part = frame[(frame.identity >= lo) & (frame.identity < hi)]
        if part.empty:
            continue
        ax.plot(
            [lo, min(hi, 1.0)],
            [part.delta.mean()] * 2,
            color=SERIES[1],
            linewidth=2.4,
            solid_capstyle="butt",
            zorder=2,
        )
        # The band count, because two of the four bands are thin and a reader has to see that
        # before reading a step as a trend.
        ax.text(
            (lo + min(hi, 1.0)) / 2,
            part.delta.mean(),
            f"n={len(part)}",
            ha="center",
            va="bottom",
            fontsize=8,
            color=SERIES[1],
        )
    rho = metrics.spearman(frame.identity.to_numpy(), frame.delta.to_numpy())
    ax.set_xlabel("identity of the nearest domain left in training", fontsize=9, color=MUTED)
    ax.set_ylabel("two-tower − 1NN, auPRC", fontsize=9, color=MUTED)
    ax.set_title(
        f"…and the margin concentrates at low identity (Spearman {rho:+.2f})",
        fontsize=10,
        color=INK,
        loc="left",
    )
    ax.set_xlim(0, 1)
    _style(ax)
    fig.suptitle(
        "The margin lives where there is nothing close to copy from",
        fontsize=11,
        color=INK,
        x=0.01,
        ha="left",
    )
    fig.tight_layout()
    fig.savefig(out / "ght_margin.png", dpi=200, facecolor=SURFACE)
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.parse_args()

    from snp2prot.config import FIGURE_DIR

    FIGURE_DIR.mkdir(parents=True, exist_ok=True)
    folds = pd.read_parquet(config.FOLD_TABLE)
    tfs = pd.read_parquet(config.TF_TABLE)
    baselines = (
        pd.read_parquet(config.BASELINE_TABLE)
        if config.BASELINE_TABLE.exists()
        else pd.DataFrame(columns=["regime", "method", "negset", "aupr", "auroc"])
    )
    for frame in (folds, tfs):
        if "protein_mode" not in frame:
            frame["protein_mode"] = "real"

    preflight_figure(FIGURE_DIR)
    filters_figure(FIGURE_DIR)
    margin_figure(tfs, baselines, FIGURE_DIR)
    settings_figure(folds, baselines, FIGURE_DIR)
    negatives_figure(folds, FIGURE_DIR)
    per_tf_figure(tfs, panel.load(), FIGURE_DIR)
    for name in sorted(FIGURE_DIR.glob("ght_*.png")):
        print(f"wrote {name}")


if __name__ == "__main__":
    main()
