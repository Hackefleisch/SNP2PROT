"""Two figures that need no projection caveat.

    .venv/bin/python docs/talk/dump_shared_space.py     # first, ~3 min
    .venv/bin/python docs/talk/make_score_figures.py    # then, ~3 min

**fig7 — the claim itself, unprojected.** Two histograms of the model's score (a cosine) over
every (protein, 8-mer) pair, split by whether the pair binds. No dimensionality reduction, one
number per pair, nothing to argue about. Split again into the proteins the model trained on and
the 173 it never saw, so the generalisation gap is on the same axis.

**fig8 — the DNA landscape, with proteins placed on it.** A UMAP of the 32,896 8-mers *alone*.
One modality, so every neighbour is a real neighbour and nothing is imposed — the problem that
forced `fig6` to add cross-modal edges cannot arise. Each protein is then placed at the 8-mer it
scores highest, which is a stated rule rather than a layout artefact: whether that marker lands
inside the protein's measured binding sites is then a genuine read-out, not a drawing choice.
"""

from __future__ import annotations

import sys
import warnings
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent / "figures"
sys.path.insert(0, str(ROOT / "src"))

from snp2prot import corpus, distances, splits  # noqa: E402

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#8f8e88"
GRID = "#e6e5e0"
DNA_BG = "#dedcd6"
BLUE = "#2a78d6"
ORANGE = "#eb6834"
RED = "#e34948"

#: The same three proteins as `fig6`, so the two figures can be read against each other.
HIGHLIGHT = [
    ("Creb1", "bZIP", "seen in training"),
    ("Tbx5", "T-box", "seen in training"),
    ("Lhx5", "Homeodomain", "held out — never seen"),
]

mpl.rcParams.update(
    {
        "figure.facecolor": SURFACE,
        "axes.facecolor": SURFACE,
        "savefig.facecolor": SURFACE,
        "font.family": "DejaVu Sans",
        "font.size": 11,
        "text.color": INK,
        "axes.labelcolor": INK2,
        "axes.edgecolor": GRID,
        "xtick.color": INK2,
        "ytick.color": INK2,
        "xtick.major.size": 0,
        "ytick.major.size": 0,
        "grid.color": GRID,
        "legend.frameon": False,
        "figure.dpi": 110,
        "savefig.dpi": 220,
        "savefig.bbox": "tight",
    }
)


def despine(ax, keep=("bottom",)):
    for side in ("top", "right", "left", "bottom"):
        ax.spines[side].set_visible(side in keep)


def load():
    space = np.load(ROOT / "data/processed/talk_shared_space.npz", allow_pickle=True)
    labels = np.load(ROOT / "data/processed/kmer_matrix.npz", allow_pickle=True)["label"]
    view = corpus.domains()
    fold = next(
        f for f in splits.all_regimes(view, distances.DomainDistances.load()) if f.label == "P3/all"
    )
    return space["protein"], space["dna"], labels, view, fold


# ----------------------------------------------------- fig 7: the scores, no projection
def fig_histograms(protein, dna, labels, fold):
    test = np.zeros(len(protein), dtype=bool)
    test[fold.test] = True

    fig, axes = plt.subplots(1, 2, figsize=(13.0, 5.2), sharex=True, sharey=True)
    bins = np.linspace(-0.7, 0.9, 130)
    panels = [
        (~test, "A   The 1,165 proteins the model trained on"),
        (test, "B   The 173 point-mutation variants it never saw"),
    ]

    for ax, (rows, title) in zip(axes, panels, strict=True):
        block = protein[rows] @ dna.T
        lab = labels[rows]
        binds = block[lab == 1]
        does_not = block[lab == 0]

        ax.hist(
            does_not,
            bins=bins,
            density=True,
            color=DNA_BG,
            lw=0,
            label=f"does not bind  (n={does_not.size:,})",
        )
        ax.hist(
            binds,
            bins=bins,
            density=True,
            color=BLUE,
            alpha=0.85,
            lw=0,
            label=f"binds  (n={binds.size:,})",
        )
        for value, colour in ((np.median(does_not), MUTED), (np.median(binds), BLUE)):
            ax.axvline(value, color=colour, lw=1.6, ls="--", zorder=4)
        ax.text(
            np.median(does_not) - 0.02,
            0.97,
            f"median {np.median(does_not):+.2f}",
            transform=ax.get_xaxis_transform(),
            fontsize=10,
            color=INK2,
            ha="right",
            va="top",
        )
        ax.text(
            np.median(binds) + 0.02,
            0.97,
            f"median {np.median(binds):+.2f}",
            transform=ax.get_xaxis_transform(),
            fontsize=10,
            color=BLUE,
            ha="left",
            va="top",
        )
        ax.set_title(title, loc="left", fontsize=12.5, color=INK, pad=12)
        ax.set_xlabel("the model's score for that (protein, 8-mer) pair — a cosine")
        ax.set_yticks([])
        despine(ax)
        ax.legend(
            loc="upper left",
            bbox_to_anchor=(0.0, -0.16),
            fontsize=10,
            labelcolor=INK2,
            ncol=1,
            handletextpad=0.5,
            labelspacing=0.5,
            borderpad=0.0,
        )
        print(
            f"  {title.split('  ')[0]}: binds {np.median(binds):+.3f}, "
            f"does not {np.median(does_not):+.3f}"
        )

    axes[0].set_ylabel("share of pairs")
    axes[1].text(
        1.0,
        -0.40,
        "arm A4 (ESM-DBP), trained on fold P3/all — wild types in training, the 173 "
        "point-mutation\nvariants held out; test AUPR 0.873.   ·   SNP2PROT @ 9390e41",
        transform=axes[1].transAxes,
        ha="right",
        va="top",
        fontsize=8.5,
        color=MUTED,
        linespacing=1.6,
    )
    fig.suptitle(
        "Does the model separate binding from non-binding? — every pair, no projection",
        x=0.005,
        ha="left",
        fontsize=13.5,
        color=INK,
        y=1.02,
    )
    fig.subplots_adjust(wspace=0.08)
    fig.savefig(OUT / "fig7_score_separation.png")
    plt.close(fig)
    print("fig7 written")


# -------------------------------------- fig 8: a DNA-only map, proteins placed by their best hit
def fig_landscape(protein, dna, labels, view):
    import umap

    reducer = umap.UMAP(n_neighbors=25, min_dist=0.25, metric="cosine", random_state=0)
    d_xy = reducer.fit_transform(dna)
    print(f"UMAP over {len(dna):,} 8-mers done")

    scores = protein @ dna.T
    best = scores.argmax(axis=1)  # each protein's highest-scoring 8-mer

    fig, axes = plt.subplots(1, 4, figsize=(18.0, 5.6))
    fig.subplots_adjust(wspace=0.08)

    for ax, (gene, family, provenance) in zip(axes[:3], HIGHLIGHT, strict=True):
        hits = view.index[(view.gene == gene) & (view.dbd_family == family)].to_numpy()
        if not len(hits):
            raise SystemExit(f"{gene} ({family}) is not in the corpus")
        row = int(hits[0])
        binds = labels[row] == 1
        landed = bool(binds[best[row]])

        ax.scatter(d_xy[~binds, 0], d_xy[~binds, 1], s=1.0, c=DNA_BG, lw=0, alpha=0.55)
        ax.scatter(
            d_xy[binds, 0],
            d_xy[binds, 1],
            s=17,
            c=BLUE,
            lw=0,
            label=f"the {int(binds.sum())} 8-mers {gene} binds",
        )
        ax.scatter(
            d_xy[best[row], 0],
            d_xy[best[row], 1],
            s=260,
            facecolor="none",
            edgecolor=RED,
            lw=2.6,
            zorder=5,
            label=f"{gene}, placed at its top-scoring 8-mer",
        )
        ax.set_title(f"{gene}   ·   {family}", loc="left", fontsize=13, color=INK, pad=12)
        ax.text(
            0.0,
            1.005,
            f"{provenance}   ·   top-scoring 8-mer {'IS' if landed else 'is NOT'} one it binds",
            transform=ax.transAxes,
            fontsize=9.5,
            color=BLUE if landed else ORANGE,
            va="bottom",
        )
        _under(ax)

    ax = axes[3]
    n_binders = (labels == 1).sum(axis=0)
    ramp = mpl.colors.LinearSegmentedColormap.from_list(
        "blues_hi", plt.get_cmap("Blues")(np.linspace(0.10, 1.0, 256))
    )
    order = np.argsort(n_binders)
    sc = ax.scatter(
        d_xy[order, 0], d_xy[order, 1], s=1.8, c=np.log1p(n_binders[order]), cmap=ramp, lw=0
    )
    cax = ax.inset_axes([0.0, -0.075, 0.62, 0.030])
    cb = fig.colorbar(sc, cax=cax, orientation="horizontal")
    ticks = [0, 1, 5, 25, 100, 400]
    cb.set_ticks(np.log1p(ticks))
    cb.set_ticklabels([str(t) for t in ticks])
    cb.outline.set_visible(False)
    cb.ax.tick_params(length=0, colors=INK2, labelsize=9, pad=2)
    cb.set_label("proteins that bind this 8-mer", color=INK2, fontsize=9.5, labelpad=4)
    ax.set_title("Every 8-mer, by how many bind it", loc="left", fontsize=13, color=INK, pad=12)
    ax.text(
        0.0,
        -0.20,
        f"{int((n_binders == 0).sum()):,} of 32,896 8-mers are bound by nothing\n"
        f"the busiest is bound by {int(n_binders.max())} of 1,338 proteins",
        transform=ax.transAxes,
        fontsize=9.5,
        color=INK2,
        va="top",
        linespacing=1.6,
    )

    for ax in axes:
        ax.set_xticks([])
        ax.set_yticks([])
        despine(ax, keep=())

    landed_all = labels[np.arange(len(protein)), best] == 1
    fig.text(
        0.5,
        0.02,
        "UMAP of the 32,896 DNA 8-mers alone, as the trained model represents them — one kind of "
        "point, so every neighbour is a real neighbour and nothing is imposed.\nEach protein is "
        "then placed at the single 8-mer it scores highest. Across all 1,338 proteins that marker "
        f"lands on an 8-mer they really bind {landed_all.mean():.0%} of the time.\n"
        "Arm A4 (ESM-DBP), trained on fold P3/all — wild types in training, the 173 "
        "point-mutation variants held out; test AUPR 0.873.   ·   SNP2PROT @ 9390e41",
        ha="center",
        fontsize=9.5,
        color=MUTED,
        linespacing=1.6,
    )
    fig.subplots_adjust(bottom=0.36)
    fig.savefig(OUT / "fig8_dna_landscape.png")
    plt.close(fig)
    print(
        f"fig8 written — top-scoring 8-mer is a true positive for "
        f"{landed_all.sum()}/{len(landed_all)} proteins"
    )


def _under(ax):
    ax.legend(
        loc="upper left",
        bbox_to_anchor=(0.0, -0.02),
        frameon=False,
        fontsize=9.5,
        markerscale=1.4,
        labelcolor=INK2,
        handletextpad=0.4,
        labelspacing=0.55,
        borderpad=0.0,
    )


if __name__ == "__main__":
    protein, dna, labels, view, fold = load()
    fig_histograms(protein, dna, labels, fold)
    fig_landscape(protein, dna, labels, view)
