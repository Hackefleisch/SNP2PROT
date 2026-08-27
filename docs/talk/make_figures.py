"""Figures for the 2026-08-27 collaborator talk.

Reads only cached artifacts already in `data/processed/` — nothing is retrained and
nothing is recomputed from the 44 M-row table. Every number plotted here traces to
`reports/nn_baseline.md`, `reports/training.md` or `docs/RESULTS.md` at commit 9390e41.

    .venv/bin/python docs/talk/make_figures.py

Palette: the data-viz reference instance, light mode. Baseline is always orange and the
model always blue, in every figure, so the eye carries one mapping across slides.
"""

from __future__ import annotations

import warnings
from pathlib import Path

import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyArrowPatch, FancyBboxPatch

warnings.filterwarnings("ignore")

ROOT = Path(__file__).resolve().parents[2]
OUT = Path(__file__).resolve().parent / "figures"
OUT.mkdir(exist_ok=True)

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#8f8e88"
GRID = "#e6e5e0"
GREY_PT = "#c9c8c2"

BLUE = "#2a78d6"
ORANGE = "#eb6834"
AQUA = "#1baf7a"
VIOLET = "#4a3aa7"
RED = "#e34948"

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
        "axes.linewidth": 1.0,
        "xtick.color": INK2,
        "ytick.color": INK2,
        "xtick.major.size": 0,
        "ytick.major.size": 0,
        "grid.color": GRID,
        "grid.linewidth": 1.0,
        "legend.frameon": False,
        "figure.dpi": 110,
        "savefig.dpi": 220,
        "savefig.bbox": "tight",
    }
)


#: Every figure carries its own provenance, because a figure travels into a slide deck alone.
#: Says which arm and which fold produced it, or says plainly that no model was involved.
COMMIT = "SNP2PROT @ 9390e41"


def stamp(ax, text, dy=-0.19):
    """A muted provenance line under an axes, right-aligned."""
    ax.text(
        1.0,
        dy,
        f"{text}   ·   {COMMIT}",
        transform=ax.transAxes,
        ha="right",
        va="top",
        fontsize=8.5,
        color=MUTED,
    )


def _under(ax):
    """Hang a frameless legend below the panel, where no data point can collide with it."""
    ax.legend(
        loc="upper left",
        bbox_to_anchor=(0.0, -0.02),
        frameon=False,
        fontsize=9.5,
        markerscale=1.5,
        labelcolor=INK2,
        handletextpad=0.4,
        labelspacing=0.55,
        borderpad=0.0,
    )


def despine(ax, keep=("left", "bottom")):
    for side in ("top", "right", "left", "bottom"):
        ax.spines[side].set_visible(side in keep)


# --------------------------------------------------------------------------- data
def load():
    emb = np.load(ROOT / "data/processed/embeddings/A4.npz", allow_pickle=True)
    domains = emb["domains"]
    vectors = emb["vectors"]

    import sys

    sys.path.insert(0, str(ROOT / "src"))
    from snp2prot import corpus

    view = corpus.domains().set_index("dbd_seq")
    view = view.loc[domains]

    base = pd.read_parquet(ROOT / "data/processed/nn_baseline_domains.parquet")
    s2 = base[(base.regime == "S2") & (base.k == 1)].set_index("domain")
    ident = s2.neighbour_identity.reindex(domains).to_numpy()

    return domains, vectors, view, ident


# ------------------------------------------------------------------- fig 1: umap
def fig_umap(domains, vectors, view, ident):
    import umap

    reducer = umap.UMAP(n_neighbors=15, min_dist=0.15, metric="cosine", random_state=0)
    xy = reducer.fit_transform(vectors)

    fam = view.dbd_family.to_numpy()
    counts = pd.Series(fam).value_counts()
    top3 = list(counts.index[:3])
    is_var = view.is_variant.to_numpy()
    in_var_cluster = view.cluster_size.to_numpy() > 1

    fig, axes = plt.subplots(1, 3, figsize=(16.2, 5.4))
    fig.subplots_adjust(wspace=0.10)

    # panel A — the three thickest families
    ax = axes[0]
    other = ~np.isin(fam, top3)
    ax.scatter(
        xy[other, 0],
        xy[other, 1],
        s=9,
        c=GREY_PT,
        lw=0,
        alpha=0.9,
        label=f"other 53 families  (n={other.sum()})",
    )
    for colour, name in zip((BLUE, ORANGE, AQUA), top3, strict=True):
        m = fam == name
        ax.scatter(xy[m, 0], xy[m, 1], s=13, c=colour, lw=0, label=f"{name}  (n={m.sum()})")
    ax.set_title("A   56 protein families", loc="left", fontsize=12.5, color=INK, pad=12)
    _under(ax)

    # panel B — how close the nearest neighbour is
    ax = axes[1]
    ramp = mpl.colors.LinearSegmentedColormap.from_list(
        "blues_hi", plt.get_cmap("Blues")(np.linspace(0.28, 1.0, 256))
    )
    sc = ax.scatter(xy[:, 0], xy[:, 1], s=13, c=ident, cmap=ramp, vmin=0.15, vmax=1.0, lw=0)
    cax = ax.inset_axes([0.0, -0.075, 0.62, 0.030])
    cb = fig.colorbar(sc, cax=cax, orientation="horizontal")
    cb.outline.set_visible(False)
    cb.ax.tick_params(length=0, colors=INK2, labelsize=9, pad=2)
    ax.set_title(
        "B   Identity to nearest relative",
        loc="left",
        fontsize=12.5,
        color=INK,
        pad=12,
    )

    # panel C — the variant clusters
    ax = axes[2]
    ax.scatter(
        xy[~in_var_cluster, 0],
        xy[~in_var_cluster, 1],
        s=9,
        c=GREY_PT,
        lw=0,
        label=f"1,057 clusters of one  (n={(~in_var_cluster).sum()})",
    )
    ax.scatter(
        xy[in_var_cluster & ~is_var, 0],
        xy[in_var_cluster & ~is_var, 1],
        s=46,
        c=VIOLET,
        lw=0,
        label="wild type of a variant cluster (108)",
    )
    ax.scatter(
        xy[is_var, 0], xy[is_var, 1], s=14, c=RED, lw=0, label="single-residue variant (173)"
    )
    ax.set_title("C   Where the mutations are", loc="left", fontsize=12.5, color=INK, pad=12)
    _under(ax)

    for ax in axes:
        ax.set_xticks([])
        ax.set_yticks([])
        despine(ax, keep=())
    fig.text(
        0.5,
        0.02,
        "UMAP of the ESM-DBP domain embedding (arm A4) — 1,338 points, identical coordinates "
        "in all three panels.  No trained model: this is the corpus, not a result.\n"
        "Panel B is identity to the nearest domain a model could have trained on, under the "
        f"S2 split.   ·   {COMMIT}",
        ha="center",
        fontsize=9.5,
        color=MUTED,
        linespacing=1.6,
    )
    fig.subplots_adjust(bottom=0.30)

    fig.savefig(OUT / "fig1_protein_axis_umap.png")
    plt.close(fig)
    np.savez(OUT / "_umap_coords.npz", xy=xy, domains=domains)
    print("fig1 written")


# --------------------------------------------------- fig 2: the identity-band result
BANDS = [
    ("< 0.30", 706, 0.032, 0.050, 0.030),
    ("0.30–0.50", 1383, 0.365, 0.316, 0.327),
    ("0.50–0.70", 315, 0.771, 0.676, 0.690),
    ("0.70–0.90", 500, 0.946, 0.837, 0.849),
    ("> 0.90", 915, 0.910, 0.881, 0.883),
]


def fig_bands():
    labels = [b[0] for b in BANDS]
    n = [b[1] for b in BANDS]
    base = [b[2] for b in BANDS]
    a4 = [b[4] for b in BANDS]

    x = np.arange(len(labels))
    w = 0.36
    fig, ax = plt.subplots(figsize=(9.6, 5.2))
    ax.set_axisbelow(True)
    ax.yaxis.grid(True)

    ax.bar(
        x - w / 2 - 0.012, base, w, color=ORANGE, label="copy the nearest known protein", zorder=3
    )
    ax.bar(x + w / 2 + 0.012, a4, w, color=BLUE, label="our two-tower model", zorder=3)

    for xi, (b, m) in enumerate(zip(base, a4, strict=True)):
        ax.text(xi - w / 2 - 0.012, b + 0.018, f"{b:.2f}", ha="center", fontsize=9.5, color=INK2)
        ax.text(xi + w / 2 + 0.012, m + 0.018, f"{m:.2f}", ha="center", fontsize=9.5, color=INK2)

    ax.set_xticks(x)
    ticks = [f"{lab}\n{c:,} proteins" for lab, c in zip(labels, n, strict=True)]
    ax.set_xticklabels(ticks, fontsize=10)
    ax.set_xlabel(
        "sequence identity between the held-out protein and its closest training protein",
        labelpad=12,
    )
    ax.set_ylabel("AUPR   (how well the binding sites are recovered)")
    ax.set_ylim(0, 1.02)
    despine(ax)
    ax.legend(loc="upper left", fontsize=10.5, labelcolor=INK2)
    ax.set_title(
        "Everything depends on one variable — and it is not the model",
        loc="left",
        fontsize=13.5,
        color=INK,
        pad=14,
    )
    stamp(
        ax,
        "every held-out domain pooled over all 19 folds · model = arm A4 (ESM-DBP) · "
        "baseline = nearest-neighbour lookup, k=1",
        dy=-0.26,
    )
    fig.savefig(OUT / "fig2_identity_bands.png")
    plt.close(fig)
    print("fig2 written")


# ------------------------------------------------------ fig 3: model vs bar by regime
REGIMES = [
    ("S1\nrandom\nproteins", 0.7690, 0.7066),
    ("S2\nnew protein\ngroups", 0.2947, 0.2799),
    ("P1\nwhole family\nunseen", 0.0244, 0.0035),
    ("P2\nmutations in\nthin families", 0.8822, 0.7499),
    ("P3\npoint\nmutations", 0.9008, 0.8652),
]


def fig_regimes():
    labels = [r[0] for r in REGIMES]
    base = [r[1] for r in REGIMES]
    model = [r[2] for r in REGIMES]
    x = np.arange(len(labels))
    w = 0.36

    fig, ax = plt.subplots(figsize=(11.0, 5.4))
    ax.set_axisbelow(True)
    ax.yaxis.grid(True)
    ax.bar(
        x - w / 2 - 0.012, base, w, color=ORANGE, label="copy the nearest known protein", zorder=3
    )
    ax.bar(x + w / 2 + 0.012, model, w, color=BLUE, label="our two-tower model", zorder=3)
    for xi, (b, m) in enumerate(zip(base, model, strict=True)):
        ax.text(xi - w / 2 - 0.012, b + 0.018, f"{b:.2f}", ha="center", fontsize=9.5, color=INK2)
        ax.text(xi + w / 2 + 0.012, m + 0.018, f"{m:.2f}", ha="center", fontsize=9.5, color=INK2)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, fontsize=10, linespacing=1.4)
    ax.set_ylabel("AUPR, averaged over held-out proteins")
    ax.set_ylim(0, 1.02)
    despine(ax)
    ax.legend(loc="upper center", fontsize=10.5, labelcolor=INK2, ncol=2)
    ax.set_title(
        "Five ways of holding proteins back — the model beats the bar on none of them",
        loc="left",
        fontsize=13.5,
        color=INK,
        pad=14,
    )
    stamp(
        ax,
        "mean over the folds of each regime · model = arm A4 (ESM-DBP) · baseline scored on "
        "the identical held-out domains",
        dy=-0.30,
    )
    fig.savefig(OUT / "fig3_regimes.png")
    plt.close(fig)
    print("fig3 written")


# ------------------------------------------------------------- fig 4: the two towers
def fig_towers():
    fig, ax = plt.subplots(figsize=(11.5, 5.0))
    ax.set_xlim(0, 11.5)
    ax.set_ylim(-0.45, 5.0)
    ax.axis("off")

    def box(x, y, w, h, text, edge, face, fs=10.5, weight="normal"):
        ax.add_patch(
            FancyBboxPatch(
                (x, y),
                w,
                h,
                boxstyle="round,pad=0.02,rounding_size=0.12",
                linewidth=1.6,
                edgecolor=edge,
                facecolor=face,
            )
        )
        ax.text(
            x + w / 2,
            y + h / 2,
            text,
            ha="center",
            va="center",
            fontsize=fs,
            color=INK,
            weight=weight,
            linespacing=1.5,
        )

    def arrow(x0, y0, x1, y1, colour):
        ax.add_patch(
            FancyArrowPatch(
                (x0, y0),
                (x1, y1),
                arrowstyle="-|>",
                mutation_scale=13,
                linewidth=1.5,
                color=colour,
                shrinkA=2,
                shrinkB=2,
            )
        )

    blue_f, orange_f = "#eaf2fd", "#fdeee7"

    box(0.15, 3.35, 2.0, 1.0, "DNA-binding\ndomain sequence", BLUE, "none")
    box(2.55, 3.35, 2.1, 1.0, "protein language\nmodel  (frozen)", BLUE, blue_f)
    box(5.05, 3.35, 1.9, 1.0, "one vector\nper protein", BLUE, blue_f)
    arrow(2.15, 3.85, 2.55, 3.85, BLUE)
    arrow(4.65, 3.85, 5.05, 3.85, BLUE)

    box(0.15, 0.65, 2.0, 1.0, "8 bp DNA site\n(all 32,896)", ORANGE, "none")
    box(2.55, 0.65, 2.1, 1.0, "small CNN\nover both strands", ORANGE, orange_f)
    box(5.05, 0.65, 1.9, 1.0, "one vector\nper 8-mer", ORANGE, orange_f)
    arrow(2.15, 1.15, 2.55, 1.15, ORANGE)
    arrow(4.65, 1.15, 5.05, 1.15, ORANGE)

    box(
        7.45,
        1.7,
        3.8,
        1.65,
        "shared 256-dimensional space\n\nbinding pairs pulled together,\n"
        "non-binding pairs pushed apart",
        INK2,
        "none",
        fs=11,
    )
    arrow(6.95, 3.85, 9.35, 3.35, BLUE)
    arrow(6.95, 1.15, 9.35, 1.70, ORANGE)

    ax.text(
        9.35,
        4.15,
        "score  =  cos(protein, DNA) / τ",
        ha="center",
        fontsize=12,
        color=INK,
        weight="bold",
    )
    ax.text(0.15, 4.62, "The model: two encoders, one space", fontsize=13.5, color=INK)
    ax.text(
        0.15,
        0.12,
        "472,705 trainable parameters — the language model is frozen and never trained; "
        "its output is read from a cache",
        fontsize=9.5,
        color=MUTED,
    )
    ax.text(
        11.35,
        -0.22,
        f"schematic — no data, no trained model   ·   {COMMIT}",
        ha="right",
        fontsize=8.5,
        color=MUTED,
    )

    fig.savefig(OUT / "fig4_two_towers.png")
    plt.close(fig)
    print("fig4 written")


# ---------------------------------------------------------- fig 5: what one row is
def fig_matrix():
    z = np.load(ROOT / "data/processed/kmer_matrix.npz", allow_pickle=True)
    e = z["escore"]
    rng = np.random.default_rng(0)

    # the 8-mers anything actually binds, and a random slice of proteins
    hot = np.argsort(-(e > 0.45).sum(axis=0))[:900]
    rows = rng.choice(e.shape[0], 200, replace=False)
    sub = e[np.ix_(rows, np.sort(hot))]

    # seriate both axes by the first principal component, so related motifs sit together
    def order_by_pc(m):
        c = m - m.mean(axis=0, keepdims=True)
        _, _, vt = np.linalg.svd(c, full_matrices=False)
        return np.argsort(c @ vt[0])

    sub = sub[order_by_pc(sub)][:, order_by_pc(sub.T)]

    fig = plt.figure(figsize=(12.8, 5.0))
    gs = fig.add_gridspec(
        2, 2, width_ratios=[2.35, 1], height_ratios=[1, 0.055], wspace=0.30, hspace=0.32
    )
    ax = fig.add_subplot(gs[0, 0])
    cax = fig.add_subplot(gs[1, 0])
    axb = fig.add_subplot(gs[0, 1])

    im = ax.imshow(sub, aspect="auto", cmap="Blues", vmin=0.0, vmax=0.5, interpolation="nearest")
    cb = fig.colorbar(im, cax=cax, orientation="horizontal")
    cb.outline.set_visible(False)
    cb.ax.tick_params(length=0, colors=INK2, labelsize=9)
    cb.set_label("E-score — how strongly this protein binds this 8-mer", color=INK2, fontsize=9.5)
    ax.set_xlabel("900 DNA 8-mers  (the ones anything binds)", color=INK2, fontsize=10)
    ax.set_ylabel("200 of the 1,338 protein domains", color=INK2, fontsize=10)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(
        "A   Every protein, measured against every 8-mer",
        loc="left",
        fontsize=12.5,
        color=INK,
        pad=12,
    )

    axb.set_axisbelow(True)
    axb.xaxis.grid(True)
    vals = [92_628, 748_240, 43_173_980]
    names = ["binds\n0.21%", "no call\n1.7%", "does not bind\n98.1%"]
    cols_ = [BLUE, GREY_PT, "#dcdbd5"]
    ypos = np.arange(3)[::-1]
    axb.barh(ypos, vals, color=cols_, height=0.5, zorder=3)
    axb.set_xscale("log")
    axb.set_yticks(ypos)
    axb.set_yticklabels(names, fontsize=10.5)
    for y, v in zip(ypos, vals, strict=True):
        axb.text(v * 1.35, y, f"{v:,}", va="center", fontsize=10, color=INK2)
    axb.set_xlim(3e4, 8e9)
    axb.set_xticks([])
    axb.set_xticks([], minor=True)
    axb.tick_params(axis="x", which="both", length=0)
    despine(axb, keep=("bottom",))
    axb.set_xlabel("rows in the merged table  (log scale)", color=INK2, fontsize=9.5)
    axb.set_title(
        "B   44,014,848 rows, 466 negatives per positive",
        loc="left",
        fontsize=12.5,
        color=INK,
        pad=12,
    )

    stamp(
        axb,
        "measured E-scores from the merged table · no model, no split, no prediction",
        dy=-0.42,
    )
    fig.savefig(OUT / "fig5_what_the_data_is.png")
    plt.close(fig)
    print("fig5 written")


if __name__ == "__main__":
    fig_bands()
    fig_regimes()
    fig_towers()
    fig_matrix()
    domains, vectors, view, ident = load()
    fig_umap(domains, vectors, view, ident)
