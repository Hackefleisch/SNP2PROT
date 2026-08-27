"""The `ML_PLAN.md` 5.2 figure — one UMAP of the shared protein–DNA space, four panels.

    .venv/bin/python docs/talk/dump_shared_space.py        # first, ~3 min
    .venv/bin/python docs/talk/make_shared_space_figure.py # then, ~2 min

Panels 1-3 are the plan's own construction: all proteins as dots, one highlighted per panel,
with the DNA behind it coloured by whether *that* protein binds it. Panel 4 colours every DNA
point by how many proteins bind it, which turns the promiscuity that made the first version of
this figure ill-defined into the subject of a panel.

**It illustrates; it does not evidence** (`ML_PLAN.md` 5.2). A contrastively-trained space
shows clusters by construction. The quantitative companion is `fig2_identity_bands.png`.
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

from snp2prot import corpus  # noqa: E402

SURFACE = "#fcfcfb"
INK = "#0b0b0b"
INK2 = "#52514e"
MUTED = "#8f8e88"
DNA_BG = "#dedcd6"
PROT_BG = "#a9a7a0"
BLUE = "#2a78d6"
ORANGE = "#eb6834"
RED = "#e34948"

#: One protein per panel. Recognisable to a biologist, three different folds, and the third
#: was **held out** of this fold's training — `P3/all` trains on wild types and tests on the
#: 173 variants, so Lhx5 is a protein the model never saw.
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
        "legend.frameon": False,
        "figure.dpi": 110,
        "savefig.dpi": 220,
        "savefig.bbox": "tight",
    }
)


def _topk(sim, k):
    """Indices of the `k` largest entries per row, and their cosine distances, sorted."""
    idx = np.argpartition(-sim, k, axis=1)[:, :k]
    got = np.take_along_axis(sim, idx, axis=1)
    order = np.argsort(-got, axis=1)
    return np.take_along_axis(idx, order, axis=1), 1.0 - np.take_along_axis(got, order, axis=1)


def joint_knn(protein, dna, k_same=15, k_cross=15, block=4096):
    """A neighbour graph over both modalities, with cross-modal edges forced in.

    **This is a deliberate choice and it has to be declared.** A plain joint k-NN does not
    produce this figure: measured on the trained space, a protein's nearest protein sits at
    cosine **+0.93** and an 8-mer's nearest 8-mer at **+0.79**, while a protein sits at
    **+0.41** from the 8-mers it actually binds. So every point's true nearest neighbours are
    its own modality, the graph has no cross-modal edges at all, and UMAP lays the two out as
    two separate islands — which says nothing except that same-modality points resemble each
    other.

    What the figure is about is the model's *ranking*: which 8-mers a protein scores highest.
    So each point keeps `k_same` neighbours of its own kind and is additionally joined to its
    `k_cross` best-scoring partners of the other kind. The edges are real cosines in the shared
    space — the score of a pair *is* their cosine — but which edges enter the graph is imposed,
    not discovered. Read the layout as an illustration of the ranking, never as evidence that
    the space is organised this way.
    """
    n_p, n_d = len(protein), len(dna)
    cross = protein @ dna.T

    pp_i, pp_d = _topk(protein @ protein.T, k_same + 1)
    pd_i, pd_d = _topk(cross, k_cross)
    dp_i, dp_d = _topk(cross.T, k_cross)

    dd_i = np.empty((n_d, k_same + 1), dtype=np.int64)
    dd_d = np.empty((n_d, k_same + 1), dtype=np.float32)
    for start in range(0, n_d, block):
        stop = min(start + block, n_d)
        dd_i[start:stop], dd_d[start:stop] = _topk(dna[start:stop] @ dna.T, k_same + 1)

    idx = np.vstack(
        [
            np.hstack([pp_i, pd_i + n_p]),
            np.hstack([dd_i + n_p, dp_i]),
        ]
    ).astype(np.int32)
    dist = np.vstack([np.hstack([pp_d, pd_d]), np.hstack([dd_d, dp_d])]).astype(np.float32)

    order = np.argsort(dist, axis=1)
    idx = np.take_along_axis(idx, order, axis=1)
    dist = np.take_along_axis(dist, order, axis=1)
    print(f"joint k-NN: {idx.shape[1]} neighbours per point, {k_cross} of them cross-modal")
    return idx, dist, None


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


def main() -> None:
    import umap

    space = np.load(ROOT / "data/processed/talk_shared_space.npz", allow_pickle=True)
    protein, dna = space["protein"], space["dna"]
    labels = np.load(ROOT / "data/processed/kmer_matrix.npz", allow_pickle=True)["label"]
    view = corpus.domains()

    joint = np.vstack([protein, dna])
    knn = joint_knn(protein, dna)
    reducer = umap.UMAP(
        n_neighbors=knn[0].shape[1],
        min_dist=0.2,
        metric="cosine",
        precomputed_knn=knn,
        random_state=0,
    )
    xy = reducer.fit_transform(joint)
    p_xy, d_xy = xy[: len(protein)], xy[len(protein) :]
    print(f"UMAP over {len(joint):,} points done")

    fig, axes = plt.subplots(1, 4, figsize=(18.0, 5.6))
    fig.subplots_adjust(wspace=0.08)

    for ax, (gene, family, provenance) in zip(axes[:3], HIGHLIGHT, strict=True):
        hits = view.index[(view.gene == gene) & (view.dbd_family == family)].to_numpy()
        if not len(hits):
            raise SystemExit(f"{gene} ({family}) is not in the corpus")
        row = int(hits[0])
        binds = labels[row] == 1

        ax.scatter(d_xy[~binds, 0], d_xy[~binds, 1], s=1.0, c=DNA_BG, lw=0, alpha=0.45)
        ax.scatter(p_xy[:, 0], p_xy[:, 1], s=9, c=PROT_BG, lw=0, alpha=0.6)
        ax.scatter(
            d_xy[binds, 0],
            d_xy[binds, 1],
            s=13,
            c=BLUE,
            lw=0,
            label=f"the {int(binds.sum())} 8-mers {gene} binds",
        )
        ax.scatter(
            p_xy[row, 0],
            p_xy[row, 1],
            s=190,
            c=RED,
            lw=1.8,
            edgecolor=SURFACE,
            zorder=5,
            label=f"{gene} — {provenance}",
        )
        ax.set_title(f"{gene}   ·   {family}", loc="left", fontsize=13, color=INK, pad=12)
        _under(ax)

    # panel 4 — the DNA landscape, coloured by promiscuity
    ax = axes[3]
    n_binders = (labels == 1).sum(axis=0)
    ramp = mpl.colors.LinearSegmentedColormap.from_list(
        "blues_hi", plt.get_cmap("Blues")(np.linspace(0.10, 1.0, 256))
    )
    order = np.argsort(n_binders)  # busy sites drawn last, so they are not buried
    sc = ax.scatter(
        d_xy[order, 0],
        d_xy[order, 1],
        s=1.8,
        c=np.log1p(n_binders[order]),
        cmap=ramp,
        lw=0,
    )
    cax = ax.inset_axes([0.0, -0.075, 0.62, 0.030])
    cb = fig.colorbar(sc, cax=cax, orientation="horizontal")
    ticks = [0, 1, 5, 25, 100, 400]
    cb.set_ticks(np.log1p(ticks))
    cb.set_ticklabels([str(t) for t in ticks])
    cb.outline.set_visible(False)
    cb.ax.tick_params(length=0, colors=INK2, labelsize=9, pad=2)
    cb.set_label("proteins that bind this 8-mer", color=INK2, fontsize=9.5, labelpad=4)
    quiet = int((n_binders == 0).sum())
    ax.set_title("Every 8-mer, by how many bind it", loc="left", fontsize=13, color=INK, pad=12)
    ax.text(
        0.0,
        -0.20,
        f"{quiet:,} of 32,896 8-mers are bound by nothing\n"
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
        for side in ("top", "right", "left", "bottom"):
            ax.spines[side].set_visible(False)

    fig.text(
        0.5,
        0.015,
        "One UMAP of the trained model's shared 256-d space — 1,338 proteins and 32,896 DNA "
        "8-mers together, identical coordinates in all four panels.\n"
        "Layout note: each point is joined to its 15 nearest neighbours of its own kind and to "
        "its 15 best-scoring partners of the other kind — the cross-modal edges are imposed, "
        "so read the layout as an illustration of the model's ranking, not as evidence.\n"
        "Arm A4 (ESM-DBP), trained on fold P3/all — wild types in training, the 173 "
        "point-mutation variants held out; test AUPR 0.873.   ·   SNP2PROT @ 9390e41",
        ha="center",
        fontsize=9.5,
        color=MUTED,
        linespacing=1.6,
    )
    fig.subplots_adjust(bottom=0.40)
    fig.savefig(OUT / "fig6_shared_space.png")
    plt.close(fig)
    print("fig6 written")


if __name__ == "__main__":
    main()
