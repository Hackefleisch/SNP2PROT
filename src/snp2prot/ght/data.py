"""The window corpus as resident arrays, and the views a fold reads.

One `(N, 301)` `uint8` block for every TF's windows plus the index columns that say which TF,
which split, which negative set and which label each row is. Built once per process.

**Positives and `shades` go to the device; `random` and `aliens` stay on the host.** The first
pair is 356 MB of tokens across both chromosome splits and is what every training step and every
validation pass touches; the second is 2.4 GB and is touched exactly once per run, at the final
evaluation, where it is streamed in chunks. On an 8 GB card that is the difference between a
training loop with no host traffic and one that pages.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd

from snp2prot.ght import config, windows

#: The negative sets kept resident. `shades` is the primary set (`GHT_PLAN.md` D4) and the only
#: one anything trains on.
RESIDENT_SETS = ("", "shades")


@dataclass
class WindowCorpus:
    """Every panel TF's windows in one set of arrays, indexed by row."""

    tfs: np.ndarray  # (n_tf,) names, the order every protein array is indexed by
    tokens: np.ndarray  # (N, 301) uint8
    tf_index: np.ndarray  # (N,) int32 into `tfs`
    label: np.ndarray  # (N,) int8
    split: np.ndarray  # (N,) 'Train' | 'Test'
    negset: np.ndarray  # (N,) '' | 'shades' | 'random' | 'aliens'
    chrom: np.ndarray  # (N,)
    ratio: np.ndarray  # (n_tf,) realised negatives-per-positive of the subsampled sets

    @property
    def n_windows(self) -> int:
        return len(self.label)

    @property
    def width(self) -> int:
        return self.tokens.shape[1]

    def rows(
        self,
        tfs: tuple[str, ...] | None = None,
        chromosomes: tuple[str, ...] | None = None,
        negset: str | None = None,
    ) -> np.ndarray:
        """Row positions matching all of the given constraints.

        `negset` selects the positives **plus** that one negative set, which is the only
        combination anything scores: a mixture of two negative sets has no chance level.
        """
        keep = np.ones(self.n_windows, dtype=bool)
        if tfs is not None:
            wanted = np.isin(self.tfs, np.asarray(tfs))
            keep &= wanted[self.tf_index]
        if chromosomes is not None:
            keep &= np.isin(self.chrom, np.asarray(chromosomes))
        if negset is not None:
            keep &= (self.label == 1) | (self.negset == ("" if negset == "positives" else negset))
        return np.flatnonzero(keep)

    def chromosomes(self, split: str) -> list[str]:
        return sorted(set(np.unique(self.chrom[self.split == split]).tolist()))

    def window_counts(self, split: str) -> dict[str, int]:
        """Windows per chromosome on one split, for the validation carve."""
        rows = (self.split == split) & np.isin(self.negset, np.asarray(RESIDENT_SETS))
        names, counts = np.unique(self.chrom[rows], return_counts=True)
        return dict(zip(names.tolist(), counts.tolist(), strict=True))

    def summary(self) -> pd.DataFrame:
        frame = pd.DataFrame(
            {
                "tf": self.tfs[self.tf_index],
                "split": self.split,
                "set": np.where(self.negset == "", "positives", self.negset),
                "label": self.label,
            }
        )
        return frame.groupby(["split", "set"]).agg(n=("label", "size")).reset_index()


def load(panel_tfs, resident_only: bool = False) -> WindowCorpus:
    """Read every TF's `.npz` and concatenate. `resident_only` drops `random` and `aliens`."""
    names = [str(t) for t in panel_tfs]
    blocks, index, label, split, negset, chrom, ratios = [], [], [], [], [], [], []
    for i, tf in enumerate(names):
        w = windows.TFWindows.load(tf)
        keep = (
            np.isin(w.negset, np.asarray(RESIDENT_SETS))
            if resident_only
            else np.ones(len(w.label), dtype=bool)
        )
        blocks.append(w.tokens[keep])
        index.append(np.full(int(keep.sum()), i, dtype=np.int32))
        label.append(w.label[keep])
        split.append(w.split[keep])
        negset.append(w.negset[keep])
        chrom.append(w.chrom[keep])
        ratios.append(w.ratio)
    return WindowCorpus(
        tfs=np.asarray(names, dtype=np.str_),
        tokens=np.concatenate(blocks),
        tf_index=np.concatenate(index),
        label=np.concatenate(label),
        split=np.concatenate(split).astype(np.str_),
        negset=np.concatenate(negset).astype(np.str_),
        chrom=np.concatenate(chrom).astype(np.str_),
        ratio=np.asarray(ratios, dtype=np.float32),
    )


def load_embeddings(arm: str, panel: pd.DataFrame) -> tuple[np.ndarray, str]:
    """`(n_tf, width)` pooled protein vectors in panel order, and the model name.

    Read through `snp2prot.embeddings.DomainEmbeddings` so the GHT arm and the PBM arm cannot
    drift apart in what a "pooled embedding" is; only the path differs. Keyed on `dbd_seq` rather
    than on the TF name, for the same reason the PBM arm is: the sequence is what was embedded,
    and a name is a label that could be attached to a different construct tomorrow.
    """
    from snp2prot.embeddings import DomainEmbeddings

    table = DomainEmbeddings.load(arm, config.embedding_table(arm))
    return table.vectors[table.rows_for(list(panel.dbd_seq))], table.model


#: How the protein tower is fed. `real` is the model; the other two are the controls that decide
#: whether the protein is a variable or an ornament.
PROTEIN_MODES = ("real", "constant", "shuffled")


def ablate(vectors: np.ndarray, mode: str = "real", seed: int = 0) -> np.ndarray:
    """The protein-side ablation, and it is the control this whole arm stands or falls on.

    A two-tower model trained on 33 TFs' peaks against their own flanks can reach a high auPRC
    **without using the protein at all**: "does this window look like a bound region" is a real and
    largely protein-independent signal — open chromatin, CpG islands, promoter composition — and a
    DNA tower alone can learn it. Setting 1 cannot distinguish that from a model that conditions on
    the protein, because every test TF was also a training TF.

    So the number that matters is not the model's auPRC but the **gap** between `real` and:

    - **`constant`** — every TF handed the same vector, so the model is structurally protein-blind
      and can only answer "is this DNA bindable by something". Anything `real` scores above this is
      what conditioning on the protein bought.
    - **`shuffled`** — the panel's vectors permuted, so the protein axis is present, carries the
      same distribution, and is *wrong*. It separates "uses the protein" from "benefits from having
      a per-TF free parameter of any kind".

    `constant` uses the panel mean rather than zeros: zeros would land the tower on a `LayerNorm`
    of a constant, which is well defined but is also a point no real embedding is near, and the
    control should differ from the model in what it knows, not in where it sits.
    """
    if mode == "real":
        return vectors
    if mode == "constant":
        return np.repeat(vectors.mean(axis=0, keepdims=True), len(vectors), axis=0)
    if mode == "shuffled":
        rng = np.random.default_rng(seed)
        order = rng.permutation(len(vectors))
        # A derangement, not a permutation: a permutation can leave a TF with its own vector, and
        # at n = 33 the expected number of such fixed points is exactly 1.
        for i in range(len(order)):
            if order[i] == i:
                j = (i + 1) % len(order)
                order[i], order[j] = order[j], order[i]
        return vectors[order]
    raise ValueError(f"unknown protein mode {mode!r}; have {PROTEIN_MODES}")


def load_residue_embeddings(arm: str, panel: pd.DataFrame) -> tuple[np.ndarray, np.ndarray, str]:
    """`(n_tf, L, width)` per-residue vectors padded to the longest domain, plus a `(n_tf, L)` mask.

    Dense rather than ragged because the panel is 33 domains of 59-219 aa: the padded block is
    37 MB and lives on the device beside the windows, where a CSR layout would buy nothing and
    cost a gather per step. The PBM arm's `ResidueEmbeddings` is ragged for the opposite reason —
    1,338 domains of 30-378 aa is two thirds padding.
    """
    from snp2prot.embeddings import ResidueEmbeddings

    table = ResidueEmbeddings.load(arm, config.residue_embedding_table(arm))
    order = {str(s): i for i, s in enumerate(table.domains)}
    rows = [order[str(s)] for s in panel.dbd_seq]
    blocks = [table.residues(i) for i in rows]
    longest = max(len(b) for b in blocks)
    dense = np.zeros((len(blocks), longest, table.width), dtype=np.float32)
    mask = np.zeros((len(blocks), longest), dtype=bool)
    for i, block in enumerate(blocks):
        dense[i, : len(block)] = block
        mask[i, : len(block)] = True
    return dense, mask, table.model
