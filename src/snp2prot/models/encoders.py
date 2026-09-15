"""The two towers: a CNN over 8-mers and a projection over pooled protein embeddings.

Both end in an L2-normalised vector of the same width `D`, which is what makes the arms
comparable — "structure wins" can never be "structure had more width" (`ML_PLAN.md` §4.2).

**The asymmetry between the towers is deliberate and is the main capacity decision**
(`docs/TRAINING.md` §5). The protein axis has 1,338 points and a 1,280-d input, so a single
linear projection to `D = 256` is already 245 parameters per training protein; the DNA axis has
32,896 points and every one of them is seen at every step. The two towers are in completely
different data regimes, so the protein tower defaults to a bare linear map and the DNA tower is
allowed to be a small convolutional network.
"""

from __future__ import annotations

import numpy as np
import torch
from torch import nn

#: Index per base. The complement is `3 - i`, which is what makes the reverse complement of a
#: token array a reversal and a subtraction rather than a table lookup.
BASES = "ACGT"
BASE_INDEX = {b: i for i, b in enumerate(BASES)}


def tokenise(kmers: np.ndarray) -> np.ndarray:
    """`(n,)` of 8-mer strings -> `(n, 8)` int64 base indices."""
    lengths = {len(k) for k in kmers}
    if len(lengths) != 1:
        raise ValueError(f"8-mers of mixed length: {sorted(lengths)}")
    try:
        return np.array([[BASE_INDEX[b] for b in k] for k in kmers], dtype=np.int64)
    except KeyError as exc:
        raise ValueError(f"non-ACGT base in the vocabulary: {exc}") from None


def reverse_complement(tokens: torch.Tensor) -> torch.Tensor:
    """Reverse the sequence and complement every base.

    Needed as a *computation* rather than an index because the stored vocabulary is
    revcomp-collapsed: 32,896 = (4⁸ + 4⁴)/2 8-mers, 256 palindromic, and **zero** pairs with
    both strands present (measured, `docs/TRAINING.md` §1). So the reverse complement of a
    stored 8-mer is, with 256 exceptions, not itself a stored 8-mer and cannot be looked up.
    """
    return 3 - tokens.flip(-1)


class DNAEncoder(nn.Module):
    """One-hot 8×4 -> a small CNN -> `D`, mean-pooled over the two strands.

    **Reverse-complement mean pooling makes the embedding strand-symmetric by construction**
    (`ML_PLAN.md` §4.1). The array measured both strands together, so the E-score is inherently
    symmetric and each stored 8-mer is one arbitrary representative of a pair. An encoder
    without this would have to learn the equivalence from data it never sees, and would carry a
    strand asymmetry the assay does not have.

    **The convolution is not followed by a global pool.** Eight positions is short and a motif
    is position-specific — which base sits at position 3 is the signal, not how often a pattern
    occurs — so the channels are flattened and projected. A global pool would discard exactly
    the information the tower exists to encode.
    """

    def __init__(self, width: int, channels: int = 64, layers: int = 2, kernel: int = 3):
        super().__init__()
        blocks: list[nn.Module] = []
        in_channels = len(BASES)
        for _ in range(layers):
            blocks += [
                nn.Conv1d(in_channels, channels, kernel, padding=kernel // 2),
                nn.GELU(),
            ]
            in_channels = channels
        self.conv = nn.Sequential(*blocks)
        self.project = nn.Linear(channels * 8, width)

    def _one_strand(self, tokens: torch.Tensor) -> torch.Tensor:
        one_hot = nn.functional.one_hot(tokens, len(BASES)).float()  # (n, 8, 4)
        hidden = self.conv(one_hot.transpose(1, 2))  # (n, channels, 8)
        return self.project(hidden.flatten(1))

    def forward(self, tokens: torch.Tensor) -> torch.Tensor:
        """`(n, 8)` token indices -> `(n, D)`, L2-normalised."""
        pooled = 0.5 * (self._one_strand(tokens) + self._one_strand(reverse_complement(tokens)))
        return nn.functional.normalize(pooled, dim=-1)


class ProteinTower(nn.Module):
    """Pooled protein-LM vector -> `D`, L2-normalised.

    Linear by default, for the reason in the module docstring: with 1,338 training proteins,
    a hidden layer is a deliberate experiment and not the obvious choice. `hidden` makes it one.

    **The input is LayerNormed, and without it the arm comparison measures the wrong thing.**
    Every arm ends in a learned projection to a common `D` so that `A1` and `A4` "enter that
    projection on equal terms" (`snp2prot.embeddings`) — but they did not. Measured 2026-08-26:
    ESM-2's pooled vectors have norms 4.83-9.85 while ESM-DBP's have 0.76-1.24, so with a bias on
    the projection `‖b‖ / ‖Wx‖` was **0.13 for `A1` and 1.16 for `A4`** — for one arm the bias was
    a correction, for the other it outweighed the signal and set the output direction.

    The cost was not a handicap but an erasure. `A4`'s embeddings are the better separated of the
    two (mean pairwise cosine 0.752 against `A1`'s 0.874), yet at initialisation the tower emitted
    **0.8930 for both arms, identical to four decimals**: the bias had flattened away exactly the
    difference the `A1` -> `A4` delta exists to measure. With the `LayerNorm` the arms come
    through as 0.891 and 0.759, which is what their own structure says.

    LayerNorm rather than an L2-normalised input, which is a trap: it sets `‖Wx‖ ≈ 0.26` beside
    `‖b‖ = 0.251` and reproduces the same bias domination for *both* arms, equalising them at
    0.94 / 0.87 by collapsing both. It also discards the embedding norm, which is worth losing —
    it correlates with domain length at +0.243 on `A4` against +0.077 on `A1`, so it is closer to
    an artefact of pooling than to anything about the protein.
    """

    def __init__(self, in_features: int, width: int, hidden: int = 0, dropout: float = 0.0):
        """**`dropout` means two different things, depending on `hidden`.** With a hidden layer it
        drops hidden units, the usual sense. Without one there is no hidden layer to drop, so it
        drops *input features* — a different regulariser (it augments the frozen embedding rather
        than thinning a representation the model built). Both are defensible and both default to
        0.0; what is not defensible is reading a swept `dropout` across the two shapes as one
        quantity.
        """
        super().__init__()
        # LayerNorm first: dropping features and then renormalising would rescale whatever
        # survived, which is not what the dropout is for.
        layers: list[nn.Module] = [nn.LayerNorm(in_features)]
        if hidden:
            layers += [
                nn.Linear(in_features, hidden),
                nn.GELU(),
                nn.Dropout(dropout),
                nn.Linear(hidden, width),
            ]
        else:
            layers += [nn.Dropout(dropout), nn.Linear(in_features, width)]
        self.net = nn.Sequential(*layers)

    def forward(self, vectors: torch.Tensor) -> torch.Tensor:
        return nn.functional.normalize(self.net(vectors), dim=-1)


#: Kernel widths of the genomic DNA tower, in base pairs. Real DBD motifs run roughly 6-20 bp —
#: E-boxes short, nuclear-receptor dimers and POU sites long — and a single-layer convolution
#: followed by a global max pool has no composition across layers, so one width caps what a
#: filter can express (`docs/GHT_PLAN.md` §5). The `w = 8` bank is kept deliberately: it is the
#: PBM arm's unit, so the option of a bridge exists at zero cost, but it does not set `w` for
#: the rest.
GENOMIC_WIDTHS = (8, 12, 16, 20)


class GenomicDNAEncoder(nn.Module):
    """One-hot `(B, 4, L)` over a 301 bp window -> `D`, max-pooled over positions and strands.

    **The PWM baseline, generalised.** A PWM best-hit score is one fixed filter followed by a
    global maximum over the sequence; this is a bank of *learned* filters at four widths followed
    by the same maximum, and ArChIPelago is a fixed ensemble of such filters feeding a random
    forest. Three methods on one conceptual line, which is the argument as much as the model.

    **Global max pooling, and it is the opposite of `DNAEncoder`'s choice — on purpose.** That
    encoder sees 8 positions where *which base sits at position 3* is the signal, so it flattens
    and projects. Here the input is 301 positions holding one site among ~293 offsets, the peak is
    a coverage pile-up whose exact offset is an artefact of random fragmentation
    (`GHT_PLAN.md` §1), and the question is whether a motif occurs **anywhere** in the window. A
    max over positions answers that and makes positional overfitting structurally impossible: the
    encoder cannot learn "position 150 matters" because it cannot see position.

    **Max over both strands, not mean.** `DNAEncoder` mean-pools the strands because it encodes a
    strand-symmetric *identity* — the array measured both strands together and each stored 8-mer
    is one arbitrary representative of a pair. Here we are *scanning*, and the semantics wanted is
    "best hit on either strand", which is what a PWM best-hit computes. So the forward and
    reverse-complement position axes are concatenated and the maximum taken over both at once.
    """

    def __init__(
        self,
        width: int,
        widths: tuple[int, ...] = GENOMIC_WIDTHS,
        channels_per_width: int = 64,
        dropout: float = 0.0,
    ):
        super().__init__()
        self.widths = tuple(int(w) for w in widths)
        self.convolutions = nn.ModuleList(
            [nn.Conv1d(len(BASES), channels_per_width, w) for w in self.widths]
        )
        pooled = channels_per_width * len(self.widths)
        self.dropout = nn.Dropout(dropout)
        self.project = nn.Linear(pooled, width)

    @staticmethod
    def one_hot(tokens: torch.Tensor) -> torch.Tensor:
        """`(B, L)` int64 base indices -> `(B, 4, L)` float."""
        return nn.functional.one_hot(tokens.long(), len(BASES)).float().transpose(1, 2)

    def pooled(self, tokens: torch.Tensor) -> torch.Tensor:
        """`(B, L)` tokens -> `(B, sum(channels))` best hit per filter, over both strands."""
        forward = self.one_hot(tokens)
        # `reverse_complement` works on token indices, which is why it is applied before the
        # one-hot rather than by flipping channels afterwards — same result, one fewer copy.
        reverse = self.one_hot(reverse_complement(tokens))
        parts = []
        for conv in self.convolutions:
            both = torch.cat([conv(forward), conv(reverse)], dim=-1)
            parts.append(nn.functional.gelu(both).amax(dim=-1))
        return torch.cat(parts, dim=-1)

    def forward(self, tokens: torch.Tensor) -> torch.Tensor:
        """`(B, L)` token indices -> `(B, D)`, L2-normalised."""
        return nn.functional.normalize(self.project(self.dropout(self.pooled(tokens))), dim=-1)


#: Kernel widths of the per-residue protein tower, in residues. A DBD holds local patterns of
#: varying length at positions not comparable across families — the homeodomain's `WFQNRR` in
#: helix 3, the bZIP basic region before its leucine heptad, the Cys/His spacing of `zf-C4`, the
#: ELK signature in Ets — running roughly 3 to 15 residues (`GHT_PLAN.md` §6.1).
RESIDUE_WIDTHS = (3, 7, 15)


class ResidueProteinTower(nn.Module):
    """Per-residue protein-LM vectors -> `D`, max-pooled over residues. **Off by default.**

    `GHT_PLAN.md` §6.2, and the reason it defaults off is capacity rather than evidence. The naive
    form is `Conv1d(1280, 64, 16)` = 1.31 M parameters in a single bank; against 33 training
    proteins that is ~40,000 parameters per protein where `docs/TRAINING.md` §5 calls the PBM
    arm's **245 per protein** the central engineering constraint. The pointwise reduction below
    brings it to ~158 k, which is still ~4.8 k per protein — an order of magnitude better and
    still an order of magnitude worse than the linear tower. So it is a deliberate experiment with
    weight decay and dropout, run against the pooled default as the control, and both towers'
    parameter counts are reported either way.

    ```
    LayerNorm(1280)              per residue
    Linear(1280 -> d_r)          pointwise, shared across positions
    Conv1d(d_r, c, k)            k in RESIDUE_WIDTHS
    amax over residues           mask-aware
    Linear(sum(c) -> D)
    ```

    **Masking is not optional.** Domains run 59-219 aa in this panel, so a batch is padded, and a
    maximum taken over a pad position silently invents a motif. Padded positions are set to
    `-inf` before the pool, which is why the mask is a required argument rather than an option.

    **Max over residues means "is this motif present anywhere".** That mirrors the DNA tower and
    is the motif-detection semantics wanted. It discards *where*, which matters for variant
    effects — attention pooling is the natural second experiment, not the first.
    """

    def __init__(
        self,
        in_features: int,
        width: int,
        reduced: int = 64,
        widths: tuple[int, ...] = RESIDUE_WIDTHS,
        channels_per_width: int = 32,
        dropout: float = 0.0,
    ):
        super().__init__()
        self.widths = tuple(int(k) for k in widths)
        self.norm = nn.LayerNorm(in_features)
        self.reduce = nn.Linear(in_features, reduced)
        self.convolutions = nn.ModuleList(
            # `padding=k // 2` so a domain shorter than the widest kernel still produces output;
            # the mask then removes whatever the padding invented at the ends.
            [nn.Conv1d(reduced, channels_per_width, k, padding=k // 2) for k in self.widths]
        )
        self.dropout = nn.Dropout(dropout)
        self.project = nn.Linear(channels_per_width * len(self.widths), width)

    def forward(self, residues: torch.Tensor, mask: torch.Tensor) -> torch.Tensor:
        """`(B, L, 1280)` per-residue vectors and a `(B, L)` bool mask -> `(B, D)`, normalised.

        `mask` is True at real residues. Pad positions are zeroed before the convolution and
        `-inf` after it, so neither the kernel's input nor the maximum ever sees one.
        """
        hidden = self.reduce(self.norm(residues)) * mask.unsqueeze(-1)
        hidden = hidden.transpose(1, 2)  # (B, reduced, L)
        parts = []
        for conv in self.convolutions:
            activated = nn.functional.gelu(conv(hidden))[..., : mask.shape[1]]
            blocked = activated.masked_fill(~mask.unsqueeze(1), float("-inf"))
            parts.append(blocked.amax(dim=-1))
        pooled = torch.cat(parts, dim=-1)
        return nn.functional.normalize(self.project(self.dropout(pooled)), dim=-1)
