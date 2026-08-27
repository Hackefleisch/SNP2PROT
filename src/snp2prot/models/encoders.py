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
