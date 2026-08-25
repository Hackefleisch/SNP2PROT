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
    """

    def __init__(self, in_features: int, width: int, hidden: int = 0, dropout: float = 0.0):
        super().__init__()
        if hidden:
            self.net: nn.Module = nn.Sequential(
                nn.Linear(in_features, hidden),
                nn.GELU(),
                nn.Dropout(dropout),
                nn.Linear(hidden, width),
            )
        else:
            self.net = nn.Sequential(nn.Dropout(dropout), nn.Linear(in_features, width))

    def forward(self, vectors: torch.Tensor) -> torch.Tensor:
        return nn.functional.normalize(self.net(vectors), dim=-1)
