"""The two towers, the null anchor and the temperature, as one module.

`docs/TRAINING.md` is the design document; this is the assembly. The model owns three things
beyond the encoders:

**The null anchor** — a learned vector in the shared space that is not any 8-mer, appended to
the DNA embedding table so that the table has `K + 1` rows. It is L2-normalised like every other
point in the space, which keeps it plottable on §5.2's UMAP and keeps every score a cosine.

**The temperature** — stored as `log(1/τ)` and clamped, the CLIP convention. Learned by default:
it is one fewer thing to sweep across 19 folds × 4 arms, and the scale it has to find depends on
`D`, which varies between experiments. **On this data the learned scale reaches its ceiling and
stays there** (`docs/TRAINING.md` §9), so the clamp is load-bearing rather than a formality, and
whether it is binding is reported with every run.

**The DNA table is computed once per step, not once per pair.** Every domain in a batch is scored
against the same 32,896 8-mers, so the encoder runs once over the vocabulary and the result is
shared across the batch — which is what makes a step cost well under a GFLOP
(`docs/TRAINING.md` §4.3).
"""

from __future__ import annotations

import math

import torch
from torch import nn

from snp2prot.models.encoders import DNAEncoder, ProteinTower


class TwoTower(nn.Module):
    """Protein tower, DNA tower, a null anchor and a temperature."""

    def __init__(
        self,
        protein_features: int = 1280,
        width: int = 256,
        dna_channels: int = 64,
        dna_layers: int = 2,
        protein_hidden: int = 0,
        protein_dropout: float = 0.0,
        temperature: float = 0.07,
        learn_temperature: bool = True,
        max_logit_scale: float = 100.0,
    ):
        super().__init__()
        self.protein = ProteinTower(protein_features, width, protein_hidden, protein_dropout)
        self.dna = DNAEncoder(width, dna_channels, dna_layers)
        self.null = nn.Parameter(torch.randn(width) * 0.01)
        logit_scale = torch.tensor(math.log(1.0 / temperature))
        self.logit_scale = nn.Parameter(logit_scale, requires_grad=learn_temperature)
        self.max_logit_scale = float(max_logit_scale)

    @property
    def temperature(self) -> float:
        scale = self.logit_scale.detach().exp().clamp(max=self.max_logit_scale)
        return float(1.0 / scale)

    @property
    def temperature_is_clamped(self) -> bool:
        """Whether the learned scale is sitting on its ceiling.

        Measured on `P3/all` (`docs/TRAINING.md` §9): it reaches the ceiling by about step
        4,000 and stays there, so the temperature is effectively a fixed hyperparameter from
        that point rather than a learned one. That is the clamp doing its job — a runaway scale
        saturates the softmax — but it should be visible rather than silent, so it is logged.
        """
        return bool(self.logit_scale.detach().exp() >= self.max_logit_scale)

    def dna_table(self, tokens: torch.Tensor) -> torch.Tensor:
        """`(K, 8)` tokens -> `(K + 1, D)`, the whole vocabulary plus the null anchor.

        Called once per step. The null is normalised alongside the real 8-mers so that every
        entry of the table is a point on the same unit sphere.
        """
        embedded = self.dna(tokens)
        null = nn.functional.normalize(self.null, dim=-1).unsqueeze(0)
        return torch.cat([embedded, null], dim=0)

    def score(self, protein_vectors: torch.Tensor, dna_table: torch.Tensor) -> torch.Tensor:
        """`(B, 1280)` and `(K + 1, D)` -> `(B, K + 1)` logits.

        The temperature is clamped rather than trusted: an unclamped learned scale can run away
        early in training and saturate the softmax, which is the failure CLIP's own clamp exists
        to prevent.
        """
        scale = self.logit_scale.clamp(max=math.log(self.max_logit_scale)).exp()
        return self.protein(protein_vectors) @ dna_table.T * scale

    def forward(self, protein_vectors: torch.Tensor, tokens: torch.Tensor) -> torch.Tensor:
        return self.score(protein_vectors, self.dna_table(tokens))

    def parameter_counts(self) -> dict[str, int]:
        """Trainable parameters per tower — reported with every result (`ML_PLAN.md` §4.2).

        Equal `D` controls the shared space but not the capacity feeding it, so a comparison
        between arms is only readable next to these numbers.
        """
        return {
            "protein_tower": sum(p.numel() for p in self.protein.parameters() if p.requires_grad),
            "dna_tower": sum(p.numel() for p in self.dna.parameters() if p.requires_grad),
            "null_anchor": self.null.numel(),
            "total": sum(p.numel() for p in self.parameters() if p.requires_grad),
        }
