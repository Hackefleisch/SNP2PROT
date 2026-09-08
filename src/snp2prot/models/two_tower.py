"""The two towers, the null anchor and the temperature, as one module.

`docs/TRAINING.md` is the design document; this is the assembly. The model owns four things
beyond the encoders:

**The null anchor** — a learned vector in the shared space that is not any 8-mer, appended to
the DNA embedding table so that the table has `K + 1` rows. It is L2-normalised like every other
point in the space, which keeps it plottable on §5.2's UMAP and keeps every score a cosine.

Its job is to give a domain with **no positive 8-mer** something to be the target of, so that
such a row produces a gradient rather than a vacuous average over an empty set
(`snp2prot.models.loss`). It was also expected to become a decision threshold and **does not** —
measured, a median of 15,159 of 32,460 8-mers outscore it. Do not read it as an operating point.

**The temperature** — stored as `log(1/τ)` and clamped, the CLIP convention. Learned by default:
it is one fewer thing to sweep across 19 folds × 4 arms, and the scale it has to find depends on
`D`, which varies between experiments.

**The clamp is load-bearing, and it is applied to the parameter rather than to the forward pass.**
Measured 2026-08-26 on `P3/all`: with the ceiling raised to 10,000 the scale never converges — it
passes 100 at about step 4,000 and is still climbing at 148 by step 6,000, which is what InfoNCE
does as the training loss falls toward zero. Clamping inside `score` instead left a **dead zone**:
above the ceiling the gradient is exactly 0, so with weight decay the parameter oscillated across
the boundary and without it froze at 4.60882 to five decimals with no force acting on it at all.
`clamp_temperature`, called after each optimiser step, is CLIP's own arrangement — the gradient
stays live and the value is simply pinned.

**No reported metric depends on any of this.** `score` is `scale x cosine` and every metric in
`snp2prot.evaluation.metrics` ranks *within* one domain, so a positive scalar cannot reorder
anything: AUPR, AUROC, precision@k and Spearman are invariant to the temperature to machine
precision, and end to end the fold scores 0.8558 clamped at 100 against 0.8554 running free to
148. What the temperature does set is the *dynamics* — how peaked the softmax is, and therefore
how much of the negative-side gradient lands on the hardest few 8-mers rather than being spread
across all 32,000 measured negatives. Whether the clamp is binding is reported with every run.

**The calibration head** — two scalars, `calibration_scale` and `calibration_bias`, that turn a
cosine into the logit of `P(binds)`. It is what gives the model a decision rule, and it is
trained by the `λ · L_bce` term of `snp2prot.models.loss.hybrid_loss` rather than by InfoNCE,
which cannot supply one: that loss is a per-row softmax and so is exactly invariant to a
per-protein offset (`T38`). At `λ = 0` the head takes no gradient and sits at its
initialisation — present, inert, and not to be read.

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

#: Where the calibration head starts. The scale is deliberately modest — a cosine lives in
#: [-1, 1], so 10 spans roughly [-10, 10] in logit space, which is wide enough to express
#: p = 5e-5 to p = 1 - 5e-5 and narrow enough not to saturate at initialisation.
CALIBRATION_SCALE_INIT = 10.0
#: The default base rate for the calibration bias, replaced per fold by
#: `TwoTower.set_calibration_prior`. The corpus rate is ~0.0021 positives per scored cell.
CALIBRATION_PRIOR_INIT = 0.002


def _logit(p: float) -> float:
    return float(math.log(p / (1.0 - p)))


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
        # The calibration head: two scalars, and its own, because `logit_scale` is saturated.
        self.calibration_scale = nn.Parameter(torch.tensor(math.log(CALIBRATION_SCALE_INIT)))
        self.calibration_bias = nn.Parameter(torch.tensor(_logit(CALIBRATION_PRIOR_INIT)))

    @torch.no_grad()
    def clamp_temperature(self) -> None:
        """Pin the learned scale to its ceiling, in place, after an optimiser step.

        CLIP's arrangement, and deliberately not a clamp inside `score`: clamping the forward
        pass makes the gradient exactly 0 above the ceiling, which strands the parameter in a
        dead zone where only weight decay can move it. Clamping the value keeps the gradient
        live and simply refuses to let it grow.
        """
        self.logit_scale.data.clamp_(max=math.log(self.max_logit_scale))

    @property
    def temperature(self) -> float:
        return float(1.0 / self.logit_scale.detach().exp())

    @property
    def temperature_is_clamped(self) -> bool:
        """Whether the learned scale is sitting on its ceiling.

        Measured on `P3/all`: it reaches the ceiling by about step 4,000 and stays there, so the
        temperature is effectively a fixed hyperparameter from that point rather than a learned
        one — and raising the ceiling to 10,000 shows it would keep climbing indefinitely. That
        is the clamp doing its job, but it should be visible rather than silent, so it is logged.
        """
        return bool(self.logit_scale.detach().exp() >= self.max_logit_scale * (1 - 1e-6))

    def dna_table(self, tokens: torch.Tensor) -> torch.Tensor:
        """`(K, 8)` tokens -> `(K + 1, D)`, the whole vocabulary plus the null anchor.

        Called once per step. The null is normalised alongside the real 8-mers so that every
        entry of the table is a point on the same unit sphere.
        """
        embedded = self.dna(tokens)
        null = nn.functional.normalize(self.null, dim=-1).unsqueeze(0)
        return torch.cat([embedded, null], dim=0)

    def cosine(self, protein_vectors: torch.Tensor, dna_table: torch.Tensor) -> torch.Tensor:
        """`(B, 1280)` and `(K + 1, D)` -> `(B, K + 1)` cosine similarities, no temperature.

        Both heads read this, so the matmul — the one real cost of a step on the protein axis —
        happens once whether or not the calibration term is in play.
        """
        return self.protein(protein_vectors) @ dna_table.T

    def score(self, protein_vectors: torch.Tensor, dna_table: torch.Tensor) -> torch.Tensor:
        """`(B, 1280)` and `(K + 1, D)` -> `(B, K + 1)` logits.

        No clamp here. The scale is bounded by `clamp_temperature` after each optimiser step,
        which keeps the gradient live rather than zeroing it past the ceiling.
        """
        return self.cosine(protein_vectors, dna_table) * self.logit_scale.exp()

    def calibrate(self, cosine: torch.Tensor) -> torch.Tensor:
        """Cosine similarities -> logits of `P(binds)`, on a scale shared across proteins.

        **Its own scale and bias, not the temperature.** Reusing `logit_scale` would look
        economical and would be wrong twice: it is clamped at `max_logit_scale` and reaches that
        ceiling by about step 4,000, so from there it is a constant and could not move to fit a
        probability; and the temperature's job is how peaked the softmax is, which is a
        statement about the *ranking* dynamics and unrelated to where the binding/non-binding
        line falls. Two jobs, two parameters.

        **The bias is initialised at the logit of the base rate.** With 466 negatives per
        positive, a head starting at `p = 0.5` spends its first thousands of steps discovering
        that almost everything is negative, and the cheapest way down is to flatten the scale —
        which destroys the ordering the other term is building. Starting at the prior is the
        standard fix for detection-scale imbalance and costs one line. `set_calibration_prior`
        replaces the default with the training pool's own rate.

        The null anchor's column is dropped: it carries no label, so it has no target.
        """
        return self.calibration_scale.exp() * cosine[..., :-1] + self.calibration_bias

    @torch.no_grad()
    def set_calibration_prior(self, positive_rate: float) -> None:
        """Move the calibration bias to `logit(positive_rate)` of the *training* pool.

        Called once, before the first step, from `snp2prot.training.Trainer.train`. It reads
        held-out labels never — the rate comes from the fold's training rows — so it is a
        property of what the model is allowed to see.
        """
        if not 0.0 < positive_rate < 1.0:
            raise ValueError(f"positive_rate must lie in (0, 1), not {positive_rate}")
        self.calibration_bias.fill_(_logit(positive_rate))

    def probabilities(self, protein_vectors: torch.Tensor, dna_table: torch.Tensor) -> torch.Tensor:
        """`(B, K)` calibrated `P(binds)` over the real 8-mers — the decision-rule output."""
        return torch.sigmoid(self.calibrate(self.cosine(protein_vectors, dna_table)))

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
            "calibration": self.calibration_scale.numel() + self.calibration_bias.numel(),
            "total": sum(p.numel() for p in self.parameters() if p.requires_grad),
        }
