"""The genomic two-tower: a protein tower, a 301 bp DNA tower, and a calibration head.

`GHT_PLAN.md` §7. **Not a subclass of `snp2prot.models.TwoTower`, and the difference is the
point.** That model owns a null anchor and a temperature because its objective is a per-protein
softmax over a complete, shared 32,896-8-mer vocabulary; the anchor gives a protein that binds
nothing a target inside that softmax. GHT windows are per-TF and essentially disjoint — GCM1 and
FLI1 share 34 of ~71,000 and ~145,000 — so there is no shared vocabulary to normalise over, no
row to anchor, and the primary objective is binary cross-entropy on an explicit label.

What **is** reused, deliberately and verbatim in form:

- `cosine(protein, dna)` -> `calibrate(cosine)` -> `BCEWithLogits`. Cosine-then-scale rather than
  concatenating the two vectors into an MLP, because that is what makes this a genuine two-tower
  with a shared space — the architectural claim of the talk. An MLP over the concatenation would
  score better and prove nothing about a shared representation.
- `calibration_scale` and `calibration_bias` under those exact names, with
  `set_calibration_prior` initialising the bias at the training pool's own base rate. The PBM
  arm's `CALIBRATION_PRIOR_INIT = 0.002` is the PBM base rate and is wrong here by two orders of
  magnitude, which is why the prior is always set per fold rather than defaulted.
- `logit_scale`, the CLIP temperature, clamped after each optimiser step — used **only** by the
  optional per-protein InfoNCE term. At `mu = 0`, the headline setting, it takes no gradient.

The null anchor is absent rather than untrained: a parameter that exists and is never used is a
thing a later reader has to rule out.
"""

from __future__ import annotations

import math

import torch
from torch import nn

from snp2prot.models.encoders import (
    GENOMIC_WIDTHS,
    RESIDUE_WIDTHS,
    GenomicDNAEncoder,
    ProteinTower,
    ResidueProteinTower,
)
from snp2prot.models.two_tower import CALIBRATION_SCALE_INIT, _logit


class GenomicTwoTower(nn.Module):
    """Protein tower, genomic DNA tower, a calibration head and a temperature."""

    def __init__(
        self,
        protein_features: int = 1280,
        width: int = 256,
        dna_widths: tuple[int, ...] = GENOMIC_WIDTHS,
        dna_channels: int = 64,
        dna_dropout: float = 0.0,
        protein_hidden: int = 0,
        protein_dropout: float = 0.0,
        protein_tower: str = "pooled",
        residue_widths: tuple[int, ...] = RESIDUE_WIDTHS,
        residue_reduced: int = 64,
        residue_channels: int = 32,
        temperature: float = 0.07,
        learn_temperature: bool = True,
        max_logit_scale: float = 100.0,
        calibration_prior: float = 0.33,
    ):
        super().__init__()
        # Two shapes of protein tower behind one name, because everything downstream — the
        # cosine, the calibration head, the checkpoint — is identical between them and only the
        # call signature differs. `GHT_PLAN.md` §6 defaults this to `pooled` on capacity grounds.
        self.protein_tower = str(protein_tower)
        if self.protein_tower == "residue":
            self.protein = ResidueProteinTower(
                protein_features,
                width,
                reduced=residue_reduced,
                widths=residue_widths,
                channels_per_width=residue_channels,
                dropout=protein_dropout,
            )
        elif self.protein_tower == "pooled":
            self.protein = ProteinTower(protein_features, width, protein_hidden, protein_dropout)
        else:
            raise ValueError(f"protein_tower must be 'pooled' or 'residue', not {protein_tower!r}")
        self.dna = GenomicDNAEncoder(width, dna_widths, dna_channels, dna_dropout)
        logit_scale = torch.tensor(math.log(1.0 / temperature))
        self.logit_scale = nn.Parameter(logit_scale, requires_grad=learn_temperature)
        self.max_logit_scale = float(max_logit_scale)
        self.calibration_scale = nn.Parameter(torch.tensor(math.log(CALIBRATION_SCALE_INIT)))
        self.calibration_bias = nn.Parameter(torch.tensor(_logit(calibration_prior)))

    @torch.no_grad()
    def clamp_temperature(self) -> None:
        """Pin the learned scale to its ceiling after an optimiser step — CLIP's arrangement."""
        self.logit_scale.data.clamp_(max=math.log(self.max_logit_scale))

    @property
    def temperature(self) -> float:
        return float(1.0 / self.logit_scale.detach().exp())

    @property
    def temperature_is_clamped(self) -> bool:
        return bool(self.logit_scale.detach().exp() >= self.max_logit_scale * (1 - 1e-6))

    def embed_protein(self, vectors: torch.Tensor, mask: torch.Tensor | None = None):
        """`(B, 1280)` pooled, or `(B, L, 1280)` per residue with its mask -> `(B, D)`."""
        if self.protein_tower == "residue":
            if mask is None:
                raise ValueError("the per-residue tower needs a mask; a max over padding is a lie")
            return self.protein(vectors, mask)
        return self.protein(vectors)

    def cosine(
        self,
        protein_vectors: torch.Tensor,
        tokens: torch.Tensor,
        protein_mask: torch.Tensor | None = None,
    ) -> torch.Tensor:
        """Row-wise cosine of matched pairs: `(B, 1280)` and `(B, L)` -> `(B,)`.

        **Row-wise, not an outer product**, and that is the shape difference from the PBM arm in
        one line. There every domain is scored against the same complete vocabulary, so the step
        is a matmul. Here a window belongs to exactly one TF, so a batch is a list of *pairs* and
        scoring an unmatched pair would be asking a question the data cannot answer.
        """
        return (self.embed_protein(protein_vectors, protein_mask) * self.dna(tokens)).sum(dim=-1)

    def cosine_grouped(
        self,
        protein_vectors: torch.Tensor,
        group: torch.Tensor,
        tokens: torch.Tensor,
        protein_mask: torch.Tensor | None = None,
    ) -> torch.Tensor:
        """Score `(B, L)` windows against `U` distinct proteins, `group` saying which is which.

        **Each distinct protein is embedded once per pass, not once per window**, and for the
        per-residue tower that is the difference between running and not running. A batch of
        8,192 windows gathered per window is an `(8192, 219, 1280)` float tensor — 8.55 GB, which
        is larger than the card. Embedded per distinct protein it is `(33, 219, 1280)`, 37 MB.
        The pooled tower does not need the fix but gets it anyway: 33 rows instead of 8,192.

        Autograd is unaffected — indexing the embedded result is a gather, and the gradient sums
        back into the protein it came from, which is what a repeated protein should contribute.
        """
        embedded = self.embed_protein(protein_vectors, protein_mask)
        return (embedded[group] * self.dna(tokens)).sum(dim=-1)

    def calibrate(self, cosine: torch.Tensor) -> torch.Tensor:
        """Cosine -> the logit of `P(binds)`, on a scale shared across proteins.

        Its own scale and bias rather than the temperature, for the reason
        `snp2prot.models.two_tower.TwoTower.calibrate` gives: the temperature governs how peaked a
        softmax is, which is a statement about ranking dynamics, and it is clamped besides.
        """
        return self.calibration_scale.exp() * cosine + self.calibration_bias

    @torch.no_grad()
    def set_calibration_prior(self, positive_rate: float) -> None:
        """Move the calibration bias to `logit(positive_rate)` of the **training** pool.

        Never the test pool's rate and never a constant: with `shades` the base rate is about
        0.33 and with `random` about 0.09, so a head defaulted to the PBM arm's 0.002 would spend
        its first thousands of steps discovering that this corpus is not 466:1.
        """
        if not 0.0 < positive_rate < 1.0:
            raise ValueError(f"positive_rate must lie in (0, 1), not {positive_rate}")
        self.calibration_bias.fill_(_logit(positive_rate))

    def forward(
        self,
        protein_vectors: torch.Tensor,
        tokens: torch.Tensor,
        protein_mask: torch.Tensor | None = None,
    ) -> torch.Tensor:
        """`(B,)` calibrated logits — the decision-rule output, and what the loss reads."""
        return self.calibrate(self.cosine(protein_vectors, tokens, protein_mask))

    def parameter_counts(self) -> dict[str, int]:
        """Trainable parameters per tower. Reported with every result (`ML_PLAN.md` §4.2)."""
        return {
            "protein_tower": sum(p.numel() for p in self.protein.parameters() if p.requires_grad),
            "dna_tower": sum(p.numel() for p in self.dna.parameters() if p.requires_grad),
            "calibration": self.calibration_scale.numel() + self.calibration_bias.numel(),
            "total": sum(p.numel() for p in self.parameters() if p.requires_grad),
        }


def genomic_loss(
    logits: torch.Tensor,
    labels: torch.Tensor,
    cosine: torch.Tensor,
    group: torch.Tensor,
    logit_scale: torch.Tensor,
    contrastive_weight: float = 0.0,
) -> tuple[torch.Tensor, dict[str, float]]:
    """`L = BCE + mu * per-protein InfoNCE`, with `mu = 0` **exactly** plain BCE.

    The knob mirrors `bce_weight` in the PBM arm and for the same reason: one line, and the
    headline setting is the one that needs no extra term. The two-tower claim is architectural
    and does not need a contrastive loss to be true (`GHT_PLAN.md` §7).

    The contrastive term is a softmax **within one protein's rows of this batch** — its own
    positives against its own negatives, which is InfoNCE with `K` = that protein's negative
    count. Across proteins it would be meaningless: another TF's window is not a negative for
    this TF in any sense the data supports, it is simply unmeasured.

    `group` is the protein index of every row. A protein contributing no positive or no negative
    to this batch is skipped rather than contributing a degenerate term.
    """
    bce = nn.functional.binary_cross_entropy_with_logits(logits, labels.float())
    parts = {"bce": float(bce.detach()), "infonce": 0.0}
    if not contrastive_weight:
        return bce, parts

    scaled = cosine * logit_scale.exp()
    terms = []
    for g in torch.unique(group):
        rows = group == g
        positive = rows & (labels == 1)
        negative = rows & (labels == 0)
        if not (positive.any() and negative.any()):
            continue
        negatives = scaled[negative]
        # One softmax per positive, over that positive against every negative of the same
        # protein — `-log p(positive)` averaged over the protein's positives.
        spread = negatives.expand(int(positive.sum()), -1)
        candidates = torch.cat([scaled[positive].unsqueeze(1), spread], dim=1)
        terms.append(-torch.log_softmax(candidates, dim=1)[:, 0].mean())
    if not terms:
        return bce, parts
    infonce = torch.stack(terms).mean()
    parts["infonce"] = float(infonce.detach())
    return bce + contrastive_weight * infonce, parts
