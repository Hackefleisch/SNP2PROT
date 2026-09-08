"""The two-tower contrastive model. `docs/TRAINING.md` is the design document."""

from snp2prot.models.encoders import DNAEncoder, ProteinTower, tokenise
from snp2prot.models.loss import calibration_bce, hybrid_loss, multi_positive_infonce
from snp2prot.models.two_tower import TwoTower

__all__ = [
    "DNAEncoder",
    "ProteinTower",
    "TwoTower",
    "calibration_bce",
    "hybrid_loss",
    "multi_positive_infonce",
    "tokenise",
]
