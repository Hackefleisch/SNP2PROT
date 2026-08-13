"""Split protocols. Phase 5.

Three regimes, always reported separately — a single pooled number hides exactly the
failure this project is designed to detect (brief §5b):

- `leave_one_variant_out`  — hold out single variants within a cluster. Tests DNA-side and
  fine protein-side generalization. Optimistic by construction.
- `leave_one_cluster_out`  — hold out whole `wt_id` clusters. Tests whether the model has
  learned anything transferable across protein sequence space.
- `leave_one_family_out`   — hold out whole `dbd_family` groups. The hardest regime.
"""

from __future__ import annotations

import pandas as pd


def leave_one_variant_out(df: pd.DataFrame, seed: int = 42):
    raise NotImplementedError("Phase 5")


def leave_one_cluster_out(df: pd.DataFrame, seed: int = 42):
    raise NotImplementedError("Phase 5")


def leave_one_family_out(df: pd.DataFrame, seed: int = 42):
    raise NotImplementedError("Phase 5")
