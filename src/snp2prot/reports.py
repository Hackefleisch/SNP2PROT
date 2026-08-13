"""Generation of the markdown reports under `reports/`.

These are review artifacts for the project owner, not logs: each is regenerated from the
tables and committed, so a diff shows how the dataset changed between phases.

- `binarization_summary_<source>.md` — n_pos / n_neg / n_gray, positive rate, and the
  `raw_score` distribution with the applied cutoffs marked (brief §4).
- `cluster_inventory.md`  — per `wt_id`: variant count and `n_mut_from_wt` distribution;
  per family: cluster count. Flags any family contributing >40% of rows (brief §5b).
- `overlap.md`            — (dbd_seq, dna_seq) pairs seen in more than one source and
  whether their labels agree; the disagreement rate reads directly as label noise.
"""

from __future__ import annotations

import pandas as pd

#: A family above this share of total rows gets flagged — C2H2 will otherwise swamp
#: everything and the model can win by memorizing one family.
FAMILY_DOMINANCE_FLAG = 0.40


def binarization_summary(df: pd.DataFrame, source: str) -> str:
    raise NotImplementedError("Phase 1")


def cluster_inventory(df: pd.DataFrame) -> str:
    raise NotImplementedError("Phase 1")


def overlap_report(frames: dict[str, pd.DataFrame]) -> str:
    raise NotImplementedError("Phase 5")
