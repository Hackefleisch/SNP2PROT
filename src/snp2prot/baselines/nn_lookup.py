"""Nearest-neighbour lookup baseline. Phase 5.

For a held-out DBD, copy the binary binding profile of the most sequence-similar DBD in
the training set (% identity over the aligned domain). Under leave-one-cluster-out this is
the number that decides whether the dataset supports learning at all: a model that cannot
beat it is a lookup table with extra steps, and we want to know that cheaply and early.
"""

from __future__ import annotations

import pandas as pd


def fit_predict(train: pd.DataFrame, test: pd.DataFrame) -> pd.Series:
    raise NotImplementedError("Phase 5")
