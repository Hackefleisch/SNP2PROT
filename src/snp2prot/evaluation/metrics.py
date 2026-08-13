"""Metrics for scoring the nearest-neighbour baseline and, later, models. Phase 5.

Always reported per split regime, never pooled: a single number across
leave-one-variant-out and leave-one-cluster-out hides exactly the generalization gap the
project exists to measure.
"""

from __future__ import annotations

from typing import Any


def compute_metrics(y_true: Any, y_pred: Any) -> dict[str, float]:
    raise NotImplementedError("Phase 5")
