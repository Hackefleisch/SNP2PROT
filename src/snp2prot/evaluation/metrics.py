"""Evaluation metrics. TODO: decide the headline metric and stick to it."""

from __future__ import annotations

from typing import Any


def compute_metrics(y_true: Any, y_pred: Any) -> dict[str, float]:
    raise NotImplementedError("TODO")
