"""Config loading and project paths.

Single source of truth for where things live, so no module hardcodes a path.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

PROJECT_ROOT = Path(__file__).resolve().parents[2]
DATA_DIR = PROJECT_ROOT / "data"
RAW_DIR = DATA_DIR / "raw"
INTERIM_DIR = DATA_DIR / "interim"
PROCESSED_DIR = DATA_DIR / "processed"
RESULTS_DIR = PROJECT_ROOT / "results"
CHECKPOINT_DIR = RESULTS_DIR / "checkpoints"
FIGURE_DIR = RESULTS_DIR / "figures"
CONFIG_DIR = PROJECT_ROOT / "configs"


@dataclass
class Config:
    """Run configuration. TODO: extend as the pipeline solidifies."""

    name: str = "default"
    seed: int = 42
    data: dict[str, Any] = field(default_factory=dict)
    model: dict[str, Any] = field(default_factory=dict)
    training: dict[str, Any] = field(default_factory=dict)


def load_config(path: str | Path) -> Config:
    """Load a YAML config file into a `Config`."""
    raise NotImplementedError("TODO: implement once the config schema is settled")
