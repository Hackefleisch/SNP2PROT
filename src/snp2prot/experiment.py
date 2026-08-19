"""Access to `configs/experiment.yaml`, the modelling config.

The mirror of `snp2prot.thresholds`, and separate from it on purpose: changing a binarization
cutoff invalidates the dataset, every report and every provenance row, while changing a fold
count or a learning rate does not (`docs/ML_PLAN.md` §9.2, decided 2026-08-19). Two things
with different blast radii do not belong in one file.

Seeds are read from here rather than passed at a call site, so that "which seed produced this
number" has one answer and it is recorded with the run.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml

from snp2prot.config import CONFIG_DIR

DEFAULT_EXPERIMENT_FILE = CONFIG_DIR / "experiment.yaml"


def load(path: str | Path | None = None) -> dict[str, Any]:
    """Load the modelling config as a plain dict."""
    p = Path(path) if path is not None else DEFAULT_EXPERIMENT_FILE
    if not p.exists():
        raise FileNotFoundError(f"experiment config not found: {p}")
    with p.open() as fh:
        return yaml.safe_load(fh)


def section(key: str, path: str | Path | None = None) -> dict[str, Any]:
    """Return one block, e.g. `section("splits")`."""
    cfg = load(path)
    if key not in cfg:
        raise KeyError(f"no {key!r} block in the experiment config; have {sorted(cfg)}")
    return cfg[key]
