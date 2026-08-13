"""Access to `configs/thresholds.yaml`.

Parsers call `load()` and read the block for their assay. A parser that writes a numeric
cutoff into its own source is a bug — the value belongs in the YAML so the owner can move
it without touching code, and so every row can record what was actually applied.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

import yaml

from snp2prot.config import CONFIG_DIR

DEFAULT_THRESHOLD_FILE = CONFIG_DIR / "thresholds.yaml"


def load(path: str | Path | None = None) -> dict[str, Any]:
    """Load the threshold config as a plain dict."""
    p = Path(path) if path is not None else DEFAULT_THRESHOLD_FILE
    if not p.exists():
        raise FileNotFoundError(f"threshold config not found: {p}")
    with p.open() as fh:
        return yaml.safe_load(fh)


def for_assay(assay_key: str, path: str | Path | None = None) -> dict[str, Any]:
    """Return one assay's block, e.g. `for_assay("pbm")`."""
    cfg = load(path)
    if assay_key not in cfg:
        raise KeyError(f"no threshold block for {assay_key!r}; have {sorted(cfg)}")
    return cfg[assay_key]
