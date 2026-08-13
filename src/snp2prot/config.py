"""Project paths.

Single source of truth for where things live. No other module builds a path by hand or
uses a relative path — import from here so that scripts, notebooks and tests all agree
regardless of the working directory they are launched from.

Directory contract (brief §1.3):
  data/raw/<source>/       append-only, never edited, never written by a parser
  data/interim/<source>/   per-source parsed output, schema-conforming Parquet
  data/processed/          merged training table(s)
  data/testsets/           Tier 4 held-out sets, physically separate from training data
"""

from __future__ import annotations

from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[2]

DATA_DIR = PROJECT_ROOT / "data"
RAW_DIR = DATA_DIR / "raw"
INTERIM_DIR = DATA_DIR / "interim"
PROCESSED_DIR = DATA_DIR / "processed"
EXTERNAL_DIR = DATA_DIR / "external"
#: Tier 4 evaluation sets. Deliberately NOT under processed/ — nothing that reads the
#: training table should be able to reach these by globbing.
TESTSET_DIR = DATA_DIR / "testsets"

CONFIG_DIR = PROJECT_ROOT / "configs"
REPORTS_DIR = PROJECT_ROOT / "reports"
DOCS_DIR = PROJECT_ROOT / "docs"

RESULTS_DIR = PROJECT_ROOT / "results"
FIGURE_DIR = RESULTS_DIR / "figures"
CHECKPOINT_DIR = RESULTS_DIR / "checkpoints"

PROVENANCE_FILE = PROJECT_ROOT / "PROVENANCE.md"


def raw_dir(source: str, create: bool = False) -> Path:
    """Raw download directory for one source, e.g. `raw_dir("BAR15A")`."""
    p = RAW_DIR / source
    if create:
        p.mkdir(parents=True, exist_ok=True)
    return p


def interim_dir(source: str, create: bool = False) -> Path:
    """Parsed-output directory for one source."""
    p = INTERIM_DIR / source
    if create:
        p.mkdir(parents=True, exist_ok=True)
    return p


def interim_table(source: str) -> Path:
    """Canonical schema-conforming Parquet for one source."""
    return interim_dir(source) / f"{source}.parquet"


def relative_to_raw(path: str | Path) -> str:
    """Render a path as the `source_file` column wants it: relative to `data/raw/`."""
    return str(Path(path).resolve().relative_to(RAW_DIR))
