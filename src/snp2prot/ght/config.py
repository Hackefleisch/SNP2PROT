"""Paths for the GHT-SELEX arm, beside `snp2prot.config` rather than inside it.

The same rule as everywhere else — nothing builds a path by hand — but kept in its own module
because this arm writes nowhere the PBM arm reads and reads nothing the PBM arm writes. The one
shared artifact is the protein language model itself, which is `snp2prot.config.MODEL_DIR`.
"""

from __future__ import annotations

from pathlib import Path

from snp2prot.config import EXTERNAL_DIR, INTERIM_DIR, PROCESSED_DIR, RAW_DIR, REPORTS_DIR

#: The Zenodo MEX-ArChIPelago benchmark release, `10.5281/zenodo.15555251` (`GHT_PLAN.md` D2).
MEX_DIR = RAW_DIR / "mex_archipelago"
DATASETS_ZIP = MEX_DIR / "datasets.zip"
PWMS_ZIP = MEX_DIR / "pwms.zip"
#: The 139 TFs with at least one approved ChIP-seq *and* one approved GHT-SELEX dataset.
PANEL_FILE = MEX_DIR / "TFs_CHS+AFS.txt"

#: Codebook metadata: which experiments were approved, and the construct amino-acid sequences.
GHT_METADATA = RAW_DIR / "codebook_ghtselex" / "GHT_SELEX_Metadata_web_2025_09_22.xlsx"
TF_METADATA = RAW_DIR / "codebook_htselex" / "TF_and_Plasmid_Metadata_web_2025_06_06.xlsx"

#: hg38, plain gzip (not BGZF), so it is read in one pass rather than indexed.
HG38 = EXTERNAL_DIR / "genomes" / "hg38.fa.gz"

GHT_INTERIM = INTERIM_DIR / "ght"
GHT_PROCESSED = PROCESSED_DIR / "ght"

#: One row per TF that entered the panel: the construct, the padded domain, the admission call.
PANEL_TABLE = GHT_INTERIM / "panel.parquet"
#: One row per panel TF that did NOT enter, with the rejection reason. Reported, never repaired.
REJECTED_TABLE = GHT_INTERIM / "rejected.parquet"

#: One row per 301 bp window: tf, split, negative set, locus, and the sequence as uint8 tokens.
#: Written per TF so a rebuild of one TF does not rewrite the corpus.
WINDOW_DIR = GHT_INTERIM / "windows"

#: Pooled ESM-2 vectors, one per panel domain. Same shape as the PBM arm's, different corpus.
EMBEDDING_DIR = GHT_PROCESSED / "embeddings"

CHECKPOINT_DIR = GHT_PROCESSED / "checkpoints"
#: Per-fold summary rows and per-TF rows, mirroring `training_folds` / `training_domains`.
FOLD_TABLE = GHT_PROCESSED / "ght_folds.parquet"
TF_TABLE = GHT_PROCESSED / "ght_tfs.parquet"
#: The §8.1 pre-flight: step against train loss, validation and **test**, for one fold.
PREFLIGHT_TABLE = GHT_PROCESSED / "ght_preflight.parquet"
BASELINE_TABLE = GHT_PROCESSED / "ght_baselines.parquet"

REPORT_DIR = REPORTS_DIR


def window_file(tf: str) -> Path:
    """One TF's windows, every split and every negative set in one file."""
    return WINDOW_DIR / f"{tf}.npz"


def embedding_table(arm: str) -> Path:
    return EMBEDDING_DIR / f"{arm}.npz"


def residue_embedding_table(arm: str) -> Path:
    return EMBEDDING_DIR / f"{arm}_residues.npz"


def checkpoint_file(arm: str, regime: str, fold: str, seed: int, tag: str = "") -> Path:
    """Trained weights for one (arm, fold, seed, configuration).

    **`tag` is not decoration.** The grid trains the same (arm, fold, seed) several times over —
    once for the model, once protein-blind, once protein-shuffled, once from the PBM arm's tower,
    once with the per-residue tower — and without a distinguishing component every one of those
    would overwrite the last, leaving a directory of checkpoints whose names say nothing about
    which configuration produced them. It is the same guard `snp2prot.config.checkpoint_file`
    applies to `bce_weight` in the PBM arm, for the same reason.

    `:` is legal in a fold name and awkward in a filename, so it becomes `-`; the fold's real name
    is stored inside the checkpoint, which is what a reader should trust.
    """
    suffix = f"_{tag}" if tag else ""
    return CHECKPOINT_DIR / f"ght_{arm}_{regime}_{fold.replace(':', '-')}_seed{seed}{suffix}.pt"


def run_tag(protein_mode: str = "real", tower: str = "pooled") -> str:
    """A filename-safe stamp for one protein configuration: `""`, `constant`, `pbm_frozen-res`.

    Empty for the default configuration, so the model's own checkpoints keep the plain names the
    rest of the arm refers to and only the variants carry a suffix.
    """
    parts = [] if protein_mode == "real" else [protein_mode]
    if tower != "pooled":
        parts.append(tower)
    return "-".join(parts)
