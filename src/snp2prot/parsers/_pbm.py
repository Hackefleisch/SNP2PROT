"""Machinery common to every universal-PBM source, whatever ships it.

A universal PBM scores a protein against all 32,896 non-redundant 8-mers and reports a
Wilcoxon-style **E-score**, a rank statistic bounded to [-0.5, 0.5]. That is true of
UniPROBE's per-accession archives and of the GEO deposits behind CIS-BP alike, so the parts
that depend only on the assay live here and the parts that depend on a distributor's file
layout stay in that distributor's module.

**The E-score column is identified, never assumed.** Reading it by position was a real defect
(open item #36): `usecols=[0, 2]` held for four UniPROBE accessions and silently returned
median intensity for two others and a Z-score for a third, briefly giving the corpus a 4:1
negative:positive ratio with `raw_score` up to 776,106. The identification below is the single
implementation of that rule; no source is allowed its own copy.
"""

from __future__ import annotations

import re

import numpy as np
import pandas as pd

from snp2prot import schema

#: A universal PBM E-score is a rank statistic bounded to this interval by construction.
#: The bound is what lets the column be identified when a file has no header.
ESCORE_MIN, ESCORE_MAX = -0.5, 0.5

#: A universal PBM scores every one of the non-redundant 8-mers, so a complete table has
#: exactly this many rows. Fewer means the depositor published only part of it.
N_NONREDUNDANT_8MERS = 32896

ASSAY = "PBM"
SCORE_TYPE = "pbm_escore"
#: A PBM 8-mer score aggregates over many flanking contexts, so no specific flank was observed.
DNA_CONTEXT = "core_only"

#: Header spellings seen so far: "E-score", "e_score", "enrichment score" (Cell08), "E-Score"
#: (GEO). A named column wins outright over the value-range rule.
ESCORE_HEADER_RE = re.compile(r"e[-_ ]?score|enrichment", re.I)


class EscoreColumnError(ValueError):
    """Raised when a table's E-score column cannot be identified unambiguously."""


def escore_column(df: pd.DataFrame, has_header: bool, label: str) -> pd.Series:
    """Return the E-score column of a raw 8-mer table, or raise saying why it cannot.

    `label` names the file or sample in the error, since a rejection is only useful if it
    says which input it is about.

    A header naming the column wins. Otherwise the column must be bounded to
    [-0.5, 0.5] **and** actually take negative values: the bound alone is not enough, because
    p-values and q-values sit in [0, 0.5] too, but an E-score over 32,896 8-mers always runs
    negative since most 8-mers are not bound, while a probability cannot.

    Zero candidates means no E-score is present; more than one means the file concatenates
    several experiments side by side. Both are errors, never a guess.
    """
    if df.shape[1] < 2:
        raise EscoreColumnError(f"{label}: only {df.shape[1]} column(s)")

    numeric = {}
    for i in range(1, df.shape[1]):
        vals = pd.to_numeric(df.iloc[:, i], errors="coerce")
        if vals.notna().sum() >= len(df) * 0.9:
            numeric[i] = vals

    if has_header:
        named = [i for i in numeric if ESCORE_HEADER_RE.search(str(df.columns[i]))]
        if len(named) == 1:
            return numeric[named[0]]

    candidates = [
        (i, v)
        for i, v in numeric.items()
        if v.min() >= ESCORE_MIN and v.max() <= ESCORE_MAX and v.min() < 0
    ]
    if not candidates:
        # Distinguish "no E-score here" from "an E-score column that never goes negative",
        # which means the depositor published only the enriched end of the table.
        bounded = [
            (i, v) for i, v in numeric.items() if v.min() >= ESCORE_MIN and v.max() <= ESCORE_MAX
        ]
        if bounded:
            i, v = bounded[0]
            raise EscoreColumnError(
                f"{label}: column {i} looks like an E-score but never goes negative "
                f"(min {v.min():.4f}, {len(df):,} rows). This is a TRUNCATED table listing "
                f"only enriched 8-mers, not the full {N_NONREDUNDANT_8MERS:,}. Labelling it "
                f"would call that protein's top hits non-binding"
            )
        raise EscoreColumnError(
            f"{label}: no column lies within [{ESCORE_MIN}, {ESCORE_MAX}] — this file "
            f"carries no E-score (columns: {list(df.columns)[:6]})"
        )
    if len(candidates) > 1:
        raise EscoreColumnError(
            f"{label}: {len(candidates)} columns look like E-scores (indices "
            f"{[i for i, _ in candidates]}); the file probably concatenates experiments"
        )
    return candidates[0][1]


def check_complete(n_rows: int, label: str) -> None:
    """Reject a partial 8-mer table, saying which kind of partial it is."""
    if n_rows == N_NONREDUNDANT_8MERS:
        return
    short = N_NONREDUNDANT_8MERS - n_rows
    kind = (
        "truncated to the enriched end"
        if n_rows < N_NONREDUNDANT_8MERS * 0.5
        else f"near-complete but {short} 8-mer(s) short"
    )
    raise EscoreColumnError(
        f"{label}: {n_rows:,} rows, expected {N_NONREDUNDANT_8MERS:,} — {kind}. The design is "
        f"fully crossed, so every protein must be scored against every 8-mer; admitting a "
        f"partial table would give that protein a different DNA axis"
    )


def binarize(escores: np.ndarray, pos_cut: float, neg_cut: float) -> np.ndarray:
    """Per-experiment binarization. Never call this on scores pooled across experiments."""
    lab = np.full(len(escores), schema.LABEL_GRAY, dtype=np.int8)
    lab[escores >= pos_cut] = schema.LABEL_BIND
    lab[escores <= neg_cut] = schema.LABEL_NONBIND
    return lab


def combine_replicates(
    tables: list[tuple[str, pd.DataFrame]], pos_cut: float, neg_cut: float
) -> tuple[pd.Series, np.ndarray, np.ndarray]:
    """Binarize each replicate independently, then combine. Returns (dna_seq, label, mean).

    Replicates can sit on different array designs whose intensity scales differ by an order of
    magnitude, so E-scores are never averaged *before* thresholding. Replicates that disagree
    on a given 8-mer fall to the gray band rather than being resolved by majority or by mean.

    Each table is `(label, frame with dna_seq + escore)`. They must cover the same 8-mers.
    """
    labels, scores, keys = [], [], None
    for label, exp in tables:
        exp = exp.sort_values("dna_seq", kind="stable").reset_index(drop=True)
        if keys is None:
            keys = exp["dna_seq"]
        elif not keys.equals(exp["dna_seq"]):
            raise ValueError(f"{label}: 8-mer set differs from the first replicate")
        e = exp["escore"].to_numpy(dtype=float)
        labels.append(binarize(e, pos_cut, neg_cut))
        scores.append(e)
    if not labels:
        raise EscoreColumnError("no readable replicate")

    stacked = np.vstack(labels)
    agreed = (stacked == stacked[0]).all(axis=0)
    label = np.where(agreed, stacked[0], schema.LABEL_GRAY).astype(np.int8)
    return keys, label, np.vstack(scores).mean(axis=0)


def build_frame(
    *,
    dna_seq: pd.Series,
    label: np.ndarray,
    raw_score: np.ndarray,
    dbd_seq: str,
    dbd_family: str,
    dbd_source: str,
    wt_id: str,
    n_mut_from_wt: int,
    mut_positions: str,
    protein_id: str,
    species: str,
    pos_cut: float,
    neg_cut: float,
    source_dataset: str,
    source_file: str,
) -> pd.DataFrame:
    """Assemble one experiment's rows. Scalars broadcast, so each string is stored once."""
    frame = pd.DataFrame({"dna_seq": dna_seq.to_numpy(), "label": label, "raw_score": raw_score})
    frame["dbd_seq"] = dbd_seq
    frame["dbd_family"] = dbd_family
    frame["dbd_source"] = dbd_source
    frame["wt_id"] = wt_id
    frame["n_mut_from_wt"] = n_mut_from_wt
    frame["mut_positions"] = mut_positions
    frame["protein_id"] = protein_id
    frame["species"] = species
    frame["dna_context"] = DNA_CONTEXT
    frame["score_type"] = SCORE_TYPE
    frame["threshold_pos"] = pos_cut
    frame["threshold_neg"] = neg_cut
    frame["assay"] = ASSAY
    frame["stringency"] = ""
    frame["neg_provenance"] = np.where(label == schema.LABEL_NONBIND, "assayed_unbound", None)
    frame["source_dataset"] = source_dataset
    frame["source_file"] = source_file
    return frame
