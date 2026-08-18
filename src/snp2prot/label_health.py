"""Which domain records carry usable positive evidence, and which only look like they do.

`E >= 0.45` is a statement about ranks, and it is comparable across experiments in the way it
claims to be: it survives any monotone change of laser power, protein concentration or
scanner. What it does not equalise is **sensitivity**. A noisy or weak array shuffles the
probe ranking slightly, every shuffle costs the statistic, and the whole distribution
compresses toward zero. Two labs measuring one protein at a fixed cutoff differ 1.6x in the
number of positives they call, and in the tail far worse — `Hoxa2` is 165 against 17
(`reports/threshold_review.md`).

At the extreme a whole source falls under the cutoff. `Cell09`'s per-domain maximum E-score
has a median of 0.428: 12 of its 17 domains never reach 0.45 at their single best 8-mer, so
they contribute rows that say only "does not bind". That is not the same claim as a variant
that has genuinely lost binding, and the difference matters because the second is the signal
this dataset exists to carry.

**Nothing here changes a label.** The cutoff stays at 0.45 (`T20`, closed 2026-08-18: it is
the best absolute cutoff measured, and the alternatives break `neg_provenance`). This module
records, per domain record, whether the record has positive evidence at all and — when it does
not — whether anything nearby does. Filtering on it is a *training* decision, made at
featurization, exactly as cluster-size restriction is (`docs/DECISIONS.md`, 2026-08-14).

## The three verdicts

A "record" is one `(dbd_seq, source_dataset)` pair: one domain as one source measured it.

- **`ok`** — the record has at least one 8-mer at or above the cutoff.
- **`dead_variant`** — no positives, but another record in the same cluster **and the same
  source** does. Same series, same lab, same array design, so its own relatives are the
  control: this is a measurement of a protein that does not bind, and it is evidence.
- **`no_evidence`** — no positives and no such control. Indistinguishable, from E-scores
  alone, from an experiment that was too weak to see anything. Keep the rows, know what they
  are, and let the training filter decide.

The discriminator is deliberately conservative: it asks only for a positive *somewhere in the
same cluster and source*, which is the narrowest control that means anything.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

from snp2prot import schema, thresholds
from snp2prot.config import LABEL_HEALTH_TABLE, PROTEIN_TABLE

OK = "ok"
DEAD_VARIANT = "dead_variant"
NO_EVIDENCE = "no_evidence"


def summarise(df: pd.DataFrame, cutoff: float | None = None) -> pd.DataFrame:
    """One row per `(dbd_seq, source_dataset)` with its label counts and best E-score.

    Takes a parsed source table (or several concatenated) and returns the per-record
    statistics. `verdict` is filled in by `classify`, which needs the whole corpus.

    `n_at_cutoff` counts 8-mers whose stored score reaches the cutoff, which is **not** the
    same as `n_pos`: labels are binarized per experiment and replicates that disagree become
    `LABEL_GRAY`, while `raw_score` is their mean. A record with `n_at_cutoff > 0` and
    `n_pos == 0` was not too weak to see — its replicates did not agree on what they saw.
    """
    cut = float(thresholds.for_assay("pbm")["positive"]) if cutoff is None else cutoff
    df = df.assign(_at_cutoff=(df["raw_score"] >= cut).astype("int64"))
    grouped = df.groupby(["dbd_seq", "source_dataset"], sort=False)
    out = grouped.agg(
        protein_id=("protein_id", "first"),
        wt_id=("wt_id", "first"),
        dbd_family=("dbd_family", "first"),
        n_rows=("label", "size"),
        max_escore=("raw_score", "max"),
        n_at_cutoff=("_at_cutoff", "sum"),
    )
    labels = grouped["label"].value_counts().unstack(fill_value=0)
    for value, name in (
        (schema.LABEL_BIND, "n_pos"),
        (schema.LABEL_NONBIND, "n_neg"),
        (schema.LABEL_GRAY, "n_gray"),
    ):
        out[name] = labels[value] if value in labels.columns else 0
    return out.reset_index()


def classify(records: pd.DataFrame) -> pd.DataFrame:
    """Add `verdict` and `has_control`, judging each record against its own cluster.

    A record with no positives is `dead_variant` only if another record of the same `wt_id`
    **from the same source** has some — the control has to have been measured by the same lab
    on the same array design, or it is not a control.
    """
    records = records.copy()
    with_positives = records[records.n_pos > 0]
    controls = set(zip(with_positives.wt_id, with_positives.source_dataset, strict=True))
    records["has_control"] = [
        (wt, src) in controls for wt, src in zip(records.wt_id, records.source_dataset, strict=True)
    ]
    records["verdict"] = OK
    silent = records.n_pos == 0
    records.loc[silent & records.has_control, "verdict"] = DEAD_VARIANT
    records.loc[silent & ~records.has_control, "verdict"] = NO_EVIDENCE

    # How deep the source that made this measurement sees at all. Carried per record because
    # `no_evidence` means two different things depending on it: one silent domain in a source
    # whose median domain reaches 0.495 is most likely a protein that binds nothing, while a
    # silent domain in a source whose median is 0.428 is a source that cannot see.
    records["source_median_max_escore"] = records.groupby("source_dataset").max_escore.transform(
        "median"
    )
    return records


def with_gene_names(records: pd.DataFrame, path: str | Path | None = None) -> pd.DataFrame:
    """Attach the protein table's `gene`, for a report a person can read.

    `protein_id` is empty for every CIS-BP record — Table S6 publishes gene and species but no
    accession (`docs/DECISIONS.md`, `T2b`) — so the accession alone cannot name a record. The
    protein table covers 1,224 of 1,335 domains (`T18`); anything it misses keeps whatever
    `protein_id` it has.
    """
    p = Path(path) if path else PROTEIN_TABLE
    records = records.copy()
    if not p.exists():
        records["gene"] = records["protein_id"]
        return records
    proteins = pd.read_parquet(p, columns=["dbd_seq", "source_dataset", "gene"])
    proteins = proteins.drop_duplicates(["dbd_seq", "source_dataset"])
    merged = records.merge(proteins, on=["dbd_seq", "source_dataset"], how="left")
    merged["gene"] = merged["gene"].fillna(merged["protein_id"]).replace("", pd.NA)
    merged["gene"] = merged["gene"].fillna("(unnamed)")
    return merged


def source_summary(records: pd.DataFrame) -> pd.DataFrame:
    """Per source: how much of it is silent, and how high it reaches when it is not.

    **`frac_silent` is the honest source-level statistic**, not `median_max_escore`. The
    median is dragged down by a source that deliberately contains dead variants — on a
    synthetic source that is half dead it falls below the cutoff, and `BAR15A` only escapes
    that because its silent fraction is 20%. `Cell09` stands out on both: 71% silent, and a
    median max E of 0.428 against 0.495 for the sources around it.

    Neither figure decides a record's verdict — that is `classify`, and it turns on whether
    the record has a control in its own cluster and source. These columns are context for
    reading the verdicts, and for noticing a source that should be looked at again.
    """
    per_source = records.groupby("source_dataset")
    out = per_source.agg(
        n_records=("dbd_seq", "size"),
        median_max_escore=("max_escore", "median"),
        n_pos=("n_pos", "sum"),
        n_rows=("n_rows", "sum"),
    )
    out["positive_rate"] = out.n_pos / out.n_rows
    for verdict in (DEAD_VARIANT, NO_EVIDENCE):
        out[f"n_{verdict}"] = per_source.verdict.apply(lambda s, v=verdict: int((s == v).sum()))
    out["n_silent"] = out[f"n_{DEAD_VARIANT}"] + out[f"n_{NO_EVIDENCE}"]
    out["frac_silent"] = out.n_silent / out.n_records
    return out.reset_index().sort_values("frac_silent", ascending=False)


def load(path: str | Path | None = None) -> pd.DataFrame:
    """The written table. Built by `scripts/build_label_health.py`."""
    p = Path(path) if path else LABEL_HEALTH_TABLE
    if not p.exists():
        raise FileNotFoundError(f"{p} not found — run scripts/build_label_health.py")
    return pd.read_parquet(p)


def usable(df: pd.DataFrame, path: str | Path | None = None, keep_dead: bool = True):
    """Drop the rows of records with no positive evidence and no control.

    The training-time filter this table exists for. `keep_dead=True` keeps `dead_variant`
    records, which is the point of distinguishing them: a variant that lost binding is the
    most informative negative in the corpus.
    """
    records = load(path)
    drop = {NO_EVIDENCE} if keep_dead else {NO_EVIDENCE, DEAD_VARIANT}
    excluded = records[records.verdict.isin(drop)]
    keys = set(zip(excluded.dbd_seq, excluded.source_dataset, strict=True))
    if not keys:
        return df
    pairs = zip(df["dbd_seq"], df["source_dataset"], strict=True)
    mask = [(seq, src) not in keys for seq, src in pairs]
    return df[mask]
