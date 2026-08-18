"""Reconciling the same domain measured by two sources, at merge time.

47 domains are stored by more than one source — 95 records, because one is stored by three.
Within a source, replicates are reconciled by agreement (`docs/METHODS.md` §6). Across
sources they are not, and cannot be: the median pair agrees on 46% of the 8-mers either
called positive, and two pairs agree on almost nothing (`C:LIN14B:NAP` 12 positives against
126, sharing none; `C:Cell08:Tlx2` 20 against 14, sharing three).

Left as they are, those records hand a sequence model **identical input with two different
labels**, and no split can separate them because they sit in one cluster by construction.

**The rule, decided by the owner on 2026-08-18: the record with more positives wins — unless
one of them belongs to a variant series, in which case the series wins.** The loser is dropped
from the merged table and stays in `data/interim/`, so nothing is destroyed and
`reports/overlap.md` still measures the noise floor from the full evidence.

The series exception covers 7 of the 48 dropped records and exists because the alternative
quietly breaks the thing the corpus is for. `BAR15A`'s `ARX_REF` has 188 positives against
`Cell08`'s 206, so more-positives would take `Cell08`'s copy — and leave ARX's five `BAR15A`
variants to be compared against a wild type measured by a different lab on a different array,
across a 46% cross-source noise floor that dwarfs any single-residue effect. Where one
candidate's source also supplies the other domains of that cluster, that source wins: a
variant series is measured end to end by one lab or it measures nothing.

Why more positives, rather than an agreement-style reconciliation:

* the failure mode of a PBM is *missing* binding, not inventing it — a weak array compresses
  its E-score distribution and calls fewer positives (`docs/METHODS.md` §5.1), so the deeper
  measurement is the more informative one;
* intersecting the two would inherit the worse array's sensitivity everywhere, and on
  `C:LIN14B:NAP` would leave a domain with no positives at all;
* it needs no new threshold, and it is auditable — every resolution is listed in
  `reports/overlap.md`.

Two properties of the corpus make this safe today, and both should be re-checked whenever a
source is added, because the rule is not safe in general:

1. **No duplicate pair disagrees about whether the protein binds at all.** All 95 records
   carry positives, so the rule never discards a measured non-binding (`label_health`'s
   `dead_variant`). If it ever would, that is a `T21` case and not a dedup case.
2. **No ties.** The margin is 1.8x at the median and 10x at the worst.
"""

from __future__ import annotations

import pandas as pd

#: Ordered, so a tie falls through to the next column and the result never depends on row
#: order. `in_series` comes first by the decision above; there are no ties in the corpus
#: today, and a rule that left them undefined would be a rebuild-to-rebuild difference
#: waiting to happen.
PRIORITY = ["in_series", "n_pos", "max_escore", "source_dataset"]


def duplicate_records(records: pd.DataFrame) -> pd.DataFrame:
    """The records whose `dbd_seq` is stored by more than one source."""
    return records[records.duplicated("dbd_seq", keep=False)]


def resolve(records: pd.DataFrame) -> pd.DataFrame:
    """Add `keep` and `in_series` to a per-record summary: one record survives per `dbd_seq`.

    `records` needs `dbd_seq`, `source_dataset`, `wt_id`, `n_pos` and `max_escore` — the
    columns `snp2prot.label_health` already writes, so the usual call is
    `merge.resolve(label_health.load())`.

    `in_series` marks a record whose source supplies more than one domain of its cluster,
    i.e. a wild type sitting beside its own variants.
    """
    out = records.copy()
    per_cluster = out.groupby(["wt_id", "source_dataset"])["dbd_seq"].transform("size")
    out["in_series"] = per_cluster > 1
    ranked = out.sort_values(
        PRIORITY, ascending=[False, False, False, True], kind="stable"
    ).drop_duplicates("dbd_seq", keep="first")
    out["keep"] = out.index.isin(ranked.index)
    return out


def deduplicate(df: pd.DataFrame, records: pd.DataFrame) -> pd.DataFrame:
    """Drop the losing records' rows from a row-level frame.

    Applied at merge, never at parse: `data/interim/` keeps every measurement, and the
    merged table keeps one per domain.
    """
    resolution = resolve(records)
    dropped = resolution[~resolution.keep]
    keys = set(zip(dropped.dbd_seq, dropped.source_dataset, strict=True))
    if not keys:
        return df
    pairs = zip(df["dbd_seq"], df["source_dataset"], strict=True)
    return df[[(seq, src) not in keys for seq, src in pairs]]
