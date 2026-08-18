"""Resolving one domain measured by two sources.

47 domains are stored twice and their labels do not agree — the median pair on 46% of the
8-mers either called positive. One record has to win, and which one is a decision with a
sharp edge: taking the deeper measurement everywhere would, in seven clusters, replace a
variant series' own wild type with another lab's copy.
"""

from __future__ import annotations

import pandas as pd

from snp2prot import merge


def rec(dbd: str, source: str, wt: str, n_pos: int, max_escore: float = 0.49, n_rows: int = 100):
    return {
        "dbd_seq": dbd,
        "source_dataset": source,
        "wt_id": wt,
        "n_pos": n_pos,
        "max_escore": max_escore,
        "n_rows": n_rows,
    }


def kept(resolution: pd.DataFrame) -> set[tuple[str, str]]:
    survivors = resolution[resolution.keep]
    return set(zip(survivors.dbd_seq, survivors.source_dataset, strict=True))


def test_the_deeper_measurement_wins():
    records = pd.DataFrame([rec("AAA", "DEEP", "C:X:1", 126), rec("AAA", "SHALLOW", "C:X:1", 12)])
    assert kept(merge.resolve(records)) == {("AAA", "DEEP")}


def test_a_domain_measured_once_always_survives():
    records = pd.DataFrame([rec("AAA", "ONLY", "C:X:1", 3)])
    assert kept(merge.resolve(records)) == {("AAA", "ONLY")}


def test_a_variant_series_keeps_its_own_wild_type():
    """`BAR15A`'s `ARX_REF` has 188 positives against `Cell08`'s 206. Taking Cell08's would
    leave ARX's five BAR15A variants compared against another lab's wild type, across a noise
    floor far larger than the substitution being measured."""
    records = pd.DataFrame(
        [
            rec("WT", "SERIES", "C:S:ARX", 188),
            rec("WT", "DEEPER", "C:S:ARX", 206),
            rec("VAR1", "SERIES", "C:S:ARX", 150),
            rec("VAR2", "SERIES", "C:S:ARX", 0),
        ]
    )
    resolution = merge.resolve(records)
    assert kept(resolution) == {("WT", "SERIES"), ("VAR1", "SERIES"), ("VAR2", "SERIES")}
    flags = resolution.set_index(["dbd_seq", "source_dataset"])["in_series"]
    assert flags[("WT", "SERIES")] and not flags[("WT", "DEEPER")]


def test_the_series_rule_only_applies_when_the_source_supplies_more_of_the_cluster():
    """One source holding one domain of a cluster is not a series — positives decide."""
    records = pd.DataFrame([rec("WT", "LONE", "C:S:X", 10), rec("WT", "DEEPER", "C:S:X", 99)])
    assert kept(merge.resolve(records)) == {("WT", "DEEPER")}


def test_ties_resolve_deterministically():
    """No tie exists in the corpus today, and one must not depend on row order if it appears."""
    a = pd.DataFrame([rec("AAA", "BBB", "C:X:1", 5, 0.47), rec("AAA", "AAA_SRC", "C:X:1", 5, 0.47)])
    assert kept(merge.resolve(a)) == kept(merge.resolve(a.iloc[::-1]))
    assert kept(merge.resolve(a)) == {("AAA", "AAA_SRC")}  # source name breaks the last tie

    by_score = pd.DataFrame(
        [rec("AAA", "DIM", "C:X:1", 5, 0.46), rec("AAA", "BRIGHT", "C:X:1", 5, 0.49)]
    )
    assert kept(merge.resolve(by_score)) == {("AAA", "BRIGHT")}


def test_deduplicate_drops_the_losing_rows_only():
    records = pd.DataFrame([rec("AAA", "DEEP", "C:X:1", 126), rec("AAA", "SHALLOW", "C:X:1", 12)])
    rows = pd.DataFrame(
        {
            "dbd_seq": ["AAA", "AAA", "BBB"],
            "source_dataset": ["DEEP", "SHALLOW", "SHALLOW"],
            "label": [1, 0, 1],
        }
    )
    out = merge.deduplicate(rows, records)
    assert list(zip(out.dbd_seq, out.source_dataset, strict=True)) == [
        ("AAA", "DEEP"),
        ("BBB", "SHALLOW"),
    ]
