"""Which records carry positive evidence.

The distinction under test is the one the table exists for: a variant that measurably lost
binding is evidence, and a domain from an experiment that saw nothing is not, even though
both look identical in the label column.
"""

from __future__ import annotations

import pandas as pd
import pytest

from snp2prot import label_health, schema

CUTOFF = 0.45


def record(dbd: str, source: str, wt: str, scores: list[float], labels: list[int]):
    """Rows for one (domain, source) record, one per 8-mer."""
    return pd.DataFrame(
        {
            "dbd_seq": dbd,
            "source_dataset": source,
            "protein_id": f"P_{dbd}",
            "wt_id": wt,
            "dbd_family": "Homeodomain",
            "raw_score": scores,
            "label": labels,
        }
    )


BINDER = record("AAA", "SRC", "C:SRC:WT", [0.48, 0.46, 0.10], [1, 1, 0])
#: Same series as BINDER, and nothing reaches the cutoff — the shape of a dead variant.
SILENT_SIBLING = record("AAB", "SRC", "C:SRC:WT", [0.31, 0.20, 0.10], [0, 0, 0])
#: Its own cluster, its own source, nothing above the cutoff anywhere near it.
LONE_SILENT = record("CCC", "WEAK", "C:WEAK:LONE", [0.43, 0.22, 0.05], [0, 0, 0])


def build(*frames):
    df = pd.concat(frames, ignore_index=True)
    return label_health.classify(label_health.summarise(df, cutoff=CUTOFF))


def test_counts_come_out_per_record():
    out = label_health.summarise(pd.concat([BINDER, SILENT_SIBLING]), cutoff=CUTOFF)
    assert len(out) == 2
    binder = out[out.dbd_seq == "AAA"].iloc[0]
    assert (binder.n_pos, binder.n_neg, binder.n_rows) == (2, 1, 3)
    assert binder.max_escore == pytest.approx(0.48)


def test_a_silent_record_with_a_binder_in_its_own_series_is_a_dead_variant():
    out = build(BINDER, SILENT_SIBLING, LONE_SILENT)
    verdicts = dict(zip(out.dbd_seq, out.verdict, strict=True))
    assert verdicts["AAA"] == label_health.OK
    assert verdicts["AAB"] == label_health.DEAD_VARIANT
    assert verdicts["CCC"] == label_health.NO_EVIDENCE


def test_the_control_has_to_come_from_the_same_source():
    """A different lab measuring the same cluster is not a control for this experiment —
    that is the whole point of `per_experiment` binarization."""
    other_lab = record("AAA", "OTHER", "C:SRC:WT", [0.49, 0.47, 0.10], [1, 1, 0])
    out = build(other_lab, SILENT_SIBLING)
    keys = zip(out.dbd_seq, out.source_dataset, strict=True)
    verdicts = dict(zip(keys, out.verdict, strict=True))
    assert verdicts[("AAB", "SRC")] == label_health.NO_EVIDENCE


def test_reaching_the_cutoff_without_a_positive_label_is_recorded_separately():
    """13 real records look silent only because their replicates disagreed: the stored score
    is the replicate mean, the label is `LABEL_GRAY`."""
    disputed = record("DDD", "SRC", "C:SRC:DDD", [0.46, 0.30], [schema.LABEL_GRAY, 0])
    out = label_health.summarise(disputed, cutoff=CUTOFF).iloc[0]
    assert out.n_pos == 0
    assert out.n_at_cutoff == 1
    assert out.n_gray == 1


def test_source_health_separates_a_weak_source_from_a_dead_variant():
    out = build(BINDER, SILENT_SIBLING, LONE_SILENT)
    per_source = label_health.source_summary(out).set_index("source_dataset")
    assert per_source.loc["SRC", f"n_{label_health.DEAD_VARIANT}"] == 1
    assert per_source.loc["WEAK", f"n_{label_health.NO_EVIDENCE}"] == 1
    assert per_source.loc["WEAK", "frac_silent"] == 1.0
    # And the caveat the docstring makes: a source that is half dead variants has a median
    # below the cutoff without being blind, which is why frac_silent is the statistic to read.
    assert per_source.loc["SRC", "median_max_escore"] < CUTOFF
    assert per_source.loc["SRC", "frac_silent"] == 0.5


def test_usable_drops_unsupported_records_and_keeps_dead_variants(tmp_path):
    out = build(BINDER, SILENT_SIBLING, LONE_SILENT)
    table = tmp_path / "label_health.parquet"
    out.to_parquet(table, index=False)
    rows = pd.concat([BINDER, SILENT_SIBLING, LONE_SILENT], ignore_index=True)

    kept = label_health.usable(rows, path=table)
    assert set(kept.dbd_seq) == {"AAA", "AAB"}

    strict = label_health.usable(rows, path=table, keep_dead=False)
    assert set(strict.dbd_seq) == {"AAA"}


def test_load_says_what_to_run_when_the_table_is_missing(tmp_path):
    with pytest.raises(FileNotFoundError, match="build_label_health"):
        label_health.load(tmp_path / "absent.parquet")
