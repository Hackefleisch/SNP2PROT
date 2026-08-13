"""The validator is the gate every parser passes through, so it gets tested hardest."""

from __future__ import annotations

import pandas as pd
import pytest

from snp2prot import schema


def test_good_frame_passes(good_frame):
    rep = schema.validate(good_frame, source="TESTSRC")
    assert rep.ok, str(rep)
    assert rep.stats["n_pos"] == 2
    assert rep.stats["n_neg"] == 2
    assert rep.stats["n_clusters"] == 1


def test_pair_id_is_deterministic_and_order_free():
    a = schema.make_pair_id("ACDE", "TAAT", "PBM", "")
    assert a == schema.make_pair_id("ACDE", "TAAT", "PBM", "")
    # Field boundaries must not be collapsible by concatenation.
    assert schema.make_pair_id("AC", "DETAAT", "PBM", "") != a


def test_stringency_separates_otherwise_identical_rows():
    lo = schema.make_pair_id("ACDE", "TAAT", "B1H", "2mM_3AT")
    hi = schema.make_pair_id("ACDE", "TAAT", "B1H", "10mM_3AT")
    assert lo != hi


def test_empty_frame_matches_schema():
    df = schema.empty_frame()
    assert list(df.columns) == list(schema.COLUMN_NAMES)


def test_coerce_rejects_invented_columns(good_frame):
    df = good_frame.assign(my_extra_column=1)
    with pytest.raises(ValueError, match="not in the schema"):
        schema.coerce(df)


def test_coerce_rejects_missing_required_column(good_frame):
    with pytest.raises(ValueError, match="required column missing"):
        schema.coerce(good_frame.drop(columns=["dbd_seq"]))


def test_negative_without_provenance_fails(good_frame):
    df = good_frame.copy()
    df.loc[df["label"] == 0, "neg_provenance"] = pd.NA
    rep = schema.validate(df, source="TESTSRC")
    assert not rep.ok
    assert any("neg_provenance" in e for e in rep.errors)


def test_unknown_neg_provenance_fails(good_frame):
    df = good_frame.copy()
    df.loc[df["label"] == 0, "neg_provenance"] = "probably_not_binding"
    rep = schema.validate(df, source="TESTSRC")
    assert not rep.ok


def test_padded_dna_is_rejected(good_frame):
    """A padded site must never reach the stored table (brief §3)."""
    df = good_frame.copy()
    df.loc[0, "dna_seq"] = "NNTAATTAGCNN"
    df.loc[0, "dna_len"] = 12
    rep = schema.validate(df, source="TESTSRC")
    assert not rep.ok
    assert any("non-ACGT" in e for e in rep.errors)


def test_dna_over_20bp_is_rejected(good_frame):
    df = good_frame.copy()
    df.loc[0, "dna_seq"] = "ACGT" * 6
    df.loc[0, "dna_len"] = 24
    rep = schema.validate(df, source="TESTSRC")
    assert any("longer than 20" in e for e in rep.errors)


def test_dna_len_must_match(good_frame):
    df = good_frame.copy()
    df.loc[0, "dna_len"] = 7
    rep = schema.validate(df, source="TESTSRC")
    assert any("dna_len" in e for e in rep.errors)


def test_mut_positions_must_agree_with_count(good_frame):
    df = good_frame.copy()
    df.loc[df["n_mut_from_wt"] == 1, "mut_positions"] = "11,12"
    rep = schema.validate(df, source="TESTSRC")
    assert any("mut_positions" in e for e in rep.errors)


def test_mut_position_beyond_domain_is_rejected(good_frame):
    """The off-by-one / isoform-offset bug the brief calls the most likely silent failure."""
    df = good_frame.copy()
    df.loc[df["n_mut_from_wt"] == 1, "mut_positions"] = "9999"
    rep = schema.validate(df, source="TESTSRC")
    assert any("outside" in e for e in rep.errors)


def test_tampered_pair_id_is_caught(good_frame):
    df = good_frame.copy()
    df.loc[0, "pair_id"] = "deadbeef"
    rep = schema.validate(df, source="TESTSRC")
    assert any("pair_id" in e for e in rep.errors)


def test_duplicate_pairs_are_caught(good_frame):
    df = pd.concat([good_frame, good_frame.iloc[[0]]], ignore_index=True)
    rep = schema.validate(df, source="TESTSRC")
    assert any("duplicate" in e for e in rep.errors)


def test_ragged_cluster_is_allowed_but_warned(good_frame):
    """Distance is alignment-based now, so differing lengths are legal — but they are also
    how padding clipped by a short construct shows up, so they are surfaced."""
    df = good_frame.copy()
    df.loc[2, "dbd_seq"] = df.loc[2, "dbd_seq"][:-5]
    df = schema.add_pair_ids(df)
    rep = schema.validate(df, source="TESTSRC")
    assert rep.ok, str(rep)
    assert any("differing length" in w for w in rep.warnings)


def test_mut_positions_are_bounded_by_the_reference_not_the_variant(good_frame):
    """A variant with a deletion is shorter than its reference, yet names reference positions."""
    df = good_frame.copy()
    ref_len = len(df.loc[0, "dbd_seq"])
    df.loc[2, "dbd_seq"] = df.loc[2, "dbd_seq"][:-5]  # variant 5 residues shorter
    df.loc[2, "mut_positions"] = str(ref_len)  # legal in the reference frame
    df.loc[2, "n_mut_from_wt"] = 1
    df.loc[3, "dbd_seq"] = df.loc[3, "dbd_seq"][:-5]
    df = schema.add_pair_ids(df)
    rep = schema.validate(df, source="TESTSRC")
    assert rep.ok, str(rep)


def test_variants_without_a_reference_row_are_rejected(good_frame):
    """mut_positions is meaningless without the reference that defines the frame."""
    df = good_frame[good_frame["n_mut_from_wt"] > 0].copy()
    rep = schema.validate(df, source="TESTSRC")
    assert any("no reference row" in e for e in rep.errors)


def test_lowercase_dbd_is_rejected(good_frame):
    df = good_frame.copy()
    df.loc[0, "dbd_seq"] = df.loc[0, "dbd_seq"].lower()
    df = schema.add_pair_ids(df)
    rep = schema.validate(df, source="TESTSRC")
    assert not rep.ok


def test_inverted_thresholds_are_caught(good_frame):
    df = good_frame.copy()
    df["threshold_pos"] = 0.1
    rep = schema.validate(df, source="TESTSRC")
    assert any("threshold" in e for e in rep.errors)


def test_bad_label_value_is_caught(good_frame):
    df = good_frame.copy()
    df.loc[0, "label"] = 2
    rep = schema.validate(df, source="TESTSRC")
    assert any("label" in e for e in rep.errors)


def test_gray_band_label_is_allowed(good_frame):
    df = good_frame.copy()
    df.loc[0, "label"] = schema.LABEL_GRAY
    df.loc[0, "raw_score"] = 0.40
    rep = schema.validate(df, source="TESTSRC")
    assert rep.ok, str(rep)
    assert rep.stats["n_gray"] == 1


def test_raise_if_failed(good_frame):
    df = good_frame.copy()
    df.loc[0, "label"] = 2
    with pytest.raises(schema.SchemaValidationError):
        schema.validate(df, source="TESTSRC").raise_if_failed()
