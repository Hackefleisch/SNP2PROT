"""The threshold file is the owner's control surface; assert its shape, not its values."""

from __future__ import annotations

import pytest

from snp2prot import thresholds


def test_loads():
    cfg = thresholds.load()
    assert cfg["version"] >= 2
    for block in ("pbm", "b1h", "snp_selex", "tier4", "domain"):
        assert block in cfg


def test_domain_block_shape():
    """dbd_seq is the PADDED Pfam envelope; see docs/DOMAIN_POLICY.md."""
    d = thresholds.for_assay("domain")
    assert d["padding_aa"] > 0, "unpadded envelopes drop variants onto their own wild type"
    # The three admission conditions default to strict.
    assert d["allow_mixed_families"] is False
    assert d["allow_repeat_arrays"] is False
    assert d["allow_no_domain"] is False


def test_pbm_band_is_ordered_and_per_experiment():
    pbm = thresholds.for_assay("pbm")
    assert pbm["positive"] > pbm["negative"], "gray band must be non-inverted"
    assert pbm["per_experiment"] is True, "pooling E-scores across array designs is a bug"


def test_tier4_is_never_binarized_at_parse_time():
    assert thresholds.for_assay("tier4")["binarize_at_parse"] is False


def test_unknown_assay_raises():
    with pytest.raises(KeyError):
        thresholds.for_assay("nope")
