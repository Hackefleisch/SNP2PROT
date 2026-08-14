"""The domain admission policy.

These tests pin the three conditions from docs/DOMAIN_POLICY.md. They use synthetic hits
rather than real HMM scans so they run without the Pfam files or the raw archives.
"""

from __future__ import annotations

from pathlib import Path

import pyhmmer
import pytest

from snp2prot import domains

CFG = {
    "padding_aa": 10,
    "allow_mixed_families": False,
    "allow_repeat_arrays": False,
    "allow_no_domain": False,
    "hmm_dir": "data/external/pfam",
}
SEQ = "".join("ACDEFGHIKL" for _ in range(10))  # 100 aa


def hit(fam, a, b):
    return domains.DomainHit(fam, a, b)


def test_padding_is_applied_on_both_sides():
    call = domains.call_domain(SEQ, [hit("Homeodomain", 30, 60)], CFG)
    assert call.ok
    assert (call.start, call.end) == (20, 70)
    assert call.sequence == SEQ[19:70]
    assert len(call.sequence) == 51  # 31 + 10 + 10


def test_padding_clips_at_sequence_bounds():
    """A domain near either terminus must not produce an out-of-range slice."""
    call = domains.call_domain(SEQ, [hit("Homeodomain", 3, 20)], CFG)
    assert (call.start, call.end) == (1, 30)
    call = domains.call_domain(SEQ, [hit("Homeodomain", 80, 98)], CFG)
    assert (call.start, call.end) == (70, 100)
    assert call.sequence == SEQ[69:100]


def test_condition_1_and_2_mixed_families_rejected():
    call = domains.call_domain(SEQ, [hit("Pou", 10, 40), hit("Homeodomain", 55, 90)], CFG)
    assert not call.ok
    assert call.rejection == domains.MIXED_FAMILIES
    assert call.sequence is None
    assert call.family == "Homeodomain, Pou"


def test_condition_2_repeat_array_rejected():
    """A C2H2 array is several small folds on flexible linkers, not one continuous unit."""
    hits = [hit("zf-C2H2", 10, 32), hit("zf-C2H2", 40, 62), hit("zf-C2H2", 68, 90)]
    call = domains.call_domain(SEQ, hits, CFG)
    assert not call.ok
    assert call.rejection == domains.REPEAT_ARRAY


def test_repeat_array_can_be_readmitted_by_config():
    hits = [hit("zf-C2H2", 20, 32), hit("zf-C2H2", 40, 62)]
    call = domains.call_domain(SEQ, hits, {**CFG, "allow_repeat_arrays": True})
    assert call.ok
    assert (call.start, call.end) == (10, 72)  # spans both fingers, padded


def test_no_domain_rejected():
    call = domains.call_domain(SEQ, [], CFG)
    assert not call.ok
    assert call.rejection == domains.NO_DOMAIN


def test_condition_3_flags_mutations_outside_the_stored_region():
    call = domains.call_domain(SEQ, [hit("Homeodomain", 30, 60)], CFG)  # stored 20-70
    assert domains.audit_variant_positions(call, [35]) == []
    assert domains.audit_variant_positions(call, [20, 70]) == []  # inclusive bounds
    assert domains.audit_variant_positions(call, [19, 71]) == [19, 71]


def test_condition_3_is_what_padding_buys():
    """VSX1_G160D: mutation 5 residues N-terminal of the Pfam start.

    Unpadded it escapes the stored region and the variant collapses onto its wild type;
    padded it is retained.
    """
    unpadded = domains.call_domain(SEQ, [hit("Homeodomain", 15, 71)], {**CFG, "padding_aa": 0})
    padded = domains.call_domain(SEQ, [hit("Homeodomain", 15, 71)], CFG)
    assert domains.audit_variant_positions(unpadded, [10]) == [10]
    assert domains.audit_variant_positions(padded, [10]) == []


def test_rejected_call_reports_every_position_as_escaped():
    call = domains.call_domain(SEQ, [], CFG)
    assert domains.audit_variant_positions(call, [5, 9]) == [5, 9]


def test_rebase_positions_are_one_based_within_dbd_seq():
    call = domains.call_domain(SEQ, [hit("Homeodomain", 30, 60)], CFG)  # stored 20-70
    assert domains.rebase_positions(call, [20]) == [1]
    assert domains.rebase_positions(call, [35]) == [16]
    # the rebased position must index the same residue in the stored slice
    assert call.sequence[domains.rebase_positions(call, [35])[0] - 1] == SEQ[34]


@pytest.mark.parametrize("pad", [0, 5, 10, 25])
def test_stored_slice_always_matches_reported_offsets(pad):
    call = domains.call_domain(SEQ, [hit("Homeodomain", 30, 60)], {**CFG, "padding_aa": pad})
    assert call.sequence == SEQ[call.start - 1 : call.end]


# --- co-located hits --------------------------------------------------------------------
# Pfam models one domain with several families in places: bZIP_1/bZIP_2 both match the same
# basic-leucine-zipper, zf-H2C2_2 overlaps zf-C2H2, Homeodomain overlaps Homeobox_KN. Before
# this was handled, mixed_families rejected most of the bZIP family.


def test_colocated_families_are_one_domain_not_two():
    hits = [hit("bZIP_1", 18, 74), hit("bZIP_2", 18, 70)]
    call = domains.call_domain(SEQ, hits, CFG)
    assert call.ok, "two Pfam models of one bZIP must not read as two domains"
    # Both models alias to one family label, so the reported family cannot depend on
    # which of them happened to score higher.
    assert call.family == "bZIP"


def test_the_higher_scoring_model_is_the_one_kept():
    """Which hit survives collapsing decides the stored boundaries, so it must be the best."""
    a = domains.DomainHit("Homeodomain", 18, 74, score=40.0)
    b = domains.DomainHit("Forkhead", 20, 70, score=95.0)
    kept = domains.collapse_colocated([a, b])
    assert len(kept) == 1
    assert kept[0] == b, "the higher-scoring model must define the domain span"


def test_separated_families_are_still_rejected():
    """A real two-domain construct (PAX + homeodomain) must still fail condition 1."""
    call = domains.call_domain(SEQ, [hit("PAX", 4, 40), hit("Homeodomain", 55, 95)], CFG)
    assert not call.ok
    assert call.rejection == domains.MIXED_FAMILIES


def test_partial_overlap_below_threshold_still_counts_as_two():
    call = domains.call_domain(SEQ, [hit("PAX", 10, 50), hit("Homeodomain", 46, 90)], CFG)
    assert not call.ok, "a small overlap is two domains touching, not one described twice"


def test_collapse_keeps_a_zinc_finger_array_an_array():
    """Fingers sit side by side, so collapsing must not turn an array into one domain."""
    hits = [hit("zf-C2H2", 10, 32), hit("zf-C2H2", 40, 62), hit("zf-C2H2", 68, 90)]
    assert len(domains.collapse_colocated(hits)) == 3
    assert domains.call_domain(SEQ, hits, CFG).rejection == domains.REPEAT_ARRAY


def test_bzip_models_report_one_family_name():
    """bZIP_1 and bZIP_2 are models of one family; dbd_family must not depend on which won."""
    a = domains.call_domain(SEQ, [hit("bZIP_1", 18, 74)], CFG)
    b = domains.call_domain(SEQ, [hit("bZIP_2", 18, 70)], CFG)
    assert a.family == b.family == "bZIP"


def test_lookalike_families_are_not_merged():
    """TF_AP-2 (mammalian) and AP2 (plant) are different domains despite similar names."""
    assert domains.canonical_family("TF_AP-2") == "TF_AP-2"
    assert domains.canonical_family("AP2") == "AP2"


# --------------------------------------------------------------------------------------
# Library loading (TODO.md T6). `scan` now loads an HMM library once per process and holds
# it, and the project's Pfam-A is pressed into HMMER's binary format. Both are speed
# changes, so what matters is that neither alters a single domain call.
# --------------------------------------------------------------------------------------

HOMEODOMAIN_HMM = Path("data/external/pfam/PF00046.hmm")
ARX = "AGSDSEEGLLKRKQRRYRTTFTSYQLEEQERAFQKTHYPDVFTREELAMRLDLTEARVQVWFQNRRAKWRKREKAGAQTHPPGLPF"


@pytest.mark.skipif(not HOMEODOMAIN_HMM.exists(), reason="Pfam HMM not downloaded")
def test_pressed_and_unpressed_libraries_agree(tmp_path):
    """Pressing is a storage format, not a scoring change.

    `scripts/press_pfam.py` cuts the library load from 19.5 s to 0.5 s, which is only safe
    because the calls come out identical. Checked here on one small family rather than the
    2.2 GB Pfam-A, but it is the same code path.
    """
    plain = tmp_path / "plain" / HOMEODOMAIN_HMM.name
    pressed = tmp_path / "pressed" / HOMEODOMAIN_HMM.name
    for dest in (plain, pressed):
        dest.parent.mkdir()
        dest.write_bytes(HOMEODOMAIN_HMM.read_bytes())

    with pyhmmer.plan7.HMMFile(pressed) as f:
        pyhmmer.hmmer.hmmpress(f, pressed)
    with pyhmmer.plan7.HMMFile(pressed) as f:
        assert f.is_pressed()
    with pyhmmer.plan7.HMMFile(plain) as f:
        assert not f.is_pressed()

    before = domains.scan({"ARX": ARX}, hmm_path=plain)
    after = domains.scan({"ARX": ARX}, hmm_path=pressed)
    assert before == after
    assert [h.family for h in after["ARX"]] == ["Homeodomain"]


@pytest.mark.skipif(not HOMEODOMAIN_HMM.exists(), reason="Pfam HMM not downloaded")
def test_library_is_loaded_once_per_path():
    """The whole saving is that a build reuses one load across all its sources."""
    a = domains._library(str(HOMEODOMAIN_HMM.resolve()))
    b = domains._library(str(HOMEODOMAIN_HMM.resolve()))
    assert a is b
