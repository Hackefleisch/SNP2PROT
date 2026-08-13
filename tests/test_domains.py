"""The domain admission policy.

These tests pin the three conditions from docs/DOMAIN_POLICY.md. They use synthetic hits
rather than real HMM scans so they run without the Pfam files or the raw archives.
"""

from __future__ import annotations

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
