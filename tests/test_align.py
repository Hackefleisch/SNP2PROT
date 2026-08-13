"""Alignment-based distance between a reference domain and its variants."""

from __future__ import annotations

import pytest

from snp2prot.align import edit_profile

REF = "ACDEFGHIKLMNPQRSTVWY"


def test_identical_sequences_have_no_edits():
    e = edit_profile(REF, REF)
    assert (e.n_edits, e.positions) == (0, ())


def test_substitution_is_located_in_the_reference_frame():
    e = edit_profile(REF, REF[:4] + "W" + REF[5:])
    assert e.n_edits == 1
    assert e.positions == (5,)
    assert (e.n_substitutions, e.n_insertions, e.n_deletions) == (1, 0, 0)


def test_insertion_is_counted_and_attributed():
    e = edit_profile(REF, REF[:5] + "W" + REF[5:])
    assert e.n_edits == 1
    assert e.n_insertions == 1
    assert 1 <= e.positions[0] <= len(REF)


def test_deletion_is_counted():
    e = edit_profile(REF, REF[:5] + REF[6:])
    assert e.n_edits == 1
    assert e.n_deletions == 1


def test_terminal_truncation_is_free():
    """dbd_seq is padded then CLIPPED by the construct, so two versions of one domain differ
    at the termini for reasons that have nothing to do with the protein."""
    assert edit_profile("XXXX" + REF + "ZZZ", REF).n_edits == 0
    assert edit_profile(REF, REF[3:]).n_edits == 0
    assert edit_profile(REF, REF[:-4]).n_edits == 0


def test_terminal_gaps_can_be_charged_if_asked():
    e = edit_profile(REF, REF[3:], free_end_gaps=False)
    assert e.n_edits == 3


def test_internal_indel_still_counts_when_ends_are_free():
    variant = REF[:8] + REF[10:]  # two residues removed from the middle
    e = edit_profile(REF, variant)
    assert e.n_edits == 2
    assert e.n_deletions == 2


def test_one_position_is_emitted_per_edited_residue():
    """The validator requires len(mut_positions) == n_mut_from_wt."""
    e = edit_profile(REF, REF[:5] + "WW" + REF[7:])
    assert len(e.positions) == e.n_edits
    assert e.positions_str.count(",") == e.n_edits - 1


@pytest.mark.parametrize("n", [1, 3, 6])
def test_substitution_count_matches_a_known_edit(n):
    variant = "".join("W" if i < n else c for i, c in enumerate(REF))
    assert edit_profile(REF, variant).n_edits == n
