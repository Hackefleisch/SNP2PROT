"""The canonical domain sequence.

The property under test is **stability**: the same domain gives the same string however
truncated the input was. Every case here is drawn from something the corpus actually did.
"""

from __future__ import annotations

from snp2prot import canonical

PAD = 10
#: A stand-in protein. The "domain" is the middle stretch; flanks are distinguishable so a
#: window that borrows from the wrong place is visible in the assertion, not just wrong.
PROTEIN = ("M" * 30) + ("ACDEFGHIKL" * 5) + ("W" * 30)
ENV_LO, ENV_HI = 31, 80  # 1-based inclusive, the 50-residue "domain" inside PROTEIN


def test_a_construct_covering_the_window_needs_no_reference():
    call = canonical.canonicalise(PROTEIN, ENV_LO, ENV_HI, pad=PAD)
    assert call.ok and call.origin == canonical.FROM_CONSTRUCT and call.full_window
    assert call.sequence == PROTEIN[ENV_LO - PAD - 1 : ENV_HI + PAD]
    assert len(call.sequence) == (ENV_HI - ENV_LO + 1) + 2 * PAD
    assert call.envelope == (PAD + 1, PAD + 50)


def test_the_same_domain_clipped_differently_gives_the_same_string():
    """This is the whole point: VENTX is 76 aa from BAR15A and 57 aa from CIS-BP."""
    generous = PROTEIN[10:100]  # plenty of flank either side
    bare = PROTEIN[ENV_LO - 1 : ENV_HI]  # zero flank, like a CIS-BP oligo-pool construct
    a = canonical.canonicalise(generous, ENV_LO - 10, ENV_HI - 10, pad=PAD)
    b = canonical.canonicalise(bare, 1, 50, reference=PROTEIN, pad=PAD)
    assert a.ok and b.ok
    assert a.origin == canonical.FROM_CONSTRUCT and b.origin == canonical.FROM_REFERENCE
    assert a.sequence == b.sequence


def test_the_construct_supplies_its_own_residues_mutations_included():
    """A reference extends; it must never overwrite an engineered mutation."""
    mutant = PROTEIN[ENV_LO - 1 : ENV_HI]
    mutant = mutant[:20] + "P" + mutant[21:]  # a point mutant inside the domain
    call = canonical.canonicalise(mutant, 1, 50, reference=PROTEIN, pad=PAD)
    assert call.ok and call.origin == canonical.FROM_REFERENCE
    assert call.sequence[PAD + 20] == "P", "the mutation was overwritten by the reference"
    assert (
        call.sequence
        != canonical.canonicalise(
            PROTEIN[ENV_LO - 1 : ENV_HI], 1, 50, reference=PROTEIN, pad=PAD
        ).sequence
    )


def test_a_construct_with_no_reference_keeps_what_it_has_flagged():
    """An engineered chimera has no natural protein to extend from, and its termini are real.

    ROG18A's FoxJ3/FoxN3 hybrids are synthetic and were the corpus's only non-homeodomain
    protein-axis depth. Nothing was cut off them, so there is nothing to repair — the window
    is short on purpose and says so.
    """
    bare = PROTEIN[ENV_LO - 1 : ENV_HI]
    call = canonical.canonicalise(bare, 1, 50, pad=PAD)
    assert call.ok and call.origin == canonical.SHORT
    assert call.full_window is False
    assert call.sequence == bare


def test_a_reference_the_construct_does_not_belong_to_is_not_used():
    """Mis-attribution is worse than absence — the PP15 rule.

    The construct is still emitted, because it is a real assayed protein, but nothing is
    borrowed from a reference it demonstrably does not match.
    """
    unrelated = "QQQQQ" + ("RSTVWY" * 9)
    call = canonical.canonicalise(unrelated, 6, 50, reference=PROTEIN, pad=PAD)
    assert call.origin == canonical.SHORT
    assert PROTEIN[:5] not in call.sequence, "residues were taken from the wrong protein"


def test_mutations_inside_the_envelope_do_not_spend_the_edit_budget():
    """ROG18A's chimeras differ from their parent by 6-8 substitutions and are still placed."""
    chimera = list(PROTEIN[ENV_LO - 1 : ENV_HI])
    for i in (3, 9, 15, 21, 27, 33, 39):  # seven substitutions, all inside the domain
        chimera[i] = "P" if chimera[i] != "P" else "G"
    call = canonical.canonicalise("".join(chimera), 1, 50, reference=PROTEIN, pad=PAD)
    assert call.ok, "in-domain mutations must not count against placement"
    assert call.placement.mismatches_inside >= 6
    assert call.placement.edits_outside == 0


def test_a_real_protein_terminus_yields_a_short_window_and_says_so():
    """Not truncation: the protein genuinely stops. Keep it, flagged."""
    short_ref = PROTEIN[ENV_LO - 3 :]  # only 2 residues before the domain exist at all
    construct = short_ref[:60]
    call = canonical.canonicalise(construct, 3, 52, reference=short_ref, pad=PAD)
    assert call.ok and not call.full_window
    assert call.sequence.startswith(short_ref[:2])


def test_padding_is_what_keeps_an_out_of_envelope_variant_distinguishable():
    """BAR15A's VSX1 G160D mutates 6 residues before the Pfam start."""
    wt = PROTEIN
    var = PROTEIN[: ENV_LO - 7] + "P" + PROTEIN[ENV_LO - 6 :]
    a = canonical.canonicalise(wt, ENV_LO, ENV_HI, pad=PAD)
    b = canonical.canonicalise(var, ENV_LO, ENV_HI, pad=PAD)
    assert a.sequence != b.sequence, "the padding must preserve a variant just outside the domain"
    trimmed_a = a.sequence[PAD:-PAD]
    trimmed_b = b.sequence[PAD:-PAD]
    assert trimmed_a == trimmed_b, "and without it the two would be identical"
