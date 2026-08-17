"""The project-internal canonical domain sequence.

`dbd_seq` used to be the Pfam envelope plus padding **clipped wherever the assayed construct
happened to end**, so one domain could be stored as two different strings depending on how much
flanking sequence a lab chose to clone. VENTX is 76 aa from BAR15A and 57 aa from CIS-BP. That
difference is cloning, not biology, and it reaches a model as though it were real.

This module produces a form whose defining property is **stability**: the same domain in gives
the same string out, however truncated the input was. Length still varies with the biology — an
insertion or deletion genuinely is a different sequence and must stay one.

The window is the **Pfam envelope +/- `padding_aa`**. The padding is part of the definition, not
an addition to it: BAR15A's VSX1 G160D mutates 6 residues before the Pfam start, and a
representation trimmed to the bare envelope collapses that variant onto its own wild type while
it carries a different label.

**A reference only ever extends; it never replaces.** Every position the construct covers is
taken from the construct, carrying whatever mutations were engineered into it; the reference
supplies only the flank the construct is missing. A short sequence caused by a real deletion, or
by a true protein terminus, therefore stays short — that is what was assayed.

Three quarters of admitted constructs already cover the window and need no reference at all.
`docs/DECISIONS.md` §11 records the measurement behind every threshold here.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any

from Bio import Align
from Bio.Align import substitution_matrices

from snp2prot import thresholds

#: Why a construct could not be canonicalised. Plain strings so they land in reports unchanged.
NO_REFERENCE = "no_reference"
UNPLACEABLE = "unplaceable_in_reference"

#: How the window was obtained.
FROM_CONSTRUCT = "construct"  # the construct already carried the whole window
FROM_REFERENCE = "reference"  # flank was borrowed
SHORT = "construct_short"  # neither was possible; what exists is emitted

#: Accept a construct->reference placement at up to this many edits OUTSIDE the Pfam envelope.
#: Inside it, mutations are the subject matter of the dataset and are never counted. Measured
#: over the 462 constructs that already have a reference: 0 edits at the median, 97.6% within
#: five, and the curve is flat past that -- ten buys 0.2% more. What five rejects is not
#: marginal (Hoxc11 at 111 edits over 57% coverage), which is the point.
MAX_EDITS_OUTSIDE = 5

#: The construct must align into its reference essentially in full. One that only partly aligns
#: is a different protein or a different isoform, not a truncation to repair.
#:
#: 0.95 rather than a rounder number because the observed distribution has a clean gap there.
#: Rejections sit at 0.987 / 0.986 / 0.986 and then at 0.842 and below, with nothing between.
#: The three above the gap are the right protein carrying a couple of vector-derived residues
#: at the termini -- `GR09:Gcn4` against `P03069`, zero edits outside the envelope -- which is
#: a known property of clone constructs. Everything at 0.842 and below is a different protein.
MIN_COVERAGE = 0.95


@dataclass(frozen=True)
class Placement:
    """Where a construct sits inside a reference, and how well it fits.

    Positions are read off the alignment rather than assumed from a constant offset, because a
    construct can carry an indel relative to its reference -- `Mlx` and `Rfx3` differ from
    theirs by 25-54 gaps at full coverage, which is an alternative isoform. A fixed offset
    would silently shift the window on the far side of any such gap.
    """

    #: 0-based reference indices of the construct's first and last aligned residue.
    ref_start: int
    ref_end: int
    #: 0-based reference indices of the Pfam envelope's first and last residue.
    env_ref_start: int
    env_ref_end: int
    coverage: float
    mismatches_outside: int
    gaps_outside: int
    mismatches_inside: int
    #: reference index -> construct index, for every aligned pair.
    ref_to_construct: dict[int, int] = field(default_factory=dict, repr=False)

    @property
    def edits_outside(self) -> int:
        return self.mismatches_outside + self.gaps_outside

    def acceptable(self, max_edits: int = MAX_EDITS_OUTSIDE, min_cov: float = MIN_COVERAGE) -> bool:
        return self.coverage >= min_cov and self.edits_outside <= max_edits


@dataclass(frozen=True)
class CanonicalCall:
    """The canonical sequence for one construct, or why there is none."""

    sequence: str | None
    #: "construct" when the construct already covered the window, "reference" when flank had to
    #: be borrowed. Recorded because it is the difference between a string we derived and one
    #: we partly imported.
    origin: str = ""
    #: 1-based inclusive envelope bounds within `sequence`.
    envelope: tuple[int, int] = (0, 0)
    #: False when the window could not be filled to the full padding on both sides. True for
    #: the great majority; False for a real protein terminus, and for an engineered construct
    #: with no natural protein to extend from.
    full_window: bool = True
    rejection: str | None = None
    placement: Placement | None = None

    @property
    def ok(self) -> bool:
        return self.rejection is None


def padding(config: dict[str, Any] | None = None) -> int:
    cfg = config or thresholds.load()["domain"]
    return int(cfg["padding_aa"])


def covers(construct: str, env_start: int, env_end: int, pad: int) -> bool:
    """Does the construct itself carry the whole window? 1-based inclusive envelope bounds."""
    return env_start - pad >= 1 and env_end + pad <= len(construct)


def _aligner() -> Align.PairwiseAligner:
    return Align.PairwiseAligner(
        mode="local",
        substitution_matrix=substitution_matrices.load("BLOSUM62"),
        open_gap_score=-11,
        extend_gap_score=-1,
    )


def place(construct: str, reference: str, env_start: int, env_end: int) -> Placement | None:
    """Locate a construct inside its reference and count how badly it fits.

    Edits are counted **outside the Pfam envelope only**. Inside it the construct is expected to
    differ -- that is where the engineered mutations live, and those residues are kept verbatim
    -- while outside is the region whose residues are about to be borrowed, so that is what has
    to match. A flat ceiling over the whole construct would reject ROG18A's chimeras, which
    differ from their parent by 6-8 substitutions yet are perfectly well placed.
    """
    if not construct or not reference:
        return None
    try:
        aln = _aligner().align(reference, construct)[0]
    except (ValueError, IndexError):
        return None

    ref_row, con_row = str(aln[0]), str(aln[1])
    ref_pos = int(aln.coordinates[0][0])
    con_pos = int(aln.coordinates[1][0])

    mm_in = mm_out = gap_out = aligned = 0
    ref_start = ref_end = env_ref_start = env_ref_end = None
    mapping: dict[int, int] = {}

    for r, c in zip(ref_row, con_row, strict=False):
        inside = env_start <= con_pos + 1 <= env_end
        if r != "-" and c != "-":
            mapping[ref_pos] = con_pos
            if ref_start is None:
                ref_start = ref_pos
            ref_end = ref_pos
            if inside:
                if env_ref_start is None:
                    env_ref_start = ref_pos
                env_ref_end = ref_pos
                if r != c:
                    mm_in += 1
            elif r != c:
                mm_out += 1
        elif not inside:
            gap_out += 1
        if c != "-":
            aligned += 1
            con_pos += 1
        if r != "-":
            ref_pos += 1

    if ref_start is None or env_ref_start is None:
        return None
    return Placement(
        ref_start=ref_start,
        ref_end=ref_end,
        env_ref_start=env_ref_start,
        env_ref_end=env_ref_end,
        coverage=aligned / len(construct),
        mismatches_outside=mm_out,
        gaps_outside=gap_out,
        mismatches_inside=mm_in,
        ref_to_construct=mapping,
    )


def conflicting_constructs(constructs: dict[str, str]) -> bool:
    """Do these constructs, which share one canonical sequence, differ OUTSIDE it?

    Two constructs collapsing to the same canonical domain are normally either byte-identical
    or differ only in how much flank each lab chose to clone -- both harmless. But a genuine
    substitution outside the window is a real difference the model will never see, attached to
    measurements that may differ because of it. `Cell09:HLH-25` and `HLH-27` are two distinct
    *C. elegans* genes with identical 76 aa domains and four substitutions just outside.

    Owner's decision, 2026-08-17: bin them. This is admission condition 3 applied symmetrically
    -- the same reasoning that rejects a variant whose mutation escapes the stored region.

    The test is a **substitution in the overlap**, not a difference in length. Constructs are
    aligned with free end gaps, so one carrying extra residues at the N terminus and another
    at the C terminus is not a conflict: they agree everywhere they overlap. Only a mismatched
    residue counts. Returns True when the group must be discarded.
    """
    seqs = sorted(set(constructs.values()))
    if len(seqs) < 2:
        return False
    al = Align.PairwiseAligner(
        mode="global",
        substitution_matrix=substitution_matrices.load("BLOSUM62"),
        open_gap_score=-11,
        extend_gap_score=-1,
        end_gap_score=0,
    )
    for i, a in enumerate(seqs):
        for b in seqs[i + 1 :]:
            try:
                aln = al.align(a, b)[0]
            except (ValueError, IndexError):
                return True
            top, bottom = str(aln[0]), str(aln[1])
            if any(x != y and x != "-" and y != "-" for x, y in zip(top, bottom, strict=False)):
                return True
    return False


def canonicalise(
    construct: str,
    env_start: int,
    env_end: int,
    reference: str | None = None,
    pad: int | None = None,
    max_edits: int = MAX_EDITS_OUTSIDE,
    min_coverage: float = MIN_COVERAGE,
) -> CanonicalCall:
    """The canonical sequence for one construct. `env_start`/`env_end` are 1-based inclusive.

    Where the construct already covers the window this is a pure slice and no reference is
    consulted; the result is provably the string the full-length protein would have given.
    """
    pad = padding() if pad is None else pad
    if covers(construct, env_start, env_end, pad):
        return CanonicalCall(
            sequence=construct[env_start - pad - 1 : env_end + pad],
            origin="construct",
            envelope=(pad + 1, pad + env_end - env_start + 1),
        )

    if reference:
        p = place(construct, reference, env_start, env_end)
        if p is not None and p.acceptable(max_edits, min_coverage):
            win_lo = max(0, p.env_ref_start - pad)
            win_hi = min(len(reference) - 1, p.env_ref_end + pad)
            seq = "".join(
                construct[p.ref_to_construct[i]] if i in p.ref_to_construct else reference[i]
                for i in range(win_lo, win_hi + 1)
            )
            return CanonicalCall(
                sequence=seq,
                origin=FROM_REFERENCE,
                envelope=(p.env_ref_start - win_lo + 1, p.env_ref_end - win_lo + 1),
                full_window=(
                    p.env_ref_start - pad >= 0 and p.env_ref_end + pad <= len(reference) - 1
                ),
                placement=p,
            )

    # Neither route filled the window. Emit what the construct actually has, flagged.
    #
    # An engineered chimera cloned as a bare domain has no natural protein to extend from --
    # ROG18A's FoxJ3/FoxN3 hybrids are synthetic, and they were the corpus's only
    # non-homeodomain protein-axis depth. Their termini are not truncations: nothing was cut
    # off, the assayed protein simply is that sequence. Discarding them would throw away real
    # measurements to satisfy a rule about repairing damage that was never done.
    #
    # The cost is that a short window can differ in length from a fully-canonical copy of the
    # same domain, which is the variation this module exists to remove. That is resolved at
    # merge, not here: a short entry subsumed by a full one is dropped as redundant, so no two
    # stored rows describe one domain at two lengths.
    lo = max(1, env_start - pad)
    hi = min(len(construct), env_end + pad)
    return CanonicalCall(
        sequence=construct[lo - 1 : hi],
        origin=SHORT,
        envelope=(env_start - lo + 1, env_end - lo + 1),
        full_window=False,
    )
