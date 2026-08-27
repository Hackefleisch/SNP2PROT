"""Distance between a reference domain and one of its variants.

Cluster members used to be required to be the same length, because `n_mut_from_wt` was a
Hamming distance. That excluded any variant carrying an insertion or deletion, and it also
excluded variants that differ only because our own padding was clipped. So distance is now
taken from a **pairwise alignment**, which admits indels.

Two properties of that alignment are load-bearing:

**Terminal gaps are free.** `dbd_seq` is the Pfam envelope plus `padding_aa` residues each
side, clipped where the assayed construct ends. Two versions of the very same domain can
therefore differ in length purely by how much padding fitted: in ROG18A every bare domain is
83-85 aa, but the padded ones run 95-105 aa because some constructs allowed the full 10
residues of N-terminal padding and others only 2. Charging for those terminal gaps would
report 8 phantom indels between a protein and its own chimera. Internal indels still count —
those are real.

**Free terminal gaps need a floor, or they stop measuring distance at all.** With end gaps
weighted zero, an arbitrarily long prefix of one sequence and suffix of the other can be
discarded for nothing, and the edits are then counted only over whatever window survives. For
two unrelated domains the optimiser exploits exactly that: a 72 aa and a 60 aa Myb domain
align on six residues, score +14 against −22 for the honest alignment, and report **3 edits
where a charged alignment reports 59**. `overlap` records how much of the shorter sequence
actually aligned so a caller can refuse that answer; `snp2prot.clusters` does, at
`cluster.min_overlap`. Legitimate clipping is unaffected — when the difference really is
padding, the shorter sequence is fully contained and `overlap` is 1.0.

**Positions are reported in the reference frame.** Every variant in a cluster is then
described in one coordinate system, which is what makes `mut_positions` comparable across a
cluster and mappable onto a single predicted structure. An inserted residue is attributed to
the reference position it precedes, and one entry is emitted per edited residue, so
`len(mut_positions) == n_mut_from_wt` always holds.
"""

from __future__ import annotations

import functools
from dataclasses import dataclass

from Bio import Align
from Bio.Align import substitution_matrices


@dataclass(frozen=True)
class EditProfile:
    """How a variant differs from its reference."""

    n_edits: int
    #: One 1-based reference-frame position per edited residue, ascending. Duplicates occur
    #: when several residues are inserted at the same point.
    positions: tuple[int, ...]
    #: Substitutions / insertions / deletions, for reporting.
    n_substitutions: int = 0
    n_insertions: int = 0
    n_deletions: int = 0
    #: Residues aligned to a residue, as a fraction of the shorter sequence. 1.0 means one
    #: sequence is wholly contained in the other, which is what differently-clipped padding
    #: looks like. A low value means the aligner discarded most of both sequences into free
    #: terminal gaps, and `n_edits` then describes a window rather than the domains.
    overlap: float = 1.0
    #: The same quantity before it was divided: residues aligned to a residue. `snp2prot.distances`
    #: stores this rather than the ratio, so a caller can change its mind about the guard without
    #: a rebuild — and so it does not have to multiply `overlap` back out and round.
    n_aligned: int = 0

    @property
    def positions_str(self) -> str:
        return ",".join(str(p) for p in self.positions)


@functools.lru_cache(maxsize=1)
def _aligner(free_end_gaps: bool = True) -> Align.PairwiseAligner:
    a = Align.PairwiseAligner(mode="global")
    a.substitution_matrix = substitution_matrices.load("BLOSUM62")
    a.open_gap_score = -11.0
    a.extend_gap_score = -1.0
    if free_end_gaps:
        # Terminal truncation is an artefact of padding clipping, not a difference between
        # the proteins; see the module docstring.
        a.open_end_gap_score = 0.0
        a.extend_end_gap_score = 0.0
    return a


def edit_profile(reference: str, variant: str, free_end_gaps: bool = True) -> EditProfile:
    """Align `variant` to `reference` and describe every edit in reference coordinates."""
    if reference == variant:
        return EditProfile(0, (), n_aligned=len(reference))

    alignment = _aligner(free_end_gaps).align(reference, variant)[0]
    ref_row, var_row = str(alignment[0]), str(alignment[1])

    # Terminal gap runs are skipped entirely when end gaps are free, so that a domain padded
    # to 10 residues and the same domain clipped to 2 are not called different.
    lo, hi = 0, len(ref_row)
    if free_end_gaps:
        while lo < hi and (ref_row[lo] == "-" or var_row[lo] == "-"):
            lo += 1
        while hi > lo and (ref_row[hi - 1] == "-" or var_row[hi - 1] == "-"):
            hi -= 1

    aligned = sum(
        1 for i in range(lo, hi) if ref_row[i] != "-" and var_row[i] != "-"
    )  # residue-to-residue, the part of the comparison that actually happened

    positions: list[int] = []
    n_sub = n_ins = n_del = 0
    ref_pos = sum(1 for c in ref_row[:lo] if c != "-")  # residues of the reference consumed

    for i in range(lo, hi):
        r, v = ref_row[i], var_row[i]
        if r != "-":
            ref_pos += 1
        if r == v:
            continue
        if r == "-":
            # Inserted residue: attribute it to the reference position it precedes.
            positions.append(min(ref_pos + 1, len(reference)))
            n_ins += 1
        elif v == "-":
            positions.append(ref_pos)
            n_del += 1
        else:
            positions.append(ref_pos)
            n_sub += 1

    return EditProfile(
        n_edits=len(positions),
        positions=tuple(sorted(positions)),
        n_substitutions=n_sub,
        n_insertions=n_ins,
        n_deletions=n_del,
        overlap=aligned / min(len(reference), len(variant)),
        n_aligned=aligned,
    )
