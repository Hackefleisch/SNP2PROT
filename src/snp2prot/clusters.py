"""Grouping canonical domains into clusters, for building splits.

A cluster answers one question: which rows must never be split across train and test, because
their sequences are close enough that a model could copy an answer rather than derive it. It
is **not** a claim that its members are interchangeable — every member is its own training row
of (sequence, 8-mer) -> label, and a cluster never collapses them.

**The algorithm is CD-HIT's greedy incremental clustering** (Li & Godzik 2006; the same as
MMseqs2 `--cluster-mode 2`), not something invented here. Sort by decreasing length, the
longest sequence becomes a cluster representative, and every remaining sequence is compared
**only to representatives** -- joining the first it is close enough to, or founding a new
cluster. Three properties earn it the job:

* **the longest member represents**, so the least-clipped form of a domain is the reference;
* **comparison is never transitive**, which is what stops chaining. Single-linkage on this
  corpus produced a 35-domain blob at a threshold of 5 edits against 8 at 1 edit -- A near B
  near C, with A and C unrelated, is not a cluster;
* **every member is within `max_edits` of the representative**, which is exactly the invariant
  `mut_positions` needs: one coordinate frame per cluster (`docs/DECISIONS.md` §2).

Its known weakness is order dependence -- a sequence joins the first representative it matches
rather than its best. Ties are broken deterministically by sequence so a rebuild reproduces the
same clusters.

Distance comes from `snp2prot.align`, which scores BLOSUM62 with affine gaps and **free
terminal gaps**. Implemented here rather than by installing MMseqs2: the aligner already has
the right semantics, and a full pairwise pass over the corpus takes seconds. The algorithm is
what is being reused; a C++ dependency would add nothing at this scale.
"""

from __future__ import annotations

from dataclasses import dataclass, field

from snp2prot import align


@dataclass
class Cluster:
    """One representative and the domains that fall within `max_edits` of it."""

    representative: str
    members: list[str] = field(default_factory=list)

    @property
    def size(self) -> int:
        return len(self.members)


def cluster(
    sequences: dict[str, str],
    max_edits: int,
    same_family_only: bool = True,
    families: dict[str, str] | None = None,
) -> list[Cluster]:
    """Greedy incremental clustering of `{key: sequence}`.

    Only sequences of the same Pfam family are compared when `same_family_only`, which is both
    correct and a large saving: two domains of different folds are never one cluster, and the
    comparison is skipped rather than computed and discarded.
    """
    fams = families or {}
    # Longest first, then by sequence, so the result does not depend on dict ordering.
    order = sorted(sequences, key=lambda k: (-len(sequences[k]), sequences[k], k))

    reps: list[tuple[str, str, str]] = []  # (key, sequence, family)
    clusters: dict[str, Cluster] = {}
    for key in order:
        seq = sequences[key]
        fam = fams.get(key, "")
        joined = None
        for rkey, rseq, rfam in reps:
            if same_family_only and fam and rfam and fam != rfam:
                continue
            # No length prefilter. Terminal gaps are free, so two sequences of very
            # different length can still be zero edits apart -- which is precisely the
            # differently-clipped case this dataset exists to reconcile. Filtering on length
            # would silently refuse to merge them. Greedy clustering only ever compares
            # against representatives, so the cost of checking honestly is small.
            if align.edit_profile(rseq, seq).n_edits <= max_edits:
                joined = rkey
                break
        if joined is None:
            reps.append((key, seq, fam))
            clusters[key] = Cluster(representative=key, members=[key])
        else:
            clusters[joined].members.append(key)
    return list(clusters.values())


def assign(clusters: list[Cluster]) -> dict[str, str]:
    """member key -> representative key, for stamping a cluster id onto rows."""
    return {m: c.representative for c in clusters for m in c.members}
