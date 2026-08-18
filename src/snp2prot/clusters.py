"""Grouping canonical domains into clusters, for building splits.

A cluster answers one question: which rows must never be split across train and test, because
their sequences are close enough that a model could copy an answer rather than derive it. It
is **not** a claim that its members are interchangeable — every member is its own training row
of (sequence, 8-mer) -> label, and a cluster never collapses them.

**The algorithm is CD-HIT's greedy incremental clustering** (Li & Godzik 2006; the same as
MMseqs2 `--cluster-mode 2`), not something invented here. Sort by decreasing length, the
longest sequence becomes a cluster **seed**, and every remaining sequence is compared
**only to seeds** -- joining the first it is close enough to, or founding a new
cluster. Three properties earn it the job:

* **the longest member seeds**, so the least-clipped form of a domain is what others are
  matched against, which is the direction the padding artefact runs;
* **comparison is never transitive**, which is what stops chaining. Single-linkage on this
  corpus produced a 35-domain blob at a threshold of 5 edits against 8 at 1 edit -- A near B
  near C, with A and C unrelated, is not a cluster;
* **every member is within `max_edits` of the seed**, which bounds how far apart a cluster
  can spread.

**The seed is not the reference.** Two jobs were being done by one field. The seed decides
*membership* and is an artefact of the algorithm — with every member of a variant series the
same padded length, it was settled by the alphabetical order of the amino-acid string, which
made `HOXD13_S316C` the frame for its own wild type and injected position 50 into all seven
siblings' `mut_positions`. The **reference** decides the *coordinate frame* and is the
**medoid**: the member with the smallest total distance to the others. For a wild type plus
k single substitutions that provably picks the wild type, since it sits one edit from each
while any variant sits two from the rest. Decided 2026-08-18 (`docs/DECISIONS.md`, `D3`).

One invariant weakens and is therefore checked rather than assumed: members are within
`max_edits` of the *seed* by construction, but only within `2 * max_edits` of the reference by
the triangle inequality. Measured over the corpus the worst case is 5 of 5, and
`build_clusters.py` reports any cluster that exceeds it.

Its known weakness is order dependence -- a sequence joins the first seed it matches
rather than its best. Ties are broken deterministically by sequence so a rebuild reproduces the
same clusters.

Distance comes from `snp2prot.align`, which scores BLOSUM62 with affine gaps and **free
terminal gaps**. Implemented here rather than by installing MMseqs2: the aligner already has
the right semantics, and a full pairwise pass over the corpus takes seconds. The algorithm is
what is being reused; a C++ dependency would add nothing at this scale.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path

import pandas as pd

from snp2prot import align
from snp2prot.config import CLUSTER_TABLE


@dataclass
class Cluster:
    """One seed and the domains that fall within `max_edits` of it."""

    #: The founding member: longest-first, and what membership was decided against.
    seed: str
    members: list[str] = field(default_factory=list)

    @property
    def size(self) -> int:
        return len(self.members)


def cluster(
    sequences: dict[str, str],
    max_edits: int,
    same_family_only: bool = True,
    families: dict[str, str] | None = None,
    min_overlap: float = 0.0,
) -> list[Cluster]:
    """Greedy incremental clustering of `{key: sequence}`.

    Only sequences of the same Pfam family are compared when `same_family_only`, which is both
    correct and a large saving: two domains of different folds are never one cluster, and the
    comparison is skipped rather than computed and discarded.

    `min_overlap` is the guard on free terminal gaps. Without it, two unrelated domains can be
    called a few edits apart because the aligner parked almost all of both sequences in gaps
    that cost nothing — 31 of 202 memberships in this corpus were formed that way, aligning on
    3-10% of the shorter sequence (`docs/DECISIONS.md`, `T25`). Set it from
    `cluster.min_overlap`; the default of 0 preserves the old behaviour for callers that
    deliberately want it, such as reproducing a historical inventory.
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
            # against seeds, so the cost of checking honestly is small.
            profile = align.edit_profile(rseq, seq)
            if profile.n_edits <= max_edits and profile.overlap >= min_overlap:
                joined = rkey
                break
        if joined is None:
            reps.append((key, seq, fam))
            clusters[key] = Cluster(seed=key, members=[key])
        else:
            clusters[joined].members.append(key)
    return list(clusters.values())


def assign(clusters: list[Cluster]) -> dict[str, str]:
    """member key -> seed key, for stamping a cluster id onto rows."""
    return {m: c.seed for c in clusters for m in c.members}


def medoid(sequences: list[str]) -> str:
    """The member with the smallest total distance to the others: a cluster's reference.

    Ties break by length then sequence, so the choice is deterministic. With one member the
    answer is itself; with two the total distances are equal and the longer one wins, which
    is a coin toss between two real proteins and costs nothing either way.
    """
    if len(sequences) == 1:
        return sequences[0]
    distance = {
        (a, b): align.edit_profile(a, b).n_edits for a in sequences for b in sequences if a != b
    }
    return min(
        sequences,
        key=lambda c: (sum(distance[(c, o)] for o in sequences if o != c), -len(c), c),
    )


def references(clusters: list[Cluster], sequences: dict[str, str] | None = None) -> dict[str, str]:
    """member key -> the key of the reference its edits should be expressed against.

    `sequences` is the same map `cluster` was given; omit it when the keys *are* the
    sequences, which is how the corpus calls this.
    """
    seqs = sequences or {}
    out: dict[str, str] = {}
    for c in clusters:
        by_seq = {seqs.get(m, m): m for m in c.members}
        ref = by_seq[medoid(list(by_seq))]
        for m in c.members:
            out[m] = ref
    return out


# ---------------------------------------------------------------------------------------
# The cluster inventory: one row per cluster, written once and read cheaply.
#
# "Give me clusters with at least 5 domains" is a training-time question (`docs/DECISIONS.md`
# section 1), and answering it from the row tables means grouping 45 million rows by `wt_id`
# and counting distinct `dbd_seq` -- a minute of work to learn 1,133 numbers. A stored
# `cluster_size` column would answer it by Parquet predicate pushdown, but it would touch
# `schema.py` and every parser's output to denormalise a per-cluster fact onto 45 million
# rows. The side table costs one file, changes no parser, and keeps the 22-column schema
# frozen.
# ---------------------------------------------------------------------------------------

#: Columns of the inventory. `n_domains` is what "cluster size" means throughout the project:
#: distinct canonical domains, NOT rows and NOT constructs -- one domain assayed by two
#: sources is one domain.
INVENTORY_COLUMNS = (
    "wt_id",
    "dbd_family",
    "seed",
    "reference",
    "n_domains",
    "n_variants",
    "max_edits",
    "n_sources",
    "n_constructs",
    "n_rows",
)


def inventory(domains: pd.DataFrame) -> pd.DataFrame:
    """Build the per-cluster inventory from per-construct rows.

    `domains` needs one row per (source, construct) with `wt_id`, `dbd_seq`, `dbd_family`,
    `n_mut_from_wt`, `source_dataset` and `n_rows`. Aggregating constructs rather than corpus
    rows keeps this cheap; `n_rows` is summed from what the caller counted.
    """
    grouped = domains.groupby("wt_id", sort=True)
    out = pd.DataFrame(
        {
            "dbd_family": grouped["dbd_family"].first(),
            # The reference is the domain everything else is described against: the one at
            # zero edits from itself. Exactly one per cluster, by construction.
            "reference": grouped.apply(
                lambda g: g.loc[g["n_mut_from_wt"].idxmin(), "dbd_seq"], include_groups=False
            ),
            # The seed is what membership was decided against, carried through so a cluster
            # can be traced back to the assignment that formed it.
            "seed": grouped["seed"].first(),
            "n_domains": grouped["dbd_seq"].nunique(),
            "n_variants": grouped.apply(
                lambda g: int(g.loc[g["n_mut_from_wt"] > 0, "dbd_seq"].nunique()),
                include_groups=False,
            ),
            "max_edits": grouped["n_mut_from_wt"].max().astype("int32"),
            "n_sources": grouped["source_dataset"].nunique(),
            "n_constructs": grouped.size(),
            "n_rows": grouped["n_rows"].sum(),
        }
    ).reset_index()
    return (
        out[list(INVENTORY_COLUMNS)]
        .sort_values(["n_domains", "wt_id"], ascending=[False, True])
        .reset_index(drop=True)
    )


def load(path: str | Path | None = None) -> pd.DataFrame:
    """Read the cluster inventory. Raises if `build_clusters.py` has not been run."""
    p = Path(path) if path else CLUSTER_TABLE
    if not p.exists():
        raise FileNotFoundError(f"no cluster inventory at {p}\nrun scripts/build_clusters.py")
    return pd.read_parquet(p)


def ids_with_at_least(n_domains: int, path: str | Path | None = None) -> set[str]:
    """`wt_id`s whose cluster holds at least `n_domains` distinct domains."""
    inv = load(path)
    return set(inv.loc[inv["n_domains"] >= n_domains, "wt_id"])


def select(df: pd.DataFrame, min_domains: int, path: str | Path | None = None) -> pd.DataFrame:
    """Rows belonging to clusters of at least `min_domains` domains.

    The point of the inventory: this is a membership test against a set of ids, not a
    group-by over the rows being filtered.
    """
    return df[df["wt_id"].isin(ids_with_at_least(min_domains, path))]
