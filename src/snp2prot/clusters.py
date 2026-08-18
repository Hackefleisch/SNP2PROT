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
from pathlib import Path

import pandas as pd

from snp2prot import align
from snp2prot.config import CLUSTER_TABLE


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
            # against representatives, so the cost of checking honestly is small.
            profile = align.edit_profile(rseq, seq)
            if profile.n_edits <= max_edits and profile.overlap >= min_overlap:
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
    "representative",
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
            # The representative is the domain the cluster id was named after: the one at
            # zero edits from itself. Ties cannot occur -- CD-HIT has exactly one.
            "representative": grouped.apply(
                lambda g: g.loc[g["n_mut_from_wt"].idxmin(), "dbd_seq"], include_groups=False
            ),
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
