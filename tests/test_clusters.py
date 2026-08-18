"""Cluster construction.

Clusters exist to keep near-identical domains out of opposite sides of a split. They never
collapse rows: every member stays its own (sequence, 8-mer) -> label training example.
"""

from __future__ import annotations

from snp2prot import clusters

BASE = "ACDEFGHIKLMNPQRSTVWY" * 4  # 80 aa


def mutate(seq: str, *positions: int) -> str:
    out = list(seq)
    for p in positions:
        out[p] = "P" if out[p] != "P" else "G"
    return "".join(out)


def test_near_identical_domains_land_in_one_cluster():
    seqs = {"wt": BASE, "v1": mutate(BASE, 5), "v2": mutate(BASE, 5, 11), "far": BASE[::-1]}
    got = clusters.cluster(seqs, max_edits=3)
    assign = clusters.assign(got)
    assert assign["wt"] == assign["v1"] == assign["v2"]
    assert assign["far"] != assign["wt"]


def test_the_longest_member_is_the_seed():
    """So the least-clipped form of a domain becomes the reference, not whichever came first."""
    seqs = {"short": BASE[5:75], "long": BASE, "mid": BASE[2:78]}
    got = clusters.cluster(seqs, max_edits=3)
    assert len(got) == 1
    assert got[0].seed == "long"


def test_clustering_is_not_transitive():
    """A near B near C, with A and C far apart, must not become one cluster.

    Single-linkage on the real corpus produced a 35-domain blob at 5 edits against 8 at 1 —
    chaining is why CD-HIT compares only against seeds.
    """
    a = BASE
    b = mutate(BASE, 1, 2, 3)
    c = mutate(BASE, 1, 2, 3, 40, 41, 42)
    got = clusters.cluster({"a": a, "b": b, "c": c}, max_edits=3)
    assign = clusters.assign(got)
    assert assign["a"] == assign["b"], "b is within 3 of a"
    assert assign["c"] != assign["a"], "c is 6 from a and must not chain in through b"


def test_different_families_are_never_compared():
    seqs = {"x": BASE, "y": BASE}
    fams = {"x": "Homeodomain", "y": "Forkhead"}
    got = clusters.cluster(seqs, max_edits=3, families=fams)
    assert len(got) == 2, "identical sequences of different families must stay apart"


def test_result_does_not_depend_on_input_order():
    seqs = {"wt": BASE, "v1": mutate(BASE, 5), "v2": mutate(BASE, 9)}
    a = clusters.assign(clusters.cluster(seqs, max_edits=2))
    b = clusters.assign(clusters.cluster(dict(reversed(list(seqs.items()))), max_edits=2))
    assert a == b


def test_every_member_is_within_the_threshold_of_its_seed():
    """The invariant mut_positions depends on: one coordinate frame per cluster."""
    from snp2prot import align

    seqs = {f"s{i}": mutate(BASE, *range(i)) for i in range(6)}
    for c in clusters.cluster(seqs, max_edits=2):
        for m in c.members:
            assert align.edit_profile(seqs[c.seed], seqs[m]).n_edits <= 2


def test_a_minimum_overlap_refuses_a_membership_built_on_free_terminal_gaps():
    """31 of the corpus's 202 memberships were formed this way, aligning on 3-10% of the
    shorter sequence (`docs/DECISIONS.md`, `T25`). Without the floor these two join; with it
    they are two clusters."""
    seqs = {"a": "WWWWWWWWWW" + "KKKKK", "b": "KKKKK" + "YYYYYYYYYY"}
    fams = {"a": "F", "b": "F"}

    unguarded = clusters.cluster(seqs, max_edits=5, families=fams)
    assert len(unguarded) == 1

    guarded = clusters.cluster(seqs, max_edits=5, families=fams, min_overlap=0.6)
    assert len(guarded) == 2


def test_the_floor_does_not_separate_a_domain_from_its_own_clipped_copy():
    """Which is the whole reason terminal gaps are free: same domain, less padding."""
    full = "MMMMMMMMMM" + "KKKKKAAAAAKKKKKAAAAA" + "WWWWWWWWWW"
    clipped = "KKKKKAAAAAKKKKKAAAAA"
    seqs = {"full": full, "clipped": clipped}
    fams = {"full": "F", "clipped": "F"}
    guarded = clusters.cluster(seqs, max_edits=5, families=fams, min_overlap=0.6)
    assert len(guarded) == 1


def test_the_reference_is_the_medoid_not_the_seed():
    """A wild type plus k single substitutions: the wild type sits one edit from each, any
    variant sits two from the rest, so the medoid is the wild type. Before this, the seed was
    settled by the alphabetical order of the sequence, which made `HOXD13_S316C` the frame
    for its own wild type and put position 50 in all seven siblings' `mut_positions`."""
    wt = "ACDEFGHIKLMNPQRSTVWY" * 2
    variants = {f"v{i}": wt[:i] + "W" + wt[i + 1 :] for i in (3, 9, 17, 25)}
    seqs = {"wt": wt, **variants}
    fams = dict.fromkeys(seqs, "F")

    cs = clusters.cluster(seqs, max_edits=2, families=fams)
    assert len(cs) == 1
    assert clusters.references(cs, seqs)["wt"] == "wt"


def test_the_reference_of_a_singleton_is_itself():
    seqs = {"only": "ACDEFGHIKL"}
    cs = clusters.cluster(seqs, max_edits=2, families={"only": "F"})
    assert clusters.references(cs, seqs) == {"only": "only"}


def test_the_medoid_is_deterministic_when_two_members_tie():
    """With two members the total distances are equal; the longer one wins, then the
    sequence. A coin toss between two real proteins, but it must not depend on row order."""
    a, b = "ACDEFGHIKLMNPQ", "ACDEFGHIKLMNPQRST"
    assert clusters.medoid([a, b]) == clusters.medoid([b, a]) == b
