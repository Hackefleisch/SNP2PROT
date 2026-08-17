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


def test_the_longest_member_is_the_representative():
    """So the least-clipped form of a domain becomes the reference, not whichever came first."""
    seqs = {"short": BASE[5:75], "long": BASE, "mid": BASE[2:78]}
    got = clusters.cluster(seqs, max_edits=3)
    assert len(got) == 1
    assert got[0].representative == "long"


def test_clustering_is_not_transitive():
    """A near B near C, with A and C far apart, must not become one cluster.

    Single-linkage on the real corpus produced a 35-domain blob at 5 edits against 8 at 1 —
    chaining is why CD-HIT compares only against representatives.
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


def test_every_member_is_within_the_threshold_of_its_representative():
    """The invariant mut_positions depends on: one coordinate frame per cluster."""
    from snp2prot import align

    seqs = {f"s{i}": mutate(BASE, *range(i)) for i in range(6)}
    for c in clusters.cluster(seqs, max_edits=2):
        for m in c.members:
            assert align.edit_profile(seqs[c.representative], seqs[m]).n_edits <= 2
