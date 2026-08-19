"""Split regimes.

The regimes exist to make one failure visible: a model that memorises per-protein profiles
scores well wherever a near-relative of the held-out domain stayed in training. Every test
here is about what a regime *removes*, because that is the only thing that makes its number
mean anything.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from snp2prot import distances, splits

BASE = "ACDEFGHIKLMNPQRSTVWY" * 4  # 80 aa


def mutate(seq: str, *positions: int) -> str:
    out = list(seq)
    for p in positions:
        out[p] = "P" if out[p] != "P" else "G"
    return "".join(out)


@pytest.fixture
def domains() -> pd.DataFrame:
    """A wild type with two variants, a near-identical domain in its own cluster, and six
    unrelated singletons across two families."""
    rows = [
        (BASE, "Homeodomain", "wt", False),
        (mutate(BASE, 5), "Homeodomain", "wt", True),
        (mutate(BASE, 5, 9), "Homeodomain", "wt", True),
        # One edit from the wild type but clustered separately: the T27 case S2 exists for.
        (mutate(BASE, 70), "Homeodomain", "neighbour", False),
    ]
    for i in range(6):
        family = "Forkhead" if i % 2 else "Homeodomain"
        rows.append((BASE[::-1][:60] + "AC" * (i + 1), family, f"single{i}", False))
    frame = pd.DataFrame(rows, columns=["dbd_seq", "dbd_family", "wt_id", "is_variant"])
    return frame.sort_values("dbd_seq").reset_index(drop=True)


def test_a_fold_holds_out_domains_and_the_two_halves_partition_the_corpus(domains):
    for fold in splits.s1_random(domains, n_folds=5, seed=1):
        assert not set(fold.test) & set(fold.train)
        assert len(fold.test) + len(fold.train) == len(domains)


def test_every_domain_is_held_out_exactly_once_across_the_folds(domains):
    seen = np.concatenate([f.test for f in splits.s1_random(domains, n_folds=5, seed=1)])
    assert sorted(seen) == list(range(len(domains)))


def test_s2_holds_out_near_identical_domains_together_even_across_clusters(domains):
    """`T27`: greedy clustering compares only to seeds, so a one-edit pair can land in two
    clusters. Holding out by cluster would leave one of them in training."""
    d = distances.build(domains.dbd_seq.to_numpy(), processes=1)
    components = distances.connected_components(d, min_identity=0.9, min_overlap=0.6)
    wt = int(np.flatnonzero(domains.dbd_seq == BASE)[0])
    near = int(np.flatnonzero(domains.dbd_seq == mutate(BASE, 70))[0])
    assert domains.wt_id[wt] != domains.wt_id[near]  # different clusters ...
    assert components[wt] == components[near]  # ... one component

    for fold in splits.s2_components(domains, d, n_folds=2, seed=1, min_identity=0.9):
        assert (wt in set(fold.test)) == (near in set(fold.test))


def test_components_do_not_cross_families_when_families_are_given(domains):
    d = distances.build(domains.dbd_seq.to_numpy(), processes=1)
    labels = distances.connected_components(
        d, min_identity=0.0, min_overlap=0.0, families=domains.dbd_family.to_numpy()
    )
    for component in set(labels):
        assert domains.dbd_family[labels == component].nunique() == 1


def test_a_looser_identity_floor_groups_more_domains_together(domains):
    """`D5`: the floor is what decides how hard S2 is, so it has to be the thing that moves.

    At 0.9 only the near-identical domains group; at 0.4 the unrelated singletons of one
    family chain into it, which is the safe direction for a split — over-grouping removes more
    from training than necessary but can never leak.
    """
    # 80 aa apiece: one edit apart, then 16 edits apart (80% identity), then 32 (60%).
    sequences = np.array(
        [BASE, mutate(BASE, 1), mutate(BASE, *range(0, 64, 4)), mutate(BASE, *range(0, 64, 2))]
    )
    d = distances.build(sequences, processes=1)
    families = np.array(["Homeodomain"] * 4)

    tight = distances.connected_components(d, 0.9, 0.6, families)
    assert len(set(tight)) == 3  # only the one-edit pair groups
    loose = distances.connected_components(d, 0.5, 0.6, families)
    assert len(set(loose)) == 1  # all four chain into one
    assert np.bincount(loose).max() > np.bincount(tight).max()


def test_p1_removes_a_whole_family_and_leaves_the_rest(domains):
    (fold,) = splits.p1_family(domains, "Homeodomain")
    assert set(domains.dbd_family.to_numpy()[fold.test]) == {"Homeodomain"}
    assert "Homeodomain" not in set(domains.dbd_family.to_numpy()[fold.train])


def test_p2_holds_out_whole_variant_clusters_but_keeps_the_singletons(domains):
    """The gentler mirror of P1: the family stays in training, the mutated members do not."""
    (fold,) = splits.p2_variant_clusters(domains, exclude_family="Forkhead")
    held = domains.iloc[fold.test]
    assert set(held.wt_id) == {"wt"}
    assert held.is_variant.sum() == 2  # the two variants and their reference
    assert (~domains.iloc[fold.train].is_variant).all()


def test_p3_holds_out_variants_only_so_every_wild_type_stays_in_training(domains):
    """Which is what makes the baseline the exact null hypothesis for C1: with no variant in
    training, a held-out variant's nearest neighbour is its own wild type."""
    fold = splits.p3_variants(domains, fraction=1.0, seed=1, name="all")
    assert domains.is_variant.to_numpy()[fold.test].all()
    assert not domains.is_variant.to_numpy()[fold.train].any()
    assert BASE in set(domains.dbd_seq.to_numpy()[fold.train])


def test_p3_draws_differ_between_repeats_and_repeat_between_runs(domains):
    a = splits.p3_variants(domains, 0.5, seed=7, name="half", repeat=0)
    b = splits.p3_variants(domains, 0.5, seed=7, name="half", repeat=1)
    again = splits.p3_variants(domains, 0.5, seed=7, name="half", repeat=0)
    assert list(a.test) == list(again.test)
    assert len(a.test) == len(b.test) == 1


def test_the_digest_identifies_the_held_out_set_not_the_fold_name(domains):
    a, b = splits.s1_random(domains, n_folds=2, seed=1)
    assert a.digest(domains) != b.digest(domains)
    same = splits.s1_random(domains, n_folds=2, seed=1)[0]
    assert same.digest(domains) == a.digest(domains)
