"""Carving a validation set out of a fold's training side.

A `Fold` is a train/test division and nothing else, so early stopping on its test set would
leak (`docs/TRAINING.md` §7a). These pin that the carve never touches the test set and that it
respects the grouping it claims to.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from snp2prot import distances, splits

BASE = "ACDEFGHIKLMNPQRSTVWY" * 4


def mutate(seq: str, *positions: int) -> str:
    out = list(seq)
    for p in positions:
        out[p] = "P" if out[p] != "P" else "G"
    return "".join(out)


@pytest.fixture
def corpus_frame() -> pd.DataFrame:
    """One near-identical pair and 18 mutually unrelated domains.

    Deliberately *not* a series of progressively-mutated copies of one sequence: single linkage
    chains those into a single component, which leaves a fold with no training set at all and
    tests the fallback rather than the grouping.
    """
    rng = np.random.default_rng(0)
    alphabet = np.array(list("ACDEFGHIKLMNPQRSTVWY"))
    rows = [(BASE, "Homeodomain", "wt", False), (mutate(BASE, 5), "Homeodomain", "wt", True)]
    for i in range(18):
        rows.append(("".join(rng.choice(alphabet, 80)), "Homeodomain", f"s{i}", False))
    frame = pd.DataFrame(rows, columns=["dbd_seq", "dbd_family", "wt_id", "is_variant"])
    return frame.drop_duplicates("dbd_seq").sort_values("dbd_seq").reset_index(drop=True)


@pytest.fixture
def dist(corpus_frame):
    return distances.build(corpus_frame.dbd_seq.to_numpy(), processes=1)


def test_the_validation_set_comes_out_of_training_and_never_out_of_test(corpus_frame, dist):
    fold = splits.s1_random(corpus_frame, n_folds=4, seed=1)[0]
    train, validation = splits.validation_split(fold, corpus_frame, dist, 0.25, seed=1)
    assert not set(validation) & set(fold.test)
    assert set(train) | set(validation) == set(fold.train)
    assert not set(train) & set(validation)


def test_a_zero_fraction_means_a_fixed_step_budget_and_no_validation(corpus_frame, dist):
    fold = splits.s1_random(corpus_frame, n_folds=4, seed=1)[0]
    train, validation = splits.validation_split(fold, corpus_frame, dist, 0.0, seed=1)
    assert list(train) == list(fold.train)
    assert len(validation) == 0


def test_component_grouping_keeps_near_identical_domains_on_one_side(corpus_frame, dist):
    """Otherwise the validation set contains a near-twin of something in training, and early
    stopping selects the checkpoint that memorised best."""
    fold = splits.s1_random(corpus_frame, n_folds=2, seed=3)[0]
    train, validation = splits.validation_split(
        fold, corpus_frame, dist, 0.3, seed=3, grouping="component"
    )
    components = distances.connected_components(dist, 0.9, 0.6, corpus_frame.dbd_family.to_numpy())
    assert not set(components[train]) & set(components[validation])


def test_s2_uses_component_grouping_without_being_told_to(corpus_frame, dist):
    """`grouping="regime"` means "the regime's own unit", and S2's unit is the component."""
    fold = splits.s2_components(corpus_frame, dist, n_folds=2, seed=1, min_identity=0.9)[0]
    train, validation = splits.validation_split(fold, corpus_frame, dist, 0.3, seed=1)
    components = distances.connected_components(dist, 0.9, 0.6, corpus_frame.dbd_family.to_numpy())
    assert len(validation)
    assert not set(components[train]) & set(components[validation])


def test_the_carve_is_reproducible_from_its_seed(corpus_frame, dist):
    fold = splits.s1_random(corpus_frame, n_folds=4, seed=1)[0]
    first = splits.validation_split(fold, corpus_frame, dist, 0.25, seed=5)[1]
    again = splits.validation_split(fold, corpus_frame, dist, 0.25, seed=5)[1]
    other = splits.validation_split(fold, corpus_frame, dist, 0.25, seed=6)[1]
    assert list(first) == list(again)
    assert list(first) != list(other)


def test_a_grouping_that_would_leave_no_training_data_falls_back(corpus_frame, dist):
    """Single linkage at a loose floor can put an entire fold in one component. Taking it as
    validation would leave nothing to train on, so the carve degrades to the ungrouped one
    rather than returning an empty training set mid-grid."""
    fold = splits.s1_random(corpus_frame, n_folds=2, seed=1)[0]
    train, validation = splits.validation_split(
        fold, corpus_frame, dist, 0.3, seed=1, grouping="component"
    )
    # A floor of 0.0 joins everything into one component, which is the degenerate case.
    degenerate = splits.validation_split(
        fold, corpus_frame, dist, 0.3, seed=1, grouping="component", min_identity=0.0
    )
    assert len(train) and len(validation)
    assert len(degenerate[0]) and len(degenerate[1])


def test_the_validation_set_gets_its_own_digest(corpus_frame, dist):
    """A run records which domains were held out, and that applies to the validation slice as
    much as to the test set (`ML_PLAN.md` §9.2)."""
    fold = splits.s1_random(corpus_frame, n_folds=4, seed=1)[0]
    _, validation = splits.validation_split(fold, corpus_frame, dist, 0.25, seed=1)
    digest = splits.digest_of(corpus_frame, validation)
    assert len(digest) == 12
    assert digest != splits.digest_of(corpus_frame, fold.test)
    assert digest == splits.digest_of(corpus_frame, np.array(list(reversed(validation))))
