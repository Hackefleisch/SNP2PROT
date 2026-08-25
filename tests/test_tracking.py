"""The experiment-tracking wrapper.

Two properties matter and neither is about MLflow: a run must be unable to exist without a
record of what it held out, and nothing in the training path may break when the tracker is
absent.
"""

from __future__ import annotations

import pytest

from snp2prot import tracking


def test_a_run_without_digests_is_refused():
    """`ML_PLAN.md` §9.2: log the split, not just the hyperparameters. A regime name and a seed
    do not reconstruct which domains were in test once the corpus grows."""
    with pytest.raises(ValueError, match="held out"):
        with tracking.start_run("x", digests={}):
            pass


def test_a_disabled_run_is_a_working_no_op():
    """Training must not depend on the tracker being installed or enabled."""
    with tracking.start_run("x", digests={"test": "abc"}, enabled=False) as run:
        assert not run.active
        run.log_metrics({"aupr": 0.5}, step=1)
        run.log_artifact(__file__)


def test_nested_config_becomes_flat_parameters():
    """MLflow params are flat, and the config is not."""
    flat = tracking._flatten({"model": {"dna": {"channels": 64}, "width": 256}, "seed": 1})
    assert flat == {"model.dna.channels": 64, "model.width": 256, "seed": 1}


def test_non_finite_metrics_are_dropped_rather_than_logged():
    """A domain with no positives has an undefined AUPR, and `nan` in a tracked metric turns a
    whole run's chart into a gap."""
    assert not tracking._finite(float("nan"))
    assert not tracking._finite(float("inf"))
    assert tracking._finite(0.0)
