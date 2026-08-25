"""Experiment tracking, behind a wrapper thin enough to swap.

`ML_PLAN.md` §9.2 chose **MLflow on a local backend** — no server, no cloud, which matches a
project whose data is git-ignored and whose provenance is deliberately self-contained. The
owner's brief was "MLflow, W&B, whatever", and this wrapper is what makes the "whatever" a
one-line change later instead of a rewrite.

**The backend is a local SQLite file, not the `mlruns/` file store the plan named.** MLflow now
refuses the filesystem store outright — it is in maintenance mode and raises unless an opt-out
environment variable is set. SQLite satisfies every requirement §9.2 actually gave (local, no
server, no cloud, one git-ignored path) and is the supported route, so the change is to the URI
and to nothing else. Both the database and the artifacts live under `mlruns/`.

**The rule this module exists to enforce: log the split, not just the hyperparameters.** It is
the project's own standard applied to modelling. Rule 5 makes every row of the dataset carry the
cutoffs actually applied, so a table describes its own binarization; a run is the same. A metric
is meaningless without knowing exactly which domains were held out, and "regime `S2`, fold 3,
seed 20260819" does not pin that down once the corpus changes — so `start_run` requires the
digests and refuses a run without them.

Tracking is optional at the edges: if MLflow is not installed, or `enabled=False`, every call
becomes a no-op and the run still produces its report and its artifacts. Nothing in the training
path may depend on the tracker being there.
"""

from __future__ import annotations

import contextlib
from collections.abc import Iterator
from pathlib import Path
from typing import Any

from snp2prot.config import PROJECT_ROOT

#: The local backend. Git-ignored already.
TRACKING_DIR = PROJECT_ROOT / "mlruns"
DATABASE = TRACKING_DIR / "mlflow.db"
ARTIFACT_DIR = TRACKING_DIR / "artifacts"


def available() -> bool:
    try:
        import mlflow  # noqa: F401
    except ImportError:
        return False
    return True


@contextlib.contextmanager
def start_run(
    name: str,
    experiment_name: str = "snp2prot",
    params: dict[str, Any] | None = None,
    digests: dict[str, str] | None = None,
    enabled: bool = True,
) -> Iterator[Run]:
    """A tracked run. `digests` names the held-out sets and is mandatory.

    Not optional because a run without it cannot be interpreted later: the fold name and the
    seed do not reconstruct which domains were in test once the corpus grows.
    """
    if not digests:
        raise ValueError(
            "a run must record which domains were held out, not just its hyperparameters "
            "(ML_PLAN.md §9.2) — pass digests={'test': ..., 'validation': ...}"
        )
    if not (enabled and available()):
        yield Run(None)
        return

    import mlflow

    TRACKING_DIR.mkdir(parents=True, exist_ok=True)
    ARTIFACT_DIR.mkdir(parents=True, exist_ok=True)
    mlflow.set_tracking_uri(f"sqlite:///{DATABASE}")
    if mlflow.get_experiment_by_name(experiment_name) is None:
        mlflow.create_experiment(experiment_name, artifact_location=ARTIFACT_DIR.as_uri())
    mlflow.set_experiment(experiment_name)
    with mlflow.start_run(run_name=name):
        mlflow.log_params(_flatten(params or {}))
        mlflow.log_params({f"digest.{k}": v for k, v in digests.items()})
        yield Run(mlflow)


class Run:
    """A handle that is either MLflow or nothing, with the same surface either way."""

    def __init__(self, backend):
        self._backend = backend

    @property
    def active(self) -> bool:
        return self._backend is not None

    def log_metrics(self, metrics: dict[str, float], step: int | None = None) -> None:
        if self._backend is not None:
            clean = {k: float(v) for k, v in metrics.items() if _finite(v)}
            self._backend.log_metrics(clean, step=step)

    def log_artifact(self, path: str | Path) -> None:
        if self._backend is not None:
            self._backend.log_artifact(str(path))


def _finite(value: Any) -> bool:
    try:
        return float(value) == float(value) and abs(float(value)) != float("inf")
    except (TypeError, ValueError):
        return False


def _flatten(params: dict[str, Any], prefix: str = "") -> dict[str, Any]:
    """Nested config dicts become dotted keys, because MLflow params are flat."""
    out: dict[str, Any] = {}
    for key, value in params.items():
        name = f"{prefix}{key}"
        if isinstance(value, dict):
            out |= _flatten(value, f"{name}.")
        else:
            out[name] = value
    return out
