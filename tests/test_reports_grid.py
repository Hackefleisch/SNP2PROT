"""`write_report` against a frame with exactly the columns `run_fold` produces.

The bug this exists to stop: a report writer reading a column that no longer exists, or a section
that is defined and never called. Both are invisible until the very end of a five-hour grid, and
both happened — `row.spearman` survived the metric's deletion, and four sections were written,
unit-tested in isolation, and never wired into `write_report`.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
run_grid = pytest.importorskip("run_grid")

from snp2prot import experiment  # noqa: E402


def _summary_columns() -> set[str]:
    """Every key `training.run_fold` puts on a summary row, read off the source.

    Deliberately not a hand-written list: a hand-written list drifts from the producer, which is
    exactly the failure mode here.
    """
    import re

    src = (Path(__file__).resolve().parents[1] / "src/snp2prot/training.py").read_text()
    body = src.split("summary |= {")[1].split("\n    }")[0]
    keys = set(re.findall(r'"(\w+)":', body))
    keys |= set(re.findall(r'summary\["(\w+)"\]\s*=', src))
    keys |= set(re.findall(r'"(\w+)": _nanmean|"(\w+)": float', src)) and keys
    # macro_average's own keys, which land on the same row
    keys |= {
        "n_domains",
        "n_scored",
        "n_undefined",
        "aupr",
        "aupr_median",
        "auroc",
        "recall_at_precision",
        "precision_at_10",
        "precision_at_50",
        "precision_at_100",
        "random_mean",
        "random_p95",
    }
    keys |= {f"n_scored_{m}" for m in ("aupr", "auroc", "recall_at_precision", "precision_at_50")}
    # joined on by `run_grid.baseline_by_fold`, not produced by `run_fold`
    keys |= {
        "baseline_aupr",
        "baseline_median",
        "delta",
        "baseline_k5_aupr",
        "baseline_k5_median",
        "delta_k5",
    }
    return keys


def _frame(n_seeds: int = 1) -> pd.DataFrame:
    rows = []
    for arm in ("A1", "A4"):
        for regime, fold in (
            ("S1", "fold-0"),
            ("S2", "fold-0"),
            ("P1", "Homeodomain"),
            ("P3", "all"),
        ):
            for i in range(n_seeds):
                row = dict.fromkeys(_summary_columns(), 0.0)
                row |= {
                    "arm": arm,
                    "regime": regime,
                    "fold": fold,
                    "model_seed": 20260819 + i,
                    "seed": 20260819,
                    "digest": "abc123",
                    "code_commit": "9390e41",
                    "code_dirty": "false",
                    "checkpoint": f"x/{arm}_{regime}_{fold}.pt",
                    "aupr": 0.5 + 0.01 * i,
                    "aupr_final": 0.49,
                    "selection_gain": 0.01,
                    "aupr_median": 0.6,
                    "auroc": 0.9,
                    "chance_aupr": 0.002,
                    "random_mean": 0.0021,
                    "random_p95": 0.0031,
                    "best_step": 14000.0,
                    "steps_run": 15000.0,
                    "suppression": 0.49 if regime == "P3" else np.nan,
                    "n_scored_suppression": 18.0 if regime == "P3" else 0.0,
                    "baseline_aupr": 0.4,
                    "baseline_median": 0.4,
                    "delta": 0.1,
                }
                rows.append(row)
    return pd.DataFrame(rows)


def _per_domain() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "arm": ["A1"] * 4,
            "regime": ["P3"] * 4,
            "fold": ["all"] * 4,
            "domain": list("WXYZ"),
            "aupr": [0.1, 0.2, np.nan, 0.4],
        }
    )


def test_the_report_renders_from_exactly_what_run_fold_produces(tmp_path, monkeypatch):
    monkeypatch.setattr(run_grid, "PROCESSED_DIR", tmp_path)
    pd.DataFrame({"domain": ["W"]}).to_parquet(tmp_path / "c1_variants.parquet")
    out = tmp_path / "training.md"
    run_grid.write_report(out, _frame(), _per_domain(), experiment.load())
    assert out.exists() and out.stat().st_size > 0


@pytest.mark.parametrize(
    "heading",
    [
        "## Per fold",
        "## Is it better than random?",
        "## Did it notice the mutation?",
        "## Does selecting on validation beat the model at the budget?",
        "## The C1 evaluation set",
        "## Configuration",
    ],
)
def test_every_section_is_actually_wired_into_the_report(heading, tmp_path, monkeypatch):
    """Four of these were written, unit-tested and never called."""
    monkeypatch.setattr(run_grid, "PROCESSED_DIR", tmp_path)
    pd.DataFrame({"domain": ["W"]}).to_parquet(tmp_path / "c1_variants.parquet")
    out = tmp_path / "training.md"
    run_grid.write_report(out, _frame(), _per_domain(), experiment.load())
    assert heading in out.read_text()


def test_the_report_says_so_when_there_is_no_error_bar(tmp_path, monkeypatch):
    monkeypatch.setattr(run_grid, "PROCESSED_DIR", tmp_path)
    pd.DataFrame({"domain": ["W"]}).to_parquet(tmp_path / "c1_variants.parquet")
    out = tmp_path / "training.md"

    run_grid.write_report(out, _frame(n_seeds=1), _per_domain(), experiment.load())
    assert "One seed per fold" in out.read_text()

    run_grid.write_report(out, _frame(n_seeds=3), _per_domain(), experiment.load())
    text = out.read_text()
    assert "## Across seeds (3 per fold)" in text
    assert "±" in text, "multiple seeds must show a spread on the per-fold table"


def test_no_deleted_metric_is_referenced_anywhere_in_the_writers():
    """Spearman was removed as a model metric; two report writers kept reading the column."""
    for name in ("run_grid.py", "run_nn_baseline.py"):
        src = (Path(__file__).resolve().parents[1] / "scripts" / name).read_text()
        assert "spearman" not in src.lower(), f"{name} still references the deleted metric"
