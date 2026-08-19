"""The C1 evaluation set.

It answers `D6`: `P3/all`'s mean AUPR is 0.928 with a median of 1.000, so the variants that
carry claim C1 are invisible in it. The set is the tail, and what these tests pin is which
variants land in it and why — a wrong membership rule would quietly redefine the claim.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from snp2prot.evaluation import c1


@pytest.fixture
def per_domain() -> pd.DataFrame:
    """One `P3/all` row per variant, plus a row from another fold that must be ignored."""
    rows = [
        # domain, fold, k, aupr, verdict, family, spearman
        ("well_predicted", "all", 1, 0.99, "ok", "Homeodomain", 0.95),
        ("just_above_cut", "all", 1, 0.75, "ok", "Homeodomain", 0.80),
        ("poorly_predicted", "all", 1, 0.31, "ok", "Homeodomain", 0.40),
        ("lost_binding", "all", 1, np.nan, "dead_variant", "Forkhead", 0.30),
        ("unverifiable", "all", 1, np.nan, "no_evidence", "Forkhead", 0.20),
        ("other_fold", "half:draw-0", 1, 0.10, "ok", "Homeodomain", 0.10),
        ("top_k_copy", "all", 5, 0.10, "ok", "Homeodomain", 0.10),
    ]
    return pd.DataFrame(
        rows, columns=["domain", "fold", "k", "aupr", "verdict", "family", "spearman"]
    )


def test_it_takes_the_poorly_predicted_and_the_ones_with_no_aupr_at_all(per_domain):
    got = c1.select(per_domain, max_baseline_aupr=0.7)
    assert set(got.domain) == {"poorly_predicted", "lost_binding"}


def test_a_variant_whose_wild_type_predicts_it_is_not_in_the_set(per_domain):
    """The point of the set: most single substitutions do not change what a domain binds, and
    those variants have no headroom to measure C1 in."""
    got = c1.select(per_domain, max_baseline_aupr=0.7)
    assert "well_predicted" not in set(got.domain)
    assert "just_above_cut" not in set(got.domain)


def test_no_evidence_records_are_excluded_because_silence_is_not_evidence(per_domain):
    """`T21`: with no control there is no telling a true non-binder from a failed assay."""
    assert "unverifiable" not in set(c1.select(per_domain).domain)
    assert "unverifiable" in set(c1.select(per_domain, keep_no_evidence=True).domain)


def test_only_the_p3_all_rows_at_k_one_are_considered(per_domain):
    """Other folds hold overlapping variants and the top-k form is a different predictor;
    either would put the same variant in the set twice."""
    got = c1.select(per_domain, max_baseline_aupr=0.7)
    assert "other_fold" not in set(got.domain)
    assert "top_k_copy" not in set(got.domain)
    assert got.domain.is_unique


def test_losing_binding_is_recorded_as_a_different_reason_from_scoring_low(per_domain):
    """They are not the same evidence: one is a profile predicted badly, the other is no
    profile to predict."""
    got = c1.select(per_domain, max_baseline_aupr=0.7).set_index("domain")
    assert got.loc["lost_binding", "c1_reason"] == c1.LOST_BINDING
    assert got.loc["poorly_predicted", "c1_reason"] == c1.POORLY_PREDICTED


def test_the_worst_predicted_variants_come_first(per_domain):
    got = c1.select(per_domain, max_baseline_aupr=0.7)
    assert got.domain.iloc[0] == "lost_binding"  # no AUPR at all ranks ahead of a low one


def test_the_summary_is_stratified_by_family(per_domain):
    """84 of 173 variants are homeodomain, so a pooled number is a homeodomain number."""
    variants = per_domain[(per_domain.fold == "all") & (per_domain.k == 1)]
    got = c1.summarise(c1.select(per_domain, max_baseline_aupr=0.7), variants).set_index("family")
    assert got.loc["Homeodomain", "n_variants"] == 3
    assert got.loc["Homeodomain", "n_in_set"] == 1
    assert got.loc["Forkhead", "n_dead"] == 1
