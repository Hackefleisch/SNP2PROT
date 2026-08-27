"""The per-domain view of the merged table.

It computes no new fact — every column comes from a companion table that is already the
authority on it — so what is tested is that the join keeps the merged table's domain set and
that `is_variant` means what the split regimes assume it means.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from snp2prot import corpus

WT = "RKRGRQTYTRYQTLELEKEFHFNRYLTRRRRIEIAHALCLTERQIKIWFQNRRMKWKKEN"
VARIANT = WT[:10] + "A" + WT[11:]
OTHER = "MGSKKPRLIWTPQLHKRFVDAVAHLGIKNAVPKTIMQLMNVEGLTRENVASHLQKYRLYL"


@pytest.fixture
def records() -> pd.DataFrame:
    """Three domain records, one of them the same domain measured by a second source."""
    rows = [
        (WT, "SRC_A", "C:SRC_A:WT", "Homeodomain", 100, "ok"),
        (VARIANT, "SRC_A", "C:SRC_A:WT", "Homeodomain", 4, "ok"),
        (OTHER, "SRC_A", "C:SRC_A:OTHER", "Myb", 0, "no_evidence"),
        (OTHER, "SRC_B", "C:SRC_A:OTHER", "Myb", 30, "ok"),
    ]
    frame = pd.DataFrame(
        rows, columns=["dbd_seq", "source_dataset", "wt_id", "dbd_family", "n_pos", "verdict"]
    )
    return frame.assign(max_escore=0.5, protein_id="P1", gene="g")


@pytest.fixture
def inventory() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "wt_id": ["C:SRC_A:WT", "C:SRC_A:OTHER"],
            "reference": [WT, OTHER],
            "n_domains": [2, 1],
        }
    )


def test_one_row_per_domain_survives_a_duplicate(records, inventory, monkeypatch):
    monkeypatch.setattr(corpus.clusters, "load", lambda *a, **k: inventory)
    got = corpus.domains(records)
    assert len(got) == 3
    assert got.dbd_seq.is_unique
    # The deeper of the two measurements of OTHER is the one that reaches the merged table.
    assert got.loc[got.dbd_seq == OTHER, "source_dataset"].item() == "SRC_B"


def test_a_domain_is_a_variant_when_it_is_not_its_reference(records, inventory, monkeypatch):
    monkeypatch.setattr(corpus.clusters, "load", lambda *a, **k: inventory)
    got = corpus.domains(records).set_index("dbd_seq")
    assert not got.loc[WT, "is_variant"]
    assert got.loc[VARIANT, "is_variant"]
    assert not got.loc[OTHER, "is_variant"]


def test_rows_are_sorted_by_sequence_because_that_is_an_identity(records, inventory, monkeypatch):
    """The 8-mer matrix and the distance matrix are both indexed by it, and a mismatched pair
    is otherwise undetectable."""
    monkeypatch.setattr(corpus.clusters, "load", lambda *a, **k: inventory)
    got = corpus.domains(records)
    assert list(got.dbd_seq) == sorted(got.dbd_seq)


def test_the_filter_drops_no_evidence_and_keeps_dead_variants(records, inventory, monkeypatch):
    """`T21`: a silent record with no control cannot be told from a failed assay, while a
    variant that measurably lost binding is the most informative negative in the corpus."""
    monkeypatch.setattr(corpus.clusters, "load", lambda *a, **k: inventory)
    got = corpus.domains(records)
    got.loc[got.dbd_seq == VARIANT, "verdict"] = "dead_variant"
    got.loc[got.dbd_seq == OTHER, "verdict"] = "no_evidence"
    mask = corpus.trainable(got)
    assert mask[got.dbd_seq == VARIANT].item()
    assert not mask[got.dbd_seq == OTHER].item()
    assert not corpus.trainable(got, keep_dead=False)[got.dbd_seq == VARIANT].item()


# --- the row-order guarantee ------------------------------------------------------------


def test_require_aligned_accepts_tables_in_the_same_order():
    corpus.require_aligned(
        ["A", "B", "C"], matrix=["A", "B", "C"], distances=np.array(["A", "B", "C"])
    )


def test_require_aligned_catches_a_reordering_that_a_set_check_would_miss():
    """The realistic failure: same domains, different order after a rebuild. A length check and
    a set check both pass; only comparing the sequences catches it."""
    reference = ["A", "B", "C"]
    reordered = ["C", "B", "A"]
    assert len(reordered) == len(reference) and set(reordered) == set(reference)
    with pytest.raises(ValueError, match="first difference at row 0"):
        corpus.require_aligned(reference, matrix=reordered)


def test_require_aligned_names_the_table_and_reports_the_lengths():
    with pytest.raises(ValueError, match=r"embeddings table .* 2 domains against 3"):
        corpus.require_aligned(["A", "B", "C"], embeddings=["A", "B"])


def test_require_aligned_uses_the_reference_name_it_is_given():
    with pytest.raises(ValueError, match="the 8-mer matrix has"):
        corpus.require_aligned(["A", "B"], "the 8-mer matrix", embeddings=["A", "Z"])


def test_require_aligned_checks_every_table_not_just_the_first():
    with pytest.raises(ValueError, match="distances table"):
        corpus.require_aligned(["A", "B"], matrix=["A", "B"], distances=["A", "Z"])
