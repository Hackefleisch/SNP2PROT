"""The per-cluster inventory — `TODO.md` T3's cheap answer to "how big is this cluster?"."""

from __future__ import annotations

import pandas as pd
import pytest

from snp2prot import clusters

WT = "RKRGRQTYTRYQTLELEKEFHFNRYLTRRRRIEIAHALCLTERQIKIWFQNRRMKWKKEN"
VAR = WT[:10] + "A" + WT[11:]
LONE = "MEVTSQSTLPPGFRFHPTDEELIVYYLRNQTMSKPCPVSIIPEVDIYKFDPWQLPEKTEF"


@pytest.fixture
def constructs() -> pd.DataFrame:
    """Two clusters: one with a reference plus a variant measured twice, one singleton.

    `seed` is the variant here, not the reference — the case `D3` is about: membership was
    decided against whichever domain the greedy pass happened to reach first, and the
    coordinate frame is the medoid.
    """
    return pd.DataFrame(
        [
            # same domain, two sources -> two constructs, ONE domain
            ("C:X:Hd", WT, "P1", "Homeodomain", 0, "SRC_A", 32896, VAR),
            ("C:X:Hd", WT, "P2", "Homeodomain", 0, "SRC_B", 32896, VAR),
            ("C:X:Hd", VAR, "P1", "Homeodomain", 1, "SRC_A", 32896, VAR),
            ("C:Y:Nac", LONE, "P3", "NAM", 0, "SRC_A", 32896, LONE),
        ],
        columns=[
            "wt_id",
            "dbd_seq",
            "protein_id",
            "dbd_family",
            "n_mut_from_wt",
            "source_dataset",
            "n_rows",
            "seed",
        ],
    )


def test_size_counts_domains_not_constructs(constructs):
    """One domain assayed by two labs is one domain. Counting rows or constructs inflates it."""
    inv = clusters.inventory(constructs).set_index("wt_id")
    assert inv.loc["C:X:Hd", "n_domains"] == 2
    assert inv.loc["C:X:Hd", "n_constructs"] == 3
    assert inv.loc["C:X:Hd", "n_rows"] == 3 * 32896
    assert inv.loc["C:X:Hd", "n_sources"] == 2
    assert inv.loc["C:Y:Nac", "n_domains"] == 1


def test_variant_and_edit_columns(constructs):
    inv = clusters.inventory(constructs).set_index("wt_id")
    assert inv.loc["C:X:Hd", "n_variants"] == 1
    assert inv.loc["C:X:Hd", "max_edits"] == 1
    assert inv.loc["C:Y:Nac", "n_variants"] == 0
    assert inv.loc["C:Y:Nac", "max_edits"] == 0


def test_reference_is_the_zero_edit_member_and_the_seed_is_kept_separately(constructs):
    inv = clusters.inventory(constructs).set_index("wt_id")
    assert inv.loc["C:X:Hd", "reference"] == WT
    assert inv.loc["C:X:Hd", "seed"] == VAR


def test_columns_and_ordering(constructs):
    inv = clusters.inventory(constructs)
    assert list(inv.columns) == list(clusters.INVENTORY_COLUMNS)
    assert inv["n_domains"].is_monotonic_decreasing


def test_selection_by_size_round_trips_through_a_file(constructs, tmp_path):
    path = tmp_path / "clusters.parquet"
    clusters.inventory(constructs).to_parquet(path, index=False)

    assert clusters.ids_with_at_least(2, path) == {"C:X:Hd"}
    assert clusters.ids_with_at_least(1, path) == {"C:X:Hd", "C:Y:Nac"}
    assert clusters.ids_with_at_least(3, path) == set()

    rows = pd.DataFrame({"wt_id": ["C:X:Hd", "C:Y:Nac", "C:X:Hd"], "label": [1, 0, 0]})
    assert list(clusters.select(rows, 2, path)["wt_id"]) == ["C:X:Hd", "C:X:Hd"]


def test_missing_inventory_says_what_to_run(tmp_path):
    with pytest.raises(FileNotFoundError, match="build_clusters.py"):
        clusters.load(tmp_path / "nope.parquet")
