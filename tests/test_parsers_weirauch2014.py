"""CIS-BP / Weirauch 2014.

The join between the measurements and the assayed sequence is what makes this source
admissible at all, so the tests concentrate on it: the plasmid ID in a GEO sample title, the
named columns of Table S6, and the cluster keys that must not fuse two orthologues.
"""

from __future__ import annotations

import gzip

import pandas as pd
import pytest

from snp2prot.parsers import weirauch2014 as w

SAMPLE = """^SAMPLE = GSM1291226
!Sample_title = pTH1294_HK_8mer_593
!sample_table_begin
ID_REF\tVALUE\tE-Score\tZ-Score
AAAAAAAA\t40152.5\t0.402\t3.30
AAAAAAAC\t32710.8\t-0.246\t1.64
!sample_table_end
^SAMPLE = GSM1291227
!Sample_title = pTH1294_ME_8mer_462
!sample_table_begin
ID_REF\tVALUE\tE-Score\tZ-Score
AAAAAAAA\t9999.0\t0.410\t3.10
AAAAAAAC\t8888.0\t-0.300\t1.20
!sample_table_end
"""


def _soft(tmp_path, text=SAMPLE):
    p = tmp_path / "family.soft.gz"
    with gzip.open(p, "wt") as fh:
        fh.write(text)
    return p


def test_title_yields_plasmid_and_array():
    """`pTH1294_HK_8mer_593` — the plasmid is the join key, the array marks a replicate."""
    m = w.TITLE_RE.match("pTH1294_HK_8mer_593")
    assert (m.group("plasmid"), m.group("array")) == ("pTH1294", "HK")
    assert w.TITLE_RE.match("pTH5539_ME_8mer_7668").group("array") == "ME"
    assert w.TITLE_RE.match("something_else") is None


def test_both_array_designs_are_read_as_separate_replicates(tmp_path):
    got = list(w.iter_sample_tables(_soft(tmp_path), {"pTH1294"}, require_complete=False))
    assert [(p, a) for p, a, _ in got] == [("pTH1294", "HK"), ("pTH1294", "ME")]
    assert got[0][2]["dna_seq"].tolist() == ["AAAAAAAA", "AAAAAAAC"]
    assert got[0][2]["escore"].tolist() == [0.402, -0.246]


def test_an_incomplete_sample_table_is_rejected(tmp_path):
    """Every protein must be scored against all 32,896 8-mers or its DNA axis differs."""
    with pytest.raises(w._pbm.EscoreColumnError, match="expected 32,896"):
        list(w.iter_sample_tables(_soft(tmp_path), {"pTH1294"}))


def test_blocks_for_unwanted_plasmids_are_skipped(tmp_path):
    """The record is ~1.7 GB; a plasmid the policy rejected must not be parsed at all."""
    assert list(w.iter_sample_tables(_soft(tmp_path), set(), require_complete=False)) == []


def test_table_s6_columns_are_required_by_name(tmp_path):
    """Reading a column by position is what made the E-score wrong for three accessions."""
    path = tmp_path / "clones.xlsx"
    pd.DataFrame({"Plasmid ID": ["pTH1"], "Insert AA": ["ACDEF"]}).to_excel(
        path, sheet_name=w.CLONES_SHEET, index=False
    )
    with pytest.raises(ValueError, match="missing column"):
        w.load_constructs(path)


def test_constructs_load_with_their_architecture(tmp_path):
    path = tmp_path / "clones.xlsx"
    pd.DataFrame(
        {
            "Plasmid ID": ["pTH1", "notaplasmid"],
            "Gene Name": ["Foxc1", "x"],
            "Species": ["Mus musculus", "x"],
            "Insert AA": ["acdef", "ACDEF"],
            "#Flanking AAs": [15, 0],
            "Gene ID": ["ENSMUSG00000050295", "x"],
        }
    ).to_excel(path, sheet_name=w.CLONES_SHEET, index=False)
    got = w.load_constructs(path)
    assert set(got) == {"pTH1"}
    assert got["pTH1"].insert_aa == "ACDEF"  # upper-cased
    assert got["pTH1"].flanking_aa == 15
    assert got["pTH1"].gene_id == "ENSMUSG00000050295"


def test_colliding_gene_names_do_not_fuse_into_one_cluster():
    """25 gene names here name a different sequence in a different organism.

    `FOXG1` is both human and medaka; `MATALPHA2` is three different yeasts. Keying a cluster
    on the gene alone would silently merge orthologues into one `wt_id`.
    """
    groups = {
        "AAAA": {"genes": {"FOXG1"}, "plasmids": ["pTH2"]},
        "CCCC": {"genes": {"FOXG1"}, "plasmids": ["pTH9"]},
        "GGGG": {"genes": {"Foxc1"}, "plasmids": ["pTH3"]},
    }
    got = w._assign_wt_ids(groups)
    assert got["GGGG"] == "weirauch2014:Foxc1"  # unique, left clean
    assert got["AAAA"] != got["CCCC"]
    assert got["AAAA"] == "weirauch2014:FOXG1_pTH2"
    assert got["CCCC"] == "weirauch2014:FOXG1_pTH9"
