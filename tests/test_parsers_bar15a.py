"""BAR15A parser logic.

Tests the pure helpers only — the real archive is 100 MB and git-ignored, so the suite must
pass on a clean checkout. What is pinned here is the reasoning that would otherwise be easy
to regress: the two filename spellings, and mutation bookkeeping taken from sequences rather
than allele names.
"""

from __future__ import annotations

import pytest

from snp2prot.parsers import bar15a

PAGE = """
<dl><dt>Species</dt><dd>Homo sapiens</dd></dl>
<dl><dt>Domain</dt><dd>Homeobox</dd></dl>
<dl><dt>Swiss-Prot</dt><dd><a href="#">Q96QS3</a></dd></dl>
<dl><dt>ARX DNA binding domain</dt><dd><kbd>1 &nbsp;&nbsp;RRYRTTFTSY</kbd></dd></dl>
<dl><dt>ARX_REF Insert Sequence</dt>
<dd>
<kbd>1 &nbsp;&nbsp;AGSDSEEGLL KRKQ<br>51 &nbsp;RRAKWRK</kbd>
</dd></dl>
<dl><dt>ARX_L343Q Insert Sequence</dt>
<dd>
<kbd>1 &nbsp;&nbsp;AGSDSEEGLL KRKQ<br>51 &nbsp;RRAKWRQ</kbd>
</dd></dl>
"""


def test_contig_8mer_regex_matches_both_spellings():
    """HOXD13's replicates use the second spelling; missing it drops three alleles."""
    assert bar15a.CONTIG_8MER_RE.search("A/A_REF/A_REF_R1/A_REF_R1_8mers.txt")
    assert bar15a.CONTIG_8MER_RE.search("H/H_REF/H_REF_R2/H_REF_R2_8mers_11111111.txt")


def test_contig_8mer_regex_rejects_other_kmer_files():
    assert not bar15a.CONTIG_8MER_RE.search("A/A_REF/A_REF_R1/A_REF_R1_pwm.txt")
    assert not bar15a.CONTIG_8MER_RE.search("A/A_REF/A_REF_R1/A_REF_R1_8mers_11101111.txt")


def test_clean_sequence_strips_numbering_and_markup():
    got = bar15a._clean_sequence("1 &nbsp;&nbsp;AGSDSEEGLL KRKQ<br>51 &nbsp;RRAKWRK")
    assert got == "AGSDSEEGLLKRKQRRAKWRK"


def test_load_metadata_reads_inserts_and_fields(tmp_path):
    (tmp_path / "ARX.html").write_text(PAGE)
    meta = bar15a.load_metadata(tmp_path)["ARX"]
    assert meta.family == "Homeobox"
    assert meta.protein_id == "Q96QS3"
    assert meta.species == "Homo sapiens"
    assert set(meta.inserts) == {"REF", "L343Q"}
    assert meta.ref.startswith("AGSDSEEGLL")


def test_load_metadata_requires_a_reference_allele(tmp_path):
    (tmp_path / "X.html").write_text(
        "<dl><dt>X_A1G Insert Sequence</dt><dd>\n<kbd>1 AAA</kbd></dd></dl>"
    )
    with pytest.raises(ValueError, match="no REF"):
        bar15a.load_metadata(tmp_path)


def test_mutation_bookkeeping_is_one_based():
    assert bar15a._mutation_bookkeeping("ACDEF", "ACDEF") == (0, "")
    assert bar15a._mutation_bookkeeping("ACDEF", "AGDEF") == (1, "2")
    assert bar15a._mutation_bookkeeping("ACDEF", "AGDEY") == (2, "2,5")


def test_mutation_bookkeeping_rejects_length_change():
    with pytest.raises(ValueError, match="length differs"):
        bar15a._mutation_bookkeeping("ACDEF", "ACDE")


def test_anomalous_alleles_are_documented_and_only_one_is_dropped():
    """The four name/sequence disagreements found in Phase 1 stay explicitly listed."""
    assert set(bar15a.ANOMALIES) == {
        "PITX2_T114P",
        "POU3F4_A237G",
        "PAX6_R26G",
        "PROP1_R112Q",
    }
    # Only the one with no difference from REF is unusable: it would otherwise duplicate the
    # reference allele's rows on pair_id.
    assert bar15a.UNUSABLE == {"PITX2_T114P"}
    assert bar15a.UNUSABLE <= set(bar15a.ANOMALIES)
