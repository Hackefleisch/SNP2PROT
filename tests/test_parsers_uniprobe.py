"""Shared UniPROBE machinery.

The header-detection test is the important one: EMBO10's 8-mer files have no header row, so
a fixed `skiprows=1` silently drops the first 8-mer from every experiment in that panel — a
loss that leaves no trace beyond a row count nothing else checks.
"""

from __future__ import annotations

import io
import zipfile

import numpy as np
import pytest

from snp2prot import schema
from snp2prot.parsers import _uniprobe

BODY = "AAAAAAAA\tTTTTTTTT\t0.49\t100.0\t3.0\nAAAAAAAC\tGTTTTTTT\t0.10\t50.0\t1.0\n"
HEADED = "8-mer\t8-mer\tE-score\tMedian\tZ-score\n" + BODY
# Cell08 layout: seven columns, descriptive header.
WIDE = (
    "8-mer\t8-mer reverse complement\tenrichment score\tmedian intensity\tzscore"
    "\tpvalue of the enrichment\tQvalue\n"
    "AAAAAAAA\tTTTTTTTT\t0.49\t100.0\t3.0\t0.001\t0.01\n"
    "AAAAAAAC\tGTTTTTTT\t0.10\t50.0\t1.0\t0.5\t0.6\n"
)
BARE = BODY


def _zip(members: dict[str, str]) -> zipfile.ZipFile:
    buf = io.BytesIO()
    with zipfile.ZipFile(buf, "w") as z:
        for name, body in members.items():
            z.writestr(name, body)
    buf.seek(0)
    return zipfile.ZipFile(buf)


@pytest.mark.parametrize(
    "body", [HEADED, BARE, WIDE], ids=["5col-header", "5col-no-header", "7col-header"]
)
def test_header_is_detected_not_assumed(body):
    z = _zip({"G/G_8mers.txt": body})
    df = _uniprobe.read_8mer_table(z, "G/G_8mers.txt")
    assert len(df) == 2, "a headerless file must not lose its first data row"
    assert df["dna_seq"].tolist() == ["AAAAAAAA", "AAAAAAAC"]
    assert df["escore"].tolist() == [0.49, 0.10]


@pytest.mark.parametrize(
    "name",
    [
        "G/G_REF/G_REF_R1/G_REF_R1_8mers.txt",
        "G/G_REF/G_REF_R2/G_REF_R2_8mers_11111111.txt",
        "Alx3/3418.2/Alx3_3418.2_contig8mers.txt",
    ],
)
def test_contig_regex_matches_every_known_spelling(name):
    assert _uniprobe.CONTIG_8MER_RE.search(name)


@pytest.mark.parametrize(
    "name", ["G/G_pwm.txt", "G/G_8mers_top_enrichment.txt", "G/G_8mers_11101111.txt"]
)
def test_contig_regex_rejects_other_files(name):
    assert not _uniprobe.CONTIG_8MER_RE.search(name)


def test_binarize_bands():
    e = np.array([0.50, 0.45, 0.40, 0.35, -0.5])
    lab = _uniprobe.binarize(e, 0.45, 0.35)
    assert lab.tolist() == [
        schema.LABEL_BIND,
        schema.LABEL_BIND,
        schema.LABEL_GRAY,
        schema.LABEL_NONBIND,
        schema.LABEL_NONBIND,
    ]


def test_replicate_disagreement_falls_to_gray():
    """Replicates may sit on different array designs; disagreement is not averaged away."""
    a = "8-mer\t8-mer\tE-score\tM\tZ\nAAAAAAAA\tT\t0.49\t1\t1\nAAAAAAAC\tG\t0.49\t1\t1\n"
    b = "8-mer\t8-mer\tE-score\tM\tZ\nAAAAAAAA\tT\t0.48\t1\t1\nAAAAAAAC\tG\t0.10\t1\t1\n"
    z = _zip({"G/R1_8mers.txt": a, "G/R2_8mers.txt": b})
    keys, label, mean = _uniprobe.reconcile_replicates(
        z, ["G/R1_8mers.txt", "G/R2_8mers.txt"], 0.45, 0.35
    )
    assert keys.tolist() == ["AAAAAAAA", "AAAAAAAC"]
    assert label.tolist() == [schema.LABEL_BIND, schema.LABEL_GRAY]
    assert mean[1] == pytest.approx(0.295)


def test_replicates_must_cover_the_same_8mers():
    a = "8-mer\t8-mer\tE-score\tM\tZ\nAAAAAAAA\tT\t0.49\t1\t1\n"
    b = "8-mer\t8-mer\tE-score\tM\tZ\nAAAAAAAC\tG\t0.49\t1\t1\n"
    z = _zip({"G/R1_8mers.txt": a, "G/R2_8mers.txt": b})
    with pytest.raises(ValueError, match="differs from the first replicate"):
        _uniprobe.reconcile_replicates(z, ["G/R1_8mers.txt", "G/R2_8mers.txt"], 0.45, 0.35)


def test_clean_sequence_strips_numbering_and_markup():
    got = _uniprobe.clean_sequence("1 &nbsp;&nbsp;AGSDSEEGLL KRKQ<br>51 &nbsp;RRAKWRK")
    assert got == "AGSDSEEGLLKRKQRRAKWRK"


def test_detail_page_accepts_both_uniprot_labels():
    """Some pages label it `Swiss-Prot`, others `Uniprot`."""
    a = _uniprobe.parse_detail_page("<dt>Swiss-Prot</dt><dd>Q96QS3</dd>", "G")
    b = _uniprobe.parse_detail_page("<dt>Uniprot</dt><dd>P12345</dd>", "G")
    assert a.protein_id == "Q96QS3"
    assert b.protein_id == "P12345"


def test_detail_page_accepts_both_insert_label_styles():
    allele_style = _uniprobe.parse_detail_page(
        "<dt>ARX_L343Q Insert Sequence</dt><dd>\n<kbd>1 AAAA</kbd></dd>", "ARX"
    )
    clone_style = _uniprobe.parse_detail_page(
        "<dt>Clone pTH3418 insert sequence</dt><dd>\n<kbd>1 CCCC</kbd></dd>", "Alx3"
    )
    assert allele_style.inserts == {"L343Q": "AAAA"}
    assert clone_style.inserts == {"REF": "CCCC"}
