"""`reports.overlap_report` — the cross-source label agreement check."""

from __future__ import annotations

import pandas as pd

from snp2prot import reports, schema

WT = "RKRGRQTYTRYQTLELEKEFHFNRYLTRRRRIEIAHALCLTERQIKIWFQNRRMKWKKEN"
OTHER = WT[:20] + "A" + WT[21:]

#: Enough 8-mers to make a Jaccard mean something without writing out 32,896.
KMERS = ["TAATTAGC", "TAATTAGG", "GGGGCCCC", "ACGTACGT", "TTTTAAAA"]


def _frame(source: str, dbd: str, labels: list[int], scores: list[float]) -> pd.DataFrame:
    rows = [
        {
            "dbd_seq": dbd,
            "dbd_family": "Homeodomain",
            "dbd_source": "pfam_hmmer_padded",
            "wt_id": "C:TEST:Hd",
            "n_mut_from_wt": 0,
            "mut_positions": "",
            "protein_id": "P00001",
            "species": "Mus musculus",
            "dna_seq": dna,
            "dna_context": "core_only",
            "label": label,
            "raw_score": score,
            "score_type": "pbm_escore",
            "threshold_pos": 0.45,
            "threshold_neg": 0.35,
            "assay": "PBM",
            "stringency": "",
            "neg_provenance": "assayed_unbound" if label == 0 else None,
            "source_dataset": source,
            "source_file": f"{source}/fake.txt",
        }
        for dna, label, score in zip(KMERS, labels, scores, strict=True)
    ]
    return schema.coerce(pd.DataFrame(rows))


def test_no_shared_domain_reports_nothing():
    frames = {
        "A": _frame("A", WT, [1, 1, 0, 0, 0], [0.49, 0.47, 0.1, 0.1, 0.1]),
        "B": _frame("B", OTHER, [1, 1, 0, 0, 0], [0.49, 0.47, 0.1, 0.1, 0.1]),
    }
    assert "No domain is measured by more than one source." in reports.overlap_report(frames)


def test_perfect_agreement():
    labels, scores = [1, 1, 0, 0, 0], [0.49, 0.47, 0.1, 0.2, 0.3]
    frames = {"A": _frame("A", WT, labels, scores), "B": _frame("B", WT, labels, scores)}
    out = reports.overlap_report(frames)
    assert "**1.000**" in out  # median positive-call Jaccard
    assert "None." in out  # nothing flagged


def test_disagreement_is_flagged_and_the_gray_zone_is_not_a_disagreement():
    # B calls a different 8-mer positive, so the positive sets are disjoint. B's third
    # 8-mer sits between the cutoffs (label -1) where A called a negative: no call was
    # made, so that must not count against the hard-label agreement.
    frames = {
        "A": _frame("A", WT, [1, 0, 0, 0, 0], [0.49, 0.10, 0.10, 0.20, 0.30]),
        "B": _frame("B", WT, [0, 1, -1, 0, 0], [0.10, 0.49, 0.40, 0.20, 0.30]),
    }
    out = reports.overlap_report(frames)
    assert "0.000" in out
    assert "C:TEST:Hd" in out
    assert "Flagged replicates" in out
    # Four hard comparisons survive; two of them (the swapped positives) differ.
    assert "2 of 4 calls differ" in out
