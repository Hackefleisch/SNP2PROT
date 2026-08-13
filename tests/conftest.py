"""Shared fixtures: a minimal schema-conforming table to mutate in tests."""

from __future__ import annotations

import pandas as pd
import pytest

from snp2prot import schema

WT = "RKRGRQTYTRYQTLELEKEFHFNRYLTRRRRIEIAHALCLTERQIKIWFQNRRMKWKKEN"
MUT = WT[:10] + "A" + WT[11:]


@pytest.fixture
def good_frame() -> pd.DataFrame:
    """Two DBDs (a WT and one point mutant) x two 8-mers, with both label classes."""
    rows = []
    for dbd, nmut, muts in ((WT, 0, ""), (MUT, 1, "11")):
        for dna, score, label in (("TAATTAGC", 0.49, 1), ("GGGGCCCC", 0.11, 0)):
            rows.append(
                {
                    "dbd_seq": dbd,
                    "dbd_family": "Homeodomain",
                    "dbd_source": "cisbp_curated",
                    "wt_id": "TEST_WT",
                    "n_mut_from_wt": nmut,
                    "mut_positions": muts,
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
                    "source_dataset": "TESTSRC",
                    "source_file": "TESTSRC/fake.txt",
                }
            )
    # dna_len and pair_id are deliberately omitted: coerce() derives them, and the
    # fixture exercising that path is how we notice if the derivation regresses.
    return schema.coerce(pd.DataFrame(rows))
