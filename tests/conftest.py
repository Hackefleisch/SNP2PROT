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


# `scripts/` is not a package, and one of its modules carries logic worth testing directly:
# the co-binding overlap search. It is imported under `scripts_<name>` rather than moved into
# `snp2prot`, because it is a one-off analysis whose home is the script — what is pinned here is
# its arithmetic, not an API. (The filter/motif comparison DID move, to `snp2prot.ght.baselines`,
# once a second caller needed it.)
def _load_script(name: str):
    import importlib.util
    import sys

    from snp2prot.config import PROJECT_ROOT

    key = f"scripts_{name}"
    if key in sys.modules:
        return sys.modules[key]
    spec = importlib.util.spec_from_file_location(key, PROJECT_ROOT / "scripts" / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    sys.modules[key] = module
    spec.loader.exec_module(module)
    return module


_load_script("check_ght_cobinding")
