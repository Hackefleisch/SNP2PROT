#!/usr/bin/env python
"""Recompute E-scores from raw scans and score them against a parsed source.

    python scripts/validate_escores.py --manifest <tsv> [--source BAR15A] [--out <tsv>]

The manifest is one row per array, tab-separated, with a header:

    gene    array   protein_gpr             cy3_gpr
    ARX_REF 159     .../..._Alexa488_1-8.gpr.gz  .../..._Cy3_1-8.gpr.gz

`gene` must match the `gene` column of the protein table for `--source`, which is how the
recomputed 8-mers find the stored ones to be judged against. `cy3_gpr` may be empty; pair the
control to the protein scan by slide and chamber, never by date, since the double-stranding
scan is often days earlier.

What the columns mean (`reports/escore_recomputation.md` for the interpretation):

  rho       Spearman over all 32,896 8-mers
  rho_bind  Spearman over the 8-mers the stored source scores above 0.3 -- the part a label
            depends on, and the number that says whether the ranking was reproduced
  ranked_J  Jaccard against the stored positives, taking as many top 8-mers as the stored
            source called positive. Compare against 0.70-0.72, the agreement of two
            replicates within one source (`reports/overlap.md`)
  n_at_45   how many 8-mers the recomputation puts at or above the corpus's own `pbm.positive`
            cutoff, against `n_stored`. These do not match, and that is the finding
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from snp2prot import rawpbm, thresholds
from snp2prot.config import PROTEIN_TABLE, interim_table


def spearman(a: pd.Series, b: pd.Series) -> float:
    return float(np.corrcoef(a.rank(), b.rank())[0, 1])


def stored_escores(source: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    proteins = pd.read_parquet(PROTEIN_TABLE, columns=["gene", "dbd_seq", "source_dataset"])
    proteins = proteins[proteins.source_dataset == source]
    rows = pd.read_parquet(interim_table(source), columns=["dbd_seq", "dna_seq", "raw_score"])
    return proteins, rows


def main() -> None:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest", required=True, type=Path)
    ap.add_argument("--source", default="BAR15A")
    ap.add_argument("--out", type=Path)
    args = ap.parse_args()

    positive = float(thresholds.for_assay("pbm")["positive"])
    manifest = pd.read_csv(args.manifest, sep="\t").fillna({"cy3_gpr": ""})
    proteins, rows = stored_escores(args.source)
    design = rawpbm.load_design()

    out = []
    for gene, group in manifest.groupby("gene", sort=False):
        match = proteins[proteins.gene == gene]
        if match.empty:
            print(f"  {gene}: no domain under source {args.source}, skipped")
            continue
        stored = rows[rows.dbd_seq == match.iloc[0].dbd_seq].set_index("dna_seq").raw_score
        arrays = [(r.protein_gpr, r.cy3_gpr or None) for r in group.itertuples()]
        recomputed, widths = rawpbm.score_allele(arrays, design)

        both = pd.concat([stored.rename("stored"), recomputed.rename("new")], axis=1).dropna()
        stored_pos = set(both[both.stored >= positive].index)
        ranked = set(both.nlargest(len(stored_pos), "new").index) if stored_pos else set()
        binding = both[both.stored > 0.3]
        out.append(
            {
                "gene": gene,
                "arrays": len(arrays),
                "tail_width": ", ".join(f"{w:.1f}" for w in widths),
                "rho": round(spearman(both.stored, both.new), 3),
                "rho_bind": round(spearman(binding.stored, binding.new), 3)
                if len(binding) > 10
                else None,
                "ranked_J": round(len(stored_pos & ranked) / len(stored_pos | ranked), 3)
                if stored_pos
                else None,
                "max": round(both.new.max(), 3),
                "n_at_45": int((both.new >= positive).sum()),
                "n_stored": len(stored_pos),
            }
        )

    table = pd.DataFrame(out)
    print(table.to_string(index=False))
    if args.out:
        table.to_csv(args.out, sep="\t", index=False)
        print(f"\nwrote {args.out}")


if __name__ == "__main__":
    main()
