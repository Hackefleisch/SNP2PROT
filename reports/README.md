# Reports

Generated markdown, **committed to git** — unlike `results/`, which is ignored. These are
the review artifacts handed back at each phase boundary, and their diffs are the record of
how the dataset changed.

| file | produced in | contents |
|---|---|---|
| `binarization_summary_<source>.md` | Phase 1-4 | n_pos / n_neg / n_gray, positive rate, score distribution vs. applied cutoffs |
| `cluster_inventory_<source>.md` | Phase 1-4 | variants per `wt_id`, `n_mut_from_wt` distribution, clusters per family, dominance flags |
| `cluster_inventory.md` | Phase 5 | the same, for the merged table |
| `overlap.md` | 2026-08-17 | the 47 domains stored by two sources, and how far their labels agree — the label noise floor |
| `clusters.md` | 2026-08-17 | corpus-wide cluster inventory from `build_clusters.py` |
| `validation_<source>.md` | Phase 1-4 | schema validator output |
| `nn_baseline.md` | Phase 7 | the bar: nearest-neighbour lookup per fold, as a ranking **and** as a decision |
| `training.md` | Phase 7 | the 19-fold x 2-arm grid, against that bar |
| `calibration.md` | 2026-08-28 | the `λ` sweep on `L = L_infonce + λ · L_bce` — does the model have a decision rule (`T38`) |
| `pooling_check.md` | Phase 7 | pre-flight: is a variant separated from its wild type in the pooled embedding at all |
| `unpooling_check.md` | 2026-09-08 | pre-flight: does the single-residue signal survive *before* the mean, and does it beat BLOSUM62 (`T36`) |
| `archive/` | — | superseded hand-maintained files, frozen; see `TODO.md` for what replaced them |

Regenerate, never hand-edit. Open tasks and decisions live in [`TODO.md`](../TODO.md) at
the repository root; resolved ones in [`docs/DECISIONS.md`](../docs/DECISIONS.md).
