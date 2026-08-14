# Reports

Generated markdown, **committed to git** — unlike `results/`, which is ignored. These are
the review artifacts handed back at each phase boundary, and their diffs are the record of
how the dataset changed.

| file | produced in | contents |
|---|---|---|
| `binarization_summary_<source>.md` | Phase 1-4 | n_pos / n_neg / n_gray, positive rate, score distribution vs. applied cutoffs |
| `cluster_inventory_<source>.md` | Phase 1-4 | variants per `wt_id`, `n_mut_from_wt` distribution, clusters per family, dominance flags |
| `cluster_inventory.md` | Phase 5 | the same, for the merged table |
| `overlap.md` | Phase 5 | cross-source duplicate pairs and label agreement |
| `validation_<source>.md` | Phase 1-4 | schema validator output |
| `archive/` | — | superseded hand-maintained files, frozen; see `TODO.md` for what replaced them |

Regenerate, never hand-edit. Open tasks and decisions live in [`TODO.md`](../TODO.md) at
the repository root; resolved ones in [`docs/DECISIONS.md`](../docs/DECISIONS.md).
