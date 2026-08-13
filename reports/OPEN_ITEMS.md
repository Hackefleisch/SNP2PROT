# Open items

Decisions explicitly reserved for the project owner, plus questions raised by the work.
Flag here; do not resolve unilaterally.

## Reserved by the owner (brief §7)

| # | question | bearing evidence found so far | status |
|---|---|---|---|
| 1 | Keep the B1H data at all, given its selection-absence negatives? | — | open |
| 2 | Restrict the final set to DBDs with ≥5 variants (sharper clusters, much smaller set)? | — | open |
| 3 | Does the padded-20 bp or the common-core variant become primary? | — | open |

## Raised during the work

| # | item | phase | status |
|---|---|---|---|
| 4 | B1H binarization rule (recovered@10 mM → 1, not-recovered@2 mM → 0) is provisional and must be checked against Persikov et al.'s own treatment before Phase 3 closes. | 0 | open — source now available: `persikov2015_..._supp-methods.pdf` §2a covers processing and filtering of the protein selection data |
| 5 | SNP-SELEX OBS cutoffs are `null` in `configs/thresholds.yaml` — to be set from the observed score distribution in Phase 4, then confirmed by the owner. | 0 | open — **no narrative supplement exists** for Yan 2021 (only a HOMER dump PDF; the rest are `.xlsx`), so the cutoff must be derived from the main text plus the score distribution rather than quoted from a methods section |
| 6 | UniPROBE 8-mer score files may sit behind the same academic-use click-through as the probe sequences. If so this becomes a manual download and blocks Phase 1. | 0 | open — first thing to check in Phase 1 |
