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
| 6 | UniPROBE 8-mer score files may sit behind the same academic-use click-through as the probe sequences. | 0 | **closed 2026-08-13** — not gated. `downloads/BAR15A/*` 302-redirects to `filehandler.php`, which serves the file directly. No licence page, no cookie, no MTA. |

## Phase 1 findings needing a decision

| # | item | phase | status |
|---|---|---|---|
| 7 | **UniPROBE deposits 120 variant alleles; the paper says 117.** The three extra are HOXD13 I297V, N298S and Q325K, whose only files are named `*_8mers_11111111.txt` — precisely the ones a `*_8mers.txt` glob drops, which is likely how the 117 figure arises downstream. All three are parsed and included. Decide whether to keep them or restrict to the paper's 117. | 1 | open |
| 8 | **Four alleles' deposited insert sequences contradict their own names.** `PITX2_T114P` shows no substitution at all and is **excluded** (identical to REF, it would collide with the reference on `pair_id`). `POU3F4_A237G` carries two substitutions — the intended one plus A312V, POU3F4's other deposited allele — and is kept with `n_mut_from_wt=2`. `PAX6_R26G` actually changes Q27→G (off by one). `PROP1_R112Q` actually gives R112M. The latter three are kept with mutation bookkeeping taken from the sequence, never the name. Decide whether to drop them, correct them against the paper, or keep as-is. | 1 | open |
| 9 | **Negative:positive ratio is 383:1**, well above the 50-100:1 the brief expects, because E-score ≥ 0.45 is strict on a 32,896-8-mer space. Lowering `pbm.positive` would move it. Only 0.26% of rows are positive. Decide whether to relax the cutoff or down-sample negatives at training time. | 1 | open |
| 10 | **`Homeobox` is 46.2% of BAR15A rows**, over the 40% dominance flag, though this is a single-source figure and will shift once C2H2 B1H lands in Phase 3. | 1 | open, revisit at Phase 5 |
| 11 | **`dbd_seq` is the deposited clone insert, not a Pfam-trimmed domain.** Only 20 of 41 detail pages carry a separate DBD field, so the insert is the only representation available for every allele — and it is what was physically assayed. Insert lengths vary by gene (86-282 aa). This makes `dbd_source` = `uniprobe_clone_insert` and will need reconciling against HMM-derived boundaries from other sources at merge time. | 1 | open, decide by Phase 5 |
