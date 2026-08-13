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

## Phase 2 findings

| # | item | phase | status |
|---|---|---|---|
| 12 | **Noyes et al. 2008 is not in UniPROBE.** The brief pairs it with Berger 2008 as a homeodomain panel, but Noyes used bacterial one-hybrid rather than PBM and appears under no UniPROBE accession. It belongs to the Phase 3 B1H work instead. Recorded in `PROVENANCE.md` under sources attempted. | 2 | **closed** |
| 13 | **The family panels the brief asked for are identified** (`docs/UNIPROBE_ACCESSIONS.md`): ETS = `EMBO10` (parsed), forkhead = `PNAS13` (parsed) and `ROG18A` (not parsed), bZIP = `MAR17A` and `SCI09` (neither parsed). Decide how much further to go — `SCI09` (104 mouse TFs, broadest family spread) and `GR09` (89 yeast TFs) are the obvious next two. | 2 | open |
| 14 | **11 pairs of Cell08 genes have byte-identical homeodomains** — Evx1/Evx2, Gbx1/Gbx2, Hoxa5/Hoxb5, Irx2/Irx5, Lhx1/Lhx5, Lhx2/Lhx9, Lmx1a/Lmx1b, Nkx2-4/Titf1, Phox2a/Phox2b, Pitx2/Pitx3, Vax1/Vax2. They are separate experiments on the same protein sequence, so they are reconciled as replicates and share one `wt_id` (`Cell08:Evx1/Evx2`). This is also the lookup-table risk in its purest form: two "different TFs" that a sequence model cannot possibly distinguish. | 2 | open, feeds the Phase 5 split design |
| 15 | **`EMBO10` and `PNAS13` detail pages carry no `DNA binding domain` field** (0/22 and 0/20), so `dbd_seq` falls back to the clone insert for those two panels while `Cell08` uses the Pfam-trimmed domain (168/168). `dbd_source` distinguishes them. Insert lengths are tight (EMBO10 103-143 aa, PNAS13 95-138 aa) with **one outlier: `Jumeau` at 719 aa**, effectively full-length. HMMER is not installed, so no HMM trimming was attempted. | 2 | open |
| 16 | **Three DBD conventions are now in play**: `uniprobe_dbd_field` (Cell08, 57-60 aa), `uniprobe_clone_insert` (BAR15A 86-282 aa, EMBO10, PNAS13). Extends open item #11 — this has to be reconciled before the Phase 5 merge, or `dbd_seq` length alone will correlate with source. | 2 | open, decide by Phase 5 |
| 17 | **The panels add breadth but no within-cluster depth.** Cell08 + EMBO10 + PNAS13 contribute ~199 clusters of exactly one DBD each, against BAR15A's 41 clusters averaging ~4 DBDs. Leave-one-variant-out can only ever be evaluated on BAR15A rows; leave-one-cluster-out is what the panels support. Worth deciding whether more point-mutant sources are needed. | 2 | open |
