# Open items

Decisions reserved for the project owner, plus questions raised by the work. Flag here; do
not resolve unilaterally. Reconciled 2026-08-13 after the domain policy and Phase 2b — items
superseded by later work were closed rather than left to rot.

## A. Open decisions

| # | decision | bearing evidence | urgency |
|---|---|---|---|
| 2 | **Restrict the set to DBDs with ≥5 variants?** (brief §7) | Only **7 of 274 clusters** hold 5 or more domains, and 24 hold more than one. A ≥5 rule would cut the protein axis to a handful of clusters. | low — leaning no |
| 3 | **Does the padded-20 bp or the common-core dataset become primary?** (brief §7) | Currently moot: every admitted row is an 8-mer, so there is nothing to reconcile. Becomes live the moment Phase 4 (19 bp) is admitted. | **before Phase 4** |
| 5 | **SNP-SELEX OBS cutoff.** `snp_selex.positive`/`negative` are `null` in `configs/thresholds.yaml`. | No narrative supplement exists for Yan 2021 — only a HOMER dump PDF, everything else `.xlsx` — so the cutoff must come from the main text plus the observed score distribution, not a quoted methods value. | **before Phase 4** |
| 7 | **Keep 120 BAR15A variant alleles or the paper's 117?** | The three extra are HOXD13 `I297V`, `N298S`, `Q325K`, whose files are named `*_8mers_11111111.txt` — exactly what a naive glob drops, which is likely where 117 comes from. HOXD13 survives the policy with all 8 domains, so this affects one cluster. | low |
| 8 | **`PROP1_R112Q`: keep, correct, or drop?** | Its name says R→Q; the deposited sequence gives R→M. Kept with bookkeeping taken from the sequence. `PITX2_T114P` remains excluded (its insert shows no substitution at all). The other two anomalies vanished when their genes were rejected as mixed-family. | low |
| 22 | **Is 67% full-length coverage enough?** | The protein table maps **243 of 361 domains** onto a canonical UniProt sequence (338 have a sequence at all). The rest are clone constructs differing from the canonical isoform, or non-model species. Domain-level embeddings work for all; full-protein for two thirds. | **before embedding work** |
| 27 | ~~Homeodomain share~~ **RETIRED as framed — see #32.** | Row share equals domain share (every domain contributes exactly 32,896 rows), so the figure said nothing extra. The brief's 40% flag was written because "C2H2 will otherwise swamp everything"; C2H2 is now 4 domains (1.2%), so the threat it guarded against does not exist. Homeodomain is also *more* internally diverse than the minority families — median pairwise identity 34%, only 1.1% of pairs ≥90% identical, against Forkhead's 20.3%. Nothing is being swamped. | closed |
| 32 | **Three families contribute zero point variants, and 54 of 66 variants are homeodomain.** Forkhead has 6, zf-C4 5, PAX 1; ETS, HLH and zf-C2H2 have none. So the project's central question — does sensitivity to single-residue changes transfer to another fold? — cannot currently be answered. This is what #27 was obscuring, and it points at point-mutant series in a second family rather than at more single-protein panels. | **the real gap.** #30 (4 forkhead chimeras) is the only lead in hand. | **closed 2026-08-13** — owner: acceptable while point variants are not the majority of the data; subsets can be built later, and the question will resurface on its own when modelling starts. |
| 33 | **Four of seven families are too small to hold out**: PAX 3 domains, zf-C2H2 4, HLH 13. Leave-one-family-out on those tests a handful of proteins and returns noise. Only Homeodomain (223), Forkhead (46), zf-C4 (28) and ETS (23) support a meaningful family holdout. | Report LOFO only for the four viable families; treat the rest as descriptive. | **closed 2026-08-13** — owner: not worth acting on now; revisit when splits are designed. |
| 28 | **Parse `GB11` (Gordan 2011, 27 bHLH proteins)?** | Would roughly triple HLH domain count, helping #33 — but adds **zero point variants**, so it does nothing for #32, the real gap. Yields are low anyway: the policy admits about half of all constructs. | low, reassessed |
| 29 | **Drop `MAR17A:Esrrb`?** | It carries 2 unresolved `X` residues in an 89 aa zf-C4 domain. Harmless to a sequence embedder, a genuine problem for structure prediction and a 3D embedder. One domain of 340. | **before structure work** |
| 30 | **Cluster ROG18A's chimeras with their parents?** | 4 engineered forkhead chimeras sit alongside their 3 parents as 7 singleton clusters. They pass all three conditions, and `FoxN3_J3` differs from `FoxN3` by only 6 substitutions in 97 aa — real protein-axis depth, which the corpus is short of. Only that pair has matching lengths; the rest would need the cluster length-homogeneity rule relaxed. | medium |
| 14 | **How should identical-domain pairs be split?** | 11 Cell08 paralogue pairs have byte-identical homeodomains (Evx1/Evx2, Lhx2/Lhx9, …), reconciled into one `wt_id`. They are the lookup-table risk in pure form: two "different TFs" no sequence model can separate. | **before Phase 5 splits** |
| 23 | **Phases 4 and 6 need per-construct screening**, not a wholesale call. | SNP-SELEX varies DNA not protein, so condition 3 is vacuous there, but many of its 270 TFs are C2H2. Tier 4 are bHLH dimers, and MAX's substitutions are "in and around" the DBD. | at each phase |

## B. Decisions already taken

| # | decision | outcome |
|---|---|---|
| 9 | Negative:positive ratio (255:1) | **Leave the stored table intact**; control the ratio when sampling training batches. `pbm.positive` stays at 0.45. |
| 10 | Family imbalance | **Parse more panels** rather than down-sample. Homeodomain 79.5% → 64.9%. |
| 18 | `dbd_seq` definition | **Pfam envelope padded by 10 aa each side.** Documented in `configs/thresholds.yaml`, `CLAUDE.md`, `README.md`, `docs/DOMAIN_POLICY.md`. |
| 19–20 | C2H2 zinc-finger arrays | **Dropped**, in Phase 3 and in already-parsed sources alike. Reversible via `domain.allow_repeat_arrays`. |
| 24 | Protein-axis depth | **Accepted as-is** — 66 point variants in 24 clusters. Tier 4 stays held out; the Phase 5 NN baseline will show whether it suffices. |

## C. Closed — superseded or answered

| # | item | resolution |
|---|---|---|
| 1 | Keep B1H given weak negatives? | Moot — Phase 3 dropped on structural grounds regardless. |
| 4 | B1H binarization rule provisional | Moot — Phase 3 dropped. |
| 6 | Are UniPROBE 8-mer files licence-gated? | No. Direct download via `filehandler.php`; no click-through, cookie or MTA. |
| 11, 15, 16 | Competing `dbd_seq` conventions; EMBO10/PNAS13 lacking a DBD field; the 719 aa `Jumeau` outlier | All resolved by the domain policy. **Every one of the 340 domains is now `pfam_hmmer_padded`**, length 43–146 aa. Jumeau is 107 aa. `dbd_seq` length no longer correlates with source. |
| 12 | Noyes 2008 accession | Not in UniPROBE — it is B1H, not PBM. |
| 13 | Which panels supply forkhead/ETS/bZIP | Identified and all parsed: EMBO10, PNAS13, ROG18A, MAR17A, SCI09, GR09, SHO18A. |
| 17 | Panels add breadth, not depth | Confirmed and accepted — see #24. |
| 21 | Cost of the domain policy | Measured: 11.81M → 9.15M rows at the time; now 11.71M after Phase 2b. |
| 25, 26 | Duplicates of #8 and #7 | Merged into those rows. |
| 31 | SHO18A detail pages 500ing | Worked around via `detailsDef.php` with the site's own ids, verified against SHO18A/Apt. No ids guessed. |

## Phase 2c — full-Pfam scanning and the E-score column bug (2026-08-13)

| # | item | status |
|---|---|---|
| 34 | **The HMM set covered only 9 families; UniPROBE spans 86 domain labels.** Constructs from bZIP, HMG_box, GATA, T-box, Zn_clus, IRF and others were rejected as `no_domain` because the HMM was missing, not because nothing was there. Fixed by scanning the full Pfam-A library (30,134 families). The admission policy is now evaluated over a 41-family DNA-binding whitelist (`data/external/pfam/dbd_families.txt`), so a DBD beside an unrelated domain is still admitted. Effect: 340 -> 468 domains, 7 -> 31 families. | fixed |
| 35 | **UniPROBE labels predate several Pfam renames** — `Homeobox`→`Homeodomain`, `Fork_head`→`Forkhead`, `MADS`→`SRF-TF`, `E2F_TDP`→`WHD_E2F_TDP`, `Zn2Cys6`→`Zn_clus`, `BRIGHT`→`ARID`, `NR`→`zf-C4`, `AP-2`→`TF_AP-2`, `PBX`→`PBC`, `RFX`→`RFX_DNA_binding`. Each resolved by Pfam accession lookup, not string guessing. Missing `Homeobox` alone would have excluded the largest family and rejected almost everything. | fixed |
| 36 | **The E-score column was read by position and was wrong for three accessions.** `usecols=[0, 2]` holds for BAR15A/EMBO10/PNAS13/Cell08 but not for GR09 and SCI09 (headerless, median intensity at index 2, E-score at index 3) or RAD13A (`8-mer / 8-mer / Median / Z-score` — no E-score at all). The corpus briefly showed a 4:1 negative:positive ratio with `raw_score` up to 776,106. Fixed: the column is now identified by a header name where present, else by being bounded to [-0.5, 0.5] **and** taking negative values (p/q-values cannot). Files with no such column, or several, are rejected with a reason. | fixed |
| 37 | **A validator guard now makes this class of bug impossible to miss**: any `pbm_escore` row outside [-0.5, 0.5] is a hard schema error. | fixed |
| 38 | **13 of 30 registered accessions yield nothing.** Eleven are rejected entirely by the domain policy (C2H2 arrays, multi-domain, or families outside the DBD whitelist); `Path10` files carry no E-score column; `STI21B` likewise yields no admissible construct. All reported per source, none silent. | informational |
