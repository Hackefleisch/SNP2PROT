# Open items

Decisions reserved for the project owner, plus questions raised by the work. Flag here; do
not resolve unilaterally. Reconciled 2026-08-13 after the domain policy and Phase 2b — items
superseded by later work were closed rather than left to rot.

## A. Open decisions

| # | decision | bearing evidence | urgency |
|---|---|---|---|
| 2 | **Restrict the set to DBDs with ≥5 variants?** (brief §7) | Only **7 of 402 clusters** hold 5 or more domains, and 24 hold more than one. A ≥5 rule would cut the protein axis to a handful of clusters. | low — leaning no |
| 3 | **Does the padded-20 bp or the common-core dataset become primary?** (brief §7) | Currently moot: every admitted row is an 8-mer, so there is nothing to reconcile. Becomes live the moment Phase 4 (19 bp) is admitted. | **before Phase 4** |
| 5 | **SNP-SELEX OBS cutoff.** `snp_selex.positive`/`negative` are `null` in `configs/thresholds.yaml`. | No narrative supplement exists for Yan 2021 — only a HOMER dump PDF, everything else `.xlsx` — so the cutoff must come from the main text plus the observed score distribution, not a quoted methods value. | **before Phase 4** |
| 7 | ~~Keep 120 BAR15A variant alleles or the paper's 117?~~ **DECIDED 2026-08-13: keep all 120 — trust what is on disk.** | The three extra are HOXD13 `I297V`, `N298S`, `Q325K`, whose files are named `*_8mers_11111111.txt` — exactly what a naive glob drops, which is likely where 117 comes from. HOXD13 survives the policy with all 8 domains. | **closed** |
| 8 | **`PROP1_R112Q`: keep, correct, or drop?** | Its name says R→Q; the deposited sequence gives R→M. Kept with bookkeeping taken from the sequence. `PITX2_T114P` remains excluded (its insert shows no substitution at all). The other two anomalies vanished when their genes were rejected as mixed-family. | low |
| 22 | **Is 67% full-length coverage enough?** | The protein table maps **342 of 478 domains** onto a canonical UniProt sequence (445 have a sequence at all). The rest are clone constructs differing from the canonical isoform, or non-model species. Domain-level embeddings work for all; full-protein for two thirds. | **before embedding work** |
| 27 | ~~Homeodomain share~~ **RETIRED as framed — see #32.** | Row share equals domain share (every domain contributes exactly 32,896 rows), so the figure said nothing extra. The brief's 40% flag was written because "C2H2 will otherwise swamp everything"; C2H2 is now 4 domains (1.2%), so the threat it guarded against does not exist. Homeodomain is also *more* internally diverse than the minority families — median pairwise identity 34%, only 1.1% of pairs ≥90% identical, against Forkhead's 20.3%. Nothing is being swamped. | closed |
| 32 | **Three families contribute zero point variants, and 54 of 66 variants are homeodomain.** Forkhead has 6, zf-C4 5, PAX 1; ETS, HLH and zf-C2H2 have none. So the project's central question — does sensitivity to single-residue changes transfer to another fold? — cannot currently be answered. This is what #27 was obscuring, and it points at point-mutant series in a second family rather than at more single-protein panels. | **the real gap.** #30 (4 forkhead chimeras) is the only lead in hand. | **closed 2026-08-13** — owner: acceptable while point variants are not the majority of the data; subsets can be built later, and the question will resurface on its own when modelling starts. |
| 33 | **Four of seven families are too small to hold out**: PAX 3 domains, zf-C2H2 4, HLH 13. Leave-one-family-out on those tests a handful of proteins and returns noise. Reassessed after the full-Pfam fix: 31 families now, with Homeodomain 234, Forkhead 46, HLH 32, zf-C4 28, ETS 22, HMG_box 21, Zn_clus 17, T-box 11 all viable; the long tail below ~10 domains is not. | Report LOFO only for the four viable families; treat the rest as descriptive. | **closed 2026-08-13** — owner: not worth acting on now; revisit when splits are designed. |
| 28 | ~~Parse `GB11`?~~ **Moot — GB11's detail pages carry no sequence at all**, so it is unparseable regardless. | Superseded: HLH now stands at **32 domains** (was 13) after the full-Pfam fix admitted families the old HMM set had missed, so the gap this was meant to close has largely closed itself. Would still add zero point variants. | **closed** |
| 29 | **Drop `MAR17A:Esrrb`?** | It carries 2 unresolved `X` residues in an 89 aa zf-C4 domain — the only such domain in the corpus. Harmless to a sequence embedder, a genuine problem for structure prediction and a 3D embedder. One domain of 468. | **before structure work** |
| 30 | ~~Cluster ROG18A's chimeras with their parents?~~ **DECIDED 2026-08-13: yes, 6 substitutions is close enough.** Implemented, and it uncovered a larger bug — see #39. | 4 engineered forkhead chimeras sit alongside their parents as singleton clusters (6 domains total after dedup). They pass all three conditions, and `FoxN3_J3` differs from `FoxN3` by only 6 substitutions in 97 aa — real protein-axis depth, which the corpus is short of. `ROG18A:FoxN3` now clusters its reference with 3 chimeras at 5, 6 and 8 substitutions. The FoxJ3 chimeras stay singletons: their padded envelopes differ in length from the parent, so Hamming distance is undefined and the length-homogeneity rule would have to be broken to force them together. | **closed** |
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

## Phase 2d — one gene folder is not one protein (2026-08-13)

| # | item | status |
|---|---|---|
| 39 | **The panel parser merged distinct engineered proteins as replicates.** UniPROBE archives encode the construct at path level 2, and detail pages can carry several insert sequences, but the parser grouped by gene and took the first insert. So ROG18A's `FoxJ3_N3/` — six different chimeras — became one protein; NAR11's `HLH-1/` — a point-mutant series L13R/L13T/L13V — became one; LIN14B's `ANAC092_DBD` and `ANAC092_FL` were averaged together; GD13's two CLAMP constructs likewise. Fixed: inserts are keyed by their full construct name and experiments are split against those names. **ROG18A 7 -> 15 domains, NAR11 1 -> 4.** | fixed |
| 40 | **Two non-homeodomain variant clusters emerged from that fix** — the corpus's first: `NAR11:HLH-1` (bHLH, 3 point variants at one position) and `ROG18A:FoxN3` (forkhead, 3 chimeras at 5-8 substitutions). Directly relevant to the concern retired as #32. | informational |
| 41 | **`PROP1_R112Q`, resolved.** The Barrera supplement lists PROP1 **R99Q** among its variants with a reported effect; **R112Q is not in that table**. The deposited R99Q sequence checks out exactly (construct 44 -> protein 99, R->Q). R112Q has the right position and the right wild-type residue but substitutes **M, not Q**. Both R->Q and R->M are single-nucleotide changes from Arg, so neither is implausible, and nothing available decides which the clone really was. Kept as deposited, consistent with the #7 decision to trust the disk: the allele name is not stored in the schema at all, so the disagreement has no effect on the data — `dbd_seq` is the sequence and `mut_positions` is derived from it. The only exposure is a future join on allele name. | see #8 |

| 42 | **A bare `Insert sequence` label collapsed 88 GR09 genes into one cluster.** Keying inserts by their full construct name (needed for #39) assumed the label always carries a construct name. GR09's does not — it reads simply "Insert sequence" — so the key came out empty and every GR09 gene shared it, producing a single `GR09:` cluster of 45 distinct Zn_clus domains. Caught by the cluster-size distribution, not by any failure: the validator passed. Fixed by falling back to the gene name for bare and plasmid-named labels, with a regression test. | fixed |
| 43 | **Format assumptions UniPROBE has broken so far**, each caught by a distribution looking wrong rather than by an error: no header row (EMBO10); seven columns (Cell08); E-score in a different column, or absent (GR09/SCI09/RAD13A); a gene folder holding several proteins (ROG18A/NAR11/LIN14B/GD13); a bare sequence label (GR09). Per-source sanity summaries after every rebuild are the reason these were found; a green validator alone would not have caught any of them. | informational |

## Phase 2e — alignment-based cluster distance (2026-08-13)

| # | item | status |
|---|---|---|
| 44 | **Cluster members no longer have to be the same length.** Distance was Hamming, which excluded any variant with an indel and any variant whose padding had been clipped differently. It is now a pairwise alignment (BLOSUM62, affine gaps) with **free terminal gaps**, because `dbd_seq` is the padded envelope clipped by the construct: ROG18A's bare domains are all 83-85 aa while its padded ones run 95-105 aa purely from how much padding fitted. Charging for terminal gaps would report 8 phantom indels between FoxJ3 and its own chimera. | done, owner-requested |
| 45 | **`mut_positions` is now in the reference's coordinate frame**, so every member of a cluster shares one coordinate system and maps onto one predicted structure. A variant with a deletion may name a position past its own length; the validator bounds positions by the cluster's reference and rejects a cluster carrying variants but no reference row. | done |
| 46 | **ROG18A gains 9 more variants**: its 15 domains now form 3 clusters (FoxJ3 with 7 variants, FoxN3 with 4, FoxN2 with 1) instead of mostly singletons. The edit positions are compositionally exact — `_6aa` gives 6 consecutive positions, `_loop` gives 8, and `_6aa+loop` gives precisely their union — which is independent evidence the alignment places gaps correctly. | informational |
| 47 | **BAR15A refuses to rebase an allele across an indel.** Its construct-frame positions are shifted onto the padded domain by a constant offset, which an indel invalidates. No current allele triggers this; a future one would be rejected with a reason rather than silently mis-positioned. | done |

## Phase 2f — systematic audit of every source (2026-08-13)

Rather than reread each accession by hand, the invariants that past bugs violated were turned
into `scripts/audit_sources.py`, which sweeps all 30 registered accessions in ~4 minutes
without re-scanning Pfam. It checks archive coverage both ways, detail-page key collisions,
E-score readability, score shape, rows-per-domain, the shared 8-mer set, and cluster-size
outliers. Findings:

| # | item | status |
|---|---|---|
| 48 | **No construct-name mismatches anywhere.** The construct-splitting change of #39 drops no data in any source — the check that would have caught it is now part of the audit. | verified clean |
| 49 | **SCI09 loses nothing to the rejected combined files.** All 104 genes keep at least one readable per-replicate file; the 105 rejected 20-column files are redundant. | verified clean |
| 50 | **22 archive entries are protein complexes** — `Myc_Max`, `Kay_Jra`, `Da_Twi`, ten `HLH-2_*` heterodimers (Cell09), nine `CSL/NOTCH/MAML` complexes (PO10). Two or three chains form one binding unit, so none can satisfy admission condition 1. UniPROBE publishes no sequence for them, so they were already excluded — but by accident. Now an explicit rejection category in `docs/DOMAIN_POLICY.md`. | recorded |
| 51 | **PNAS08's archive is flat** — files sit at the top level with no gene folder, so the gene came out as a filename and nothing matched. Fixed by deriving the gene from the filename stem when there is no directory. **PNAS08 goes from 0 to 2 domains (65,792 rows).** | fixed |
| 52 | **PP15 must stay unmatched.** Its archive has both *Arabidopsis thaliana* and *A. lyrata* experiments, but its detail pages carry only the *thaliana* sequence. A prefix match would have assigned the thaliana sequence to the lyrata ortholog — a mis-attribution of exactly the kind the policy exists to prevent. Left unparsed deliberately. | closed, deliberate |
| 53 | **10 genes have a detail page with no sequence** (Cell08 4, MAR17A 5, GR09 1) — Hoxa3, Nkx3-1, Six6, Arid5b, Cebpa, E4f1, Xbp1 among them. UniPROBE simply does not publish those. Already reported per source by the parser. | informational |
