# Decisions and resolved items

The record of what was **decided**, what was **fixed**, and what was **excluded** — the
counterpart to [`TODO.md`](../TODO.md), which holds only what is still open. An item leaves
`TODO.md` by arriving here.

This is a log, not a specification. How the dataset is built is
[`docs/METHODS.md`](METHODS.md); what `dbd_seq` is and what the admission policy excludes is
[`docs/DOMAIN_POLICY.md`](DOMAIN_POLICY.md); the owner's brief is
[`docs/TFDNA_MERGE_BRIEF.md`](TFDNA_MERGE_BRIEF.md).

Entries before 2026-08-14 carry their original `#N` from the archived open-items list
([`reports/archive/OPEN_ITEMS_2026-08-14.md`](../reports/archive/OPEN_ITEMS_2026-08-14.md)) so
older commit messages and code comments still resolve.

---

## 1. Scope

### 2026-08-14 — The dataset is PBM only
**Decided by the owner.** B1H and SNP-SELEX are out; the corpus grows by extending PBM
coverage instead. Then merge, then modelling. Tier 4 held-out sets are still wanted.

Consequences, all of them simplifications:
- every stored row is an 8 bp site, so the assay-identity leak through `dna_len` cannot occur
  while this holds (`TODO.md` `N5`);
- **padded-20 bp vs. common-core stops being a question** (was `#3`) — there is nothing of a
  different length to reconcile;
- **the SNP-SELEX OBS cutoff stops being a question** (was `#5`) — it was the last blocker on
  a phase that no longer exists;
- `b1h:` and `snp_selex:` in `configs/thresholds.yaml` become dead config (`TODO.md` `T8`).

### 2026-08-14 — upbm Q-values may enter as their own `score_type`
If the Kock deposit turns out to publish only upbm affinity/contrast/specificity Q-values and
no probe-level or E-score data, it is still admitted — with a new `score_type` and its own
cutoffs, rather than being skipped to protect the single-E-score scale.

**The cost is real and was accepted knowingly.** Every stored row is currently `pbm_escore` on
one scale with one pair of cutoffs, and the validator makes any value outside [-0.5, 0.5] a
hard error precisely because reading the wrong column once produced a 4:1 ratio (§5, `#36`).
A second score type means a second cutoff pair in `configs/thresholds.yaml`, and by this
project's own convention that invalidates the dataset, every report and every provenance row.
Weighed against 122 alleles in 30 designed series — the single largest addition of protein-axis
depth identified — and the depth won. Tracked as `TODO.md` `T11`.

### 2026-08-14 — Homeodomain may go to ~58% of domains; depth beats balance
Admitting Kock takes homeodomain from 47.9% to roughly 58% of domains and from 67% to ~84% of
point variants. Accepted as a straight consequence of the earlier call that family imbalance
is fixed by parsing more rather than down-sampling (`#10`), and because multi-member clusters
roughly double — 28 to ~58 — which is the corpus's scarcest property.

It does make the transfer question harder to answer, not easier: whether single-residue
sensitivity generalises to another fold is what `TODO.md` `N2` is about, and this deepens the
fold the corpus is already richest in. That was the trade, made with the numbers in view.

### 2026-08-14 — Clusters are formed by sequence distance, not only by construct lineage
Two natural paralogues one substitution apart previously sat in separate singleton clusters,
because `wt_id` came from construct naming — a gene folder and its insert names. They now
merge. Worth roughly a dozen new multi-member clusters, at no acquisition cost, and it
subsumes the distance-0 case that `TODO.md` `T4` describes.

`wt_id` thereby stops meaning "one reference and the variants engineered from it" and starts
meaning "domains within *k* edits of each other". Two parameters follow that this decision did
not settle, and neither has a defensible default — the threshold, and which member of a merged
paralogue pair is the reference that `mut_positions` is expressed against. Open as `D3`.

### 2026-08-14 — Construct architecture is recorded beside the table, not in it
Flank length, affinity tag and expression system are perfectly confounded with source study,
so they go in the protein-side companion table keyed by construct — not as columns in the
22-column row schema. Keeps the schema frozen and rewrites no parser, while still letting a
modeller condition on the covariate or hold it out.

Noted at the time: the confound is already partly live. `dbd_seq` is the padded envelope
clipped where the construct ends, so flank length is readable off `dbd_seq` length today —
ROG18A's bare domains run 83-85 aa against 95-105 aa for the padded ones.

**Closed 2026-08-17 as recorded rather than as fully built.** `data/interim/proteins/` carries
`construct_seq` with `dbd_start` / `dbd_end` marking the padded envelope inside it, so the
flank actually present on each side is exact and derivable per construct — the part of the
covariate that varies within a source and is therefore the part worth having. Affinity tag and
expression system are *not* recorded and are not planned: both are constant within a source
and would be read off `source_dataset` anyway, so a column would restate the source name with
extra steps. If a source ever varies its tag across constructs, this reopens.

### 2026-08-14 — Cluster-size restriction is a training decision, not a dataset one
Was `#2`: whether to admit only DBDs with ≥5 variants. **The dataset stores everything; the
modeller filters.** Only 9 of 425 clusters hold 5 or more domains, so the rule would have cut
the protein axis to almost nothing — but that is beside the point, which is that a dataset
should not bake in a training-set choice. What remains is the convenience of filtering by
cluster size cheaply: `TODO.md` `T3`.

### 2026-08-17 — PBM acquisition closes after Kock; `T2` and `T2b` dropped
**Decided by the owner.** Extending PBM coverage ends with `T11` (Kock et al. 2024). Two
tasks go with it:

- **`T2` — survey other PBM deposits.** Individual GEO / ArrayExpress submissions and paper
  supplements, always subject to the same gate as everything else: does the deposit publish
  the assayed construct sequence? Dropped as not worth the yield now that CIS-BP has landed
  868 domains and the corpus spans 56 families.
- **`T2b` — resolve UniProt accessions for the CIS-BP domains.** Table S6 publishes gene name
  and species but no accession, so 884 domains have no full-length sequence. Resolving them
  from gene plus species across 124 organisms would not be clean, and a wrong mapping is worse
  than a missing one — the reasoning that left `PP15` unparsed. Dropped; the consequence is
  that full-length coverage stays at about a quarter of the corpus, which is what `D1` is
  about. Domain-level and construct-level representations are unaffected and complete.

### Earlier — B1H and C2H2 arrays dropped on structural grounds
`#1`, `#4`, `#19`, `#20`. Persikov's B1H varies a different subunit than the one that binds,
and a C2H2 zinc-finger array is 2-6 separate ~23-residue folds on flexible linkers, each with
its own Zn(2+) — not one continuous unit, so it fails admission condition 2. ~8,000 domains
excluded. Reversible via `domain.allow_repeat_arrays`, and the weak-negatives worry about B1H
(`#1`) became moot when the phase went.

---

## 2. What `dbd_seq` is, and how clusters are formed

### `#18` — `dbd_seq` is the Pfam envelope padded by 10 residues each side
Not the bare envelope, not what was on the array. Pfam's models start at the structural core;
for homeodomains that clips the N-terminal arm, which reaches into the DNA minor groove and
makes real base contacts. Three BAR15A variants (KLF11 R402Q, VSX1 G160D, SNAI2 D119E) mutate
residues 2-9 N-terminal of the Pfam start — trimmed to the bare envelope each collapses onto
its own wild type while carrying a different label, invisible to any model. 10 covers the
observed maximum of 9 with margin, and is set from the biology rather than fitted to those
three cases.

### `#44` — Cluster distance is alignment-based, with free terminal gaps
Distance was Hamming, which excluded every variant carrying an indel and every variant whose
padding had been clipped differently. It is now a pairwise alignment (BLOSUM62, affine gaps).
**Terminal gaps cost nothing**, because `dbd_seq` is the padded envelope *clipped where the
construct ends*: in ROG18A every bare domain is 83-85 aa while the padded ones run 95-105 aa,
purely from how much padding fitted. Charging for that reported 8 phantom indels between
FoxJ3 and its own chimera. Internal indels still count.

### `#45` — `mut_positions` is in the reference's coordinate frame
So every member of a cluster shares one coordinate system and maps onto one predicted
structure. A variant with a deletion may therefore name a position past its own length. The
validator bounds positions by the cluster's reference and rejects a cluster carrying variants
but no reference row. One entry per edited residue, so `len(mut_positions) == n_mut_from_wt`
always holds.

### `#47` — BAR15A refuses to rebase an allele across an indel
Its construct-frame positions are shifted onto the padded domain by a constant offset, which
an indel invalidates. No current allele triggers this; a future one is rejected with a reason
rather than silently mis-positioned.

### `#30` — ROG18A's engineered chimeras cluster with their parents
Six substitutions is close enough. `FoxN3_J3` differs from `FoxN3` by 6 substitutions in
97 aa — real protein-axis depth, which the corpus is short of. Implementing it exposed a much
larger defect, `#39` below.

---

## 3. Composition decisions

### 2026-08-17 — `T4` closed: the label noise floor is measured, and it is not small
Forty-seven domains are stored by more than one source — 48 duplicate constructs, 49 source
pairs of one protein against one 32,896-8-mer set. The task had two jobs and both are now
discharged.

**Job 1, a shared split group, fell out of clustering.** These pairs are distance 0, so CD-HIT
puts them in one `wt_id` by construction and no split can separate them. Nothing further to do.

**Job 2, label agreement, is measured** by `reports.overlap_report`, regenerated into
[`reports/overlap.md`](../reports/overlap.md) by `scripts/make_overlap_report.py`. Two passes
over the interim tables, because only ~3M of 45M rows are replicated.

The number that matters is not the one that looks reassuring:

| measure | value |
|---|---|
| pooled hard-label agreement | 99.98% (300 of 1,572,985 calls differ) |
| **median positive-call Jaccard** | **0.458** |
| pooled positives agreed on | 2,863 of 5,816 (49.2%) |
| median E-score rank correlation | 0.608 |

Positives are under 0.5% of an 8-mer table, so agreeing on the negatives holds pooled
agreement near 100% however badly two labs agree about binding. Against the within-source
figures already in `METHODS.md` §6 — 70.4% and 72.2% positive-set Jaccard for byte-identical
distinct genes and for technical replicates — cross-laboratory agreement is some 25 points
worse. **A model that reproduces held-out positives much past this is reproducing a source.**

Nine pairs are flagged below Jaccard 0.20 or rho 0.30, reported and not repaired. One is not
noise on any reading: `C:LIN14B:NAP` — *Arabidopsis* ANAC092 / O49255, an identical 144 aa
stored domain in both `LIN14B` and `weirauch2014` — shares **zero** positive calls between the
two deposits, at rho 0.068. Twelve positives against 126, disjoint. Open as `D4`.

A second, smaller inconsistency surfaced with it: the domain shared by `Cell08` and
`weirauch2014` under `C:Cell08:Tlx2` is byte-identical but carries `species` as *Mus musculus*
in one source and *Homo sapiens* in the other. Metadata only, no label affected, folded into
`T5`.

### `#9` — The 318:1 negative:positive ratio stays in the stored table
Controlled when sampling training batches, not by discarding rows. A PBM low-E-score negative
is real evidence of non-binding and throwing it away destroys information. `pbm.positive`
stays at 0.45.

### `#10` — Family imbalance is fixed by parsing more panels, not by down-sampling
Homeodomain has fallen from 79.5% to 47.9% of domains this way.

### `#27`, `#33` — Family-share and family-size worries, closed
The brief's 40% homeodomain flag was written because "C2H2 will otherwise swamp everything";
C2H2 is now 4 domains (0.8%), so the threat it guarded against does not exist. Homeodomain is
also *more* internally diverse than the minority families — median pairwise identity 34%, only
1.1% of pairs ≥90% identical, against Forkhead's 20.3%. Nothing is being swamped. Which
families are large enough to hold out is now a note, `TODO.md` `N4`, not a decision.

### `#24`, `#32` — Protein-axis depth accepted as-is
81 point variants in 28 clusters. `#32` framed this as the project's central gap — three
families with zero point variants and 54 of 66 variants homeodomain. The `#39`/`#46` fixes
moved Forkhead from 6 variants to 18 and HLH from 0 to 3, so the premise no longer holds as
stated. **Owner: acceptable while point variants are not the majority of the data**; subsets
can be built later, and the question resurfaces on its own when modelling starts. Carried
forward as `TODO.md` `N2`.

---

## 4. Individual constructs

### `#7` — Keep all 120 BAR15A variant alleles, not the paper's 117
Trust what is on disk. The three extra are HOXD13 `I297V`, `N298S`, `Q325K`, whose files are
named `*_8mers_11111111.txt` — exactly what a naive glob drops, which is the likely origin of
117. HOXD13 survives the policy with all 8 domains.

### `#8`, `#41` — `PROP1_R112Q` kept as deposited
Its name says residue 112 Arg→Gln; the deposited sequence gives Arg→**Met** at that same
position. Position and wild-type residue are both right; only the replacement disagrees.
PROP1's other allele `R99Q` matches its name exactly, and the Barrera supplement lists R99Q
but **not** R112Q. Both Arg→Gln and Arg→Met are single-nucleotide changes from an Arg codon,
so nothing available settles which the clone really was.

It matters more than most: construct position 57 maps to **homeodomain position 43, the first
residue of helix 3 — the recognition helix** that reads bases in the major groove.

Kept, because the allele name is stored nowhere in the schema. `dbd_seq` is the sequence and
`mut_positions` (53, within the padded domain) is derived by diffing against the reference,
never parsed from the name. The stored claim — "this domain, with M at position 53, binds
these 8-mers" — is true whichever name is right. **The one exposure is a future join on allele
name `R112Q` that assumes Gln; anyone doing that must read this entry first.**

`PITX2_T114P` remains excluded: its insert shows no substitution at all. Two further anomalies
vanished when their genes were rejected as mixed-family.

---

## 5. Defects found and fixed

Each of these was caught by a distribution looking wrong, **never by a failing validator**.
That is why per-source sanity summaries exist, and why `scripts/audit_sources.py` now sweeps
the invariants they violated. The live checklist is `TODO.md` `N3`.

### `#34` — The HMM set covered 9 families; UniPROBE spans 86 domain labels
Constructs from bZIP, HMG_box, GATA, T-box, Zn_clus, IRF and others were rejected as
`no_domain` because the HMM was missing, not because nothing was there. Fixed by scanning the
full Pfam-A library (30,134 families), with the admission policy evaluated over a DNA-binding
whitelist (`data/external/pfam/dbd_families.txt`, 41 families at the time, 65 after `T15`
below) so a DBD beside an unrelated
domain is still admitted. **340 → 468 domains, 7 → 31 families.**

### `#35` — UniPROBE labels predate several Pfam renames
`Homeobox`→`Homeodomain`, `Fork_head`→`Forkhead`, `MADS`→`SRF-TF`, `E2F_TDP`→`WHD_E2F_TDP`,
`Zn2Cys6`→`Zn_clus`, `BRIGHT`→`ARID`, `NR`→`zf-C4`, `AP-2`→`TF_AP-2`, `PBX`→`PBC`,
`RFX`→`RFX_DNA_binding`. Each resolved by Pfam accession lookup, never by string guessing.
Missing `Homeobox` alone would have excluded the largest family.

### `#36`, `#37` — The E-score column was read by position and was wrong for three accessions
`usecols=[0, 2]` holds for BAR15A/EMBO10/PNAS13/Cell08 but not for GR09 and SCI09 (headerless,
median intensity at index 2, E-score at index 3) or RAD13A (`8-mer / 8-mer / Median / Z-score`
— no E-score at all). The corpus briefly showed a 4:1 negative:positive ratio with `raw_score`
up to 776,106. The column is now identified by header name where present, else by being
bounded to [-0.5, 0.5] **and** taking negative values (p/q-values cannot). Files with no such
column, or several, are rejected with a reason. A validator guard makes the class of bug
impossible to miss: any `pbm_escore` outside [-0.5, 0.5] is a hard schema error.

### `#39`, `#40` — One gene folder is not one protein
The costliest defect so far. UniPROBE archives encode the construct at path level 2 and detail
pages can carry several insert sequences, but the parser grouped by gene and took the first
insert — **merging distinct engineered proteins as replicates**. ROG18A's `FoxJ3_N3/`, six
different chimeras, became one protein; NAR11's `HLH-1/`, a point-mutant series L13R/L13T/L13V,
became one; LIN14B's `ANAC092_DBD` and `ANAC092_FL` were averaged together; GD13's two CLAMP
constructs likewise. Fixed by keying inserts on their full construct name and splitting
experiments against those names. **ROG18A 7 → 15 domains, NAR11 1 → 4**, and it produced the
corpus's first two non-homeodomain variant clusters.

### `#42` — A bare `Insert sequence` label collapsed 88 GR09 genes into one cluster
Keying inserts by construct name (needed for `#39`) assumed the label always carries one.
GR09's reads simply "Insert sequence", so the key came out empty and every GR09 gene shared
it, producing a single `GR09:` cluster of 45 distinct Zn_clus domains. Caught by the
cluster-size distribution; the validator passed. Fixed by falling back to the gene name for
bare and plasmid-named labels, with a regression test.

### `#46` — Evidence the alignment places gaps correctly
After `#44`, ROG18A's 15 domains form 3 clusters (FoxJ3 with 7 variants, FoxN3 with 4, FoxN2
with 1) instead of mostly singletons. The edit positions are compositionally exact: `_6aa`
gives 6 consecutive positions, `_loop` gives 8, and `_6aa+loop` gives precisely their union.

### `#51` — PNAS08's archive is flat
Files sit at the top level with no gene folder, so the gene came out as a filename and nothing
matched. Fixed by deriving the gene from the filename stem when there is no directory.
**PNAS08 goes from 0 to 2 domains (65,792 rows).**

---

## 6. Sources and constructs excluded

### `#58` — `Path10`: tables truncated to the enriched end
Its files carry a valid E-score column but only 341-1,391 rows each, cut at E ≥ 0.25, against
the 32,896 of a full universal PBM. It therefore holds **no non-binding evidence at all**, and
applying our cutoffs would label its enriched 8-mers between 0.25 and 0.35 as non-binding —
the exact inverse of the truth. 18 *Plasmodium* AP2 proteins excluded; deregistered, reason in
`PROVENANCE.md`.

### `#59` — `GD09`: excluded twice over
Its tables hold 32,895 of the 32,896 non-redundant 8-mers (the palindrome `CCTTAAGG` is
missing), and the fully crossed design means a protein measured against a different 8-mer set
cannot join the others. Independently, its single protein NSY-7 — deposited twice, as `NSY-7`
and `Nsy-7` — returns no Pfam hit anywhere in its 294-residue construct. The completeness
error now distinguishes "truncated to the enriched end" (Path10, ~4% of the set) from
"near-complete but *n* short" (GD09, one short).

### `#52` — `PP15` must stay unmatched
Its archive holds both *Arabidopsis thaliana* and *A. lyrata* experiments, but its detail pages
carry only the *thaliana* sequence. A prefix match would have assigned the thaliana sequence to
the lyrata ortholog — a mis-attribution of exactly the kind the policy exists to prevent. Left
unparsed deliberately.

### `#50` — 22 archive entries are protein complexes
`Myc_Max`, `Kay_Jra`, `Da_Twi`, ten `HLH-2_*` heterodimers (Cell09), nine `CSL/NOTCH/MAML`
complexes (PO10). Two or three chains form one binding unit, so none can satisfy admission
condition 1. UniPROBE publishes no sequence for them, so they were already excluded — but by
accident. Now an explicit rejection category in [`DOMAIN_POLICY.md`](DOMAIN_POLICY.md).
**Directly relevant to Tier 4**, which is bHLH dimers: `TODO.md` `T7`.

### `#53` — 10 genes have a detail page with no sequence
Cell08 4, MAR17A 5, GR09 1 — Hoxa3, Nkx3-1, Six6, Arid5b, Cebpa, E4f1, Xbp1 among them.
UniPROBE simply does not publish those. Reported per source by the parser.

### `#38` — 12 of 30 registered accessions contribute nothing
Eleven rejected entirely by the domain policy (C2H2 arrays, multi-domain, or families outside
the DBD whitelist); the rest for the data reasons above. All reported per source, none silent.
Attrition is tabulated in [`METHODS.md`](METHODS.md) §9.1.

---

## 7. Build performance

### 2026-08-14 — The build was profiled, then made ~6x faster without changing a byte
Was `TODO.md` `T6`. `build_dataset.py --all` took ~25 minutes. Measured on real source
builds, two things dominated and neither was doing useful work.

**The Pfam-A load was 100% fixed cost.** Scanning 3 sequences took 19.5 s; scanning 100 took
19.4 s. The library is 2.2 GB of *text*, reopened once per accession, so a 30-source build
spent about ten minutes re-parsing one unchanging file to do ~0.1 s of scanning.

**Validation re-checked constant columns 32,896 times per construct.** `dbd_seq` has 489
distinct values in the corpus but every check ran per row, and
`df.apply(_positions_in_range, axis=1)` built a pandas Series for each of 16.6 M rows to test
a handful of integers — 50 s of a 139 s corpus validation on its own.

Four changes, in descending order of payoff:

| # | change | effect |
|---|---|---|
| 1 | **Press the library** (`scripts/press_pfam.py`) | 19.5 s → 3.5 s per source |
| 2 | **Cache the loaded library per process** (`domains._library`) | → 0.5 s once, then ~0.6 s per source |
| 3 | **Validator axis checks on distinct values**, weighted back to row counts | validate 162.7 s → 27.0 s over the corpus |
| 4 | **Vectorized `add_pair_ids` payload** | 2.8x; it runs twice per build, in `coerce` and again in `validate` |

**Nothing about the dataset changed, and that was checked rather than assumed:**
- all 18 sources revalidated — **every report string, every stat and every `pair_id`
  identical** to the previous implementation;
- six sources fully rebuilt and compared against the committed Parquet — **byte-identical**;
- the pressed library checked against the unpressed one by moving the `.h3*` files aside:
  identical calls, 169 hits over 150 constructs;
- 8 regression tests added. They pin *behaviour*, not speed — specifically that the
  row-weighted counts in every validator message still equal what the per-row implementation
  produced, since that weighting is the one thing a future edit could silently drop.

**Pressing is a per-machine setup step**, not a repo artifact: the `.h3*` files are 2.6 GB of
git-ignored derived data beside the library, and re-downloading Pfam-A means pressing again.
`domains._library` says so out loud when it is handed a large unpressed library rather than
silently running six times slower.

**Parallelism was deliberately left out** and is `TODO.md` `T9`. The serial fixes were worth
more at a fraction of the risk: parallelism means restructuring `build_dataset.py` and
`reconcile_replicates`, and process pools make failures harder to attribute — which matters
for a build whose per-source rejection reasons are load-bearing. Threads were measured at
**1.0x** (the CSV parser holds the GIL); processes at **7.6x**.

---

## 8. Questions closed without action

| was | resolution |
|---|---|
| `#6` Are UniPROBE 8-mer files licence-gated? | No. Direct download via `filehandler.php`; no click-through, cookie or MTA. |
| `#11`, `#15`, `#16` Competing `dbd_seq` conventions; EMBO10/PNAS13 lacking a DBD field; the 719 aa `Jumeau` outlier | All dissolved by the domain policy. Every domain is `pfam_hmmer_padded`; `dbd_seq` length no longer correlates with source. Jumeau is 107 aa. |
| `#12` Noyes 2008 accession | Not in UniPROBE — it is B1H, not PBM. |
| `#13` Which panels supply forkhead / ETS / bZIP | Identified and all parsed: EMBO10, PNAS13, ROG18A, MAR17A, SCI09, GR09, SHO18A. |
| `#14` How should identical-domain pairs be split? | Reframed, not closed. The 11 same-source Cell08 pairs it described are no longer byte-identical — the 10 aa padding separated them (Evx1/Evx2 now differ at 12 of 77 positions, Lhx2/Lhx9 at 3). The real case is 17 **cross-source** pairs: `TODO.md` `T4`. |
| `#17` Panels add breadth, not depth | Confirmed and accepted — see `#24`. |
| `#21` Cost of the domain policy | Measured: 11.81M → 9.15M rows at the time. |
| `#28` Parse `GB11`? | Moot — its detail pages carry no sequence at all, so it is unparseable regardless. HLH reached 35 domains by other means. |
| `#31` SHO18A detail pages 500ing | Worked around via `detailsDef.php` with the site's own ids, verified against SHO18A/Apt. No ids guessed. |
| `#48` Did construct-splitting drop data anywhere? | No — verified clean across all 30 accessions, and the check is now part of `audit_sources.py`. |
| `#49` Does SCI09 lose anything to its rejected combined files? | No. All 104 genes keep at least one readable per-replicate file; the 105 rejected 20-column files are redundant. |

---

## 9. Source candidates screened and rejected

### 2026-08-14 — a literature sweep for PBM sources with designed protein variation
A separate research session surveyed the literature for PBM datasets carrying reference plus
mutant alleles of the same DBD. Screening its results against the admission policy and
against what is already parsed:

| candidate | verdict |
|---|---|
| **Kock et al. 2024, *Nat Commun* 15:3110** — 30 HD allelic series, 122 alleles | **pursue** — `TODO.md` `T11`. The one large find. |
| Rogers et al. 2019, *Mol Cell* 74:245 | **already parsed as `ROG18A`** — 15 domains in 3 clusters, including the `N3(J3-6aa)` swap the sweep highlighted. |
| Liu et al. 2018, *eLife* 7:e34594 | **already parsed as `LIU18B`**, but only 2 of its constructs. The lead is the Dryad deposit, not the paper: `T14`. |
| Ibrahim et al. 2013, *Genome Res* 23:2091 — HOXD13 wt/Q325R/Q325K | **mostly held already.** `Q325K` is in BAR15A, which carries HOXD13 as 8 domains / 7 variants. Only `Q325R` would be new — one domain. |
| **Siggers et al. 2014, *Mol Cell* 55:640** — Msn2/Msn4/Com2/Usv1/Rgm1 | **rejected.** All are C2H2 zinc-finger arrays: separate ~23-residue folds on flexible linkers, failing admission condition 2. These are the exact grounds on which Phase 3 was dropped and ~8,000 domains excluded (§1). Reversible only via `domain.allow_repeat_arrays`. |
| "Mine Berger 2008 for pairs at ≤3 substitutions" | **measured, and smaller than proposed.** Over all 489 domains: 31 pairs at ≤3 edits, but 23 are cross-source and largely re-find the 17 `T4` duplicates plus their variant halos. Genuinely new same-source paralogues at ≤3 edits: 8 pairs. Now handled generally by `T12`. |
| Chu 2012; Noyes 2008; Aditham 2021; gcPBM paralogue panels | **out by assay**, as the sweep itself noted — B1H, microfluidic affinity, or an array design that does not yield the 32,896 8-mers. |

**Identifiers from that sweep are unverified** and are recorded in `TODO.md` as such. Rule 1
applies: none of them may be turned into a download URL by pattern, and each must be confirmed
against the source before use.

---

## 10. The CIS-BP gate, answered

### 2026-08-14 — CIS-BP publishes the assayed construct sequence
`TODO.md` `T1` rested on one question: does CIS-BP publish the sequence that was physically on
the array, or only its own DBD annotation of the full-length protein? Attributability — the
first admission condition — depends on the answer, and it decided whether CIS-BP was the
largest available expansion or a dead end.

**It publishes the construct.** Weirauch et al. 2014's own methods: *"All inserts were
sequence verified in full. Insert sequences and other information are available in Table
S6."* The project site describes Table S6 as *"DBD clone source material. This spreadsheet
provides information on the clone source material, the experimental construct sequences, and
the clone source contributors."*

Verified at every link rather than taken on the paper's word:

| link | contents |
|---|---|
| `TabS6_DBD_clone_information.xlsx`, sheet *Experimental constructs* | 1,032 rows with `Plasmid ID`, `Insert AA`, `#Flanking AAs`, `Backbone`, `Tag location` |
| GEO `GSE53348` | 2,064 samples, titled `pTH####_<HK\|ME>_8mer_<n>` — the plasmid ID is the join key |
| `GSM1291226` (`pTH1294_HK_8mer_593`) | `ID_REF / VALUE / E-Score / Z-Score`, **32,896 rows** |

The construct architecture is stated per row, so the confound `T13` exists to record is
published rather than inferred: 671 constructs carry 50 flanking residues, 96 carry 15, and
**265 carry none at all** — bare DBDs from an Agilent oligo pool, described in the methods as
*"These constructs did not contain any AAs flanking the DBD."*

**Measured yield through our own admission policy**, by scanning all 1,032 inserts and calling
each with `call_domain`: **731 admitted (71%)** — 240 rejected `no_domain`, 49 `repeat_array`,
12 `mixed_families` — giving 722 distinct sequences of which **677 are new**. That takes the
protein axis from 489 domains to roughly 1,166, across 106 species instead of 24, and fills
the families the corpus is thinnest in (Myb_DNA-binding 58 against our 2, AP2 39 against 3,
GATA 27 against 2, WRKY 22 against 1).

It also brings the depth the corpus lacks, because the paper's first selection strategy
deliberately populated nine DBD identity bins up to **90-99.99%**: within the new domains
alone there are 61 pairs within 3 edits and 103 within 5, led by Myb_DNA-binding — a family
with no variant depth here today.

**Nothing about the existing configuration has to change**, which is what separates this from
Kock: the same E-score statistic on the same scale, with the paper quoting `E > 0.45` as its
own significance cutoff — identical to `pbm.positive` — and HK/ME dual arrays per plasmid,
which `per_experiment: true` and `reconcile_replicates` already handle.

Two limits recorded honestly. The 240 `no_domain` rejections come from a stricter bar than
CIS-BP's own (gathering thresholds over full Pfam-A, against their 81 models at Eval < 0.01)
and should be spot-checked rather than assumed correct. And this is Weirauch 2014's own 1,032
constructs, **not** the ~2,294 TFs with PBM data the database aggregates: for aggregated
entries the construct sequence belongs to the contributing study, and much of that is UniPROBE
already held.

### 2026-08-14 — The whitelist is built from every source database's vocabulary, not one
`#34` fixed the HMM library but left the *whitelist* derived from a single database:
`build_dbd_family_list.py` read UniPROBE's `Domain` labels, so it could only ever name
families UniPROBE happened to publish. Screening CIS-BP exposed the same failure a second
time — of 240 constructs rejected as `no_domain`, **only 33 had no Pfam hit at all**; the
other 207 had a perfectly good hit belonging to a family we had no word for. Plant TCP, Dof,
NAC and SBP, Doublesex DM, the Rel homology domain of NF-kB.

The script now reads **both** UniPROBE's `Domain` labels and CIS-BP's `Pfam ID` column, each
being that database's own curation of what binds DNA, and each resolved against Pfam-A `NAME`
fields rather than trusted verbatim. **41 -> 65 families**, 24 added, none removed. The output
records which database vouched for each family.

Six CIS-BP labels predate Pfam renames and could not be resolved by string match. Rather than
guess, each was resolved **by the sequence**: scan the constructs CIS-BP gave that label to
against full Pfam-A and read off which family actually hits them. The method reproduces the
three existing hand mappings (`Homeobox`, `Fork_head`, `E2F_TDP`) exactly, which is why it was
trusted for the rest, and every result is confirmed by the target's own Pfam description —
`DUF573` -> `GeBP-like_DBD` ("DBD domain"), `RHD` -> `RHD_DNA_bind` and not `RHD_dimer`,
`EIN3` -> `EIN3_DNA-bd` and not `EIN3_N`.

Effect, measured: CIS-BP admission **731 -> 896**, and on sources already parsed **477 -> 490**
constructs with nothing lost — `LIN14B` alone goes from 1 admitted construct to 12, because
NAC is a family the corpus previously could not see. A rebuild follows (`TODO.md` `T16`).

Pinned by `tests/test_domains.py::test_whitelist_covers_the_families_recovered_from_cisbp`,
because this failure mode is silent: a missing family looks exactly like an absent domain.

### 2026-08-14 — Corpus rebuilt on the widened whitelist
`T16`. `build_dataset.py --all` re-run after `T15`; same 18 of 30 accessions contribute, all
validator-clean, and the twelve that yield nothing do so for their existing recorded reasons.

| | before | after |
|---|---:|---:|
| rows | 16,645,376 | **17,040,128** |
| distinct domains | 489 | **501** |
| clusters | 425 | **437** |
| Pfam families | 30 | **32** |
| binding / non-binding / excluded | 51,387 / 16,323,689 / 270,300 | 51,642 / 16,713,577 / 274,909 |
| admitted constructs | 568 (72.9%) | 581 (74.6%) |
| rejected `no_domain` | 55 | 42 |

Every added domain traces to the whitelist: `LIN14B` 1 → 12 constructs (11 `NAM`) and `SCI09`
68 → 69 (`MH1`). Only those two sources' reports changed, which is the cleanest available
evidence that nothing else moved.

The row-count invariant still holds exactly — **17,040,128 = 518 × 32,896**, where 518 is 501
distinct domains plus the 17 measured by two studies each. Every `dna_len` is 8. Cluster
structure is otherwise untouched: still 81 variants in 28 clusters, and the new NAC domains are
all singletons, so this bought breadth rather than depth.

The protein table was rebuilt with it: 501 domains, canonical UniProt for 462 (92%), domain
located in the full-length sequence for 354 (71%).

### 2026-08-14 — `weirauch2014` parsed; the corpus nearly triples
`T1` closed. 29,080,064 rows from 884 domains, validator-clean with no warnings, joined to
Table S6 on plasmid ID. The corpus goes from 501 domains to **1,364**, 32 families to **56**,
24 organisms to **138**, and 17,040,128 rows to **46,120,192**.

**What it bought, and what it did not.** Families holding ten or more domains went from 8 to
**27**, which turns leave-one-family-out from a handful of proteins into a real test and
largely settles the concern recorded as `#33`. Cross-source duplicate domains went from 17 to
38, so the label-agreement check of `T4` now has twice the evidence and, for the first time,
compares two different laboratories rather than two deposits from one. But **protein-axis depth
is untouched**: all 81 variants still sit in 28 clusters, and 1,293 of 1,321 clusters hold a
single domain. This was breadth, exactly as predicted before parsing.

**Two costs, both recorded rather than absorbed.** Full-length UniProt coverage fell from 71%
to 26%, because Table S6 publishes gene and species but no accession (`T2b`). And the source
carries three construct architectures — 671 with 50 flanking residues, 96 with 15, **265 with
none** — so its stored domains are systematically shorter than a UniPROBE construct of the
same protein. That is `T13`'s confound, now present in the data rather than anticipated.

Three pieces of structure came out of the work:

- **`parsers/_pbm.py`.** Everything that depends on the assay rather than the distributor —
  E-score identification, completeness, binarization, replicate combining, frame assembly —
  now has one home. `_uniprobe.py` keeps only UniPROBE's file layout and re-exports the rest,
  and five sources were re-parsed to confirm the extraction is byte-identical.
- **A latent bug, found by the extraction.** A header naming the E-score column used to
  return early, skipping the row-count check entirely, so a *headed* but truncated table would
  have been admitted silently; only headerless sources like `Path10` were ever caught. The
  check is now uniform. It matters here: GEO's tables are headed.
- **`build_protein_table.py` is no longer UniPROBE-shaped.** It derives its sources from the
  parser registry, which was right, but assumed every source publishes HTML detail pages —
  registering this parser made it raise `FileNotFoundError`. Construct lookup is now per
  source.

---

## 11. The canonical domain sequence

### 2026-08-17 — `dbd_seq` becomes construct-independent
Until now `dbd_seq` was the Pfam envelope plus 10 residues, **clipped wherever the assayed
construct happened to end**. So the same domain from two labs could be stored as two different
strings — VENTX is 76 aa from BAR15A and 57 aa from CIS-BP's zero-flank construct. That is an
artefact of cloning, and it reaches the model as if it were biology.

**Decided: the stored sequence is a project-internal canonical form.** Not "canonical" in
UniProt's sense — a definition local to this dataset, whose defining property is *stability*:
the same domain in gives the same string out, however truncated the input was. Length still
varies with the biology, because an insertion or deletion genuinely is a different sequence.

The window is **the Pfam envelope ± 10 residues**, and the padding stays: it is part of the
canonical definition rather than an addition to it, and now that it is taken from a reference
it is always fully available instead of surviving by luck. VSX1 G160D mutates 6 residues
before the Pfam start and currently survives only because BAR15A's construct happens to be
long enough.

Why the stored residues, not a fixed-width alignment: a match-state-only representation
**collapses VSX1 G160D onto its own reference** — identical sequence, different label — which
is the exact failure the padding was introduced to prevent (§2, `#18`).

### 2026-08-17 — A reference is consulted only where the construct falls short
Measured over 1,477 admitted constructs: **1,120 (76%) already cover envelope ± 10**, and for
those the window extracted from the construct is provably the string the full-length protein
would give. 37 canonical windows already come out byte-identical across sources
(`BAR15A:ARX_REF` = `Cell08:Arx`, 77 aa).

The remaining 357 need an external reference, and every one of them carries an identifier — 68
a UniProt accession, 289 a CIS-BP Gene ID. Resolution through UniProt cross-references was
verified working for Ensembl, FlyBase and Araport. Roughly 121 of the CIS-BP identifiers are
standard cross-references; the rest are assembly-specific and will partly fail.

**Expect to lose on the order of 10%.** That is accepted: a construct that cannot be placed is
discarded rather than force-fitted, which is the `PP15` rule (§6) — mis-attribution is worse
than absence.

The reference **extends**, it never replaces. The construct's own residues are kept, carrying
whatever mutations were engineered into them; the reference supplies only flanking residues the
construct is missing. So a shorter sequence caused by a real deletion or a true protein
terminus stays shorter — that is the thing that was assayed.

### 2026-08-17 — Placement is accepted at ≤5 edits outside the envelope
100% identity is impossible: the corpus is full of deliberate point mutants. Measured on the
462 constructs where a reference is already held, edits between construct and reference are
0 at the median and:

| ceiling | accepted | with ≥99% of the construct aligned |
|---|---:|---:|
| ≤1 | 94.6% | 92.4% |
| ≤2 | 96.3% | 94.2% |
| **≤5** | **97.6%** | **95.2%** |
| ≤10 | 97.8% | 95.2% |

Five sits at the knee — ten buys 0.2% more. The 11 it rejects are genuinely wrong, not
marginal: `Hoxc11` at 111 edits over 57% coverage, `Etv4` at 36% coverage.

**Edits are counted outside the Pfam envelope only.** Inside it, mutations are the subject of
the dataset and are kept unconditionally; outside is where residues are borrowed, so that is
what must match. A flat ceiling would reject ROG18A's chimeras, which differ from their parent
by 6-8 substitutions while being perfectly well placed. Mismatches and gaps are counted
separately: `Mlx` and `Rfx3` show 0-1 mismatches with 25-54 gaps, which is an alternative
isoform rather than a wrong protein.

### 2026-08-17 — Constructs differing outside the canonical window are discarded, not merged
Of 54 stored domains produced by more than one construct, 24 come from byte-identical
constructs and 26 differ only in how much flank each lab cloned — both benign. **Four carry a
genuine substitution outside the stored window** and were being silently reconciled as
replicates of one protein:

| substitutions outside | stored | constructs |
|---:|---:|---|
| 4 | 76 aa | `Cell09:HLH-25` / `Cell09:HLH-27` — two distinct *C. elegans* genes |
| 2 | 70 aa | `weirauch2014:pTH8163` / `pTH9718` |
| 1 | 106 aa | `MAR17A:Foxc1` / `BAR15A:FOXC1_REF` / `weirauch2014:pTH2673` |
| 1 | 77 aa | `Cell08:Msx3` / `MAR17A:Msx3` |

**Owner's decision: bin them.** A difference the model cannot see, attached to measurements
that may differ because of it, is a confound — the same reasoning as admission condition 3,
applied symmetrically. Until now condition 3 was enforced only in `bar15a.py`; in the panel
and CIS-BP parsers it was vacuous, because `mut_positions` is computed *from* the stored
sequences and so nothing could fall outside by construction. This makes it a real check
everywhere.

### 2026-08-17 — Clusters are built on canonical domains, never on reference proteins
Tempting, since references are construct-independent and a cluster representative need not be
a dataset row. But DNA-binding domains are conserved while the rest of the protein diverges.
Among pairs whose domains are within 5 edits, median domain identity is **96%** while median
full-length identity is **68%**, and the extremes are stark: **Irx3 and Irx4 have 94% identical
domains and 32% identical proteins.** Any protein-level clustering separates them, after which
a model trained on Irx3 predicts Irx4 for free — precisely the leakage clusters exist to
prevent. Also Hoxa10/Hoxd10 at 45%, Vax1/VAX2 at 46%.

The reference's job is to make the domain canonical and to establish protein identity for
deduplication. It is not the clustering unit.

### 2026-08-17 — CD-HIT greedy incremental is the clustering algorithm
Not invented here: the published greedy incremental algorithm (Li & Godzik 2006; MMseqs2
`--cluster-mode 2`). Sort by decreasing length, longest becomes the representative, and each
remaining sequence is compared **only to representatives**. Three properties earn it the job:

* the longest member becomes the representative, so the least-clipped form is the reference;
* comparison is never transitive, which kills chaining — single-linkage produced 35-domain
  blobs at 5 edits, against 8 at 1 edit;
* every member is within *k* of the representative, which is the invariant `mut_positions`
  needs (`#45`).

Implemented on `snp2prot.align` rather than by installing MMseqs2: the aligner already has the
right semantics, including free terminal gaps, and a full pairwise scan of the corpus takes
16 seconds. Adding a C++ dependency to a 16-second job is not warranted; the algorithm is
reused, which is the part that matters.

### 2026-08-17 — Built, and the corpus rebuilt on it
`snp2prot.canonical`, `snp2prot.references` and `snp2prot.clusters`, wired into all three
parsers, with `scripts/build_clusters.py` assigning `wt_id` corpus-wide.

| | before | after |
|---|---:|---:|
| rows | 46,120,192 | **45,495,168** |
| distinct domains | 1,364 | **1,335** |
| clusters | 1,321 | **1,133** |
| clusters holding >1 domain | 28 | **122** |
| variants | 81 | **202** |

Protein-axis depth roughly tripled, and stopped being a homeodomain story: 87 of 202 variants
are homeodomain against 54 of 81 before, with Myb 16, bHLH 15, forkhead 13, AP2 9, zf-C4 7,
SAND 5, RFX 5, TCP 5. That is the corpus finally able to ask whether sensitivity to a single
residue transfers between folds.

**Engineered constructs keep a short window** (owner, 2026-08-17). ROG18A's FoxJ3/FoxN3
chimeras are synthetic and cloned as bare domains: no natural protein exists to extend them
from, and their termini are not truncations. Discarding them cost the corpus its only
non-homeodomain variant depth, which was the wrong trade. A short entry that duplicates a
fully-canonical copy of the same domain is dropped at clustering instead, so no two rows
describe one domain at two lengths.

**Clustering threshold: 5 edits** (owner). 1,133 clusters, largest 8 — against the 35-domain
blob single-linkage produced at the same threshold, which is CD-HIT's non-transitivity earning
its place.

Four defects were found and fixed while doing this, each recorded because each was silent:

* **Stale interim output.** `build_dataset.py` left a source's previous Parquet on disk when
  that source stopped yielding, and it was still being read as part of the corpus. Two
  accessions were affected. The build now deletes it.
* **Cluster ids collided.** Naming a cluster after its representative's lineage `wt_id` is
  readable but not unique — `ROG18A:FoxJ3` names six different representatives, so six
  clusters collapsed into one, undoing the clustering for exactly the chimeras that had just
  been rescued. Collisions now take a suffix.
* **The validator assumed clusters are source-local.** They are not: a `Cell08` variant
  routinely belongs to a cluster represented by a `BAR15A` domain, which a per-source
  validator cannot see. It reported every such cluster as a missing reference. Now a warning,
  with the bounds check done corpus-wide — where it passes for all 202 variants.
* **An O(rows) alignment.** `build_clusters.py` computed an edit profile per row rather than
  per distinct domain: 45 million alignments instead of 1,335.

**Known gap:** the protein-side table covers 1,224 of 1,335 domains, because it locates a
domain by substring within its construct and a reference-extended canonical sequence is no
longer a substring of it. Tracked as `T18`.
