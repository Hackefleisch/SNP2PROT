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

### 2026-08-14 — Cluster-size restriction is a training decision, not a dataset one
Was `#2`: whether to admit only DBDs with ≥5 variants. **The dataset stores everything; the
modeller filters.** Only 9 of 425 clusters hold 5 or more domains, so the rule would have cut
the protein axis to almost nothing — but that is beside the point, which is that a dataset
should not bake in a training-set choice. What remains is the convenience of filtering by
cluster size cheaply: `TODO.md` `T3`.

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
full Pfam-A library (30,134 families), with the admission policy evaluated over a 41-family
DNA-binding whitelist (`data/external/pfam/dbd_families.txt`) so a DBD beside an unrelated
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
