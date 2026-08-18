# TODO

The working file for this project: **open tasks, open decisions and notes**. One entry per
thing, each short enough to read in ten seconds.

- Something **decided** or **finished** moves out of here into [`docs/DECISIONS.md`](docs/DECISIONS.md).
- Something **measured** about the dataset goes in [`docs/METHODS.md`](docs/METHODS.md), not here.
- Nothing in here is resolved unilaterally — tasks are executed, decisions are the owner's.

Entry ids (`T1`, `D1`, `N1`) are stable. The pre-2026-08-14 list used a different `#N`
scheme; it is frozen at [`reports/archive/OPEN_ITEMS_2026-08-14.md`](reports/archive/OPEN_ITEMS_2026-08-14.md).

---

## Where things stand

**This is a pure PBM dataset.** Decided 2026-08-14. B1H, SNP-SELEX and everything that would
have mixed assay types are out — which removes the `dna_len` leakage problem entirely, since
every row is an 8-mer and always will be.

19 PBM sources parsed and validator-clean: **45,495,168 rows, 1,335 domains in 1,133
clusters, 56 Pfam families, 137 organisms**, every protein scored against the same 32,896
8-mers.

`dbd_seq` is the **canonical domain** — envelope ± 10, construct-independent
(`docs/DECISIONS.md` §11) — and `wt_id` comes from CD-HIT clustering at 5 edits, not construct
lineage.

**The protein axis is no longer the bottleneck it was.** 122 clusters hold 202 variants,
against 28 and 81 before, and only 87 of the 202 are homeodomain (was 54 of 81): Myb 16,
bHLH 15, forkhead 13, AP2 9, zf-C4 7. `T11` and `T14` would deepen it further.

That is what makes the project's central question askable at all — **does sensitivity to
single-residue change transfer across folds?** With variants in homeodomain only it was
unanswerable; five non-homeodomain families now carry them. Still thin, and it resurfaces on
its own when modelling starts.

**The label noise floor is measured and it is not small** (`reports/overlap.md`, 2026-08-17).
47 domains are stored by two sources; the median pair agrees on **45.8%** of the 8-mers either
called positive, against 70-72% for replicates within one source. Any model that reproduces
held-out positives much past that is reproducing a laboratory. One pair — `C:LIN14B:NAP` —
agrees on nothing at all, which is `D4`.

**And the cutoff behind those labels has now been swept** (`reports/threshold_review.md`,
2026-08-18). `E >= 0.45` is the best *absolute* cutoff there is — every alternative agrees
less — but a rank-matched rule at the same stringency agrees 8-12 points better, and the
fixed cutoff silently switches `Cell09` off entirely. `T20` and `T21`.

**Timings are measured, from the full rebuild of 2026-08-17.** `build_dataset.py --all` is
**7 min 28 s** — the "4-5 min" figure carried here since 2026-08-14 was optimistic — and the
whole pipeline (build, cluster, protein table, reports, overlap, audit) is **13 minutes**.
Per-step timings and the order to run them in are `docs/METHODS.md` §10.2. **The library must
be pressed once per machine — `python scripts/press_pfam.py`** — or every source pays 19.5 s to
re-read 2.2 GB of text.

**The plan, in order:**

| # | step | state |
|---|---|---|
| 1 | **Extend PBM coverage** beyond UniPROBE | **`T11` screened 2026-08-18, decision open** — the route is `T19`, now built and measured: rescoring reproduces the ranking but needs one recalibrated cutoff, which makes it the same decision as `T20`. Plus `T14`; `T2`/`T2b` dropped 2026-08-17 |
| 2 | **Deepen the protein axis in what we already hold** | **done 2026-08-17** — clustering by sequence distance, and the cluster inventory (`T3`) with it |
| 3 | **Merge**: single table, splits, NN baseline | after step 1 — the overlap report is already done (`T4`) |
| 4 | **Modelling** | after step 3 |

A research sweep on 2026-08-14 surveyed the literature for PBM sources with designed protein
variation. Its verdict, after screening against the admission policy: **Kock 2024 is the one
large find** (`T11`), two of its other leads were already parsed, and one fails condition 2
outright. See [`docs/DECISIONS.md`](docs/DECISIONS.md) §9.

---

## Open tasks

### T18 — Reconnect the protein-side table to canonical sequences
It covers **1,224 of 1,335 domains**. `build_protein_table.py` locates a domain by substring
within its construct, and a canonical sequence extended from a reference is no longer a
substring of the construct it came from. Use `canonical.place` instead of `dbd in construct`.
Note the construct-architecture covariate (`T13`, closed 2026-08-17) rides on the same
columns: flank length is `dbd_start` and `len(construct_seq) - dbd_end`, so whatever fixes the
placement fixes that too.

### T11 — Kock et al. 2024: screened 2026-08-18, **decision open**
Steps 1-3 are done and written up in
[`reports/kock2024_screening.md`](reports/kock2024_screening.md), with the measurements behind
every claim. Nothing was downloaded into `data/raw/` and no parser was written: **whether and
how this source enters is the owner's call.** In brief —

- **Identifiers verified.** Nat Commun 15:3110, `10.1038/s41467-024-47396-0`, PMC11006913;
  data at **GEO `GSE233827`** (missing from the old note, and the important one) plus Harvard
  Dataverse `10.7910/DVN/FDQHCF`, CC0. The reported Zenodo DOI is **not data** — it is the
  `upbm` R package.
- **Yield: 87 of 93 clones admitted, 67 new domains.** The whole PAX4 series fails
  `mixed_families` (PAX + homeodomain in one construct). Variants 202 → ~269, and 13 singleton
  clusters become variant-bearing. Homeodomain share 29.7% → 33.3%, not the ~58% projected
  before CIS-BP landed.
- **DNA axis is exact** — 32,896 8-mers per allele, matching our stored `dna_seq` strings
  directly. Step 3 closed.
- **The published scores do not fit our scale.** No E-scores exist anywhere in the deposit.
  `affinityEstimate` ranks agree with our E-scores (Spearman 0.83-0.93) but have **no
  transferable scale**: the value at our `E = 0.45` boundary spans 10.82-12.14 across four
  proteins while the `E = 0.35` boundary spans 10.08-11.16. A rate-calibrated label agrees
  with our BAR15A labels at Jaccard 0.51-0.79 — inside the measured noise floor. Four options
  costed in §8 of the report; **the recommendation is `T19`**, which makes the scale question
  disappear rather than answering it.
- **37 of its 131 records are Barrera 2016 reprocessed** — our own `BAR15A` arrays, not a
  second experiment. Parse `This study` only, or `reports/overlap.md` measures a source
  agreeing with itself.

Two things fall out regardless of the decision: `BAR15A:SIX6` carries the wrong full-length
accession (`Q6P051`, a 305 aa fragment entry, not `O95475`), which is why Kock's SIX6
numbering sits 59 residues off ours; and the same-raw-data-two-pipelines agreement in §6 is an
independent read on the noise floor that belongs in `docs/METHODS.md`.

### T19 — E-scores from raw arrays: built and measured, **one decision left**
Done 2026-08-18 and written up in
[`reports/escore_recomputation.md`](reports/escore_recomputation.md). The pipeline is
[`src/snp2prot/rawpbm.py`](src/snp2prot/rawpbm.py) with
[`scripts/validate_escores.py`](scripts/validate_escores.py) and 9 tests; ~2 s per array.
Validated on **32 arrays / 14 `BAR15A` alleles** rescored from Barrera's own raw scans in
`GSE233827` and compared with the E-scores we already hold.

**The ranking reproduces; the values do not.** Spearman 0.78-0.91 over the 8-mers a protein
binds, rank-matched positive sets at Jaccard 0.57-0.77 (median 0.64) — the band between the
corpus's cross-source floor (0.458) and its within-source replicates (0.70-0.72). But
`E >= 0.45` selects 0-32 8-mers where the stored values select 129-204. The cutoff that
matches the stored positive **count** is 0.312-0.343 across all 11 live alleles, median 0.321.

So `T19`'s premise was half right: no new `score_type`, no new scale, no new
`neg_provenance` — but **one recalibrated cutoff**, near 0.32-0.35, derived from the `BAR15A`
overlap rather than chosen. Still far better than the `affinityQ` route (`T11`), and it
resolves `T20`'s dilemma: the three dead `ARX` variants top out at 0.33-0.43 and yield 0-2
positives, so an absolute floor on a recomputed E-score does what a rank rule cannot.

Three findings worth keeping: **keep saturated spots** (260 of one array's top 1,000 probes
are saturated; masking them decapitates the statistic), **normalise Cy3 as a sequence-model
residual, not a ratio** (raw division drops agreement 0.489 → 0.123), and join by
`(Column, Row)`, never `probeID`. Spatial detrending was tried and is not worth a knob.

Two things do not work and bound the method: **more replicates do not restore the range** (five
`HOXD13_REF` arrays give 11 positives against 144), most likely because Barrera's masliner
multi-power scan stitching cannot be reproduced from one scan per array; and **array QC is
unsolved** — upper-tail width separates good arrays from bad, but a dead variant is
indistinguishable from a bad array by that measure, so nothing filters on it yet.

### T20 — Fixed cutoff or rank-matched rule?
**Measured, not speculative** — [`reports/threshold_review.md`](reports/threshold_review.md).
`E >= 0.45` is the **best absolute cutoff available**: median cross-source Jaccard 0.520,
against 0.509 at 0.40 and 0.495 at 0.35. Loosening it does not recover agreement.

But the absolute-cutoff *rule* is what costs the agreement. At the same stringency, taking
each experiment's **top 100** gives 0.600 and top 50 gives 0.639 — 8 to 12 points better —
and removes the 1.6x median positive-count asymmetry between two labs measuring one protein
(worst case `Hoxa2`, 165 against 17).

The catch is `T21`'s other half: a naive rank rule would manufacture ~100 positives for each
of `BAR15A`'s 18 dead variants, which is precisely the signal the corpus exists to carry. A
rule that survives both failure modes is a **rank cap under an absolute floor**, which is two
parameters, which by convention invalidates the dataset and every report. Not to be decided
in passing — but note `T19` removes the pressure to decide it at all.

### T21 — `Cell09` and `LIU18B` are label-dead at the current cutoff
`Cell09`'s per-domain **maximum** E-score has a median of **0.428**: 12 of its 17 domains
never reach 0.45 at their single best 8-mer, so the source contributes 559,232 rows that say
only "does not bind". Its positive rate is 0.0048% against 0.14-0.49% for healthy sources.
`LIU18B` has **one** positive 8-mer in the entire source.

This is a compressed E-score distribution meeting a cutoff calibrated on other arrays — not
17 proteins that bind nothing. Distinct from `BAR15A`'s 18 zero-positive domains, which have a
healthy distribution with a dead tail and are probably real. Options: exclude the two sources,
re-binarize them at their own scale (`T20`), or accept them as negatives-only and say so in
`METHODS`. 54 of 1,383 domain records (3.9%) have no positive at all; these two sources are
most of the unexplained part.

### T22 — `BAR15A:SIX6` carries the wrong full-length protein
`protein_id` is `Q6P051`, a 305 aa TrEMBL *"SIX6 protein (Fragment)"* from a cDNA clone, not
reviewed `SIX6_HUMAN` `O95475` (246 aa). Consequence: our `full_seq` frame runs **59 residues**
off the literature's numbering, which is how it was caught — Kock's `SIX6-H141N` lands at 200
in our frame. Every other gene checked exactly (17 variants across CRX, HESX1, MSX2, NKX2-5,
PITX2, PROP1, VENTX, all delta 0). Affects the protein table only; `dbd_seq` is unchanged.
Fix the accession, rebuild the protein table, and check whether `references.py` resolved it
from a stale mapping — if so, other TrEMBL fragment entries may be doing the same thing.

### T23 — Do stored domains match the reference proteome inside the padding?
Kock's `HOXD13` clone carries `T` where `P35453` carries `D`, at protein position 261 — six
residues N-terminal of the Pfam start, so **inside** our ±10 padding. Our `BAR15A` domain
carries the same `T`, so both labs cloned from one source and the stored `dbd_seq` differs
from the reference proteome at that residue. Probably a clone polymorphism and probably
harmless, but it has never been audited: for every domain with a `full_seq`, count residues
where the padded window disagrees with the reference. A systematic offset would be an
isoform bug; scattered singletons are clone reality. Cheap, and it protects `mut_positions`.

### T14 — Check the Liu 2018 Dryad deposit
`LIU18B` is parsed, but UniPROBE published only **2 constructs** — `AncBcd` and `AncZB`, both
singletons. The paper (Liu, Onal, Datta, Rogers et al. 2018, *eLife* 7:e34594) reports a Q50K
substitution accounting for the evolution of Bicoid specificity, plus epistatically
interacting mutations, which implies more assayed constructs than UniPROBE holds. A deposit
is reported at Dryad `10.5061/dryad.pm3g4r3` with a companion GitHub repo — unverified.

Small, but it is ancestral reconstruction, so every sequence is stated explicitly and there is
no accession chasing. Reconstructed-ancestor series are also the one place where designed
protein-axis depth exists outside homeodomain point mutants.

### T15 — One domain is two different proteins, depending on the source
Left over from `T5b`, which fixed how organism names are *written* and could not fix this. The
domain under `C:Cell08:Tlx2` is byte-identical in `Cell08` and `weirauch2014`, and is stored
as *Mus musculus* by one and *Homo sapiens* by the other. A byte-identical homeodomain across
mouse and human is entirely possible, so this may be two correct records — but their labels
disagree (Jaccard 0.097, `reports/overlap.md`), and one deposit having the wrong protein would
explain both facts at once. Related: `D4`.

**`T11` would multiply this by ten and add a worse variant of it.** 9 of Kock's human domains
are byte-identical to `Cell08` mouse domains, so each becomes another species disagreement.
One is not a species question at all: **`HOXC9-K195R`, a human disease variant, is
byte-identical to `Cell08`'s mouse `Hoxc9` wild type.** One `dbd_seq`, two identities, two
labels from two assays — a sequence model would see identical input with different labels,
which is the failure condition 3 of the admission policy exists to prevent. Worth deciding the
general rule here before it arrives: does the corpus store one row per sequence or one per
(sequence, protein)?

### T24 — Parse-time gotchas held for `T11`, if it is admitted
Found during screening, all verified, none of them blocking on their own:
- **PAX4's whole series is rejected** (`PAX` + `Homeodomain` in one 234 aa construct,
  6 alleles). Expected and correct; worth stating in the parser's rejection report rather
  than discovering later.
- **3 HOXB9 clones carry only 8 aa of C-flank**, so the ±10 window cannot be filled from the
  construct — `canonical` must extend from the UniProt reference for those.
- **`NKX2-6-REF` appears twice** in Supplementary Data 4 (`NKX2-6` and `NKX2-6series2`),
  one clone assayed in two batches. Two records, one sequence.
- **SD4 is already replicate-collapsed** — one value per (allele, 8-mer). `per_experiment`
  binarization and `reconcile_replicates` have nothing to act on, and no within-source
  replicate agreement can be computed. 4 alleles rest on a single array. `T19` removes this
  too: the raw `.gpr`s are per-array.
- **Supplementary Data 3 is legacy `.xls`** and needs `xlrd`; the project reads `.xlsx` via
  `openpyxl` and has no `.xls` reader. One new dev dependency for one 84 KB file.
- **The Barrera-set rows must not enter the training table** — 37 of 131 records are our own
  `BAR15A` arrays reprocessed, 1,217,152 rows of a source agreeing with itself, which would
  corrupt the one report that measures cross-source agreement.

### T25 — Free terminal gaps can call two unrelated domains near-identical
`align.edit_profile` with `free_end_gaps: true` reports **3 edits** between a 60 aa
homeodomain-like window from `GD09`'s NSY-7 and `C:weirauch2014:lim-4` — because the global
aligner pushed everything into free terminal gaps and scored a 5-residue overlap. The setting
is right and load-bearing (`ROG18A`'s clipped padding, `docs/DECISIONS.md`), but it has no
floor: nothing requires the aligned region to be a meaningful fraction of either sequence.

Nothing in the corpus is known to be affected — clustering compares within a family and real
domains share a full core — but a heavily clipped short domain meeting a long representative
of its own family is the shape that would trigger it. Fix is a minimum-overlap guard
(reject the alignment if the ungapped overlap is below, say, 60% of the shorter sequence),
plus a sweep of existing clusters for any member joined on a short overlap.

### T26 — `.gitignore` silently drops the docs it promises to keep
`data/raw/*` excludes the *directory*, so git never descends into it and the
`!data/**/README.md` / `!data/**/HOWTO.md` negations below can never fire — a directory
excluded at one level cannot have its contents re-included. Two files are affected today and
neither is in the repo: `data/raw/weirauch2014/README.md` (has been missing all along) and
`data/external/pbm_design/README.md` (written 2026-08-18, and the only place the
join-by-position finding is recorded outside `reports/`).

Fix is to ignore files rather than directories, per data root:

```gitignore
data/raw/**
!data/raw/**/
!data/**/.gitkeep
!data/**/README.md
!data/**/HOWTO.md
```

Not applied unilaterally: it changes what the next commit picks up, and that is the owner's
to see. Rule 8 also assumes `data/raw/<source>/HOWTO.md` is tracked evidence — right now it
would not be.

### T9 — Parallelise the build
Deferred deliberately after the T6 work (see [`docs/DECISIONS.md`](docs/DECISIONS.md) §7):
the serial fixes landed first because they were worth more at a fraction of the risk. This is
the remaining headroom, roughly **4 min → 1.5-2 min**.

Measured, so the design is not guesswork:
- **Threads are useless here — 1.0x at 8, 16 and 24 workers.** pandas' CSV parser holds the
  GIL. Anyone reaching for `ThreadPoolExecutor` on this workload gets nothing.
- **Processes gave 7.6x** on 48 SCI09 files (9.7 s → 1.3 s at 16 workers), saturating there.

Two levels are needed, not one:
1. **Across sources.** Architecturally free — rule 7 already guarantees parsers share no
   state and each writes its own file.
2. **Within a source.** Unavoidable, because **SCI09 alone is 2.3 GB of the 3.2 GB read**.
   Parallelising only across sources leaves SCI09 as the critical path and the build cannot
   go below what SCI09 costs on its own. `reconcile_replicates` reads members serially.

Constraints to respect:
- **Memory, not cores, is the limit.** Cell08's frame is ~2 GB in memory and each worker
  holding the Pfam profile block adds ~1.6 GB. At 16 workers that heads past 64 GB. Size at
  **8 workers** and measure before going wider — that is already near the 7.6x ceiling.
- `zipfile.ZipFile` handles cannot be shared; each worker opens its own.
- **Per-source rejection reasons are load-bearing** (`N3` exists because of them). A process
  pool must not blur which source failed and why — `build_dataset.py` currently reports that
  per source and must continue to.

### T10 — Give `Pfam-A.hmm` a PROVENANCE row
The nine superseded per-family HMMs each have one; the 2.2 GB full library that the entire
admission policy now rests on has none. Rule 3 makes this publication evidence, not a
convenience log. Needs the download URL, release, size and sha256. The pressed `.h3*` files
are derived and need no row of their own, but the row should say the library gets pressed.

---

## Open decisions

### D1 — Is 70% full-length coverage enough?
**Before embedding work, and CIS-BP made it sharper.** The protein table maps **333 of 1,335
domains (25%)** onto a canonical UniProt sequence; **439 (33%)** have a full sequence at all.
It was 71% before CIS-BP: Table S6 publishes no UniProt accession, so none of its 884 domains
resolve to one today. Resolving them from gene plus species is possible but is a lookup this
project has not yet had to do, and for 124 organisms it will not be clean. The remainder are clone constructs
differing from the canonical isoform, or non-model species with no clean mapping.

Domain-level and construct-level embeddings work for all 1,335. Full-protein works for about
a quarter, and the missing three quarters are almost entirely one source.
The decision is whether that asymmetry is acceptable or whether full-length becomes a filter.

### D2 — Drop `MAR17A:Esrrb`?
**Before structure work.** It carries 2 unresolved `X` residues in an 89 aa zf-C4 domain — the
only such domain in the corpus, 1 of 1,335. Harmless to a sequence embedder; a genuine problem
for structure prediction and any 3D embedder. Drop it, or carry it and exclude it at
structure-generation time.

### D3 — Which member of a merged paralogue pair is the reference?
Referenced by [`docs/DECISIONS.md`](docs/DECISIONS.md) §2 but never written down here.
Clustering by sequence distance merges two natural paralogues that construct lineage kept
apart, and `mut_positions` is expressed in the reference's frame — so the merge has to pick
one, and neither has a claim. The **threshold** half of this decision is settled (5 edits,
owner); the reference-choice half is not.

### D4 — `C:LIN14B:NAP` — two deposits of ANAC092 that share no positive call
**Surfaced by `T4`'s agreement check, 2026-08-17.** *Arabidopsis* ANAC092 / O49255 is stored
by both `LIN14B` and `weirauch2014` with a byte-identical 144 aa domain. Their positive sets
are **disjoint** — 12 against 126, zero shared — at a rank correlation of 0.068 across all
32,896 8-mers. Every other replicate pair in the corpus reaches at least 0.097 Jaccard and
0.200 rho.

This is outside the range noise explains, so one of the two measurements is probably attached
to the wrong construct. It is one domain of 1,335 and no rule in the admission policy catches
it, since both deposits are internally consistent. Options: exclude the domain, keep one
deposit and say which, or keep both and let the disagreement stand as measured. Not resolved
unilaterally — but note that whatever the answer, both rows currently sit in **one cluster**,
so a split cannot separate them and a model sees one protein with two contradictory profiles.


---

## Notes

### N7 — Cluster size is a one-file lookup now
`data/interim/clusters/clusters.parquet`, one row per cluster, written by `build_clusters.py`.
Read it with `snp2prot.clusters.load` / `ids_with_at_least(n)` / `select(df, n)` rather than
grouping the row tables — and note that **size means distinct canonical domains**, not rows and
not constructs, so a domain assayed by two labs counts once. Today: 122 clusters hold more than
one domain, 34 hold three or more, 14 hold five or more.

### N1 — The row count is an exact invariant
`45,495,168 = 1,383 x 32,896`, where 1,383 is 1,335 distinct domains plus 48 measured by more
than one source. Canonicalisation raised that overlap: domains that used to differ only by how
much flank a lab cloned are now one string, so they are recognised as the same measurement
made twice — which is what `T4`'s agreement check needs. Every construct contributes exactly one full 8-mer table. If a rebuild's row count is
not a clean multiple of 32,896, something dropped or duplicated rows — the fastest single
sanity check available.

### N6 — A single-source rebuild must be followed by the cluster pass
`build_dataset.py --source X` writes `wt_id` in its lineage form (`Cell09:HLH-1`), while every
other source on disk carries the cluster form (`C:NAR11:HLH-1`) that `build_clusters.py`
rewrote. Rebuilding one source therefore desynchronises it from the corpus silently — the
table validates, the row count is right, and the domain simply leaves its cluster.

**Always run `scripts/build_clusters.py` after any single-source rebuild.** It is corpus-wide
and idempotent (it strips the `C:` prefix to recover the lineage name underneath), takes
**79-84 s**, and reports the domain and cluster counts so a mismatch is visible: 1,335 / 1,133
/ 122 multi-member / largest 8 as of 2026-08-17.

**That is necessary but it is not sufficient, and the ROG18A rebuild of 2026-08-17 showed
why.** Rebuilding a source replaces what is on disk with what the *current* parser produces,
and the rest of the corpus still holds what the parser produced whenever it was last built. If
the two differ, the mixture validates cleanly and the row-count invariant still holds — the
only visible symptom is a report diff. Rebuilding `ROG18A` regrouped its 15 chimeras from
`{8, 5, 1, 1}` clusters to three pairs and nine singletons. **The new grouping is the correct
one**: the pairwise edit distances were checked directly against `align.edit_profile`, and
under the 5-edit rule exactly three pairs qualify. The committed table simply predated a
parser change. `docs/DECISIONS.md` `#30` — "ROG18A's engineered chimeras cluster with their
parents" — describes the old grouping and needs revisiting.

The corpus is currently a mixture: `Cell09`, `NAR11`, `ROG18A`, `LIU18B` and `weirauch2014`
were rebuilt on 2026-08-17, the other 14 sources were not. **`build_dataset.py --all` followed
by `build_clusters.py` is the way to make it self-consistent again**, and it is the owner's to
run (rule 10).

### N3 — Checklist for any new PBM source
Every one of these was caught by a distribution looking wrong, never by an error or a failing
validator. Check each before trusting a new accession:

- no header row on the 8-mer tables;
- a column count that differs from the last source (7, 20, 4 have all occurred);
- the E-score in a different column, or absent entirely;
- **one gene folder holding several distinct engineered proteins** — the costliest defect so
  far, it silently averaged different proteins as replicates;
- a bare or plasmid-named insert label, which collapses every gene in the source into one
  cluster;
- a flat archive with no gene directories;
- a table truncated to the enriched end, or *n* rows short of 32,896.

Full account in [`docs/DECISIONS.md`](docs/DECISIONS.md). Run
`scripts/audit_sources.py` (~4 min) after any new source lands — it sweeps every one of these.

### N4 — 27 families are large enough to hold out
**Largely fixed by CIS-BP: 27 families now hold 10 or more domains, against 8 before.**
Homeodomain 427, HLH 100, bZIP 79, Zn_clus 74, Forkhead 69, Myb_DNA-binding 59, zf-C4 58,
AP2 41, GATA 29, Ets 26, HMG_box 26, NAM 26, WRKY 23, TCP 20, Zn_ribbon_Dof 17, ARID 17, and
eleven more. Leave-one-family-out is now a real test rather than a handful of proteins.

### N5 — Where the `dna_len` warning went
The brief warns never to pool sources without checking `dna_len`, because PBM's 8 bp, B1H's
9 bp and SNP-SELEX's 19 bp each carry a different positive rate, so length alone leaks assay
identity and with it the label prior. **The pure-PBM decision dissolves this**: every stored
row is 8 bp. It returns the moment any non-PBM assay is admitted, so the warning stays in
`CLAUDE.md` rather than being deleted.
