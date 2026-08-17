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

`build_dataset.py --all` was ~25 minutes and should now be roughly **4-5**, after the T6 work
of 2026-08-14 (`docs/DECISIONS.md` §7). **The library must be pressed once per machine —
`python scripts/press_pfam.py`** — or every source pays 19.5 s to re-read 2.2 GB of text.

**The plan, in order:**

| # | step | state |
|---|---|---|
| 1 | **Extend PBM coverage** beyond UniPROBE | active — `T1` (CIS-BP), `T11` (Kock), `T14`, `T2` |
| 2 | **Deepen the protein axis in what we already hold** | `T12` — cluster by sequence distance |
| 3 | **Merge**: single table, overlap report, splits, NN baseline | after steps 1-2 |
| 4 | **Tier 4 held-out sets** in `data/testsets/` | wanted — `T7` |
| 5 | **Modelling** | after step 3 |

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
Also the natural home for `T13`'s architecture covariates, since both touch the same loop.

### T2b — Resolve UniProt accessions for the CIS-BP domains
Spun out of `T1`, which is otherwise done. Table S6 gives gene name and species but no UniProt
accession, so 884 of 1,364 domains have no full-length sequence and `D1` went from 71%
coverage to 26%. Everything domain-level and construct-level works; only full-protein does
not. Resolving 124 organisms' gene names against UniProt will not be clean, and a wrong
mapping is worse than a missing one — the same reasoning that left `PP15` unparsed.

### T2 — Survey other PBM deposits
Individual GEO / ArrayExpress submissions and paper supplements. Lower yield per unit effort
than `T1` and subject to the same construct-sequence caveat. Do after `T1` resolves.

### T11 — Kock et al. 2024: acquire, screen, parse
The best variant source identified so far, and the reason the corpus can double its
multi-member clusters. **All identifiers below come from a research session and are
unverified — confirm each against the source before use, and do not construct URLs from
them** (rule 1).

Reported as: Kock et al. 2024, *Nat Commun* 15:3110 — 30 homeodomain allelic series, 92 HD
missense variants plus their 30 reference alleles, 122 alleles total, on the same "all
10-mer" universal array design as ROG18A (AMADID #030236). Construct provenance is stated as
DBD plus 15 aa flanks, gene-synthesized, GST fusion, Sanger-verified, with sequences in
Supplementary Data 3. Deposits reported at Harvard Dataverse `10.7910/DVN/FDQHCF` and Zenodo
`10.5281/zenodo.10460649`.

**Step 1, before anything else: what is actually in the deposit?** The reported deliverable is
upbm affinity/contrast/specificity *Q-values*, not Wilcoxon E-scores. If probe-level or
E-score data is there, this is an ordinary parse. If only Q-values, it enters as a new
`score_type` (decided 2026-08-14, `docs/DECISIONS.md` §1) with its own cutoffs, and that is a
`configs/thresholds.yaml` change with the blast radius that implies.

**Step 2: the overlap is larger than it looks.** Three of the four genes in its epistasis
probe already carry Barrera variant series here — `BAR15A:PROP1` (2 variants), `BAR15A:VENTX`
(2), `BAR15A:SIX6` (2); only PAX4 is absent. Two consequences:
- Kock's Arg→Gln at canonical HD position 31 in PROP1 and VENTX may be **the same
  substitution Barrera assayed**, on a different array architecture — a cross-lab replicate
  of a *variant*, which nothing in the corpus currently provides. Our note puts Barrera's
  PROP1 `R99Q` at HD position 30, one off, which is what a numbering-convention difference
  looks like. It may also bear on the `R112Q` question closed as unresolvable
  (`docs/DECISIONS.md` §4).
- **Coordinate frames must be reconciled, and this is the real work.** `mut_positions` is in
  the reference's frame, and the reference is the padded envelope *clipped by the construct*.
  Kock's DBD+15aa construct and Barrera's are different constructs of one protein, so merging
  them into one cluster means two frames for one `wt_id`. Rebasing already refuses to cross
  an indel by design.

**Step 3: verify the DNA axis is exact.** ROG18A used the same array design and parsed clean
at the full 32,896 8-mers, so there is precedent — but GD09 was excluded for being *one*
8-mer short, and that bar does not move.

### T13 — Record construct architecture in the protein-side table
**Decided 2026-08-14: side table, not schema columns.** Flank length, affinity tag and
expression system are perfectly confounded with source study, so a model trained across
sources can learn "which lab made this" instead of "what does this residue do".

This is not only a future risk. `dbd_seq` is the padded envelope **clipped where the construct
ends**, so flank length is already partly readable off `dbd_seq` length — ROG18A's bare domains
run 83-85 aa against 95-105 aa for the padded ones, purely from how much padding fitted. A
DBD+15aa construct admits the full 10 aa of padding; a bare-DBD construct does not.

Goes in `data/interim/proteins/` via `scripts/build_protein_table.py`, keyed by construct, so
the 22-column schema stays frozen and no parser is rewritten. Today the covariate barely
varies — all 18 sources are UniPROBE PBM — which is exactly why it is cheap to add now and
expensive once Kock's DBD+15aa GST constructs and CIS-BP's DBD+50aa T7-GST are both in.
Note `dbd_source` does not cover this: it records the annotation *method*
(`pfam_hmmer_padded`), not what was on the array.

### T14 — Check the Liu 2018 Dryad deposit
`LIU18B` is parsed, but UniPROBE published only **2 constructs** — `AncBcd` and `AncZB`, both
singletons. The paper (Liu, Onal, Datta, Rogers et al. 2018, *eLife* 7:e34594) reports a Q50K
substitution accounting for the evolution of Bicoid specificity, plus epistatically
interacting mutations, which implies more assayed constructs than UniPROBE holds. A deposit
is reported at Dryad `10.5061/dryad.pm3g4r3` with a companion GitHub repo — unverified.

Small, but it is ancestral reconstruction, so every sequence is stated explicitly and there is
no accession chasing. Reconstructed-ancestor series are also the one place where designed
protein-axis depth exists outside homeodomain point mutants.

### T3 — Make cluster size cheap to filter on
Whether to require ≥5 variants per DBD is a **training-time** choice, not a dataset one
(`D-2026-08-14-clusters`). What the dataset owes the modeller is the ability to select on it
without a full scan: today it means grouping 16.6M rows by `wt_id` and counting distinct
`dbd_seq`.

Design question inside the task: a stored `cluster_size` column is a Parquet predicate
pushdown but touches `schema.py` and every parser's output; a small `wt_id -> size` side table
or a helper in `snp2prot.splits` costs nothing and stays out of the 22-column schema. Second
is likely right — flag before implementing.

### T4 — Reconcile the 17 cross-source identical domains
Seventeen `dbd_seq` values are byte-identical across two accessions, stored once per source —
**559,232 rows (17 x 32,896) of exact protein-axis duplication**.

```
BAR15A:PITX2  = Cell08:Pitx2      EMBO10:Spdef  = SCI09:Spdef
BAR15A:PHOX2B = Cell08:Phox2b     EMBO10:Elf3   = SCI09:Elf3
BAR15A:VAX2   = Cell08:Vax2       EMBO10:Sfpi1  = SCI09:Sfpi1
BAR15A:HOXD13 = Cell08:Hoxd13     EMBO10:Ehf    = SCI09:Ehf
BAR15A:ARX    = Cell08:Arx        EMBO10:Gabpa  = SCI09:Gabpa
BAR15A:MSX2   = Cell08:Msx2       MAR17A:Tfap2a = SCI09:Tcfap2a
BAR15A:SIX6   = SCI09:Six6        Cell08:Msx3   = MAR17A:Msx3
BAR15A:FOXC1  = MAR17A:Foxc1      MAR17A:Foxm1  = PNAS13:Foxm1
PNAS13:FOXN2  = ROG18A:FoxN2
```

Two jobs at merge time, neither optional:
1. **They must share a split group.** Identical protein, identical 8-mers, two sources — any
   split that separates them leaks the test set into training. This is the lookup-table
   failure in its purest available form.
2. **Check label agreement.** Two labs measured the same domain against the same 32,896
   8-mers. Disagreement is a free measurement of the label noise floor, and the single most
   informative number the corpus can produce about its own reliability. Belongs in
   `reports/overlap.md`.

Since the 2026-08-14 decision to cluster by sequence distance (`T12`), job 1 falls out for
free: these 17 pairs are distance 0, so any distance rule merges them. What stays specific to
this task is job 2, the agreement check — merging them makes it *necessary*, because a
cluster that silently contains two labs' readings of one protein must be known to.

### T5 — Resolve the 18 `$species` placeholders
Eighteen domains carry `species` as the literal string `$species` — an unsubstituted template
placeholder on some UniPROBE detail pages, not a fault on our side. Metadata only, affects no
label. Resolve from UniProt (which carries the organism) or normalise to null. Related: the
field mixes `C. elegans` with full binomials, so organism names are not consistently
formatted.

### T7 — Screen Tier 4 per construct
Tier 4 is bHLH dimers, and it needs a per-construct call rather than a wholesale one:
- a heterodimer is two chains forming one binding unit, which is exactly what admission
  condition 1 rejects (see the 22 complexes already excluded on those grounds);
- MAX's substitutions are described as "in and around" the DBD, so condition 3 — the variation
  lies inside the stored region — has to be checked for **every variant**, not assumed.

Kd / dG stay continuous at parse time and are thresholded only at evaluation, with the cutoff
swept (`tier4.binarize_at_parse: false`).

### T8 — Retire the dead config blocks
`configs/thresholds.yaml` still carries `b1h:` and `snp_selex:` blocks with `null` cutoffs and
TODOs for phases that no longer exist. Harmless but misleading. Remove or comment as dropped —
check `snp2prot.thresholds` does not require the keys before deleting.

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
**Before embedding work, and CIS-BP made it sharper.** The protein table maps **354 of 1,364
domains (26%)** onto a canonical UniProt sequence; **462 (34%)** have a full sequence at all.
It was 71% before CIS-BP: Table S6 publishes no UniProt accession, so none of its 884 domains
resolve to one today. Resolving them from gene plus species is possible but is a lookup this
project has not yet had to do, and for 124 organisms it will not be clean. The remainder are clone constructs
differing from the canonical isoform, or non-model species with no clean mapping.

Domain-level and construct-level embeddings work for all 1,364. Full-protein works for about
a quarter, and the missing three quarters are almost entirely one source.
The decision is whether that asymmetry is acceptable or whether full-length becomes a filter.

### D2 — Drop `MAR17A:Esrrb`?
**Before structure work.** It carries 2 unresolved `X` residues in an 89 aa zf-C4 domain — the
only such domain in the corpus, 1 of 1,364. Harmless to a sequence embedder; a genuine problem
for structure prediction and any 3D embedder. Drop it, or carry it and exclude it at
structure-generation time.


---

## Notes

### N1 — The row count is an exact invariant
`45,495,168 = 1,383 x 32,896`, where 1,383 is 1,335 distinct domains plus 48 measured by more
than one source. Canonicalisation raised that overlap: domains that used to differ only by how
much flank a lab cloned are now one string, so they are recognised as the same measurement
made twice — which is what `T4`'s agreement check needs. Every construct contributes exactly one full 8-mer table. If a rebuild's row count is
not a clean multiple of 32,896, something dropped or duplicated rows — the fastest single
sanity check available.

### N2 — Protein-axis depth, and what it can now answer
81 point variants across 28 clusters; **1,293 of 1,321 clusters (98%) hold a single domain**.
By family: Homeodomain 54, Forkhead 18, zf-C4 5, HLH 3, PAX 1.

This matters for the project's central question — does sensitivity to single-residue change
transfer across folds? It used to be unanswerable, with variants in homeodomain only. Three
non-homeodomain families now carry variants, and Forkhead's 18 came from fixing a parser bug
rather than from new data. Still thin, but no longer empty. Accepted as-is; it resurfaces on
its own when modelling starts.

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

### N4 — Only 8 families are large enough to hold out
**Largely fixed by CIS-BP: 27 families now hold 10 or more domains, against 8 before.**
Homeodomain 440, HLH 103, bZIP 79, Zn_clus 74, Forkhead 72, zf-C4 62, Myb_DNA-binding 60,
AP2 44, GATA 29, Ets 26, HMG_box 26, NAM 26, WRKY 23, TCP 21, Zn_ribbon_Dof 18, ARID 17, and
eleven more. Leave-one-family-out is now a real test rather than a handful of proteins.

### N5 — Where the `dna_len` warning went
The brief warns never to pool sources without checking `dna_len`, because PBM's 8 bp, B1H's
9 bp and SNP-SELEX's 19 bp each carry a different positive rate, so length alone leaks assay
identity and with it the label prior. **The pure-PBM decision dissolves this**: every stored
row is 8 bp. It returns the moment any non-PBM assay is admitted, so the warning stays in
`CLAUDE.md` rather than being deleted.
