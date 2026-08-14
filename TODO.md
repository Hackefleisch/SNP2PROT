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

18 UniPROBE accessions parsed and validator-clean: 16,645,376 rows, 489 domains in 425
clusters, 30 Pfam families, 24 organisms, every protein scored against the same 32,896 8-mers.

`build_dataset.py --all` was ~25 minutes and should now be roughly **4-5**, after the T6 work
of 2026-08-14 (`docs/DECISIONS.md` §7). **The library must be pressed once per machine —
`python scripts/press_pfam.py`** — or every source pays 19.5 s to re-read 2.2 GB of text.

**The plan, in order:**

| # | step | state |
|---|---|---|
| 1 | **Extend PBM coverage** beyond UniPROBE | active — `T1`, `T2` |
| 2 | **Merge**: single table, overlap report, splits, NN baseline | after step 1 |
| 3 | **Tier 4 held-out sets** in `data/testsets/` | wanted — `T7` |
| 4 | **Modelling** | after step 2 |

---

## Open tasks

### T1 — Survey CIS-BP for admissible PBM data
The largest identified expansion, and now the critical path.

CIS-BP build 3.10 reports **2,294 TFs with PBM data** (plus 97 under
`PBM,_CSA_and_or_DIP-chip`) against 4,989 TFs with any directly determined motif, across 741
species. Its bulk download page serves **E-scores, Z-scores and probe intensities** per
species — the same data type we use, not just PWMs — through a POST form to
`bulk_archive.php`. Overlap with our 779 UniPROBE constructs is unknown and probably large,
since CIS-BP aggregates published data including Bulyk-lab deposits.

**Answer this before anything else: does CIS-BP publish the assayed construct sequence?**
Our admission policy rests on annotating what was physically on the array. If CIS-BP supplies
only its own DBD annotation of the full-length protein, attributability fails condition 1 and
the data cannot be admitted as things stand. That one question decides between a large
expansion and a dead end.

Then work the checklist in `N3` — every format assumption UniPROBE has broken so far.

### T2 — Survey other PBM deposits
Individual GEO / ArrayExpress submissions and paper supplements. Lower yield per unit effort
than `T1` and subject to the same construct-sequence caveat. Do after `T1` resolves.

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
**Before embedding work.** The protein table maps **342 of 489 domains** onto a canonical
UniProt sequence; **450** have a full sequence at all. The remainder are clone constructs
differing from the canonical isoform, or non-model species with no clean mapping.

Domain-level embeddings work for all 489. Full-protein embeddings work for roughly two thirds.
The decision is whether that asymmetry is acceptable or whether full-length becomes a filter.

### D2 — Drop `MAR17A:Esrrb`?
**Before structure work.** It carries 2 unresolved `X` residues in an 89 aa zf-C4 domain — the
only such domain in the corpus, 1 of 489. Harmless to a sequence embedder; a genuine problem
for structure prediction and any 3D embedder. Drop it, or carry it and exclude it at
structure-generation time.

---

## Notes

### N1 — The row count is an exact invariant
`16,645,376 = 506 x 32,896`, where 506 = 489 distinct domains + the 17 cross-source duplicates
of `T4`. Every construct contributes exactly one full 8-mer table. If a rebuild's row count is
not a clean multiple of 32,896, something dropped or duplicated rows — the fastest single
sanity check available.

### N2 — Protein-axis depth, and what it can now answer
81 point variants across 28 clusters; **397 of 425 clusters (93%) hold a single domain**.
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
Homeodomain 234, Forkhead 54, HLH 35, zf-C4 28, Ets 23, HMG_box 21, Zn_clus 17, T-box 11. The
long tail below ~10 domains gives noise under leave-one-family-out. Report LOFO for the viable
eight and treat the rest as descriptive. Revisit when splits are designed.

### N5 — Where the `dna_len` warning went
The brief warns never to pool sources without checking `dna_len`, because PBM's 8 bp, B1H's
9 bp and SNP-SELEX's 19 bp each carry a different positive rate, so length alone leaks assay
identity and with it the label prior. **The pure-PBM decision dissolves this**: every stored
row is 8 bp. It returns the moment any non-PBM assay is admitted, so the warning stays in
`CLAUDE.md` rather than being deleted.
