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

`build_dataset.py --all` was ~25 minutes and should now be roughly **4-5**, after the T6 work
of 2026-08-14 (`docs/DECISIONS.md` §7). **The library must be pressed once per machine —
`python scripts/press_pfam.py`** — or every source pays 19.5 s to re-read 2.2 GB of text.

**The plan, in order:**

| # | step | state |
|---|---|---|
| 1 | **Extend PBM coverage** beyond UniPROBE | **closes with `T11`** (Kock), plus `T14`; `T2`/`T2b` dropped 2026-08-17 |
| 2 | **Deepen the protein axis in what we already hold** | **done 2026-08-17** — clustering by sequence distance shipped |
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

### T3 — Make cluster size cheap to filter on
Whether to require ≥5 variants per DBD is a **training-time** choice, not a dataset one
(`D-2026-08-14-clusters`). What the dataset owes the modeller is the ability to select on it
without a full scan: today it means grouping 16.6M rows by `wt_id` and counting distinct
`dbd_seq`.

**Smaller since 2026-08-17:** `build_clusters.py` already computes the whole size distribution
and writes it to `reports/clusters.md`, but only as prose — nothing machine-readable, and
`snp2prot.splits` is still three `NotImplementedError` stubs.

Design question inside the task: a stored `cluster_size` column is a Parquet predicate
pushdown but touches `schema.py` and every parser's output; a small `wt_id -> size` side table
or a helper in `snp2prot.splits` costs nothing and stays out of the 22-column schema. Second
is likely right, and the cluster rebuild is now the obvious place to emit it — flag before
implementing.

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

### N4 — Only 8 families are large enough to hold out
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
