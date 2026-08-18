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

19 PBM sources parsed and validator-clean: **45,593,856 rows, 1,338 domains in 1,165
clusters, 56 Pfam families, 131 organisms**, every protein scored against the same 32,896
8-mers.

`dbd_seq` is the **canonical domain** — envelope ± 10, construct-independent
(`docs/DECISIONS.md` §11) — and `wt_id` comes from CD-HIT clustering at 5 edits, not construct
lineage.

**The protein axis is thinner than it looked.** 108 clusters hold **173 variants**, and 84 of
them are homeodomain: forkhead 13, HLH 9, zf-C4 7, then MADF, RFX, TCP and SAND at 5 each,
Myb 4, bZIP 4, DM 4. The earlier figure of 202 across 122 clusters included 31 memberships
formed on a degenerate alignment (`T25`, closed 2026-08-18) — the correction fell hardest on
Myb (16 to 4) and AP2 (9 to 1). `T14` is the only open lead that would deepen it.

That is what makes the project's central question askable at all — **does sensitivity to
single-residue change transfer across folds?** With variants in homeodomain only it was
unanswerable; five non-homeodomain families now carry them. Still thin, and it resurfaces on
its own when modelling starts.

**The label noise floor is measured and it is not small** (`reports/overlap.md`, 2026-08-17).
47 domains are stored by two sources; the median pair agrees on **45.8%** of the 8-mers either
called positive, against 70-72% for replicates within one source. Any model that reproduces
held-out positives much past that is reproducing a laboratory.

**Only one of each duplicate reaches the merged table** (`T15`/`D4` closed 2026-08-18): the
record with more positives, unless a variant series is involved, in which case the series
wins. 48 of 95 records drop at merge, 3.5% of rows — `snp2prot.merge`, applied at merge, never
at parse.

**The cutoff was swept and kept** (`reports/threshold_review.md`, `T20` closed 2026-08-18).
`E >= 0.45` is the best absolute cutoff there is, and the rank-matched alternative was
rejected: it would turn negatives into a ranking artifact and invent ~100 positives for each
`BAR15A` variant that lost binding. The sensitivity asymmetry it would have hidden — two labs
differing 1.6x in positive count on one protein — is now recorded in `docs/METHODS.md` §5.1
instead.

**The merged table exists** (`reports/merge.md`, `reports/RESULTS.md`): 44,014,848 rows, one
record per domain, 48 duplicate records dropped. `RESULTS.md` is generated from the tables and
is the thing to read first.

**Records with no positive evidence are flagged, not re-binarized** (`reports/label_health.md`,
`T21` closed 2026-08-18). 54 of 1,386 records have no positive 8-mer: 20 are variants that
measurably lost binding, with a control in their own series, and 34 have nothing to compare
against — 2.5% of the corpus, dropped at training time with `label_health.usable(df)`, not at
parse time.

**Timings are measured, from the full rebuild of 2026-08-18.** `build_dataset.py --all` is
about **7.5 minutes** and the whole pipeline — parse, cluster, protein table, label health,
reports, overlap, merge, results, audit — is **~15 minutes**.
Per-step timings and the order to run them in are `docs/METHODS.md` §10.2. **The library must
be pressed once per machine — `python scripts/press_pfam.py`** — or every source pays 19.5 s to
re-read 2.2 GB of text.

**The plan, in order:**

| # | step | state |
|---|---|---|
| 1 | **Extend PBM coverage** beyond UniPROBE | **closed 2026-08-18** — Kock 2024 screened and excluded ([`reports/kock2024_excluded.md`](reports/kock2024_excluded.md)); `T2`/`T2b` dropped 2026-08-17. `T14` is the one small lead left open |
| 2 | **Deepen the protein axis in what we already hold** | **done 2026-08-17** — clustering by sequence distance, and the cluster inventory (`T3`) with it |
| 3 | **Merge**: single table, splits, NN baseline | **merge done 2026-08-18** — `data/processed/training.parquet`, 44,014,848 rows = 1,338 x 32,896. Splits and the NN baseline are still docstring stubs |
| 4 | **Modelling** | after step 3 |

A research sweep on 2026-08-14 surveyed the literature for PBM sources with designed protein
variation. Its one large find, Kock et al. 2024, was acquired, screened and **excluded on
2026-08-18** — the deposit publishes no E-scores and its own statistic has no transferable
scale ([`reports/kock2024_excluded.md`](reports/kock2024_excluded.md)). Two of its other leads
were already parsed and one fails condition 2 outright. See
[`docs/DECISIONS.md`](docs/DECISIONS.md) §9.

---

## Open tasks

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

### T27 — Near-identical domains can sit in different clusters
Surfaced by the `T25` sweep, and **not caused by it** — this is the greedy algorithm's known
order dependence, now measured for the first time. A domain joins the *first* seed
within 5 edits, not the best, so two paralogues that are 1-5 edits apart can end up in separate
clusters if one of them matched an earlier seed first.

**13 such pairs exist** across the four families swept (Homeodomain, Myb, AP2, HLH), all at
full overlap: `BAR15A:VAX2`/`Cell08:Vax1` at 1 edit, `BAR15A:HOXC4`/`Cell08:Hoxb4` at 1-2,
`weirauch2014:six3`/`Optix` at 2-3, `Cell08:Hoxc10`/`Hoxa10` at 4, `BAR15A:VSX1`/`Cell08:Vsx1`
and `Cell08:Irx4`/`Irx6` and `BAR15A:PITX2`/`Cell08:Pitx3` at 5.

**This is leakage in leave-one-cluster-out**: hold out one cluster and a near-identical
sequence remains in training. Options are a post-pass that merges clusters whose
seeds are within the threshold (reintroducing some chaining, which single-linkage
was rejected for), a stricter split regime that groups clusters by connected component at
evaluation time only, or accepting and reporting it. The last is cheapest and honest, and the
split code is not written yet — decide it there rather than in the clustering.

### T26 — `.gitignore` silently drops the docs it promises to keep
`data/raw/*` excludes the *directory*, so git never descends into it and the
`!data/**/README.md` / `!data/**/HOWTO.md` negations below can never fire — a directory
excluded at one level cannot have its contents re-included. Two files are affected today and
neither is in the repo: `data/raw/weirauch2014/README.md` has been missing all along, and
any future `data/**/README.md` or `HOWTO.md` would go the same way.

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

_None open._

---

## Notes

### N8 — Kock et al. 2024 is excluded, but the work is not lost
Screened, rescored and **excluded 2026-08-18** (`docs/DECISIONS.md` §1). Everything is in
[`reports/kock2024_excluded.md`](reports/kock2024_excluded.md): verified identifiers (Nat
Commun 15:3110, GEO `GSE233827`, Dataverse `10.7910/DVN/FDQHCF`), the yield through our own
admission policy (**67 new domains**, 13 singleton clusters would become variant-bearing), and
the measurements that decided it — no E-scores in the deposit, and recomputed ones reproduce
the ranking but not the scale.

**It can be revived without redoing the investigation.** The report's §3.2 records the three
things that made rescoring work (keep saturated spots; Cy3 as a sequence-model residual, not a
ratio; join by `(Column, Row)`, never by probe ID) and §5 the conditions that would make it
worth it — chiefly `T20` settling on a rank-matched rule, which makes a recomputed source's
compressed scale stop mattering. The code was deleted with the decision; it is ~200 lines and
the report specifies it.

### N7 — Cluster size is a one-file lookup now
`data/interim/clusters/clusters.parquet`, one row per cluster, written by `build_clusters.py`.
Read it with `snp2prot.clusters.load` / `ids_with_at_least(n)` / `select(df, n)` rather than
grouping the row tables — and note that **size means distinct canonical domains**, not rows and
not constructs, so a domain assayed by two labs counts once. Today: 108 clusters hold more than
one domain, 28 hold three or more, 11 hold five or more.

Two of its columns name domains and they are not the same one: `seed` is what membership was
decided against, `reference` is the medoid and the frame `mut_positions` uses (`D3`).

### N1 — The row count is an exact invariant
`45,593,856 = 1,386 x 32,896`, where 1,386 is 1,338 distinct domains plus 48 measured by more
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
**79-84 s**, and reports the domain and cluster counts so a mismatch is visible: 1,338 / 1,165
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
were rebuilt on 2026-08-17 and `MAR17A` on 2026-08-18, the other 13 sources were not — and
`MAR17A` showed what that costs: re-running it dropped `Esrrb` as intended **and added
`MAR17A:Mlx`**, a domain the current parser admits and the committed table never had. **`build_dataset.py --all` followed
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
