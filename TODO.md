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

**The merged table exists** (`reports/merge.md`, [`docs/RESULTS.md`](docs/RESULTS.md)):
44,014,848 rows, one record per domain, 48 duplicate records dropped. `RESULTS.md` is
generated from the tables and is the thing to read first.

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

**The nearest-neighbour baseline is built and measured** ([`reports/nn_baseline.md`](reports/nn_baseline.md),
2026-08-19) — step 0 of the modelling plan, and it did both jobs §8.1 gave it.

**The bar**, mean per-protein AUPR, `k = 1`:

| regime | AUPR | median NN identity |
|---|---:|---:|
| `S1` random domains | **0.769** | 0.742 |
| `S2` components at >= 0.5 identity | **0.296** | 0.386 |
| `P1` homeodomain holdout | **0.024** | 0.268 |
| `P2` non-homeodomain variant clusters | **0.882** | 0.824 |
| `P3/all` every variant, wild types kept | **0.928** | 0.986 |

The pure `k = 1` lookup beats the identity-weighted top-5 average almost everywhere, so the
primary form stays the primary form.

**And it checked the splits, which is where it earned its place.** Both findings went to the
owner and both are now settled (`docs/DECISIONS.md` §2): **`S2`'s edge threshold moved from 5
edits to 0.5 identity** (`D5`), which took it from 0.751 — barely below `S1` — to 0.296; and
**C1 is evaluated on the 29 variants where the wild-type copy fails** (`D6`), because `P3/all`'s
median is 1.000 and the claim is invisible in the mean.

The whole lookup curve is one variable, which is the most useful thing the report says: AUPR by
nearest-neighbour identity runs 0.032 / 0.365 / 0.771 / 0.946 / 0.910 across the bands below 0.3,
0.3-0.5, 0.5-0.7, 0.7-0.9 and above 0.9. Any model's number has to be read against its own band.

**The sequence arms are built and the pooling pre-flight passes**
([`reports/pooling_check.md`](reports/pooling_check.md), 2026-08-19). `A1` (ESM-2 650M) and `A4`
(ESM-DBP) are cached as pooled 1280-d vectors per domain, 20 s each on the GPU.

`ML_PLAN.md` §3.1 required this measurement before any training run, because mean pooling could
in principle dilute a single-residue change to nothing:

| arm | variant vs family scale | own reference nearest | displacement ~ edits | C1 signal (AUC) |
|---|---:|---:|---:|---:|
| `A1` ESM-2 650M | 0.028 | 73% | 0.408 | 0.622 (1.3σ) |
| `A4` ESM-DBP | 0.026 | 79% | 0.410 | **0.709 (2.4σ)** |

**Pooling stays** — variants are separated, the separation scales with the number of edits, and
96% of variants sit within the nearest 5% of their family. Attention pooling (`T32`) is not to be
reached for, per §3.1's own instruction.

**And the `A1` -> `A4` comparison is already visible before any training.** The domain-adapted
model is better on every column, and on the one that matters most for C1 — does the vector move
*further* for the mutations that actually changed binding — it is 2.4 standard errors from chance
against A1's 1.3. Weak evidence at n = 28 vs 62, and a property of the embeddings rather than a
result about binding, but it is the first evidence in the project bearing on §4.2's question.

**Modelling has a plan of record now: [`docs/ML_PLAN.md`](docs/ML_PLAN.md)** (2026-08-19).
It adds a second dataset — Codebook SELEX, Jolma et al. 2026 *Nature*
`10.1038/s41586-026-10798-9` — and a contrastive two-tower model over 1 DNA and 3 protein
embeddings, for a talk. §1-§9 is the owner's plan; §10 is a review of it, unacted on.
`T28`-`T31` below are the items that review raised.

**The plan, in order:**

| # | step | state |
|---|---|---|
| 1 | **Extend PBM coverage** beyond UniPROBE | **closed 2026-08-18** — Kock 2024 screened and excluded ([`reports/kock2024_excluded.md`](reports/kock2024_excluded.md)); `T2`/`T2b` dropped 2026-08-17. `T14` is the one small lead left open |
| 2 | **Deepen the protein axis in what we already hold** | **done 2026-08-17** — clustering by sequence distance, and the cluster inventory (`T3`) with it |
| 3 | **Merge**: single table, splits, NN baseline | **done 2026-08-19** — `data/processed/training.parquet`, 44,014,848 rows = 1,338 x 32,896; five split regimes in `snp2prot.splits`; the baseline in [`reports/nn_baseline.md`](reports/nn_baseline.md) |
| 4 | **Modelling** | **step 0 and the sequence embeddings done 2026-08-19** — [`docs/ML_PLAN.md`](docs/ML_PLAN.md). `A1` and `A4` are cached and the §3.1 pre-flight passes ([`reports/pooling_check.md`](reports/pooling_check.md)); next is the two-tower harness itself |

A research sweep on 2026-08-14 surveyed the literature for PBM sources with designed protein
variation. Its one large find, Kock et al. 2024, was acquired, screened and **excluded on
2026-08-18** — the deposit publishes no E-scores and its own statistic has no transferable
scale ([`reports/kock2024_excluded.md`](reports/kock2024_excluded.md)). Two of its other leads
were already parsed and one fails condition 2 outright. See
[`docs/DECISIONS.md`](docs/DECISIONS.md) §9.

---

## Open tasks

### T33 — Codebook publishes PBM data too; check whether we already hold it
Unhurried, and **not a proposal to reopen acquisition** — that closed 2026-08-18
(`reports/kock2024_excluded.md`). But <https://codebook.ccbr.utoronto.ca/v2/index_v2.php> lists a
**PBM** panel alongside the SELEX data, with raw scanner data, Z-scores, **E-scores and 8-mer
statistics** — the same assay, the same score type and, if it is a universal PBM, the same 32,896
8-mer space our 19 sources share.

So the usual cost of a new source does not apply: no new binarization rule, no new negative
semantics, no `dna_len` question. The only real questions are whether it overlaps what we already
parsed and whether it clears the three-condition admission policy.

An hour of checking, whenever. Run `scripts/audit_sources.py` and the `N3` checklist if it turns
out to be new.

### T32 — Encoder variants, to test once a first result exists
The encoder is fixed for now at **plain one-hot with reverse-complement mean pooling**
([`ML_PLAN.md`](docs/ML_PLAN.md) §4.1, decided 2026-08-19). These four are the deferred
alternatives. All are **one-encoder swaps against a fixed pipeline**, so once the training
harness exists each is cheap — the whole label matrix fits in GPU memory (`ML_PLAN.md` §10.10).

**Protein side — attention pooling instead of mean pooling.** Protein embeddings are pooled to one
1280-d vector per domain ([`ML_PLAN.md`](docs/ML_PLAN.md) §3.1), which is what makes the four arms
comparable and matches TransBind. The risk it carries is on claim C1: pooling dilutes a
single-residue change, and §3.1 specifies a pre-flight check that measures how much. **If that
check comes back bad, attention pooling is the fallback** — still one fixed-width vector per
domain, so every reason for pooling survives, but learned rather than uniform weighting over
residues. Do not reach for it before the measurement says it is needed.

**DNA side — the four below.**

**a. Other symmetry pooling.** Mean is the chosen default; **max** is the DeepBind/Basset
convention and is not obviously worse — mean averages the two strands' evidence, max takes the
better-matching strand, which is arguably closer to what a TF does. `logsumexp` sits between
them. Cheapest of the four.

**b. Chemical feature channels.** Three binary channels alongside the four one-hot ones:
**S/W** (`GC`/`AT`, 3 vs 2 hydrogen bonds), **R/Y** (`AG`/`CT`, purine vs pyrimidine ring),
**M/K** (`AC`/`GT`, amino vs keto). These are the groupings a TF actually reads in the major and
minor groove. Note they **uniquely determine the base** — A=WRM, C=SYM, G=SRK, T=WYK — so they
are a complete encoding on their own, not redundancy; the point is to offer the network a second,
chemically grouped view. Seven channels, ~5 lines, no new dependency.

**c. DNAshape** (Rohs lab: minor groove width, propeller twist, roll, helical twist, from
pentamer lookups). Genuinely **the established biophysical DNA encoding in the TF-binding
field** — a far better answer than a genomic LM if anyone asks what "the accepted method" is.
Costs a real build: a pentamer table, and 8-mers have edge effects because a pentamer window only
covers 4 of the 8 positions.

**d. Brute-force symmetry by augmentation.** Instead of pooling architecturally, present **both
reverse complements during training** as separate examples with the same label, and drop the
pooling. Tests whether the inductive bias is needed or whether data augmentation buys the same
thing — the classic equivariance-vs-augmentation trade, and worth knowing for this task.

Note the comparison in **d** has a clean quantitative readout beyond AUPR: augmentation makes the
embedding only *approximately* symmetric, so measure `||f(s) - f(revcomp(s))||` over held-out
8-mers. Architectural pooling drives it to exactly zero by construction; augmentation does not,
and how close it gets is the answer.

### T31 — Two things to agree with the structure colleague
Everything else about the ensembles is hers to decide or ours to discuss — she is in the lab and
knows the project. Only these two need flagging, because both are cheaper to settle before the run
than after.

**a. Fix the ensemble size.** Arm A3 represents an ensemble by its **spread**
([`ML_PLAN.md`](docs/ML_PLAN.md) §4.2), and a standard deviation estimated from 10 structures is
both noisier and differently biased than one from 50. If `N` varies per domain, A3 silently
encodes *how many structures that domain got* — which correlates with nothing biological. A fixed
`N` for every domain is the clean fix; variable `N` subsampled down to the smallest before
computing spread also works. Mixed sizes used as-is is the one option that is wrong.

**b. We may only need the variant-bearing clusters.** The full ask is 1,338 domains; the
**108 variant-bearing clusters hold 281 domains**, so scoping to those cuts the run to ~21% of the
structures. That is where the ensemble earns its cost — the singletons answer no question that
needs a spread.

**But it has a consequence worth deciding with her, not around her.** Structures for only 281
domains means arms **A2 and A3 cannot run the full-axis `S1`/`S2` regimes** (§5), which use all
1,338. The comparison would then have to be either (i) restricted so *every* arm — A1 included —
sees only the 281-domain subset, keeping it apples-to-apples on a smaller axis, or (ii) full
structures after all. Option (i) is cheap and still valid; it just makes the headline
sequence-vs-structure number a subset result, and that has to be said on the slide.

### T30 — What to do with the 54 zero-positive records
`label_health.usable()` drops all 54. They are not one group, and the choice differs per group
([`ML_PLAN.md`](docs/ML_PLAN.md) §3.1):

| verdict | n | what it is |
|---|---:|---|
| `dead_variant` | 20 | a variant that measurably **lost** binding, with a control in its own series |
| `no_evidence` | 34 | no positives **and** no control to read that against |

**The 20 are the sharpest evidence in the corpus for claim C1.** A model that returns the
wild-type profile for a dead variant has failed in exactly the way the project exists to detect —
and under §6.1's `P3` regime the NN baseline predicts precisely that. Proposal: hold them out as a
dedicated C1 evaluation set, never trained on, reported separately.

**The 34 are a different question, raised 2026-08-19 by the InfoNCE decision.** An all-negative row
is genuine training signal — it pushes every 8-mer away from that protein — and §3.1 specifies a
`null` anchor as the mechanism that lets InfoNCE express it. The catch is that **without a control
you cannot tell a true non-binder from a failed assay**, so training on the 34 risks teaching the
model that a perfectly good TF binds nothing. Options: exclude (status quo), include via the null
anchor, or include as an ablation and report the difference.

Both decisions are the owner's, and both want settling before the split code is written. The
mechanism is worth building either way — it will matter more on SELEX, where all-negative proteins
may be commoner.

### T29 — Measure how reliable the SELEX negatives are, once the data lands
**Mostly resolved before it started.** The collaborators can score **every** 8-mer, positives and
negatives, so a SELEX negative is *measured and scored below threshold* — `assayed_unbound`, not
`selection_absent` — and rule 4 is satisfied ([`ML_PLAN.md`](docs/ML_PLAN.md) §2.3).

What is left is a measurement, not a decision. A SELEX negative still comes from **depletion in a
selection** rather than a direct spot reading, and enrichment is noisy for 8-mers poorly covered
in the input library, so SELEX may be **less reliable on the negative side** than PBM. Do not
assume it either way — check:

- the enrichment-score distribution and per-8-mer read coverage, when the table arrives;
- **cross-assay agreement on the shared proteins**, which the §7 transfer experiment's *easy* set
  gives for free. Compare agreement on positives against agreement on negatives; if the negative
  side is much worse, that is a property of the assay pair and gets reported as one rather than
  absorbed into a model score.

Related and already folded into `ML_PLAN.md` §7: the measured **45.8%** median cross-lab
agreement (`reports/overlap.md`) goes on the transfer figure as an explicit ceiling. Two labs
running the *same* assay agree that much; two different assays will not do better.

### T28 — Send the SELEX collaborators the k-mer spec
**The decision is made — 8-mers only** ([`ML_PLAN.md`](docs/ML_PLAN.md) §2.3, 2026-08-19). What
is left is the task of sending it, and it wants sending **before** they start work, not after.

The spec, in four lines:

- **k = 8**, nothing else;
- over the **same 32,896** non-redundant 8-mers the PBM tables already use;
- canonicalised as the **lexicographically smaller** of `(s, revcomp(s))`;
- with a **continuous score** per 8-mer, so we set our own cutoff rather than inherit theirs.

**And decide the provenance treatment when they arrive.** The k-mers are neither a download nor
our own parse, so rule 3 has no obvious slot for them. Recommended: **treat them as raw for us** —
`data/raw/codebook_selex/`, append-only, with a `PROVENANCE.md` row naming the collaborator, the
date, a one-line method description and the file's sha256, plus their code archived alongside. The
paper's own accessions (`ML_PLAN.md` §2.1, read from its data-availability statement, never
guessed) go in the same row so the chain from publication to our table is unbroken. Rule 3 makes
this publication evidence, not a convenience log.

The third line is the one that fails silently if it is missed. Verified against the tables
2026-08-19: all 19 sources carry one identical 8-mer set, and every entry is the lexicographically
smaller member of its reverse-complement pair (the 256 palindromes trivially). A table
canonicalised the other way validates, looks the right size, and joins on almost nothing.

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

_None open._ `D5` (the `S2` identity floor) and `D6` (the C1 evaluation set) were settled by the
owner on 2026-08-19 and are recorded in [`docs/DECISIONS.md`](docs/DECISIONS.md) §2.

---

## Notes

### N9 — `gene` and `wt_id` disagree for ~100 domains, `weirauch2014` mostly
Noticed 2026-08-19 while building `snp2prot.corpus`. Of the 1,057 singleton clusters, **129 have
a `gene` that is not the gene their `wt_id` is named after** — 114 of them in `weirauch2014`,
and about 15 of those are only `(unnamed)`. A singleton cluster is one domain, so the two names
describe the same record and one of them is wrong: `C:weirauch2014:hoxd13a` is a zebrafish
homeodomain and the protein table calls it `TCP20`, which is a different family entirely.

**Nothing downstream depends on it.** `gene` is a display column; labels, sequences, clusters and
every split are keyed on `dbd_seq`. It matters when a person reads a report and when a variant
has to be checked against the paper that describes it, which is `T23`'s job.

The likely site is the `weirauch2014` join — GEO SOFT and Table S6 are joined on plasmid ID
(`docs/METHODS.md`), and `wt_id` and `gene` come from different sides of it. Worth an hour before
any figure carries gene names.

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
