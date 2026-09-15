# SNP2PROT — working notes for Claude

## What this project is

A **training set** for a model that predicts whether a transcription factor binds a DNA
site, built so that it varies on **both axes at once**:

- **protein axis** — DNA-binding domain (DBD) sequence, with many single-residue variants of
  the same WT domain forming dense clusters;
- **DNA axis** — binding sites ≤20 bp, with many single-base variants.

The label is **binary**. Negatives are wanted in bulk (50–100:1 neg:pos is expected, not a
class-imbalance bug to fix). Full-length protein is not needed; DBD only.

The full specification is [docs/TFDNA_MERGE_BRIEF.md](docs/TFDNA_MERGE_BRIEF.md) — the
owner's brief, stored verbatim. **Read it before touching anything under `src/snp2prot/`.**
Where it and this file disagree about *paths*, this file wins (see the mapping below); where
they disagree about *intent*, the brief wins.

**Current state: 19 PBM sources parsed and validator-clean — 45,593,856 rows, 1,338 domains
in 1,165 clusters, 56 Pfam families, 131 organisms, every protein scored against the same
32,896 8-mers. `dbd_seq` is the canonical domain (`snp2prot.canonical`) and `wt_id` comes from
CD-HIT clustering (`scripts/build_clusters.py`), not from construct lineage. Two source
families and no more: 18 UniPROBE accessions plus CIS-BP / Weirauch 2014 (`weirauch2014`, GEO
GSE53348). Acquisition closed on 2026-08-18 when Kock et al. 2024 was screened and excluded
([reports/kock2024_excluded.md](reports/kock2024_excluded.md)); the next step is phase 5,
merge. See [TODO.md](TODO.md).**

**There is now a second arm.** [docs/GHT_PLAN.md](docs/GHT_PLAN.md) is the plan of record for
**human Codebook GHT-SELEX** on the MEX-ArChIPelago benchmark — 33 TFs, 9.2 M 301 bp genomic
windows, the benchmark's own positives, negatives and chromosome splits. It is a separate corpus
with a separate loss and **nothing in it merges into `data/processed/training.parquet`**; the two
arms share the *protein tower*, which is the scientific point. [docs/GHT_RESULTS.md](docs/GHT_RESULTS.md)
is what it measured. Code lives under `src/snp2prot/ght/`.

**Read [docs/METHODS.md](docs/METHODS.md) first** — it records how the dataset was built and
why, in enough detail to reimplement. [docs/DOMAIN_POLICY.md](docs/DOMAIN_POLICY.md) states
the admission rules; [TODO.md](TODO.md) holds open tasks and decisions, [docs/DECISIONS.md](docs/DECISIONS.md) the resolved ones.

> ### `dbd_seq` is the PADDED Pfam domain
> Not the bare Pfam envelope, and not the sequence that was on the array. It is the Pfam
> envelope plus 10 residues each side. Three BAR15A variants mutate residues 2-9 residues
> N-terminal of the Pfam start, and without the padding each collapses onto its own wild type
> with a different label. **Read [docs/DOMAIN_POLICY.md](docs/DOMAIN_POLICY.md) before using
> the column.** Settings live in the `domain:` block of `configs/thresholds.yaml`.

## The admission policy — three conditions

A construct enters the dataset only if the stored `dbd_seq` satisfies all three. They follow
from what the data is for: predict a structure, find the centre of the DNA-contacting
residues, embed the residues around it.

1. **Sole responsibility** — that subunit alone produced the measured interaction.
2. **One continuous region** — one contiguous stretch, not fragments scattered through the
   protein.
3. **The variation lies inside it** — every change that alters binding falls within the
   stored region, or a sequence model sees identical inputs with different labels.

Enforced in `snp2prot.domains`, which rejects mixed-family constructs (PAX+homeodomain,
POU+homeodomain), C2H2 zinc-finger arrays, constructs with no Pfam hit, and any variant whose
mutation escapes the padded window. Rejections are reported, never silently repaired.

## The failure this project exists to avoid

If the DBDs in the set are all distant from each other, a model can memorize per-protein
profiles and score well without learning anything — a lookup table with extra steps. Two
concrete consequences that shape the code:

1. **Never pool sources without checking `dna_len`.** PBM gives 8 bp, B1H 9 bp, SNP-SELEX
   19 bp, and each has a different positive rate, so length alone leaks assay identity and
   with it the label prior. Mitigated at *featurization* time, never in the stored table.
2. **The nearest-neighbour baseline is the bar, and it copies the neighbour's BINARY calls.**
   `snp2prot.baselines.nn_lookup` takes the most identical training DBD under the overlap guard
   and copies what it *binds*, not its E-scores. That is the whole point: the model trains on
   thresholded labels, so a baseline copying the continuous profile would be scored on
   information the model never sees — measured, that was worth 29-76% of its AUPR (`D7`,
   2026-08-27). A PBM E-score is a rank-enrichment statistic read at a cutoff, not a graded
   affinity; **the whole PBM corpus is binary, baseline and evaluation included.**
   Under `P3/all`, copying the wild type *is* the hypothesis that a mutation does nothing.
   [reports/nn_baseline.md](reports/nn_baseline.md) bands the folds by nearest-neighbour
   identity, which is what drives every one of them.

   **Compare a model against `k = 5`, not `k = 1`.** `k = 1` copies one neighbour's calls — ~50
   tied 1s over ~32,000 tied 0s — so its AUPR measures set overlap, while a model emitting 32,896
   distinct scores is judged on a ranking; the gap is partly output format. `k = 5` averages five
   neighbours weighted by identity, giving a ranking from the same binary labels. Measured
   2026-08-28, `k=1` / `k=5`: `S1` 0.475 / 0.680, `S2` 0.134 / 0.256, `P1` 0.006 / 0.007,
   `P2` 0.574 / 0.752, `P3/all` 0.660 / 0.840. The model beats `k=1` by +0.16 to +0.25 and `k=5`
   by **+0.005 to +0.052** — the second is the real margin. **Claim C1 is directionally
   supported but an order of magnitude short**, and **`TODO.md` `T38` questions whether the
   approach is usable at all**: the model emits a ranking, the task needs a call, and no
   defensible threshold exists.
   [docs/ML_RESULTS.md](docs/ML_RESULTS.md) is the reading of all of it.

4. **The ranking is not the deliverable, and as of 2026-08-28 the loss says so.** The InfoNCE
   objective is a per-row softmax and is therefore *exactly invariant to adding a constant to a
   protein's row* — it never says where a row sits, only how it is ordered, which is why the null
   anchor never became a threshold and why no post-hoc calibration could have. `L = L_infonce +
   λ · L_bce` adds the term that is not shift-invariant, and `λ = 0` is the earlier objective
   exactly. (It does not reproduce earlier *numbers*: GPU training is not reproducible across
   processes — measured at ~0.003 AUPR on `S1/fold-2` — so the sweep carries its own `λ = 0`
   control.)

   **Measured 2026-08-30 over 312 runs, and the result splits in two.** The decision rule works:
   `λ = 30` costs nothing in ranking, the per-protein call count now tracks the truth, and
   Spearman(predicted count, true count) is 0.48-0.52 wherever the model generalises. **The
   capability it was built for is absent**: dead-variant detection is 0.43-0.52 against the
   baseline's **0.492** — chance, and the literal null hypothesis, since under `P3/all` the
   baseline copies the variant's own wild type. The model returns the wild type's answer for a
   single-residue variant (Spearman with the wild type's true count **+0.77 to +0.91**, rising
   with λ), so **this was never a threshold problem — it is a representation problem**, and it
   belongs to `T36`. Scored as decisions rather than rankings the model also loses to `k = 1` by
   0.006 call F1 while beating `k = 5` by 0.035 AUPR: the margin is a ranking margin.
   `docs/ML_RESULTS.md` §9 is the reading, `docs/DECISIONS.md` §13, `docs/TRAINING.md` §2.5,
   [reports/calibration.md](reports/calibration.md).

3. **Read a per-protein AUPR against its chance level, never raw.** It is anchored to that
   protein's own positive rate, which ranges 0.0015-0.0034 across the folds — so 0.05 on `P1` and
   0.05 on `S2` are not the same statement. Every fold reports `chance_aupr` and a sampled
   random-ranking 95th percentile, and a result at or below it is flagged **at chance**. One row
   of the first grid was.

## Layout

```
TODO.md                     open tasks, open decisions, notes  <- THE WORKING FILE
docs/DECISIONS.md           decisions taken, defects fixed, sources excluded
docs/TFDNA_MERGE_BRIEF.md   the owner's spec, verbatim — source of record
docs/REFERENCES.md          dataset -> paper -> Crossref-verified DOI, plus reading order
docs/UNIPROBE_ACCESSIONS.md all 36 UniPROBE accessions, citations, family survey
docs/DOMAIN_POLICY.md       what dbd_seq is, the padding, and what gets excluded  <- READ THIS
docs/METHODS.md             publication-quality account of how the dataset was built
docs/ML_PLAN.md             the modelling plan for the talk; §10 is a review, not the plan
docs/ML_RESULTS.md          what phase 7 measured and what it rules out  <- READ WITH ML_PLAN
docs/GHT_PLAN.md            the GHT-SELEX arm: plan of record  <- READ BEFORE src/snp2prot/ght/
docs/GHT_RESULTS.md         what the GHT arm measured, and what it does not show
docs/TRAINING.md            how the two-tower model is trained: loss, null anchor, batch shape,
                            and section 2.5 the calibration term that gives it a decision rule
                            <- READ BEFORE TOUCHING src/snp2prot/models/ OR training.py
docs/papers/                paper PDFs (git-ignored); README.md there is the manifest
docs/papers_inbox/          the owner drops papers here; Claude identifies and files them
PROVENANCE.md               one row per raw file: URL, accession, timestamp, size, sha256
configs/thresholds.yaml     the ONLY config file; all binarization cutoffs live here
reports/                    generated markdown, committed; the phase-boundary deliverables
docs/RESULTS.md             what the finished dataset contains, generated from the tables

src/snp2prot/
  schema.py       unified 22-column row schema + validate(); the gate every parser passes
  thresholds.py   reads configs/thresholds.yaml
  experiment.py   reads configs/experiment.yaml -- the MODELLING config, deliberately separate
  config.py       every path in the project; nothing builds a path by hand
  domains.py      Pfam/HMMER annotation + the three-condition admission policy
  canonical.py    the project-internal canonical domain sequence  <- READ BEFORE dbd_seq
  references.py   identifier -> reference protein, cached; only used when a construct is short
  clusters.py     CD-HIT greedy incremental clustering + the cluster inventory (size lookup)
  label_health.py which records carry positive evidence; the training-time filter (T21)
  corpus.py       one row per domain of the merged table: the view splits and baselines read
  distances.py    all-vs-all domain identity, cached; feeds S2, the NN baseline and §5.2
  embeddings.py   pooled protein-LM vectors, one per domain, per arm (A1, A4)
  models/         the two towers, the null anchor and the InfoNCE loss  <- docs/TRAINING.md
  training.py     the training loop and the per-fold run driver
  tracking.py     MLflow behind a wrapper; a run cannot exist without its split digests
  merge.py        one record per domain: which of two sources' measurements survives (T15/D4)
  proteins.py     the protein-side companion table (bare/padded domain, construct, full-length)
  parsers/        one module per source, each exposing parse() -> pd.DataFrame
    _pbm.py       assay-level machinery shared by ALL universal-PBM sources
    _uniprobe.py  UniPROBE file layout on top of _pbm (formats differ per accession)
    weirauch2014.py  CIS-BP / Weirauch 2014: GEO SOFT + Table S6, joined on plasmid ID
  metadata/       CIS-BP / Pfam+HMMER / UniProt lookups shared by all parsers
  reports.py      binarization summary, cluster inventory, overlap report
  splits.py       the five regimes: S1, S2 (components at >= 0.5 identity), P1, P2, P3
  baselines/      nn_lookup — the number a model must beat
  evaluation/     metrics.py — per-protein AUPR / precision@k / R@P, macro-averaged, each with
                  its chance level and a sampled random-ranking band; `suppression` for the dead
                  variants, where every ranking metric is undefined
                  calibration.py — the NON-rank metrics (T38): ECE, the called set's F1/Jaccard,
                  the parameter-free `expected_count_rule`, and `interaction_power` /
                  `detection_auroc` — has this protein lost binding altogether
                  ceilings.py — rank_ceiling only; the ridge and RBF probes were deleted (D10)
                  c1.py — the C1 evaluation set: variants the wild-type copy fails on (D6)
  data/matrix.py  the merged table as a dense domain x 8-mer array pair, for modelling
  ght/            THE SECOND ARM -- human Codebook GHT-SELEX on the MEX-ArChIPelago benchmark.
                  Separate corpus, separate loss, separate everything except the PROTEIN TOWER,
                  which is the point. Nothing here merges into data/processed/training.parquet.
                  <- docs/GHT_PLAN.md is the plan, docs/GHT_RESULTS.md the reading
    panel.py      139 benchmark TFs -> 33, through the SAME admission policy as every source
    windows.py    the benchmark's BEDs -> 301 bp base tokens cut out of hg38
    data.py       the window corpus as resident arrays, and the protein-side ABLATIONS
    model.py      GenomicTwoTower: cosine -> calibrate -> BCE, no null anchor (GHT_PLAN 7)
    splits.py     C1 (chromosomes), G1/G2 (TFs) -- two holdout axes, both in the digest
    training.py   the loop, the per-fold driver, and the PBM protein-tower transfer
    baselines.py  PWM best-hit / sum-occupancy, and motif transfer by DBD identity

data/raw/<source>/       append-only, never edited          (git-ignored)
data/interim/<source>/   per-source parsed Parquet          (git-ignored)
data/interim/proteins/   protein table: domain at 4 levels + UniProt full-length
data/interim/clusters/   one row per cluster: size, family, variants  <- read via snp2prot.clusters
data/interim/label_health/ one row per (domain, source): verdict, label counts, best E-score
data/external/pfam/      Pfam HMMs for boundary annotation
data/external/uniprot/   cached canonical sequences
data/interim/ght/         the GHT arm: panel.parquet (33 TFs) and windows/<tf>.npz (9.2M
                         301 bp windows, ~3.3 GB)                       (git-ignored)
data/processed/ght/      the GHT arm's embeddings, checkpoints, fold/TF tables and baselines
data/external/genomes/   hg38.fa.gz, verified against UCSC's own md5   (git-ignored)
data/processed/checkpoints/ trained weights, one file per (arm, fold, seed), both the
                         validation-selected model and the model at the step budget
data/processed/          merged training table: 1,338 domains x 32,896 8-mers  (git-ignored)
                         plus lambda_sweep_folds.parquet / lambda_sweep_domains.parquet (the
                         T38 sweep, resumable -- deleting them restarts it from zero)
                         plus the cached modelling artifacts: kmer_matrix.npz (the same
                         table as arrays), distances.npz (all-vs-all domain identity),
                         nn_baseline_domains.parquet (per-domain baseline scores),
                         c1_variants.parquet (the C1 evaluation set) and
                         embeddings/<arm>.npz (pooled protein-LM vectors) and
                         training_folds.parquet / training_domains.parquet (grid results)
mlruns/                  local MLflow backend: mlflow.db + artifacts  (git-ignored)
data/external/models/    language-model checkpoints; README.md there says how they arrive
data/testsets/           held-out sets, kept physically apart; empty  (git-ignored)
```

### Path mapping vs. the brief

The brief was written against a generic layout. These are the same things under this repo's
names — **do not create the left-hand paths**, that would fork the structure in two.

| brief says | this repo uses |
|---|---|
| `config/thresholds.yaml` | `configs/thresholds.yaml` |
| `src/parsers/<source>.py` | `src/snp2prot/parsers/<source>.py` |
| `src/splits.py` | `src/snp2prot/splits.py` |
| `src/baselines/nn_lookup.py` | `src/snp2prot/baselines/nn_lookup.py` |
| Tier 4 "separate directory" | `data/testsets/` (Tier 4 dropped; the path stands) |

## Rules that are not negotiable

1. **Never invent an accession, URL, or file name — or a DOI.** Several in the brief are
   unverified, and one plausible-looking DOI for Persikov 2015 is a different paper. If
   a download 404s or a file's layout differs from the description, **stop and report exactly
   what was found** — do not substitute a "similar" dataset, and do not guess a URL pattern.
2. **Raw stays raw.** `data/raw/` is append-only. Parsers read from it and write to
   `data/interim/`. Nothing edits a raw file, ever.
3. **Every raw file gets a `PROVENANCE.md` row** — URL, accession, UTC timestamp, bytes,
   sha256, one-line description. Use `scripts/record_provenance.py`. This is publication
   evidence, not a convenience log.
4. **`neg_provenance` is mandatory on every `label == 0` row.** A PBM low-E-score negative is
   real evidence of non-binding; a sequence merely absent from a B1H selection is not. If
   these are ever pooled indistinguishably, the dataset is dead. The validator enforces it.
5. **No cutoff literals in parser code.** Read `snp2prot.thresholds`; write the applied values
   into each row's `threshold_pos` / `threshold_neg` so a table is self-describing.
6. **Never pad `dna_seq` in the stored table.** Padding is a modeling decision. The validator
   rejects non-ACGT characters for this reason.
7. **One parser per source, no cross-source logic inside it.** Merging, dedup and
   reconciliation live in the merge step, not in a parser.
8. **Blocked ≠ stuck.** A source behind a licence click-through or MTA gets a
   `data/raw/<source>/HOWTO.md` with exact manual steps, a row in the PROVENANCE manual queue,
   and then you move to the next source. Do not burn effort brute-forcing a download.
9. **Stop at each phase boundary** and hand back the reports. Do not run ahead to modeling.
10. **Never run a command that takes more than a few minutes.** Measured 2026-08-17:
    `build_dataset.py --all` is 7 min 28 s and the whole pipeline is 13 min
    (`docs/METHODS.md` §10.2). Hand over the exact command
    instead, with how long it takes and what to check in the output. Those remain the owner's to run; single small
    sources and test subsets are fine to execute here. **A single-source rebuild must be
    followed by `scripts/build_clusters.py`** or that source silently leaves its clusters
    (`TODO.md` N6).
11. **Never `git commit` unsolicited.** The owner reviews work as an uncommitted diff;
    committing removes the review surface. Finish, run tests and lint, leave the tree dirty.

## Conventions

- **Model weights live in `data/external/models/` and carry a `PROVENANCE.md` row**, like the
  Pfam HMMs — they are the same kind of thing, a third-party artifact the build depends on and
  cannot regenerate. `scripts/record_provenance.py` accepts `data/external/` paths for this.
  `build_embeddings.py` points `torch.hub.set_dir` at that directory so nothing lands in a
  home-directory cache.
- **`configs/thresholds.yaml` holds only live blocks** — `pbm:`, `domain:`, `cluster:`.
  `b1h:`, `snp_selex:` and `tier4:` were removed on 2026-08-17 with the phases they configured.
- **There are exactly two config files, and the split between them is the point.**
  `configs/thresholds.yaml` configures the *dataset* (read via `snp2prot.thresholds`);
  `configs/experiment.yaml` configures *modelling* — split regimes, folds, seeds, baseline and
  metric settings — read via `snp2prot.experiment`. Added 2026-08-19 with the NN baseline, as
  [docs/ML_PLAN.md](docs/ML_PLAN.md) §9.2 specified. They are not merged because changing a
  threshold invalidates the dataset, every report and every provenance row, while changing a fold
  count or a seed does not. Two things with different blast radii do not belong in one file.
  Runs are to be tracked with MLflow on a local backend, and a run records **which domains were
  held out**, not just its hyperparameters — `snp2prot.splits.Fold.digest` is that record.
- **Paper PDFs go in `docs/papers/`**, git-ignored, named `<firstauthor><year>_<slug>.pdf`
  (`_supp`, `_supp-<what>`, `_fig<N>` for the rest). The manifest is `docs/papers/README.md`;
  `scripts/check_papers.py` reports what is missing, unlisted, or waiting in the inbox.
- **The owner adds papers by dropping them in `docs/papers_inbox/`** under whatever name they
  downloaded with. Filing them is Claude's job, and it is done by **reading each file's
  embedded metadata and first page — never by trusting the download name.** Elsevier `mmc`
  numbers have twice turned out to be main articles rather than supplements, and the two 2008
  Cell homeodomain papers have PIIs that run opposite to their DOI order. After filing, add a
  manifest row and leave the inbox empty.
- **Parquet, not CSV**, for anything over ~10⁵ rows. Checkpoint often.
- Parser output goes through `schema.coerce(df)`, then `schema.validate(df, source)`, and the
  report's `raise_if_failed()` before any write to `data/interim/`.
- `pair_id` and `dna_len` are **derived** — `coerce()` computes them. Never hand-set them.
- **Cluster distance is alignment-based, so variants may carry insertions and deletions.**
  `snp2prot.align.edit_profile` aligns a variant to its reference (BLOSUM62, affine gaps) and
  counts every edited residue. **Terminal gaps are free**: `dbd_seq` is the padded envelope
  *clipped* where the construct ends, so two versions of one domain differ at the termini for
  reasons that have nothing to do with the protein — in ROG18A every bare domain is 83-85 aa
  while the padded ones run 95-105 aa. Internal indels still count. Settings live in the
  `cluster:` block of `configs/thresholds.yaml`.
- `mut_positions` is 1-based **in the reference's frame**, not the variant's, so every member
  of a cluster shares one coordinate system and maps onto one predicted structure. A variant
  with a deletion can therefore name a position past its own length. One entry is emitted per
  edited residue, so `len(mut_positions) == n_mut_from_wt` always holds. Isoform-offset
  off-by-ones remain the most likely silent bug — every variant must land on the residue the
  paper claims.
- Type hints on public functions; `from __future__ import annotations` at module top.
- Line length 100; `ruff check` + `ruff format` clean before a commit.

## Commands

```bash
uv venv --python 3.11 .venv && uv pip install --python .venv -e ".[dev]"
.venv/bin/python -m pytest -q
.venv/bin/python scripts/press_pfam.py                     # ONCE per machine, before any build
.venv/bin/python scripts/build_dataset.py --all            # parse -> validate -> interim
.venv/bin/python scripts/make_reports.py --source BAR15A    # regenerate reports/
.venv/bin/python scripts/audit_sources.py                  # invariant sweep, ~4 min
.venv/bin/python scripts/build_protein_table.py            # protein-side companion table
.venv/bin/python scripts/build_clusters.py                 # wt_id corpus-wide + cluster inventory, ~1.5 min
.venv/bin/python scripts/make_overlap_report.py            # cross-source label agreement, ~5 s
.venv/bin/python scripts/build_label_health.py             # which records carry positive evidence, ~12 s
.venv/bin/python scripts/build_merged.py                   # one record per domain -> data/processed/, ~1 min
.venv/bin/python scripts/make_results.py                   # regenerate docs/RESULTS.md, ~5 s
.venv/bin/python scripts/build_matrix.py                  # domain x 8-mer arrays for modelling, ~5 s
.venv/bin/python scripts/build_distances.py               # all-vs-all domain identity, ~15 s
.venv/bin/python scripts/run_nn_baseline.py [--top-k]     # the bar -> reports/nn_baseline.md, ~10 s
.venv/bin/python scripts/build_embeddings.py --arm A1     # pooled ESM-2 vectors, ~20 s on the GPU
.venv/bin/python scripts/check_pooling.py                 # the ML_PLAN 3.1 pre-flight, ~10 s
.venv/bin/python scripts/train.py --arm A1 --fold P3/all  # one fold, ~8 min -- the fast check
.venv/bin/python scripts/run_grid.py [--seeds N]          # 19 folds x 2 arms -> reports/training.md, ~5 h
.venv/bin/python scripts/run_lambda_sweep.py --budget-hours 55  # T38 sweep, 42 h; resumable
.venv/bin/python scripts/measure_ceilings.py              # is the shared space wide enough, ~15 s
.venv/bin/python -m ruff check . && .venv/bin/python -m ruff format .
.venv/bin/python scripts/record_provenance.py data/raw/<source>/<file> --url ... --desc ...
```

### The GHT-SELEX arm (`docs/GHT_PLAN.md`) — a separate corpus, a separate pipeline

```bash
.venv/bin/python scripts/build_ght_panel.py                # 139 benchmark TFs -> 33 admitted, ~40 s
.venv/bin/python scripts/build_ght_windows.py              # 9.2M 301 bp windows out of hg38, ~5 min
.venv/bin/python scripts/build_ght_embeddings.py --arm A1  # 33 domains through ESM-2, ~10 s
.venv/bin/python scripts/check_ght_cobinding.py            # what `aliens` actually asks, ~2 min
.venv/bin/python scripts/check_ght_filters.py              # are the learned filters motifs, ~20 s
.venv/bin/python scripts/preflight_ght.py --steps 20000 --eval-every 250   # THE STEP BUDGET, ~25 min
.venv/bin/python scripts/run_ght_baselines.py              # PWM + motif transfer, ~7 min
.venv/bin/python scripts/run_ght_grid.py --seeds 3 --resume                # the grid, ~2.5 h
.venv/bin/python scripts/run_ght_grid.py --seeds 1 --protein-mode constant --resume   # the control
.venv/bin/python scripts/make_ght_report.py && .venv/bin/python scripts/make_ght_figures.py
.venv/bin/python scripts/make_ght_story.py                 # the talk page -> results/ght_story.html
```

Note: system Python 3.10 has no `ensurepip`, so `python -m venv` produces a venv without pip.
Use `uv` (already installed at `~/.local/bin/uv`).

## Phases

| phase | scope | status |
|---|---|---|
| 0 | scaffold, `PROVENANCE.md`, thresholds config, schema + validator | **done** |
| 1 | UniPROBE / Barrera `BAR15A` end-to-end, cluster inventory for it alone | **done** |
| 2 | remaining UniPROBE family panels | **done** — Cell08, EMBO10, PNAS13, then SCI09, GR09, MAR17A, SHO18A, ROG18A to rebuild breadth after the policy. Survey in `docs/UNIPROBE_ACCESSIONS.md`; `GB11` (27 bHLH) is the best remaining candidate. |
| 3 | ~~Persikov B1H + Najafabadi C2H2~~ | **DROPPED** — C2H2 arrays fail condition 2; Persikov varies a different subunit than the one that binds. ~8,000 domains excluded. |
| 4 | ~~SNP-SELEX, trimmed to a 19 bp window~~ | **DROPPED 2026-08-14** — the dataset is PBM only, so every row is an 8-mer and the `dna_len` leak cannot occur |
| 5 | merge, overlap report, splits, NN baseline | **done 2026-08-19** — merge and overlap on 2026-08-18 (`data/processed/training.parquet`, `reports/merge.md`, `docs/RESULTS.md`); splits and the NN baseline on 2026-08-19 (`reports/nn_baseline.md`) |
| — | **extend PBM coverage beyond UniPROBE** | **CLOSED 2026-08-18** — CIS-BP landed as `weirauch2014`; Kock 2024 screened and excluded (`reports/kock2024_excluded.md`). The corpus is UniPROBE + CIS-BP and grows no further |
| 6 | ~~Tier 4 test sets, in `data/testsets/`~~ | **DROPPED 2026-08-17** — the owner no longer wants the bHLH dimer sets. `data/testsets/` stays as empty scaffolding for any future held-out set |
| 7 | **modelling** — contrastive two-tower over Codebook SELEX + PBM | **the sequence arms are finished, 2026-08-30.** [docs/TRAINING.md](docs/TRAINING.md) is the design; [docs/ML_RESULTS.md](docs/ML_RESULTS.md) is what they measured and **§9 is the last word: the model beats `k = 5` by +0.035 AUPR, loses to `k = 1` by 0.006 call F1, and cannot detect a lost interaction.** The 312-run λ sweep closed the decision-rule question (`T38`) and relocated the open one to the protein representation (`T36`). What remains on this arm is incremental; the structure arms (`A2`/`A3`) are the next real step and are blocked on `T31`. Steps 0 and the sequence arms done 2026-08-19 ([reports/nn_baseline.md](reports/nn_baseline.md), [reports/pooling_check.md](reports/pooling_check.md)). Plan: [docs/ML_PLAN.md](docs/ML_PLAN.md). Build order is PBM + sequence embeddings first (no collaborator dependency); SELEX and structure ensembles are blocked on others. `S2` holds out **connected components at >= 0.5 identity**, not clusters (`T27`, then `D5` on 2026-08-19 — `docs/DECISIONS.md` §2) |
| 8 | **the GHT-SELEX arm** — protein-conditioned binding on human genomic windows | **built and measured 2026-09-15.** [docs/GHT_PLAN.md](docs/GHT_PLAN.md) is the plan, [docs/GHT_RESULTS.md](docs/GHT_RESULTS.md) the reading. Two settings: held-out chromosomes (`C1`, where a per-TF PWM competes) and held-out TF (`G1`/`G2`, where it cannot run at all). Every number is reported beside a **protein-blind control**, because 38% of a TF's peaks are co-bound by another panel TF and a model with no protein input scores well on those. `docs/DECISIONS.md` §14 |

## Deferred to the owner — flag, do not resolve

Tracked in [TODO.md](TODO.md), the working file for open tasks, open decisions and notes.
Currently open: whether 70% full-length UniProt coverage is enough (before embedding work),
and whether to drop `MAR17A:Esrrb` for its 2 unresolved `X` residues (before structure work).
Phase 7 adds four: the SELEX k-mer format to request (`T28`), how to score the transfer
experiment against weak negatives (`T29`), whether the 20 `dead_variant` records become the
point-mutation evaluation set (`T30`), and the structure-ensemble spec (`T31`).
**Both of the older two are now live** — Phase 7 *is* the embedding and structure work.
If parsing turns up evidence bearing on either, add it to that file rather than acting on it.

Phase 8 adds three: three Codebook raw files have **no source URL** and the owner is the only one
who can supply them (`T39`), whether to widen the 33-TF GHT panel (`T40`), and the
`FAMILY_ALIASES` defect — the alias is applied *after* the family filter, so `Homeobox_KN`,
`SOXp` and `zf-H2C2_2` never fire and the fix re-decides all 1,338 PBM domains (`T41`).

Anything decided or finished moves to [docs/DECISIONS.md](docs/DECISIONS.md) — read it before
reopening a question, since several obvious-looking ones have already been settled with
reasons.
