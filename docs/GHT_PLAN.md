# The GHT-SELEX arm — plan of record and implementation spec

**Written 2026-09-15.** Two audiences: the project record, and a Claude Code session that will
implement this. Everything needed to write the code should be here. Where it is not, the gap is
marked **OPEN** rather than guessed.

**Why this exists.** The consortium criticised us for not collaborating and for building around a
proteoform limitation instead of working with partners' data. The project premise was plant
proteoform SELEX from a partner project that has not delivered. So we demonstrate the technology
now on human Codebook SELEX, so it can fire when the real data lands. The scientific question is
unchanged: **make the protein a variable rather than an index, so learned interactions transfer.**
Sequences now, structures next, dockings and complexes after.

Separate from [`ML_PLAN.md`](ML_PLAN.md), which is the PBM two-tower plan of record.
**Nothing here merges into `data/processed/training.parquet`.**

> **Implemented and measured 2026-09-15.** What it produced is
> [`GHT_RESULTS.md`](GHT_RESULTS.md); the decisions taken along the way are
> [`DECISIONS.md`](DECISIONS.md) §14; §9 of `GHT_RESULTS.md` lists what this plan asked
> for that was **not** built. This document is left as written, as the plan of record.

---

## 1. Background

**GHT-SELEX** is HT-SELEX on fragmented human genomic DNA instead of a random library. HEK293
genomic DNA was cut enzymatically to a **median ~64 bp**, size-selected below 200 bp, then run
through three or four cycles of bind, wash, elute, amplify, sequencing every cycle.

**A peak is a coverage pile-up and its shape is an artefact of random fragmentation.** A fragment
is retained if it spans the site. Fragments covering a site at X start at random offsets in the
~64 bases before X, so coverage is maximal at X and decays away. The apex is the **summit**.

**MAGIX** (Model-based Analysis of Genomic Intervals with eXponential enrichment, Najafabadi lab)
turns that into a number. It models abundance as multiplying by `e^b` per cycle, where `b` is the
log fold-change per cycle, conceptually tied to binding energy. Order of operations: pool all
reads into one coverage track, threshold, take the highest-coverage base per region as the summit,
define a **200 bp candidate peak** centred there, then fit `b` from per-cycle per-replicate counts
inside it. `b` is the `coefficient.br` column. `full_LL`/`reduced_LL` are the likelihood-ratio
test giving `pvalue`, corrected to `fdr`. `coefficient.ar` is **undocumented** in the README, both
papers and the MAGIX repo. Do not use it without asking the authors.

**MAGIX never scans the genome.** It only tests windows that already showed a pile-up, about half a
percent of hg38. A non-significant MAGIX row is a **near-miss, not background**. Measured across 51
peak files, 45% of rows pass FDR < 0.05 at the median and 97% at the extreme. **Never use MAGIX
non-significant rows as negatives.**

---

## 2. Decisions

| # | decision | rationale |
|---|---|---|
| D1 | Train on **Codebook GHT-SELEX** | The partners' platform, and the political requirement. PBM is the talk's second half |
| D2 | Use the **MEX-ArChIPelago benchmark release**, not raw MAGIX peaks | Ships positives, three negative sets, chromosome splits, PWMs and trained baselines. Comparable **by construction** |
| D3 | Positives = `positives.bed` | MEX peaks pooled across replicates, **301 bp centred on summits**, hg38 |
| D4 | Negatives = **`shades` primary, `random` secondary** | `shades` is locally composition-matched, which distribution matching cannot achieve. `random` exposes a composition shortcut |
| D5 | Splits = **theirs**, disjoint chromosome subsets | Stops memorising a repeat family or locus straddling the split, and keeps us comparable |
| D6 | Label = **binary** | Puts every baseline on one axis. Regressing `coefficient.br` is more elegant but no baseline occupies that axis |
| D7 | DNA encoder = **multi-width conv + global max pool over 301 bp** | See §5. Global max pooling is translation-invariant, so positional overfitting is structurally impossible |
| D8 | Protein encoder = **pooled ESM-2 (`A1`) as the default, plus an optional per-residue conv tower behind a config flag** | See §6. The conv is conceptually motivated; it is off by default on capacity grounds, not on the strength of the un-pooling report, which is contested |
| D9 | **Drop the PBM family-coverage filter** | Leftover from an abandoned zero-shot plan. Training on GHT makes coverage internal to the GHT panel |
| D10 | **Drop C2H2 arrays** (decided 2026-09-15) | They fail admission condition 2, one contiguous region. Costs most of the panel but keeps the policy honest |
| D11 | **Cite ArChIPelago, do not rerun** (decided 2026-09-15) | We follow their protocol faithfully — their positives, negatives and splits — so their published numbers apply. Rerunning buys nothing |

---

## 3. What `shades` is

For each positive summit, a fake upstream summit is drawn uniformly from **[−750, −450] bp** and a
fake downstream summit from **[+450, +750] bp**. A 301 bp window is centred on each, identically to
positives. Two per positive gives a nominal **2:1** balance; realised is often nearer **1:1**
because one flank is frequently rejected.

All negative types exclude positive regions, ENCODE blacklist regions, and any region containing
`N`. `random` and `aliens` are GC-matched to positives. `shades` is not, and does not need to be,
because local adjacency is a stronger control than distribution matching.

**The label is 0 by construction, not by measurement.** No experiment says the protein failed to
bind there. This is *absence of evidence* in the sense of rule 4 and categorically different from a
PBM negative, where a low E-score is a measurement against that specific 8-mer. If these ever share
a table with PBM rows they need distinct `neg_provenance`.

Fragments are ~64 bp and shades sit ≥450 bp away, so no fragment covering the true site can reach a
shade. The residual failure mode is a genuine second site in the flank.

> **Never compare auPRC across negative sets.** Chance auPRC is the positive fraction: ~0.33–0.50
> for `shades`, ~0.0099 for `random`/`aliens`. Their own paper notes `shades` gives the highest
> auPRC *because of the milder imbalance*. Same rule this project already applies per protein.

---

## 4. The panel

`TFs_CHS+AFS.txt` lists **139 TFs** with at least one approved ChIP-seq *and* one approved
GHT-SELEX dataset. After D10:

| | count |
|---|---|
| benchmark TFs | 139 |
| C2H2-containing, dropped by D10 | 92 |
| **usable** | **47** |
| of those, with construct AA sequence | 47 / 47 |

**Why the list is filtered on ChIP-seq at all, and why we accept it.** The intersection exists
because the benchmark was built for cross-platform transfer (train ChIP-seq, test GHT-SELEX and
back). We do not need transfer. Measured: there are **59** non-C2H2 GHT-approved TFs, **47 inside**
the benchmark and **12 outside**. All 12 have construct sequences and MAGIX peak files, but two of
them (`TIGD4`, `ZBED5`) are among the three withdrawn for mislabelled samples, leaving **10 clean
extras**, so 47 → 57, about +21%.

Taking them means building positives, three negative sets and chromosome splits ourselves for those
10, which reintroduces exactly the work D2 avoids and makes those 10 non-comparable to every
baseline. **Recommendation: stay at 47 for the talk**, and record the 10 as a known, cheap
extension. 47 is thin; §11 says how to show that honestly.

---

## 5. The DNA encoder

**Answering the shape question.** Input is one-hot `(B, 4, 301)`. A `nn.Conv1d(in_channels=4,
out_channels=C, kernel_size=w)` filter automatically spans all four channels — that is what
`in_channels` means — so the "4" is not a design choice. **`w`, the kernel length along the
sequence, is the design choice, and it is the motif width.**

**Use several widths, not one.** Real DBD motifs run roughly 6–20 bp; nuclear-receptor dimers and
POU sites are long, E-boxes short. A single-layer conv followed by global max pooling has no
composition across layers, so one width caps what a filter can express.

```
widths  w ∈ {8, 12, 16, 20}, 64 filters each  →  256 channels total
```

**On forcing `w = 8` for a PBM bridge.** Keep a `w = 8` bank so the option exists at zero cost, but
do not let it set `w` for the rest. Two reasons it should not be the plan. A PBM 8-mer is the
*entire measured unit* while a 301 bp window holds one site among ~293 positions, so the two are not
the same scale. And the existing `DNAEncoder` in `src/snp2prot/models/encoders.py` **deliberately
does not pool** — its docstring says eight positions is short and which base sits at position 3 is
the signal, so it flattens and projects. Global max pooling is the opposite choice. The honest
bridge between the arms is the **protein tower**, which is the scientific subject anyway.

**Pool with max over both strands, not mean.** The existing 8-mer encoder mean-pools the two strands
because it encodes a strand-symmetric *identity*. Here we are *scanning*, and the semantics we want
is "best hit on either strand", which is exactly what a PWM best-hit computes. So take the max over
the concatenation of forward and reverse-complement positions.

```python
# new class, e.g. GenomicDNAEncoder, in src/snp2prot/models/encoders.py
# reuse reverse_complement() already in that module
one_hot: (B, 4, L)  # L = 301
h_fwd = [conv_w(one_hot) for w in WIDTHS]  # each (B, 64, L - w + 1)
h_rev = [conv_w(rc(one_hot)) for w in WIDTHS]
h = [GELU(torch.cat([f, r], dim=-1)) for f, r in zip(h_fwd, h_rev)]
pooled = torch.cat([x.amax(dim=-1) for x in h], dim=-1)  # (B, 256) best hit per filter
z = F.normalize(self.project(pooled), dim=-1)  # (B, D), D = 256
```

**Why this is the right story as well as the right model.** A PWM best-hit score is exactly one
fixed filter plus global max pooling. So this encoder is the *learned generalisation* of the PWM
baseline, and ArChIPelago is a *fixed ensemble* of such filters feeding a random forest. Three
methods on one conceptual line.

---

## 6. The protein tower

**A convolution over residues makes conceptual sense, and the dimensionality worry is correct.**
Both are addressed here. The pooled embedding stays the default on **capacity** grounds, with the
conv tower available behind a config flag.

### 6.1 Why a conv is conceptually right

A DNA-binding domain contains local sequence patterns of *varying length*, at positions that are
not comparable across families. That is precisely what a bank of kernels of different widths
detects and what a single pooled vector cannot express. Concrete examples: the homeodomain's
invariant `WFQNRR` in helix 3, the bZIP basic region followed by a leucine-zipper heptad, the
Cys/His spacing that coordinates zinc in `zf-C4`, the ELK signature in Ets. These run roughly 3 to
15 residues. A width-3 kernel and a width-15 kernel are looking for genuinely different objects.

### 6.2 The dimensionality problem, and the fix

The naive form is badly unbalanced, exactly as suspected:

```
Conv1d(1280, 64, kernel_size=16)  =  1280 × 64 × 16  =  1,310,720 parameters   # ONE bank
```

At 47 training proteins that is ~28,000 parameters per protein. [`TRAINING.md`](TRAINING.md) §5
calls the PBM arm's **245 per protein** "the central engineering constraint". So the naive form is
two orders of magnitude out.

**Pointwise reduction first, as proposed:**

```
LayerNorm(1280)                          # per residue
Linear(1280 → d_r),  d_r = 64            # pointwise, shared across positions   81,920
Conv1d(64, 32, k) for k ∈ {3, 7, 15}     # multi-width banks                    51,200
amax over residues, mask-aware           # → (B, 96)
Linear(96 → D = 256)                                                            24,576
F.normalize
```

About **158k parameters**, an order of magnitude below naive, though still ~3.4k per protein at
n = 47. That is why it defaults off rather than on.

Four implementation notes:

- **Masking is required.** Domains vary in length, so pad and mask so `amax` never sees a pad
  position. A max over padding silently invents motifs.
- **Max over residues means "is this motif present anywhere".** That is the motif-detection
  semantics wanted, and it mirrors the DNA tower. It discards *where*, which matters for variant
  effects — attention pooling is the natural second experiment, not the first.
- **Needs per-residue embeddings.** `scripts/build_embeddings.py --per-residue` already produces
  them, roughly 600 MB per arm.
- **Treat it as a deliberate experiment with weight decay and dropout**, the same language
  `TRAINING.md` §5 applies to `protein_hidden`, and report both towers' parameter counts.

### 6.3 The un-pooling report is evidence, not a verdict

[`reports/unpooling_check.md`](../reports/unpooling_check.md) measured how far a single-residue
change moves the representation under four views, correlated against whether binding actually
changed:

| view | `A1` ρ change | `A4` ρ change | `A1` AUROC dead |
|---|---:|---:|---:|
| pooled | 0.094 | **0.209** | 0.526 |
| profile | **0.171** | 0.194 | 0.608 |
| window | 0.146 | 0.145 | 0.533 |
| site | −0.019 | −0.064 | 0.430 |
| BLOSUM62 control | 0.277 | 0.277 | 0.785 |

**It is one measurement and it is under active challenge.** n = 86 pairs, 18 of them dead, and
**56 of 86 are homeodomains**, so it is close to a statement about homeodomains rather than about
protein language models; the control's AUROC carries a standard error near 0.10. It is recorded
here for completeness and **is not the reason the conv tower defaults off** — capacity at n = 47
is. Do not cite it as settled, and do not let it decide whether the tower gets built.

The BLOSUM62 severity feature is worth adding either way, because it is two lines.

### 6.4 This is not `T36`

`T36` is a *non-linear* protein tower (`model.protein.hidden > 0`) on `A4`, raised because an RBF
kernel ridge beat both the baseline and the two-tower on three `S2` folds. Adjacent, cheap to try
here, and distinct from per-residue processing.

---

## 7. The loss

**The PBM InfoNCE does not transfer, and forcing it would be wrong.** That objective is a
per-protein softmax over the *complete, shared* 32,896-8-mer vocabulary, exact because the label
matrix is dense. GHT windows are **per-TF and essentially disjoint** — GCM1 and FLI1 share 34 of
~71,000 and ~145,000 — so there is no shared vocabulary to softmax over.

**Primary loss: binary cross-entropy on the calibrated logit.**

```
s(p, d)  =  cosine( protein_tower(p), dna_encoder(d) )      # both L2-normalised
logit    =  s * exp(calibration_scale) + calibration_bias   # TwoTower.calibrate
L_bce    =  BCEWithLogits(logit, y)                         # y ∈ {0, 1}
```

Reuse `TwoTower.cosine` and `TwoTower.calibrate` unchanged. Keep the **cosine-then-scale** form
rather than concatenating the two vectors into an MLP: that is what makes it a genuine two-tower
with a shared space, which is the architectural claim of the talk.

**Class weighting.** With `shades` at ~1:1 to 2:1, none needed. With `random` at 1:100, use
`pos_weight` in `BCEWithLogitsLoss` or subsample. **Set the calibration bias to the fold's actual
base rate** via a `set_calibration_prior`-style call, exactly as the PBM arm does — the existing
default `CALIBRATION_PRIOR_INIT = 0.002` is the PBM rate and is wrong here.

**Optional contrastive term, for continuity with the PBM arm.** Per protein, softmax over one
positive against that protein's own negatives is InfoNCE with `K` = the negative count:

```
L = L_bce + μ · L_infonce_per_protein          # μ = 0 by default
```

`μ = 0` must be *exactly* plain BCE so the knob is one line, mirroring how `bce_weight: 0` is
exactly the pure-ranking objective in the PBM arm. **Report `μ = 0` as the headline.** The
two-tower claim is architectural and does not need a contrastive loss to be true.

**The null anchor is not used here.** It exists to give a protein that binds nothing somewhere to
sit in a per-row softmax. With explicit labelled negatives and a BCE objective there is no row to
normalise, so leave `TwoTower.null` untrained and say so.

---

## 8. Training

Mirror `configs/experiment.yaml` conventions; every value carries its reasoning in a comment.

| setting | value | note |
|---|---|---|
| optimiser | AdamW, lr `1e-3`, weight decay `0.01`, warmup 100 | as the PBM arm |
| shared width `D` | 256 | equal across arms by `ML_PLAN.md` §4.2 |
| protein input | 1,280 | `esm2_t33_650M_UR50D`, frozen, precomputed |
| batch | mixed proteins and windows, **balanced across proteins** | unlike the PBM arm, which shuffles the protein axis and keeps the DNA axis complete |
| steps | measured, not guessed — see §8.1 | the PBM budget of 15,000 came from a curve |
| eval | every 100 steps | |
| early stopping | off, `patience: 0` | selection by validation, termination by budget |
| validation | 0.15, carved from **training chromosomes** | never from test chromosomes |
| seeds | `splits.seed` and `model.seed` separate | ≥3 seeds; GPU training is not reproducible across processes, ~0.003–0.012 AUPR spread |

Temperature: `learn_temperature: true`, `max_logit_scale: 100`, call `clamp_temperature()` after each
optimiser step. Log `temperature_is_clamped`.

### 8.1 Testing the time budget — do this first

**Two different budgets, measured in one pre-flight run.** Take one fold, set `patience: 0` so
nothing terminates early, run it well past where it looks done, and log validation *and* test at
every evaluation.

**(a) The step budget.** Read it off the **held-out test curve, never the training loss.** This is
`D11`'s lesson and the PBM arm got it backwards first: the training loss there falls to ~0.005 by
about step 6,000, which was read as "everything after this is memorisation" — yet held-out AUPR kept
improving for another 9,000 steps, worth +0.018 to +0.028. Memorising the training proteins and
generalising to unseen ones are simply not coupled. Record a table of step against train loss,
validation metric, test metric and temperature, and set the budget where the *test* curve flattens.
If it is still climbing at the end of the pre-flight, the pre-flight was too short.

**(b) The wall-clock budget**, so the grid is predictable before it is launched. Record, following
the shape of [`TRAINING.md`](TRAINING.md) §8.1:

| quantity | why |
|---|---|
| ms per step, steady state | multiplies out to the whole grid |
| cost of one evaluation pass | at `eval_every: 100` this is not negligible; in the PBM arm 97% of it was a numpy metric loop, not the forward pass |
| one fold, end to end | the unit the owner actually schedules |
| peak GPU memory | decides whether folds can run concurrently |

Then state the total for `folds × seeds` explicitly in the report. Under **rule 10** anything past a
few minutes is the owner's to launch, so this number is what gets handed over, not discovered
afterwards.

### 8.2 What a run records

The §6.2 convention of [`TRAINING.md`](TRAINING.md) applies unchanged: **log the split, not just the
hyperparameters.** Reuse the existing machinery rather than writing new.

- **`snp2prot.tracking.start_run(name, experiment_name, params, digests, enabled)`** — MLflow on a
  local SQLite backend under `mlruns/`, behind a wrapper thin enough to swap. It **raises without
  `digests`**, by design: a fold name and a seed do not reconstruct which proteins and which
  chromosomes were held out once the panel changes. Pass
  `digests={"test": ..., "validation": ...}` built from the held-out protein set *and* the held-out
  chromosome set — this arm has two holdout axes where the PBM arm has one, and both must be in the
  digest.
- **`tracking.code_version()`** — the commit and whether the tree was dirty, into every run and every
  report row. A digest pins the data, not the procedure.
- **Weights, both of them.** `TrainingResult` carries `best_state` (validation-selected) and
  `final_state` (at the step budget); they are different models and which is better is a per-fold
  empirical question. Write with `Trainer.save(path, result, meta)`, which stores both state dicts
  plus `model_config`, `protein_features` and `parameter_counts`, with `meta` carrying the digests and
  the code stamp. Path from **`config.checkpoint_file(arm, regime, fold, seed, tag)`**; use
  `weight_tag` for the variant suffix. Log the file as the run's MLflow artifact.
- **Metrics**: per-protein auROC and auPRC with their chance levels and null bands, logged per
  evaluation step via `Run.log_metrics`, plus the **distribution** and not only the macro-average —
  the point is which proteins fail, not that some do. Log both towers' trainable parameter counts and
  the matching nearest-neighbour number for the same fold, so the gap is computable without a join.
- **Artifacts**: a per-fold summary row and a per-domain frame, mirroring
  `training_folds.parquet` / `training_domains.parquet`, written under `data/processed/`. `run_fold`
  already returns exactly `(summary_dict, per_domain_df)`; follow it.
- **A generated report** under `reports/`, committed, in the style of `reports/training.md`, carrying
  the code stamp so a stale number announces itself.

Tracking must stay optional at the edges. If MLflow is absent the calls become no-ops and the run
still produces its report and its artifacts; nothing in the training path may depend on the tracker.

---

## 9. Evaluation — two settings, and the second is the point

**Setting 1 — held-out chromosomes, same TF.** Everyone competes: PWM, ArChIPelago (cited), 1NN,
5NN, us. **The goal is parity, not victory.** This tests the DNA encoder.

**Setting 2 — held-out TF.** PWMs and ArChIPelago are *fitted per TF* and cannot run on a protein
they have never seen. Only we and the nearest-neighbour lookup can play. This tests the protein
encoder and is the capability the talk exists to demonstrate.

Protein-level folds for Setting 2, mirroring the PBM regimes: a random-TF split (`S1` analogue) and
an identity-component split at ≥0.5 identity via `snp2prot.distances` (`S2` analogue).

Metrics per TF: auROC and auPRC, macro-averaged, **each reported against its own chance level** plus
a sampled random-ranking 95th percentile, as `snp2prot.evaluation.metrics` already does. Report both
negative sets separately and never compare auPRC across them.

**Baselines.** PWM sum-occupancy from `pwms.zip` (top-20 log-odds, already platform-separated against
leakage). ArChIPelago by citation (D11). 1NN and 5NN adapted: score a window by whether the nearest
training domain, by `snp2prot.distances` identity, has an overlapping peak.

---

## 10. Data inventory

**On disk, verified.**

| path | size | state |
|---|---|---|
| `data/external/genomes/hg38.fa.gz` | 983,659,424 B | gzip OK, `>chr1`, UCSC naming matches the BEDs |
| `data/raw/codebook_ghtselex/Peaks_MAGIX_McGill.tar.gz` | 300,569,175 B | 207 BEDs, one per TF |
| `data/raw/codebook_ghtselex/GHT_SELEX_Metadata_…xlsx` | 749,020 B | sha256 `743d7940…` |
| `data/raw/codebook_htselex/TF_and_Plasmid_Metadata_…xlsx` | 805,835 B | construct AA sequences |
| `data/raw/codebook_htselex/` reads | 9.2 GB | 684 FASTQs, complete, **unused by this plan** |

**To fetch**, Zenodo `10.5281/zenodo.10515307`, record 15555251, into `data/raw/codebook_ghtselex/`:

| file | size |
|---|---|
| `datasets.zip` | 7,368,594,592 B — positives, three negative sets, chromosome splits |
| `pwms.zip` | 12,766,984 B — the PWM baseline |
| `TFs_CHS+AFS.txt`, `datasets_approved.txt` | 895 B, 44,137 B — in scratchpad, move in |

`ARCHIPELAGO_models.zip` (71 MB) is **not needed** under D11.

---

## 11. Build order

**Phase A, data.** Unpack `datasets.zip`. Extract 301 bp sequences from hg38 for positives, `shades`
and `random`. Run construct amino-acid sequences through `snp2prot.domains` to get a padded Pfam
domain per TF — **the same admission path every other source uses. Do not string-match family
labels; that produced a wrong panel count once already.** Build `A1` embeddings.

**Phase B, model and baselines.** Add `GenomicDNAEncoder`. Wire a `ght:` block into
`configs/experiment.yaml`. Train on their chromosome splits, then the protein-held-out folds. Run
PWM, 1NN, 5NN.

**Phase C, figures.** One panel per setting, plus the PBM half, which needs no new computation.

Under **rule 10** the training runs are the owner's to launch; hand over exact commands and timings.

**Fallback ladder** — each rung presentable alone.

| rung | claim | state |
|---|---|---|
| 0 | PBM results only | **exists today** |
| 1 | We ingested their data, their benchmark, their splits | cheap |
| 2 | Works on held-out chromosomes | proves the pipeline |
| 3 | **Parity with PWM and ArChIPelago in Setting 1** | the credibility result |
| 4 | **Generalises to a held-out TF, where their methods cannot run** | the real result |
| 5 | Protein tower shared with the PBM arm | bonus |

---

## 12. Risks, open items, obligations

- **The step and wall-clock budgets are unmeasured.** §8.1 is the pre-flight that settles both and
  it is the **first** thing to run. Do not inherit 15,000 steps from the PBM arm, and read the test
  curve rather than the training loss.
- **The conv protein tower is a capacity risk, which is why it defaults off.** ~158k parameters
  against 47 training proteins is ~3.4k per protein, where the PBM arm's linear tower sits at 245
  and `TRAINING.md` §5 already calls that the central constraint. Turn it on as a deliberate
  experiment with weight decay and dropout, against the pooled default as the control, and report
  both towers' parameter counts either way.
- **47 proteins is thin** for Setting 2. Show per-fold variance rather than a single mean. The 10
  clean extras in §4 are the cheapest way to widen it if needed.
- **Family transfer fails** in our PBM results: holding out every homeodomain gives 0.0125 AUPR
  against 0.0031 chance. Expect Setting 2 to be carried by families with training depth.
- **Composition shortcut** is the realistic failure mode since `shades` is not GC-matched. Reporting
  `random` alongside catches it.
- **Construct mismatch**: many peak files come from full-length constructs while the model is fed the
  isolated domain.
- `coefficient.ar` is undocumented. Do not use it.
- **`PROVENANCE.md` has no rows** for `codebook_htselex`, `codebook_ghtselex` or the hg38 build. Rule
  3 requires them. The HT-SELEX row needs per-file sha256 across 684 files; follow the
  `BAR15A/details/` precedent of a `_manifest.tsv` plus one summary row.
