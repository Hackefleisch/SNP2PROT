# Training design — the two-tower contrastive model

How the model in [`ML_PLAN.md`](ML_PLAN.md) §3 is actually trained: the objective, the batch
shape, the null anchor, and what "the 19-fold grid" means. Written 2026-08-19, before any of it
is implemented, so that the design can be argued with while it is still cheap to change.

Nothing here overrides the plan. Where this document makes a choice the plan left open, it says
so and gives the reason; where it disagrees with the plan it says that too.

---

## 1. What is being trained

Two towers and one shared space (`ML_PLAN.md` §3):

| tower | input | shape |
|---|---|---|
| **protein** | the cached pooled PLM vector for a domain — arm `A1` or `A4` | 1280 → `D` |
| **DNA** | one-hot `8 × 4`, both strands, mean-pooled | 8×4 → `D` |

Both outputs are L2-normalised, and the score of a (domain, 8-mer) pair is their inner product
over a temperature:

```
s(p, k) = ⟨ f(p), g(k) ⟩ / τ
```

**The PLM is not part of the model.** `data/processed/embeddings/A1.npz` and `A4.npz` are already
built (`reports/pooling_check.md`); the protein tower reads a 1280-d vector and never runs a
language model. That is what makes the whole grid an afternoon rather than a week, and it is why
the arms are interchangeable: swapping `A1` for `A4` swaps one `.npz` and changes nothing else.

**The DNA tower is a small CNN, not a lookup table** (`ML_PLAN.md` §4.1). Reverse-complement
mean pooling is applied because the 8-mer vocabulary is revcomp-collapsed — verified again here:
32,896 = (4⁸ + 4⁴)/2 stored 8-mers, 256 of them palindromic, and **zero** revcomp pairs with both
strands present. Each stored 8-mer is one arbitrary representative of a pair the array measured
together, so an encoder that treated `AAAAAAAC` and `GTTTTTTT` as different inputs would be
carrying a strand asymmetry the assay does not have.

---

## 2. The objective

### 2.1 What we are actually optimising

The evaluation metric is fixed (`ML_PLAN.md` §5.3): for each held-out domain, rank all 32,896
8-mers, score that ranking with AUPR, macro-average over domains. So the training objective
should be *the same shape*: a per-protein ranking of the complete 8-mer space, with every protein
weighted equally.

Multi-positive InfoNCE over the full 8-mer axis is exactly that, which is the argument for it —
not that contrastive learning is fashionable.

### 2.2 The loss, stated exactly

For domain `p`, write

- `P_p = { k : Y[p,k] = 1 }` — its positives, median 49, mean 69, max 236, 5th percentile **1**
- `N_p = { k : Y[p,k] = 0 }` — its negatives, median 32,264
- gray cells (`Y = -1`, median 552 per domain) are in **neither**
- `ν` — the null anchor, §3 below

```
                 1                    exp s(p, k⁺)
L_p  =  −  ─────────  ·   Σ     log  ──────────────────────────────────────
              |P_p|     k⁺ ∈ P_p      exp s(p, k⁺)  +  Σ      exp s(p, k)
                                                    k ∈ N_p ∪ {ν}
```

and the batch loss is the mean of `L_p` over the domains in the batch.

### 2.3 Three choices inside that formula, and why

**(a) The denominator is exact, not sampled.** Every InfoNCE implementation you will read
approximates the negative set, because in the usual setting the negatives are the rest of the
corpus and there are too many. Here we *hold the complete label matrix* — 1,338 × 32,896, every
cell present — so there is nothing to approximate. Sampling a subset of negatives we already know
buys nothing and costs bias.

**(b) Sibling positives are NOT in the denominator.** This is the one place the design departs
from the textbook (SupCon's `L_out`, which keeps all positives in the denominator). The reason is
concrete: a domain with 236 positives, under `L_out`, has a loss that floors at `log|P_p|` and a
gradient that pushes its own positives *apart from each other*. That is right for representation
learning and wrong for retrieval. With sibling positives excluded, each positive competes only
against the negatives, the loss reaches zero exactly when every positive outranks every negative,
and that condition is AUPR = 1. The objective and the metric then agree at their optimum.

**(c) The `1/|P_p|` and the batch mean are the same normalisation the metric uses.** Positives per
domain span 1 to 236. Without per-domain normalisation a single 236-positive domain outweighs 236
one-positive domains, and the model would be trained on an objective the reported metric does not
measure. `ML_PLAN.md` §3.1 requires this.

### 2.4 How it is computed

Written as above it looks like `O(|P_p| · K)`. It is not — the denominator's negative part is
shared across all of a domain's positives, so with

```
logC_p  =  logsumexp   s(p, k)          (over k ∈ N_p ∪ {ν})
```

the whole thing collapses to

```
L_p  =  mean over k⁺ ∈ P_p  of   softplus( logC_p − s(p, k⁺) )
```

One `logsumexp` over the masked row, one `softplus` over the positives. Numerically stable, `O(K)`
per domain, and it is about six lines of PyTorch.

Gray cells are excluded by setting their logits to `−inf` before the `logsumexp`, so they enter
neither the numerator nor the denominator. **A mask, not a filter** — the matrix stays dense and
rectangular, exactly as `ML_PLAN.md` §3.1 specifies, and different domains masking different
columns costs nothing.

### 2.5 The calibration term — added 2026-08-28, and why it had to be

**Everything above optimises a ranking, and the deliverable is a call.** `T38` put it plainly:
every metric in this project is rank-based, so nothing ever required the model to commit, and
handing a biologist a model that says *these 8-mers score highest* rather than *these 8-mers
bind* is not the same product. Worse, the capability that matters most — noticing that a
variant has lost binding **entirely**, which is what flags a pathogenic mutation — is a
statement about a protein, and no per-protein ranking can express it. A protein that binds
nothing still has a first-ranked 8-mer.

**The diagnosis is one line, and it is structural.** §2.2's loss is a softmax over one row, so
`L_p` depends only on `logC_p − s(p, k⁺)` — differences *within* a row. **Add a constant to every
entry of a row and the loss does not change by a float.** The objective is exactly invariant to a
per-protein offset, so it never says where a row should sit, only how it should be ordered.

Three consequences, and the third is the one that closes off the easy fixes:

- the null anchor could not have become a threshold, whatever the arithmetic of §3.3 — there is
  no force in the objective that would place it between the classes rather than anywhere else;
- a *global* Platt scaling fitted on the validation slice cannot work either, and measured, it
  did not: one cut at precision ≥ 0.5 called between 56 and 408 8-mers per domain for proteins
  whose true counts ran 10 to 206;
- **no post-hoc step can recover an offset the training signal never constrained.** The fix has
  to be in the loss.

**The fix.** A second term that is *not* shift-invariant, and the cheapest one is binary
cross-entropy:

```
L  =  L_infonce  +  λ · L_bce
```

with `L_bce` the masked, **unweighted** BCE over the same cells the InfoNCE term uses, read
through a two-scalar calibration head — `p(p,k) = σ(α · cos(p,k) + β)`. `snp2prot.models.loss`
holds both terms and `hybrid_loss` combines them.

Four choices inside that, each with a reason:

**(a) Its own `α` and `β`, not the temperature.** Reusing `logit_scale` looks economical and is
wrong twice: it is clamped at `max_logit_scale` and reaches that ceiling by about step 4,000, so
from there it is a constant and could not move to fit a probability; and its job is how peaked
the softmax is, which is a statement about ranking dynamics and unrelated to where the
binding/non-binding line falls. Two jobs, two parameters. Neither is weight-decayed —
`calibration_bias` starts at the logit of the base rate, about −6.2, and decaying it toward 0 is
decaying it toward `p = 0.5`.

**(b) Unweighted, deliberately.** At 466:1 the reflex is to class-weight or reach for focal loss.
Both re-bias the output away from the empirical rate, which is exactly the thing being asked for
— a weighted BCE is calibrated to a class balance that does not exist. The ranking is already
supplied by the InfoNCE term; this term's only job is to be honest about the base rate.

**(c) `β` is initialised at the *training pool's* positive rate.** With 466 negatives per
positive a head starting at `p = 0.5` spends thousands of steps discovering that almost
everything is negative, and the cheapest way down is to flatten `α`, which destroys the ordering
the other term is building. Starting at the prior is the standard fix for detection-scale
imbalance. It is read off the fold's **training** rows, never the test ones.

**(d) Per-domain mean, then mean over domains** — the same normalisation as §2.3(c), so every
protein counts once in both terms and in the reported metric alike.

**`λ = 0` is §2.2's objective exactly.** The calibration head takes no gradient, and
`hybrid_loss` returns the InfoNCE value itself — verified on the last bit, on the CPU, by
`tests/test_loss.py::test_lambda_zero_is_the_pure_ranking_loss_to_the_last_bit`. That is what
makes this one knob rather than a redesign.

**It does not follow that a `λ = 0` run reproduces an old run's numbers, and it does not.**
Measured 2026-08-28: two processes running identical code with identical seeds diverge at the
first evaluation — losses `8.96644592285` and `8.96644687653` at step 50 — and the gap compounds.
That is CUDA's own non-determinism (the DNA tower's `Conv1d` backward accumulates with atomics)
and it predates this work entirely; end to end on `A1`/`S1/fold-2` it is worth about **0.003
AUPR**, and the sweep's three seeds put the within-cell spread at **0.005 (max 0.012)** — larger
than several deltas this project has reported. `configs/experiment.yaml`'s "two identical runs
agree to 4e-6" held within one process and not across invocations; it is corrected there.
**This is the reason the λ sweep carries its own `λ = 0` control on every fold** rather than
comparing against `reports/training.md`.

**The scale of `λ` is not 1, and the grid has to be logarithmic.** Measured at initialisation on
`A1`/`P3/all`: the InfoNCE term is **10.44** and the BCE term **0.0166** — the entropy of a
0.0021 Bernoulli — a ratio of **631**. And the ratio does not hold still: InfoNCE falls to ~0.005
by step 6,000 while the BCE term plateaus, so a `λ` that is negligible at step 1 can dominate by
step 10,000. `scripts/run_lambda_sweep.py` sweeps `λ ∈ {0, 1, 3, 10, 30, 100, 300, 1000, 3000,
10000}` for that reason, and `reports/calibration.md` is what it writes.

**What the term buys, and how it is measured.** `snp2prot.evaluation.calibration`:

| question | metric |
|---|---|
| is the probability honest | `expected_calibration_error`, quantile-binned |
| is the call good | `score_calls` — precision / recall / F1 / Jaccard of the named set |
| how many should we name | `expected_count_rule` — keep the top `round(Σ p)`, no free parameter |
| **did the interaction survive** | `interaction_power` = `Σ p`, and `detection_auroc` over it |
| the paired form of the last | `power_ratio` — the variant's `Σ p` over its wild type's |

The fourth row is the one `T38` was raised to obtain, and the nearest-neighbour baseline now
carries all of them too (`reports/nn_baseline.md`, "As a decision, not a ranking") — `k = 1`
emits a *set* natively, so this is the first comparison in the project where the baseline and
the model are scored on the same kind of output.

---

## 3. The null anchor

### 3.1 What it is

**A learned vector `ν ∈ R^D` on the DNA side that is not any 8-mer.** In implementation terms the
DNA embedding table is 32,897 rows rather than 32,896: the extra row is a free parameter, trained
by the same gradient as everything else, and present in **every** domain's denominator.

### 3.2 The problem it solves

A domain with no positives at all should be training signal. It says *"this protein binds none of
these 32,896 8-mers"*, which is a real statement about the shared space — it should push every
8-mer away from that protein.

Plain multi-positive InfoNCE cannot express it. Look at §2.2: the loss is an average over `P_p`.
If `P_p` is empty the average is over nothing, the row contributes no gradient, and it is silently
dropped. The row is not useless; the standard loss simply has nowhere to put it.

With the null anchor, an all-negative row is handled by **treating `ν` as its positive**:

```
                       exp s(p, ν)
L_p  =  − log  ────────────────────────────────
                exp s(p, ν) + Σ       exp s(p, k)
                               k ∈ N_p
```

which drives every real 8-mer's score down relative to a fixed reference point. And this is not a
special case bolted on — it is the same expression. Define

```
P̃_p  =  P_p        if P_p is non-empty
         { ν }      otherwise
```

and §2.2 written over `P̃_p` covers both, with the negative pool being `N_p ∪ {ν} \ P̃_p`. One
uniform loss, one code path, no branch on "does this domain have positives".

### 3.3 Its second job — which was the plan, and which does not work

The intention was this. Even for a domain with positives, `ν` sits in the denominator as one more
competitor, so the model must rank a real positive **above the null**, and `ν` would settle at a
position meaning *"good enough to call binding"* — a **learned decision threshold**, global rather
than per protein, obtained for free and calibrated jointly with everything else.

**Measured 2026-08-27 on a trained checkpoint, it does not happen** (`D9`). The anchor sits at
−24.3 where the negatives' median is −26.4 and the positives' median is +9.6, and a median of
**15,159 of 32,460 8-mers score above it**. It lands inside the negative cloud, roughly at its
middle, not between the two classes.

The reason is structural rather than a tuning failure. For the ~1,318 domains that *have*
positives, `ν` is simply another negative being pushed down; the only force raising it comes from
the domains that have none, where it is the target. §3.4 counts those: **at most two rows per
fold.** The anchor loses 1,318 to 2, and settles where that arithmetic puts it.

**Nothing may read `s(p,k) > s(p,ν)` as a binding call.** The metric for the domains that needed a
threshold is `snp2prot.evaluation.metrics.suppression`, which compares a dead variant with its own
wild type down the protein axis and needs no threshold at all.

**And the arithmetic above was not the whole reason.** §2.5 records the deeper one, found
2026-08-28: the loss is *exactly invariant to a per-protein offset*, so even an anchor that won
the 1,318-to-2 vote would have had no defined place to settle. The threshold now comes from a
second term rather than from `ν`, and `ν` keeps only its first job.

### 3.4 The honest statement: on PBM today it is nearly inert

`ML_PLAN.md` §3.1 says the mechanism applies to "exactly 54 records". Measured against the actual
split regimes, with `no_evidence` dropped by `label_health.usable` and the C1 set held out, the
number of all-negative rows **remaining in training** is:

| regime | training domains | all-negative rows among them |
|---|---:|---:|
| `S1`, per fold | 1,017 – 1,023 | 1 – 2 |
| `S2`, per fold | 1,016 – 1,023 | 0 – 2 |
| `P1` | 870 | 2 |
| `P2` | 1,133 | 0 |
| `P3`, per fold | 1,132 – 1,243 | 2 |

**At most two rows.** The 34 `no_evidence` records are excluded because a silent record with no
control cannot be told from a failed assay (`T21`). The 18 dead variants are held out under
`P3/all` by construction, which is where C1 is scored; the training-pool mask that used to remove
them everywhere was dropped on 2026-08-27 because the baseline never applied it (`D7`).

So: build the null anchor, but **do not expect it to move a PBM number**, and do not report it as
if it had. Two of the three reasons for building it survive — it removes a branch from the loss,
and it is the mechanism SELEX will need, where all-negative proteins are expected to be common
(`ML_PLAN.md` §2.2). The third, §3.3's threshold, turned out not to exist. If `T30` is later
decided in favour of training on the 34 `no_evidence` records, the mechanism is already there and
that decision becomes a config flag.

---

## 4. Batch structure — the question, answered

Three ways to structure this, and the choice matters more than usual because of the shape of the
data.

### 4.1 The three options

**(a) Row-wise over (domain, 8-mer) pairs.** Shuffle the 44,014,848 rows, take minibatches of
pairs. This is what a default dataloader does and it is wrong here on every count: it needs an
in-batch negative-sampling scheme to reconstruct a denominator we already hold exactly, it makes
the per-domain normalisation of §2.3(c) impossible because a batch contains a random assortment of
partial rows, and it does 44 million rows of I/O for something that fits in 42 MiB of RAM.

**(b) One protein per step.** Take one domain, score it against all 32,896, step. The objective is
right but the optimisation is not: batch size 1 on the protein axis gives a very noisy gradient to
the protein tower, and the DNA tower is updated from a single protein's perspective at every step.

**(c) Mixed batches on the protein axis, complete on the 8-mer axis.** ← **this is the design.**

### 4.2 The design

**Each step draws a random batch of `B` domains (64 by default) and scores every one of them
against all 32,896 8-mers plus the null.** So:

- **shuffling is on the protein axis only** — that is where the stochasticity comes from;
- **the 8-mer axis is never subsampled** — the denominator stays exact, per §2.3(a);
- the step produces a `64 × 32,897` logit matrix, and §2.4's loss is one masked `logsumexp` along
  its rows.

So the answer to "row-wise per protein or mixed batches" is: **mixed batches of proteins, each
protein complete**. Not mixed batches of pairs, and not one protein at a time.

The reason this is the right unit is that it makes the training objective **structurally identical
to the evaluation metric**. §5.3 scores a model by ranking one protein's 32,896 8-mers and
macro-averaging over proteins; §2.2's loss is a smooth surrogate for exactly that ranking, averaged
over exactly those proteins. There is no train/eval mismatch to reason about.

### 4.3 What it costs — measured, not estimated

| tensor | shape | dtype | size |
|---|---|---|---:|
| label matrix, resident | 1,338 × 32,896 | `int8` | **42.0 MiB** |
| DNA embedding table, per step | 32,896 × 256 | `fp32` | 32.1 MiB |
| logits, one step | 64 × 32,897 | `fp32` | **8.0 MiB** |
| one conv layer's activations, 128ch, both strands | — | `fp32` | ~257 MiB |

**The logits are not the memory driver — the DNA tower's activations are.** The whole 32,896-8-mer
vocabulary goes through the CNN on every step, twice (both strands), and its intermediate
activations are held for the backward pass. Three conv layers at 128 channels is roughly 800 MiB
of activations before the backward graph, which fits the RTX 3070's 8 GB but is the number to
watch when the encoder is sized. Channel count, not batch size, is the knob.

Compute is negligible: the `64 × 256 × 32,896` score matmul is 0.5 GFLOP, and the DNA tower is a
few hundred MFLOP. A step is well under a GFLOP against a card that does ~20 TFLOP/s, so **the
training loop will be bound by kernel launch overhead, not arithmetic**. There is no dataloader,
no I/O, and nothing to prefetch — everything is resident (`ML_PLAN.md` §3.1).

### 4.4 An "epoch" is not a useful unit here

1,338 domains at `B = 64` is **21 steps per epoch**. A 1,000-step run is ~48 passes over the
protein axis. Runs should therefore be specified in **steps**, and the thing that limits them is
overfitting, not compute.

### 4.5 The one real weakness, and the fallback

Because the 8-mer axis is complete at every step, **every DNA embedding takes a gradient at every
step**. There is no stochasticity at all on the DNA side — the DNA tower sees the identical 32,896
inputs every time, and only the protein side varies.

`ML_PLAN.md` §3.1's fallback is **sampled softmax**: keep all of a domain's positives, subsample
its negatives (say 4,096 of ~32,264) per step. That restores DNA-side stochasticity at a cost
worth almost nothing, since the negatives are 99.8% of the row. Two notes if it is used:

- it makes the denominator an estimate, so it re-introduces the bias §2.3(a) avoided. For a
  ranking objective the bias is generally tolerable, and a `log Q` correction removes it if it
  turns out not to be;
- it is **a knob, not a redesign** — same loss, same batch shape, one mask.

Start without it. Reach for it if training is unstable or the DNA tower overfits.

---

## 5. Capacity — the part most likely to go wrong

The protein axis has **1,338 points**. The input to the protein tower is 1,280-dimensional. A
single linear projection to `D = 256` is **327,680 parameters — 245 per training protein.**

**The tower begins with a `LayerNorm` over that input, added 2026-08-27** (`D12`). It is not a
capacity decision — it costs 2,560 parameters and changes no shape — but without it the arm
comparison measured the wrong thing. ESM-2's pooled vectors have norms 4.83–9.85 and ESM-DBP's
0.76–1.24, so on the projection's bias `‖b‖/‖Wx‖` was **0.13 for `A1` and 1.16 for `A4`**: for
one arm the bias was a correction, for the other it outweighed the signal and set the output
direction. `A4`'s embeddings are the better separated of the two (mean pairwise cosine 0.752
against `A1`'s 0.874), yet the untrained tower emitted **0.8930 for both, identical to four
decimals** — it had flattened away exactly the difference the `A1` → `A4` delta exists to
measure. With the `LayerNorm` they arrive at 0.891 and 0.759, and the tower is invariant to a
rescale of its input to ~1e-5. `ML_PLAN.md` §4.2's "equal terms" is now true rather than
assumed.

That asymmetry is the central engineering constraint and it points in one direction:

- **the protein tower must be small.** A linear projection is the default. A one-hidden-layer MLP
  should be treated as a deliberate experiment with weight decay, not as the obvious choice;
- **the DNA tower can be comparatively larger.** Its axis has 32,896 points and every one is seen
  at every step, so it is in a completely different data regime from the tower it shares a space
  with;
- `D` is a regulariser as much as a capacity choice, and it is shared by construction — `ML_PLAN.md`
  §4.2 requires every arm to project to the same `D` so that "structure wins" can never be
  "structure had more width".

**Report each tower's trainable parameter count with every result**, which §4.2 asks for explicitly:
equal `D` controls the shared space but not the capacity feeding it.

---

## 6. The 19-fold grid

### 6.1 What the 19 are

Exactly the folds `snp2prot.splits.all_regimes` already produces, and exactly the ones the
nearest-neighbour baseline was scored on (`reports/nn_baseline.md`):

| regime | folds | what is held out |
|---|---:|---|
| `S1` | 5 | random domains — the ceiling |
| `S2` | 5 | connected components at ≥ 0.5 identity — genuinely unseen proteins (`D5`) |
| `P1` | 1 | the whole Homeodomain family |
| `P2` | 1 | the non-homeodomain variant-bearing clusters |
| `P3` | 7 | all variants (×1), half (×3 seeded draws), a quarter (×3) |
| | **19** | |

**19 runs per arm.** Two arms are buildable today (`A1`, `A4`) → **38 runs**; four when the
structure ensembles land → 76, plus the TransBind re-implementation baseline (`ML_PLAN.md` §8.2).
`half` and `quarter` are repeated over seeded draws because which variants get dropped matters at
n = 86 and n = 43, and a single draw is a sample of one.

At the cost established in §4.3 the whole 38-run grid is a small number of GPU-hours. The grid is
*enumerated*, not searched — which is why `ML_PLAN.md` §9.2 rejected Hydra and Optuna and adopted
a plain loop over configs.

### 6.2 What a run records

`ML_PLAN.md` §9.2: **log the split, not just the hyperparameters.** A metric is meaningless
without knowing which domains were held out, and "regime `S2`, fold 3, seed 20260819" does not pin
that down once the corpus changes. So every run logs:

- the regime, the fold name, **both seeds** — `splits.seed` decides which domains are held out,
  `model.seed` the initialisation and the batch order, and they were one number until 2026-08-27
  — and **`Fold.digest`, the hash of the held-out domain sequences**, which `snp2prot.splits`
  already computes and the baseline report already prints;
- **the commit it was produced from, and whether the tree was dirty** (`tracking.code_version`).
  A digest pins the domains but not the procedure: the `OVERSHOOT` guard landed in the same
  commit as the first version of `reports/training.md`, so one row of that table came from a
  `validation_split` that no longer existed and the only way to notice was to recompute the
  carve and compare row counts;
- **the trained weights themselves** — both of them, the validation-selected checkpoint and the
  model at the step budget, to `data/processed/checkpoints/` and as the run's MLflow artifact.
  Without them a grid is 38 numbers and no models: nothing can be probed, re-scored on a new
  metric, or asked what it actually learned;
- the arm and the embedding file it read;
- every hyperparameter, from `configs/experiment.yaml`;
- both towers' trainable parameter counts (§5);
- the per-protein metrics **and their distribution**, not only the macro-average — `ML_PLAN.md`
  §5.3 wants the violin, because it shows *which* proteins fail rather than that some do;
- the corresponding NN-baseline number for the same fold, so the model-minus-baseline gap is
  computable per fold without a join against another artifact.

MLflow on a local `mlruns/` backend, behind a thin `snp2prot.tracking` wrapper so that "MLflow,
W&B, whatever" stays a one-line swap (§9.2). `mlruns/` is already in `.gitignore`.

### 6.3 Evaluation is the baseline's shape exactly

For a held-out domain: one protein-tower forward, score against all 32,896 DNA embeddings, rank.
That is the same `(n_test, 32896)` array `nn_lookup.fit_predict` returns, so it goes through the
same `snp2prot.evaluation.metrics.score_domain` with no special-casing, and the two are directly
comparable on every figure. The C1 evaluation set (`D6`) is scored the same way, on the same 29
variants, with its own absolute number reported next to the complement's — never as
"model beats baseline", which on that set is vacuous by construction.

---

## 7. Two things that need deciding before the grid runs

**(a) There is no validation split inside a fold.** Each of the 19 folds is a train/test division
and nothing else, so any early stopping would be early stopping *on the test set*, which leaks.
Two ways out:

1. **a fixed step budget, no early stopping** — simplest, defensible given how cheap a run is, and
   it removes a per-fold decision that could differ between arms;
2. **carve a validation slice out of the training domains**, using the same grouping rule as the
   regime (whole components for `S2`, so the validation set is as separated from training as the
   test set is), seeded, and logged with its own digest.

**Settled 2026-08-27: both, and neither alone.** The slice is carved as in (2) and still selects
a checkpoint, but it no longer *stops* the run — `patience: 0`, every fold trained to the fixed
budget of (1). The reason is §8.2: stopping fired on 5 of 38 runs, all `S2`, cutting the regime
the project hangs on off at step 2,750 with its best at 750. And because both the selected model
and the model at the budget are now kept and scored, whether the slice is a *useful* selector is
a reported number per fold instead of a premise — which matters most exactly where §7a said it
would, on `P1` and `P3`, where no slice of the training pool can imitate the test task.

**One measured cost of the slice** is recorded here rather than left to be rediscovered: on
`S2/fold-3` validation selected step 14,500 (test 0.2473) over the true peak at 12,000 (test
0.2633), a 0.016 selection loss. On `S1/fold-2` the best validation step *was* the best test
step. The new "Does selecting on validation beat the model at the budget?" section of
`reports/training.md` is where that is now read off rather than argued about.

**(b) `T30`'s 34 `no_evidence` records.** Still open, and §3.4 above is the concrete reason it now
matters less than it looked: with them excluded, the null anchor has almost nothing to do on PBM.
Including them would give it real work — and would risk teaching the model that a perfectly good TF
binds nothing. No action needed to start; the flag exists either way.

---

## 8. What the first runs measured

Built and run 2026-08-20. The numbers below are from `A1` on `P3/all` — chosen first because
its nearest-neighbour baseline was then **0.928** with a median of 1.000, so a wrong
implementation is obvious immediately rather than after the grid. (That figure is the old
continuous-profile baseline; the matched binary bar on the same fold is **0.660** — `D7`. The
fold is still the right first check for the same reason.)

### 8.1 Cost, measured

| | |
|---|---|
| step, steady state | **26.7 ms** (batch 64, complete 8-mer axis, forward + backward) |
| of which the DNA table forward | 7.3 ms — so the backward through it is the bulk |
| validation pass, 160 domains | 0.36 s |
| one fold at 15,000 steps | **~7.8 min** (was ~2.5 min at the old 6,000) |
| peak GPU | **825 MiB** |

So the 19-fold grid is about 2.5 hours per arm, and 38 runs across `A1` and `A4` is about
**5.1 hours** — 4.5 h of training plus 0.6 h of evaluation at `eval_every: 100`. 97% of that
evaluation cost is the per-domain AUPR loop in numpy rather than the forward pass (388 ms of
402 ms); batching the argsort across domains would cut it to roughly 0.1 h if it ever matters.
§4.3's prediction that activations rather than logits would dominate memory held: 825 MiB
against 8 MiB of logits.

### 8.2 The step budget, set from the curve rather than guessed

15,000 steps on `P3/all`, evaluating both validation and test every 500:

| step | train loss | validation AUPR | test AUPR | τ |
|---:|---:|---:|---:|---:|
| 500 | 4.10 | 0.457 | 0.636 | 0.041 |
| 1,500 | 2.33 | 0.549 | 0.793 | 0.020 |
| 3,500 | 0.49 | 0.636 | 0.854 | 0.011 |
| 6,000 | 0.57 | 0.635 | 0.856 | 0.010 |
| 10,000 | 0.06 | 0.640 | 0.872 | 0.010 |
| 15,000 | 0.007 | 0.656 | 0.881 | 0.010 |

**This table was read wrongly, and the correction is `D11`.** The original reading was
*"validation is flat from about step 3,500 while the training loss keeps falling to 0.007, so
everything past that is fitting the training proteins"*, and the budget was set to 6,000. But the
**test** column in the same table climbs from 0.854 at 3,500 to 0.881 at 15,000 — the conclusion
was drawn from the validation column while the answer sat next to it.

Re-measured 2026-08-26 on three folds with early stopping disabled, comparing windows of three
evaluations *within* a single run so no cross-run variation enters:

| fold | mean test AUPR @5.5–6.5k | @13.5–15k | gain |
|---|---:|---:|---:|
| `P3/all` | 0.8563 | 0.8798 | +0.024 |
| `S2/fold-3` | 0.2275 | 0.2457 | +0.018 |
| `S1/fold-2` | 0.6584 | 0.6863 | +0.028 |

All three were **still climbing at 15,000**. The training loss does reach ~0.005, so the training
proteins are memorised by about step 6,000 — but held-out AUPR improves for another 9,000 steps
after that. **Memorising the training set and generalising to unseen proteins are not coupled
here**, so "the loss is near zero" was never evidence for stopping. §5's capacity argument
predicts overfitting and none is visible in the test column at any budget tried.

Cosine LR decay was measured as the cheaper alternative — if the wobble at the end of training
were the limit, decay would fix it without more compute. It is slightly *worse* at equal budget
(0.8438 against 0.8458 on `P3/all`), so the budget is the lever and not the schedule.

**The default is therefore 15,000 steps with `patience: 0`.** Early stopping fired on 5 of 38
runs in the first grid and all five were `S2` — `A1 S2/fold-4` stopped at step 2,750 with its
best at 750, while `S2/fold-3` was allowed to run and was still gaining at 12,000. Patience
against a signal that swings 0.63→0.66 between adjacent evaluations was cutting runs off rather
than protecting them. The validation slice still selects a checkpoint; it no longer ends the run,
and **the model at the budget is kept and scored beside the selected one** so the selection is a
reported number rather than an assumption.

### 8.3 The temperature reaches its ceiling and stays there

τ falls from 0.07 to **0.0100 by about step 4,000 and pins**. That is the CLIP clamp at
`1/τ ≤ 100` binding — the model wants a sharper softmax than the clamp allows, which is what a
separable training set looks like. Raised to 10,000 the scale never converges: it passes 100 at
about step 4,000 and is still climbing at 148 by 15,000. The clamp is a real stopping device.

**No reported number depends on it, and that was measured rather than assumed.** `score` is
`scale × cosine` and every metric in `snp2prot.evaluation.metrics` ranks *within* one domain, so
a positive scalar cannot reorder anything — verified to 8 decimals, and end to end the fold
scores 0.8558 clamped at 100 against 0.8554 running free to 148.

**What it does set is which negatives the model learns from** (`D13`). A negative's share of the
gradient is its softmax weight in the shared denominator, so at the converged scale a median of
**2 of 32,460 measured negatives carry half the negative-side gradient**; the top 100 carry
99.7%. At the initial τ = 0.07 it is ~1,500. That is not evidence being discarded — a negative's
share is proportional to how wrongly close it still sits, so the other 32,458 are satisfied
constraints, and a low temperature is implicit hard-negative mining. It is kept low deliberately.

Two mechanical corrections came with this. The clamp is applied to the **parameter after each
optimiser step**, not inside `score`: clamping the forward pass makes the gradient exactly 0
above the ceiling, so the parameter either oscillated across the boundary on weight decay alone
or, with the decay removed, froze at 4.60882 to five decimals with no force acting on it. And
`logit_scale` and `null` are excluded from weight decay — neither is a weight, and CLIP excludes
the first for the same reason.

### 8.4 Where it stands against the bar

**The comparison this section made is withdrawn** (`D7`, 2026-08-27). It set the model's 0.881 on
`P3/all` against a baseline of 0.928 — but that baseline copied the neighbour's *continuous
E-score profile*, and the model trains on thresholded labels. The two were never scored on equal
information: between 29% and 76% of that baseline came from the ordering inside the copied
profile. Against the matched binary bar, `P3/all` is **0.660** and `S2/fold-0` is **0.145**.

A universal-PBM E-score is a rank-enrichment statistic read at a cutoff, not a graded affinity,
so its ordering is not a quantity to predict — which is why the corpus stores a binary label and
why the baseline now copies binary calls.

**No conclusion should be drawn until the grid has rerun.** What can be said from the single-fold
checks: `A1` reaches 0.8873 on `P3/all` at the 15,000-step budget, and on the C1 set — the 18
dead variants, scored by `suppression` — it sits at **0.4921 against a chance level of 0.5 and a
baseline of exactly 0.0**. On the claim the project exists to make, the model is at chance.

### 8.5 One deviation from the plan, for a reason outside it

`ML_PLAN.md` §9.2 specified MLflow with "a plain `mlruns/` directory". **MLflow now refuses the
filesystem store** — it is in maintenance mode and raises unless an opt-out environment variable
is set. The backend is therefore a local **SQLite** file at `mlruns/mlflow.db`, with artifacts
under `mlruns/artifacts`. That satisfies every requirement §9.2 actually gave — local, no server,
no cloud, one git-ignored path — and is the supported route rather than a deprecated one. The
change is to the URI and to nothing else.

---

## 9. Build order

All **done 2026-08-20** except the last, which is the owner's to run (`CLAUDE.md` rule 10).

1. ✅ **`configs/experiment.yaml`** — `model:` and `training:` blocks beside the existing `splits:`.
2. ✅ **`snp2prot/tracking.py`** — the wrapper, which refuses a run that does not record what it
   held out.
3. ✅ **`snp2prot/models/`** — `encoders.py` (CNN + revcomp pooling, projection), `two_tower.py`,
   `loss.py`.
4. ✅ **`snp2prot/training.py`** and **`scripts/train.py`** — the loop and the single-fold entry
   point, verified on `P3/all`.
5. ✅ **`scripts/run_grid.py`** — the 19 folds × the built arms, joined against the baseline on
   the same folds. Verified end to end on `P1` and `P2`; the fold digests match
   `reports/nn_baseline.md` exactly, so model and baseline are scored on identical held-out sets.
6. ⬜ **`reports/training.md`** — produced by the full grid run, which is about two hours:

   ```bash
   python scripts/run_grid.py          # 19 folds x 2 arms
   ```
