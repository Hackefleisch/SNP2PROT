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

### 3.3 Its second job, which is the one that matters on PBM

Even for a domain with positives, `ν` sits in the denominator as one more competitor. So the model
must rank a real positive **above the null**, and the null learns a position in the shared space
that means *"good enough to call binding"*. That is a **learned decision threshold**, obtained for
free and calibrated jointly with everything else — `s(p,k) > s(p,ν)` is a binding call without a
separately fitted cutoff.

It is deliberately **global, not per protein**. A per-protein threshold would fit the sensitivity
differences the corpus is known to have (`docs/METHODS.md` §5.1: two labs differ 1.6× in positive
count on one protein), but it could not be estimated for a held-out protein, which is the only
kind we score. A global anchor is the version that transfers.

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
control cannot be told from a failed assay (`T21`), and the 18 dead variants that would otherwise
qualify are the core of the C1 evaluation set (`D6`) and must never be trained on.

So: build the null anchor, but **do not expect it to move a PBM number**, and do not report it as
if it had. It is worth building for three reasons that survive that fact — it removes a branch
from the loss, it provides §3.3's threshold for all 1,338 rows, and it is the mechanism SELEX will
need, where all-negative proteins are expected to be common (`ML_PLAN.md` §2.2). If `T30` is later
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

- the regime, the fold name, the seed, and **`Fold.digest` — the hash of the held-out domain
  sequences**, which `snp2prot.splits` already computes and the baseline report already prints;
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

Recommendation: **(2)**, because a fixed budget tuned once on `S1` will be wrong for `P1`, whose
training set is 870 domains and whose task is much harder. But it costs a real slice of an already
small protein axis, so it is worth a decision rather than an assumption.

**(b) `T30`'s 34 `no_evidence` records.** Still open, and §3.4 above is the concrete reason it now
matters less than it looked: with them excluded, the null anchor has almost nothing to do on PBM.
Including them would give it real work — and would risk teaching the model that a perfectly good TF
binds nothing. No action needed to start; the flag exists either way.

---

## 8. What the first runs measured

Built and run 2026-08-20. The numbers below are from `A1` on `P3/all` — chosen first because
its nearest-neighbour baseline is **0.928** with a median of 1.000, so a wrong implementation is
obvious immediately rather than after the grid.

### 8.1 Cost, measured

| | |
|---|---|
| step, steady state | **26.7 ms** (batch 64, complete 8-mer axis, forward + backward) |
| of which the DNA table forward | 7.3 ms — so the backward through it is the bulk |
| validation pass, 160 domains | 0.36 s |
| one fold at 6,000 steps | **~2.5 min** |
| peak GPU | **825 MiB** |

So the 19-fold grid is about 50 minutes per arm, and 38 runs across `A1` and `A4` is under two
hours. §4.3's prediction that activations rather than logits would dominate memory held: 825 MiB
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

**Validation is flat from about step 3,500 while the training loss keeps falling to 0.007.**
Everything past that is fitting the training proteins — which is exactly what §5 predicted from
1,338 points against a 1,280-d input. The default is therefore **6,000 steps**, comfortably past
the plateau and short of the region where only the training loss moves.

### 8.3 The temperature reaches its ceiling and stays there

τ falls from 0.07 to **0.0100 by about step 4,000 and pins**. That is the CLIP clamp at
`1/τ ≤ 100` binding, not a coincidence — the model wants a sharper softmax than the clamp allows,
which is what a separable training set looks like.

The clamp is doing its job (a runaway scale saturates the softmax early and kills the gradient),
so the default is unchanged. But it means **the temperature is a fixed hyperparameter from step
4,000 onward, not a learned one**, and that should not be silent: `max_logit_scale` is now in
`configs/experiment.yaml` and every run reports whether the clamp is binding. Raising it is a
one-line experiment, and the expectation is that it increases overfitting rather than the score.

### 8.4 Where it stands against the bar

At 15,000 steps `A1` reaches **0.881** on `P3/all` against the baseline's **0.928** — so on the
fold where the baseline is strongest, the model does not yet beat copying the wild type. That is
one fold, one arm, at a default configuration nothing has been tuned on, and it is the expected
place to be hardest: `P3/all`'s baseline has a *median* of 1.000, meaning most single
substitutions genuinely do not change what a domain binds. `D6` exists precisely because the
claim does not live in that mean.

**No conclusion should be drawn until the grid has run.** The fold that matters for the headline
is `S2`, where the baseline is 0.296.

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
