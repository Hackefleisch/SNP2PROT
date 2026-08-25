# ML results — what phase 7 measured, and what it means

The counterpart to [`ML_PLAN.md`](ML_PLAN.md). The plan says what we intended to do; this says
what happened when we did it, what it rules out, and what the evidence now points at.

Generated artifacts stay where they are — [`reports/nn_baseline.md`](../reports/nn_baseline.md),
[`reports/pooling_check.md`](../reports/pooling_check.md),
[`reports/training.md`](../reports/training.md) carry the tables. This document carries the
*reading* of them, which none of them can: they are each written by one script that sees one
experiment.

**Status as of 2026-08-25: the two-tower model does not beat the nearest-neighbour baseline
(3 wins in 38 runs), and the largest single reason found so far is that its protein tower is
linear** — a kernel ridge on the same `A4` embedding beats both the model and the baseline on
every `S2` fold measured. The pooled representation remains the deeper suspect but is no longer
established as the binding constraint. Everything below is the evidence.

---

## 1. The bar (step 0, 2026-08-19)

`ML_PLAN.md` §8.1 required the nearest-neighbour lookup before any model was trained, for two
reasons: to set the level every later number is read against, and to check that the split
regimes hold out what they claim. It did both, and the second job is where it earned its place.

Mean per-protein AUPR, `k = 1`:

| regime | AUPR | median NN identity |
|---|---:|---:|
| `S1` random domains | 0.769 | 0.742 |
| `S2` components at ≥ 0.5 identity | 0.296 | 0.386 |
| `P1` homeodomain holdout | 0.024 | 0.268 |
| `P2` non-homeodomain variant clusters | 0.882 | 0.824 |
| `P3/all` every variant, wild types kept | 0.928 | 0.986 |

**The whole curve is one variable.** Banding every held-out domain by how identical its nearest
training neighbour was gives AUPR 0.032 / 0.365 / 0.771 / 0.946 / 0.910 across the bands below
0.3, 0.3–0.5, 0.5–0.7, 0.7–0.9 and above 0.9. That single fact is the most useful thing the
baseline produced, and every later number has to be read against its own band rather than
against a pooled mean.

Two decisions came out of it, both the owner's, both recorded in
[`DECISIONS.md`](DECISIONS.md) §2:

- **`D5`** — `S2`'s edge threshold moved from 5 edits to **0.5 identity**. At 5 edits the
  grouping cost the lookup only 0.018 AUPR, because copying a 70–90% identical neighbour already
  scores 0.927; at 0.5 the lookup falls to 0.296 and the regime is genuinely hard.
- **`D6`** — claim **C1** is evaluated on the **29 variants where the wild-type copy fails**,
  because `P3/all` is otherwise saturated: mean 0.928 with a *median of 1.000*, so most single
  substitutions genuinely do not change what a domain binds and there is no headroom in the mean.

## 2. The pre-flight (2026-08-19)

`ML_PLAN.md` §3.1 required a measurement before any training: does mean-pooling a protein
language model to one vector per domain still leave a single-residue change visible?

| arm | variant vs family scale | own reference nearest | displacement ~ edits | C1 signal (AUC) |
|---|---:|---:|---:|---:|
| `A1` ESM-2 650M | 0.028 | 73% | 0.408 | 0.622 (1.3σ) |
| `A4` ESM-DBP | 0.026 | 79% | 0.410 | 0.709 (2.4σ) |

It passed on §3.1's own criterion — variants are separated from their wild types, the separation
scales with the number of edits, and 96% of variants sit within the nearest 5% of their family —
so pooling stayed and attention pooling was not reached for.

**In hindsight the number to have weighted more heavily is the first column: 0.028.** The signal
that distinguishes a variant from its wild type is under 3% of the scale that distinguishes
family members. §6 below is what that costs.

## 3. The grid (2026-08-25)

19 folds × 2 sequence arms = **38 runs, 105 minutes**. Design in
[`TRAINING.md`](TRAINING.md), per-fold tables in [`reports/training.md`](../reports/training.md).

**The model beats the baseline on 3 of 38 runs.** Mean delta −0.044.

| regime | A1 | A4 |
|---|---:|---:|
| `S1` | −0.074 | −0.062 |
| `S2` | −0.024 | −0.015 |
| `P1` | **+0.023** | −0.021 |
| `P2` | −0.161 | −0.132 |
| `P3` | −0.038 | −0.036 |

### 3.1 The identity bands are where the claim dies

Pooling every held-out domain across all folds and banding by nearest-neighbour identity:

| NN identity | n | A1 model | A4 model | baseline | A1 wins |
|---|---:|---:|---:|---:|---:|
| 0.0–0.3 | 706 | 0.050 | 0.030 | 0.032 | 64% |
| 0.3–0.5 | 1383 | 0.316 | 0.327 | 0.365 | 43% |
| 0.5–0.7 | 315 | 0.676 | 0.690 | 0.771 | 26% |
| 0.7–0.9 | 500 | 0.837 | 0.849 | 0.946 | 11% |
| 0.9–1.0 | 915 | 0.881 | 0.883 | 0.910 | 13% |

An early read of `S2/fold-0` alone suggested "the model wins where lookup fails". **Pooled, it
does not.** `A1` wins only the bottom band, by 0.018, at absolute values where both methods are
failing; `A4` does not win it at all.

**Losing the top bands is structural and expected.** A linear protein projection into a shared
256-d cosine space cannot reproduce "copy your nearest neighbour" — that needs per-protein
profiles memorised, and the shared space has nowhere to put 1,338 × 32,896 of them. The question
was always whether it wins where lookup is hard. It does not.

### 3.2 The one fold it won, and why

| `S2` fold | test set | median same-family training domains | delta |
|---|---|---:|---:|
| fold-0 | 273, **100% homeodomain** | **155** | **+0.198** |
| fold-2 | 266, 26 families | 54 | −0.033 |
| fold-3 | 266, 26 families | 15 | −0.107 |
| fold-1 | 267, 37 families | 11 | −0.112 |
| fold-4 | 266, 28 families | 9 | −0.066 |

fold-0 holds out the 273-domain homeodomain component and leaves 155 homeodomains in training at
< 50% identity. A model can learn a homeodomain-general mapping from 155 examples; copying a
< 50%-identical relative cannot. `P1` is the endpoint of the same axis — zero same-family
training data, AUPR 0.048.

**Stated as a hypothesis, not a finding.** The rank correlation between same-family depth and
delta across the twelve `S1`/`S2`/`P1`/`P2` folds is only **0.21**: fold-0 is an outlier whose
distinguishing feature is 155, while eleven folds spanning 9–59 same-family domains sit flat at
−0.02 to −0.15 with no trend. The honest claim is that 155 is 2.6× more than any other fold has,
and that this is the obvious candidate explanation.

### 3.3 C1 is the shifted prior `D6` predicted

| | A1 model | A4 model | baseline |
|---|---:|---:|---:|
| C1 set, 29 variants | 0.301 | 0.277 | 0.244 |
| the other 144 variants | 0.918 | 0.922 | 0.981 |

The model gains **+0.057** on the variants where the wild-type copy fails and loses **−0.063** on
the ones where it works — near-symmetric. `D6` wrote the test in advance: *"a model that had
merely learned to distrust wild types everywhere would gain here and lose there, which is a
shifted prior and not C1."*

**Claim C1 is not supported by this run.** Building the evaluation set before seeing a model is
what makes that statement possible rather than arguable, and it is the clearest return the
project has had on writing a caveat down early.

### 3.4 Two anomalies

**`A4` collapses on `P1`**: AUPR 0.003, AUROC **0.452**, Spearman −0.116, with **61%** of
held-out homeodomains scoring *below chance* against `A1`'s 8%. Systematically anti-correlated,
not merely uninformative — something specific goes wrong in ESM-DBP's space when every
homeodomain is removed from training. Unexplained.

**`A4` beats `A1` on 15 of 19 folds** but by **+0.006** on average. §4.2 built `A4` to separate
modality from pretraining corpus; on this evidence the answer is "the corpus barely matters",
with the `P1` blow-up as the one large exception.

## 4. What was ruled out, by measurement

Three candidate bottlenecks. The diagnostics discriminate between them.

### 4.1 The shared space is not too small

An **oracle** rank-D factorisation of the E-score matrix — the best any bilinear model could do
choosing protein factors freely, with sight of test labels:

| rank | macro AUPR |
|---:|---:|
| 16 | 0.422 |
| 64 | 0.723 |
| **256** | **0.917** |
| 512 | 0.977 |

`D = 256` can represent 0.917 and the model delivers 0.27 on `S2`. **Not the constraint.**

### 4.2 The tower is not too weak — it is already at the linear ceiling

A ridge regression straight from the pooled ESM vector to the full 32,896-dim E-score profile:
unconstrained output, no bottleneck, no DNA tower, `λ` tuned **on test**. A strict upper bound on
any linear-in-ESM model, our tower included.

| fold | ridge ceiling | our model | baseline |
|---|---:|---:|---:|
| `S2/fold-0` | 0.422 | **0.589** | 0.391 |
| `S2/fold-2` | 0.294 | 0.261 | 0.294 |
| `S2/fold-4` | 0.150 | 0.116 | 0.182 |
| `P1` | 0.029 | **0.048** | 0.024 |

**The two-tower model is at or above that ceiling**, beating an optimally-regularised
unconstrained probe on two of four folds. The 256-d bottleneck and the DNA tower are not costing
anything — they are regularising.

### 4.3 Non-linearity does help — on `A4`, and by a lot

> **Correction, 2026-08-25.** This section first said non-linearity bought ~0.02 and changed
> nothing. That was measured on `A1` only. Running the same probe on `A4`, as
> `scripts/measure_ceilings.py` now does for every arm, reverses the conclusion. The claim was
> wrong because the experiment was half-run, not because the numbers were misread.

RBF kernel ridge against linear ridge. `ceiling` picks `λ` and bandwidth on test; `honest` picks
them on the fold's own validation slice — and here the two agree to three decimals, so this is
**not** an oracle artefact:

| arm | fold | linear | RBF | RBF (honest) | two-tower | NN baseline |
|---|---|---:|---:|---:|---:|---:|
| `A1` | `S2/fold-0` | 0.475 | 0.486 | 0.486 | **0.589** | 0.391 |
| `A1` | `S2/fold-2` | 0.287 | 0.290 | 0.287 | 0.261 | 0.294 |
| `A1` | `S2/fold-4` | 0.130 | 0.148 | 0.148 | 0.116 | 0.182 |
| `A4` | `S2/fold-0` | 0.419 | **0.568** | 0.568 | 0.537 | 0.391 |
| `A4` | `S2/fold-2` | 0.237 | **0.342** | 0.342 | 0.278 | 0.294 |
| `A4` | `S2/fold-4` | 0.108 | **0.196** | 0.196 | 0.153 | 0.182 |
| `A1` | `P1` | 0.028 | 0.012 | 0.012 | 0.048 | 0.024 |
| `A4` | `P1` | 0.028 | 0.011 | 0.010 | 0.003 | 0.024 |

**On `A4`, a kernel ridge beats the nearest-neighbour baseline on all three `S2` folds** —
+0.177, +0.048, +0.014 — and beats the two-tower model on all three as well. `A1` shows almost
no non-linear gain. On `P1`, where no same-family training data exists at all, non-linearity
hurts both arms.

So the two arms are not interchangeable in the way §3.4's +0.006 average suggested. ESM-DBP's
space carries information a **linear** map cannot reach, and ESM-2's does not — which is
consistent with their geometry: `A4`'s cross-family cosine distances run 0.28 against `A1`'s
0.087 (`reports/pooling_check.md`), a far more spread-out space for a kernel to work in.

### 4.4 What the three diagnostics leave standing

- **The shared space is not the constraint.** Rank-256 can represent 0.917; nothing here is near
  it.
- **The protein tower being *linear* is a real limitation, for `A4`.** A non-linear function of
  the same vector gains 0.09–0.15 on `S2` and overtakes the baseline. Our tower does not.
- **The representation is still the deeper suspect, but it is no longer proven to be the
  binding constraint.** Even the best probe measured — `A4` + RBF at 0.568 / 0.342 / 0.196 —
  sits far below the 0.917 the space could hold, and collapses to 0.011 on `P1` where
  same-family data runs out.

§3.3's shifted prior and §2's 0.028 still point at the pooled representation, and §5's
data-density reading is untouched. What has changed is the **order of the experiments**: the
cheapest unexplored lever is now a non-linear protein tower on `A4`, which is one config flag,
and it must be tried before concluding anything about pooling.

**Things that still will not help:** more steps, larger `D`, more DNA channels, longer patience.
Those address links that measure as non-binding.

## 5. The alternative reading, which may be the right one

`ML_PLAN.md` §5.1 predicted this outcome and named this cause before anything was run:

> The expected cause is a **lack of overlap between proteins** — the panel is broad and shallow,
> so holding out a family removes anything the model could have transferred from.

1,338 domains in 1,165 clusters, 1,057 of them singletons. The single fold where the model won is
the single fold with 155 same-family training domains. On that reading the result is a
**data-density** result rather than a method failure, and no encoder change fixes it.

The two readings are not exclusive and they are distinguishable: a better protein representation
should lift the folds with 9–54 same-family domains toward what fold-0 achieves with 155. If it
lifts nothing, §5 is the answer and the honest headline for the talk is about what this corpus
can and cannot support.

## 6. What follows

**Tier 0 — a non-linear protein tower on `A4`.** §4.3's correction makes this the cheapest
unexplored lever and it must come first: `model.protein.hidden` is already in
`configs/experiment.yaml` and currently 0. A kernel ridge on the same embedding beats both our
model and the baseline on every `S2` fold measured, so the two-tower model is demonstrably
leaving reachable signal behind. Hours, not days, and it decides whether the pooling work is
even the next question.

**Tier 1 — pool over the DNA-contacting residues.** The protein table already stores each
domain's Pfam envelope and padding offsets, so the recognition positions are addressable without
new data or a new model. Replacing "mean over all ~77 residues" with "mean over the contacting
subset" is a biological prior doing attention's job for free. It tests §4.4 directly, costs a
rebuild of `build_embeddings.py` and one grid re-run, and needs nobody else.

**Tier 2 — attention pooling** (`TODO.md` `T32`). §3.1 named it as the fallback and said not to
reach for it until a measurement demanded it. The measurement now demands it. Still one
fixed-width vector per domain, so every reason for pooling survives and the arms stay comparable.

**Tier 3 — structure.** A predicted structure knows which residues face the DNA, which is the
same information Tier 1 approximates by hand and Tier 2 tries to learn. Arms `A2`/`A3` are
blocked on a colleague (`T31`), and this result raises their priority: they are no longer a
"does modality matter" curiosity but the most direct attack on the measured bottleneck.

---

## Appendix — provenance of every number here

| claim | source |
|---|---|
| baseline per regime, identity bands | `reports/nn_baseline.md`, `data/processed/nn_baseline_domains.parquet` |
| pooling pre-flight | `reports/pooling_check.md` |
| 38-run grid, per fold and per domain | `reports/training.md`, `data/processed/training_folds.parquet`, `training_domains.parquet`, and the MLflow store at `mlruns/mlflow.db` |
| rank ceiling, ridge probe, kernel ridge | `reports/representation_ceiling.md`, from `scripts/measure_ceilings.py` (`T35`) |
| same-family depth per fold | `snp2prot.splits` + `snp2prot.corpus` |

Every fold in the grid carries a `digest` of its held-out domains, and those digests match
`reports/nn_baseline.md` exactly — so model and baseline are scored on identical held-out sets
throughout, and the deltas above are like-for-like rather than a join on a fold name.
