# ML results — what phase 7 measured, and what it means

The counterpart to [`ML_PLAN.md`](ML_PLAN.md). The plan says what we intended to do; this says
what happened when we did it, what it rules out, and what the evidence now points at.

Generated artifacts stay where they are — [`reports/nn_baseline.md`](../reports/nn_baseline.md),
[`reports/pooling_check.md`](../reports/pooling_check.md),
[`reports/training.md`](../reports/training.md),
[`reports/representation_ceiling.md`](../reports/representation_ceiling.md) carry the tables.
This document carries the *reading* of them, which none of them can: they are each written by one
script that sees one experiment.

**Every number here is from the grid of 2026-08-28** — 38 runs, commit `38b7aa1`, clean tree,
15,000 steps, 4.9 hours. It replaces the grid of 2026-08-25 in full, which was invalidated by the
modelling audit ([`DECISIONS.md` §12](DECISIONS.md)); nothing from that run survives here.

> **Those runs were all `model.bce_weight = 0`, which is now a setting rather than the only
> behaviour.** §8's limitation was acted on the same day: the objective gained a calibration term
> and a 312-run sweep of `λ` is measuring what it costs and what it buys
> ([`DECISIONS.md` §13](DECISIONS.md), [`reports/calibration.md`](../reports/calibration.md)).
> **`λ = 0` is the objective every run below used**, so nothing here is invalidated — but §8 is no
> longer the last word on it: **§9 is the result**, and it relocates the problem rather than
> solving it.

---

## The headline, in three sentences

**The two-tower model beats a nearest-neighbour lookup on 19 of 19 folds, but by +0.005 to +0.052
— not the +0.16 to +0.25 it appears to win by against a single-neighbour lookup.** Both baselines
are now in [`reports/nn_baseline.md`](../reports/nn_baseline.md) and
[`reports/training.md`](../reports/training.md); `k = 5` is the column to read. Four fifths of
that apparent margin is the model being asked for a ranking while the baseline is asked for a set;
against a five-neighbour lookup built from the same binary labels, the advantage is modest and on
`S2` it wins three folds and loses two. **Claim C1 is directionally supported and quantitatively far
short**: on variants that still bind, asking the model for the *variant* is no better than asking
for its wild type (−0.011 mean, 34% improve), and the two predictions are 98% rank-correlated;
on the 18 that lost binding it does call fewer actives than their wild types in 78% of cases
against a 48% base rate, but only ~8% fewer where the right answer is ~100%.

---

## 1. The bar

`snp2prot.baselines.nn_lookup` copies the most identical training domain's **binary calls**, under
the alignment overlap guard. That it copies calls rather than E-scores is the correction `D7` made
and it is the whole basis of the comparison: the model trains on thresholded labels, so a baseline
copying the continuous profile is scored on information the model is never given — worth 29–76% of
its AUPR. See [`DECISIONS.md` §12.1](DECISIONS.md).

**Two forms, both binary.** `k = 1` copies the nearest neighbour's calls — a *set*. `k = 5`
averages the calls of the five nearest, weighted by identity — a *ranking*, built from nothing but
binary labels. The second is the one that matters, and §2 explains why.

| regime | `k = 1` (a set) | `k = 5` (a ranking) | chance |
|---|---:|---:|---:|
| `S1` random domains | 0.4753 | 0.6797 | 0.0022 |
| `S2` components at ≥ 0.5 identity | 0.1345 | 0.2560 | 0.0022 |
| `P1` homeodomain holdout | 0.0058 | 0.0071 | 0.0031 |
| `P2` non-homeodomain variant clusters | 0.5741 | 0.7521 | 0.0018 |
| `P3/all` every variant, wild types kept | 0.6603 | 0.8398 | 0.0026 |

## 2. The margin is real but small — and most of the apparent win is output format

| arm | regime | model | bar | delta |
|---|---|---:|---:|---:|
| `A1` | S1 | **0.7017** | 0.4753 | **+0.2264** |
| `A1` | S2 | **0.2938** | 0.1345 | **+0.1593** |
| `A1` | P1 | **0.0125** | 0.0058 | +0.0067 |
| `A1` | P2 | **0.7649** | 0.5741 | **+0.1908** |
| `A1` | P3 | **0.8852** | 0.6387 | **+0.2466** |
| `A4` | S1 | **0.7130** | 0.4753 | **+0.2377** |
| `A4` | S2 | **0.2865** | 0.1345 | **+0.1520** |
| `A4` | P1 | 0.0047 | 0.0058 | −0.0010 |
| `A4` | P2 | **0.7946** | 0.5741 | **+0.2204** |
| `A4` | P3 | **0.8859** | 0.6387 | **+0.2472** |

37 of 38 runs beat that bar; the exception is `A4` on `P1`, where both numbers are near chance.
**Every fold also scored above its own random-ranking 95th percentile** — the first grid had one
row that did not.

**But `k = 1` outputs a set and the model outputs a ranking, and most of that gap is the
difference between the two, not learned biophysics.** A `k = 1` prediction is ~50 tied positives
above ~32,000 tied negatives; AUPR then measures set overlap and nothing else. Controlling for it
costs nothing — `k = 5` averages five neighbours' calls by identity, producing a graded score from
**binary labels only**, the same information the model has and the same shape of output:

| regime | `k = 1` | `k = 5` | model `A1` | vs `k = 5` |
|---|---:|---:|---:|---:|
| `S1` | 0.4753 | 0.6797 | 0.7017 | **+0.0220** |
| `S2` | 0.1345 | 0.2560 | 0.2938 | **+0.0378** |
| `P1` | 0.0058 | 0.0071 | 0.0125 | +0.0054 |
| `P2` | 0.5741 | 0.7521 | 0.7649 | +0.0128 |
| `P3` | 0.6387 | 0.8332 | 0.8852 | **+0.0520** |

**The honest margin is +0.005 to +0.052, not +0.16 to +0.25.** Roughly four fifths of the apparent
advantage over `k = 1` was the model being asked for a ranking while the baseline was asked for a
set.

And per fold it is not uniform. Under `S2` the model wins three folds and **loses two** —
`fold-1` by 0.003 and `fold-3` by 0.044 — with the regime mean carried almost entirely by
`fold-0`, where it wins by 0.215. On one seed, with the seed spread unmeasured (§6), `S2` is
better described as *the model and a five-neighbour lookup are close, with one fold where the
model is clearly ahead* than as a win.

## 3. What replaced this section

An earlier draft compared the model against a lookup copying the neighbour's **continuous E-score
profile**, called it "the incumbent method", and read the model as losing to it on 15 of 19 folds.
That comparison is removed and will not return. A universal-PBM E-score is a rank-enrichment
statistic that the field reads at a cutoff — the binarisation is done by the researchers who
produce it, not by us — so its ordering is not a quantity to predict, and a baseline built on it
is scored with information that exists nowhere else in this project (`D7`,
[`DECISIONS.md` §12.1](DECISIONS.md)). **The whole corpus is binary: training, baseline,
evaluation.** §2's `k = 5` is the ranked comparison, and it is built from labels alone.

## 4. Claim C1: a weak signal in the right direction, an order of magnitude short

C1 is that the model can infer which 8-mers a *mutated* TF engages, without having seen that
mutation. The evaluation set is the 41 held-out variants whose own wild type does not predict them
(`D6`), and it splits in two — the two halves behave completely differently.

| group | n | bar | `A1` | `A4` |
|---|---:|---:|---:|---:|
| C1, **poorly predicted** (AUPR defined) | 23 | 0.0787 | **0.5934** | **0.6090** |
| C1, **dead** — scored by `suppression` | 18 | 0.0000 | **0.4785** | 0.4567 |
| the other 132 variants | 132 | 0.7624 | 0.9266 | 0.9377 |

**The first half works.** On 23 variants where copying the wild type scores 0.079, the model
scores 0.593. It is recovering binding profiles that a lookup gets badly wrong, which is real.

**The second half does not, and it is the half that matters.** `suppression` asks: of the sites
the wild type binds, what fraction does the model rank *lower* in the variant? **1.0** means it
saw the mutation abolish binding; **0.5** means the sites moved at random; **0.0** is what copying
the wild type gives by construction. The model scores **0.4785** and **0.4567** — at chance, on
both arms, over 18 variants whose binding measurably vanished.

So the model has *not* learned that a point mutation destroyed binding. It rearranges those
profiles no better than chance. This is the sharpest evidence the corpus holds for C1 and it is
negative, and it is not softened by the first half: predicting a *changed* profile better than a
lookup is a weaker claim than detecting that binding was *abolished*.

**Note what this is not.** It is not the "shifted prior" the old grid reported. A model that had
merely learned to distrust wild types everywhere would gain on the C1 set and lose on the
complement; this one is better on the complement (0.927) than on C1 (0.593), which is what a model
that finds C1 genuinely harder looks like.

## 5. Family transfer fails

`P1` holds out all 428 homeodomains. Chance is 0.0031 and the random-ranking 95th percentile is
0.0036.

| | AUPR | × chance | above the null? |
|---|---:|---:|---|
| `A1` | 0.0125 | 4.0× | yes |
| `A4` | 0.0047 | 1.5× | barely |
| lookup, `k = 1` | 0.0058 | 1.8× | — |
| lookup, `k = 5` | 0.0071 | 2.3× | — |

`A1` is measurably above chance and roughly twice the lookup, but 0.0125 is not a working model of
an unseen fold; `A4` is indistinguishable from guessing. **Nothing here transfers across a Pfam
family.** The `LayerNorm` fix (`D12`) was tried against exactly this and moved `A4` from
0.0035 to 0.0047 — the norm mismatch was real and was not the cause.

## 6. `A1` versus `A4` cannot be read

| regime | `A1` | `A4` | `A4 − A1` |
|---|---:|---:|---:|
| S1 | 0.7017 | 0.7130 | +0.0113 |
| S2 | 0.2938 | 0.2865 | −0.0073 |
| P1 | 0.0125 | 0.0047 | −0.0078 |
| P2 | 0.7649 | 0.7946 | +0.0296 |
| P3 | 0.8852 | 0.8859 | +0.0006 |

ESM-DBP is ahead on `S1`, `P2` and `P3`, behind on `S2` and `P1`. **Every one of those differences
is within the run-to-run noise, and that noise is unmeasured.** The grid ran at one seed;
`model.seed` now exists separately from `splits.seed` and two seeds of one fold differed by 0.016
at a 300-step budget — larger than four of the five numbers above. `T37` measures it with
`run_grid.py --seeds 3`. Until then **no claim about the pretraining corpus is supported either
way**, which is the honest state of the `C2` comparison.

What *is* visible is a consistent pattern rather than a magnitude: `A4` is better where training
contains same-family relatives and worse where it does not, matching its embeddings being more
family-clustered (within-Homeodomain cosine 0.963 vs 0.731 across, against `A1`'s 0.933 and 0.874).

## 7. What is ruled out

**The shared space is not the constraint.** A rank-256 factorisation of the **binary label**
matrix reaches **0.9379** in sample, against 0.294 achieved on `S2` and 0.885 on `P3`. A representable solution
far better than anything trained exists at `D = 256`. This is a *lower* bound on what rank 256 can
express, and a lower bound licenses a conclusion when it comes out high — which it does.
([`representation_ceiling.md`](../reports/representation_ceiling.md).)

**The step budget is no longer obviously the constraint, but it may still bind.** 15,000 steps,
raised from 6,000 for a measured +0.02–0.03 (`D11`). Median best step is **13,900**, and **18 of
38 runs selected a checkpoint at or past 14,000** — so the curve has flattened but has not clearly
stopped. A further increase is cheap to test and has not been tested.

**Validation selection is close to a wash.** Selecting on the validation slice beats taking the
model at the budget by **+0.0046** overall, and loses on 15 of 38 folds. It is worth keeping —
it costs nothing now that both models are stored — but it is not doing much work.

**Two things this grid cannot rule out**, because the probes that claimed to were deleted for being
lower bounds read as upper ones (`D10`): whether a stronger or non-linear protein tower would help,
and whether the pooled representation is the ceiling. Those are `T36` and `T34`, and both are
settled by ablation — an achievable number on the same folds — rather than by a probe.

## 8. The limitation that qualifies everything above

**The model produces a ranking. The dataset, the task and any biologist using it need a call, and
the model has no rule for making one.** Raised by the owner on 2026-08-28 and open as
[`TODO.md`](../TODO.md) `T38`.

Every metric here — AUPR, AUROC, precision@k, R@P0.5 — is rank-based, so none of them ever asks
the model to commit to *this 8-mer binds, that one does not*. Turning the output into the binary
answer the assay produces needs a threshold, and the two candidates both fail:

- the **null anchor**, which was designed to be exactly this, calls a median of **15,396 of
  32,896** 8-mers positive on the validation domains, against a true median of 44 (§ and
  [`DECISIONS.md` §12.2](DECISIONS.md));
- a **global cut calibrated on validation** at precision ≥ 0.5 averages sensibly — ~65 calls per
  domain against a true median of 44 — but per domain it calls between **56 and 408**, an order of
  magnitude of spread, for proteins whose real counts run 10 to 206.

At 466:1 negative-to-positive the capability that matters is discriminating true non-binders, and
that is precisely what an uncalibrated threshold fails to deliver. **So a per-protein AUPR of 0.70
says the model orders one protein's 8-mers well relative to each other; it does not establish that
its scores mean the same thing across proteins, and the measurements above suggest they do not.**
Read every number in this document with that attached.

**Why neither candidate could have worked, found 2026-08-28.** The loss is a per-row softmax, so
it depends only on differences *within* a protein's row and is **exactly invariant to adding a
constant to that row**. Nothing in the objective ever said where a row should sit. That rules out
the null anchor on principle rather than on arithmetic, and it rules out the global cut and every
other post-hoc calibration with it — none of them can recover a quantity training never
constrained. The fix is a second term in the loss that is not shift-invariant, and it is running:
[`DECISIONS.md` §13](DECISIONS.md), [`reports/calibration.md`](../reports/calibration.md).

**One number from that work belongs here already**, because it sets the bar the rest of this
document is silent about. Scored as a *decision* rather than a ranking, the nearest-neighbour
baseline detects a variant that binds nothing at AUROC **0.492** on `P3/all` — chance, and
correctly so, since under that regime it copies the variant's own wild type and is therefore the
literal hypothesis *the mutation does nothing*. **No model in this document has been shown to
beat it**, because until now nothing here could be asked the question.

## 9. The decision rule, measured — and the capability it was meant to unlock is not there

The λ sweep `T38` called for finished 2026-08-30: **312 runs (280 distinct configurations), 42.0
hours, no failures**, over `λ ∈ {0, 1, 3, 10, 30, 100, 300, 1000, 3000, 10000}`, both arms, 3 seeds
on the diagnostic folds. [`reports/calibration.md`](../reports/calibration.md) is the table; this
is the reading. **The noise floor is 0.005 AUPR** — the mean within-cell spread across the three
seeds, max 0.012 — and nothing below it is read here.

`T38` set three conditions. **Two are met and the third, which is the one the work was for, is
not.**

### 9.1 The ranking survives — `λ ≤ 100` is free

Paired against each run's own `λ = 0` control on the same arm, fold and seed, so the
fold-to-fold variance cancels:

| λ | pairs | mean Δ AUPR | on `S2` alone |
|---:|---:|---:|---:|
| 1 – 10 | 12 each | −0.000 to +0.002 | −0.001 to +0.007 |
| **30** | 52 | **−0.0005** | **+0.003** |
| 100 | 52 | −0.0019 | −0.003 |
| 300 | 52 | −0.0065 | −0.011 |
| 1000 – 10000 | 12 each | −0.005 to −0.010 | −0.009 to −0.011 |

So calibration is genuinely free up to `λ = 100` and costs a real amount from `λ = 300` up. On
`S2` — the regime the project hangs on — `λ = 30` is the only value that is not negative.

### 9.2 The calibration arrives, and transfers exactly as far as the ranking does

`T38`'s original complaint was a global cut calling 56 to 408 8-mers for proteins whose true
counts ran 10 to 206. That is fixed: under `expected_count_rule` — keep the top `round(Σ p)`, no
free parameter — the median call count now tracks the truth (`S1` 49 against 49, `P2` 34 against
35, `P3` 70 against 55).

And it is per-protein informative, not a global average dressed up: Spearman between the
predicted count and the true count is **0.48–0.52** on `P3`, `S1` and `P2`. **But 0.11–0.18 on
`S2` and −0.09 to −0.05 on `P1`.** The model knows how much a protein binds when it has seen its
relatives and does not otherwise — the same generalisation boundary §5 found for the ranking, and
calibration does not cross it either.

### 9.3 It cannot see a lost interaction. It never could, and the loss was not the reason

The capability the work was for. On `P3/all`, where all 18 dead variants are held out:

| | model, best λ | nearest-neighbour lookup |
|---|---:|---:|
| dead-variant detection AUROC | **0.43 – 0.52** | **0.492** |
| median 8-mers called for a protein that binds nothing | **89** | — |
| median 8-mers called for a protein that does bind (true 79) | 75 | — |
| dead variants for which it called **zero** 8-mers | **0 of 18**, at every λ | — |

**The model calls more 8-mers for proteins that bind nothing than for proteins that bind.**
Across the full grid the detection AUROC runs 0.33–0.56 and is at or below chance nearly
everywhere; the one region above it (`S2/fold-3`, 0.76–0.82) rests on **3** variants.

The mechanism is unambiguous. For those 18 variants, Spearman between the model's predicted
interaction power and **its wild type's true binding count** is **+0.77 to +0.91 — and it rises
with λ.** The model is predicting *this protein binds about as much as its wild type does*, which
is the null hypothesis, and better calibration only makes it a more confident null hypothesis.
`suppression`, the rank-based metric that predates this work, agrees and is unmoved: 0.46–0.54 at
every λ, which is random.

**So the failure is upstream of the objective.** The protein tower reads a pooled embedding of a
domain that differs from its wild type at one residue, and returns the wild type's answer. No
weighting of the two loss terms can fix that, and the sweep is the evidence that none does. This
is the same wall §4's C1 result hit from the other side — and consistent with it, the C1 set's
AUPR is unchanged at every λ (0.59–0.61), with call F1 0.17–0.23 against 0.66–0.79 on the rest.

### 9.4 The new comparison, and it is not flattering

Scoring calls rather than rankings finally makes the baseline and the model commensurable, since
`k = 1` always emitted a set. Pooled over all 19 folds:

| | model | lookup | delta |
|---|---:|---:|---:|
| AUPR (vs `k = 5`) | 0.628 | 0.593 | **+0.035** |
| **call F1 (vs `k = 1`)** | **0.451** | **0.457** | **−0.006** |

**The +0.035 AUPR margin does not survive being asked for a call.** Per regime the model wins as
a decision only on `S2` (+0.013 to +0.016) and on `P3` at `λ = 300` (+0.013); it loses on `P1`,
`P2` and `S1`. This is `D7`'s lesson a second time — there the baseline's apparent strength came
from ordering information the model never had, here the model's apparent strength comes from
emitting an ordering the task never asked for. **Nothing in §2's margin is invalidated; §2 is
simply not a statement about the deliverable.**

### 9.5 What to set, and what to stop claiming

**`model.bce_weight: 30`** if calibrated output is wanted: free on every regime, the only λ that
is positive on `S2`, and count-tracking as good as any. Nothing above 100 is worth its cost.

**Drop the pathogenic-variant claim.** On the evidence here the model cannot detect a lost
interaction, and it fails in the specific way that is worst for that use — confidently calling
~90 sites for a protein that binds none. That is not a threshold problem, which is what `T38`
took it to be; it is a representation problem, and it puts `T36` (a protein tower that can
resolve a point mutation) on the critical path rather than in the backlog.

## 10. What follows

1. **`T36` — a protein tower that can resolve a point mutation.** §9.3 relocated the problem:
   the model returns the wild type's answer for a single-residue variant, and no objective fixes
   that. This was a backlog item and is now the critical path for the project's central claim.
2. **`T38` — measured, and only half answered.** The decision rule exists and is honest (§9.1,
   §9.2); the capability it was built for is absent (§9.3). What is still open is what to claim:
   §9.4 says the model's margin is a ranking margin and not a decision margin.
3. **`T37` — seeds. Partly answered, and the news is bad.** The sweep's three seeds put the
   within-cell spread at **0.005 mean, 0.012 max**, on top of a cross-process non-determinism
   floor of the same size (`DECISIONS.md` §13.5). Several deltas in §3 and §6 are inside it.
4. **C1 is where the project's claim lives and it is currently negative.** §4's second half is the
   number to move, and §9.3 says why it has not moved. It is worth asking whether 18 variants can
   support the claim at all, and what evidence would.
5. **A longer budget**, given §7 — cheap, and the curve has not clearly stopped.

---

## Appendix — provenance of every number here

| number | source |
|---|---|
| model AUPR, delta, suppression, selection, seeds | [`reports/training.md`](../reports/training.md), from `scripts/run_grid.py`, commit `38b7aa1` |
| matched bar, chance, random-ranking band, C1 set | [`reports/nn_baseline.md`](../reports/nn_baseline.md), from `scripts/run_nn_baseline.py` |
| rank ceiling | [`reports/representation_ceiling.md`](../reports/representation_ceiling.md) |
| `k = 5` ranked binary bar | [`reports/nn_baseline.md`](../reports/nn_baseline.md), `run_nn_baseline.py --top-k` |
| embedding geometry in §6 | `scripts/check_pooling.py` and the audit measurements in [`DECISIONS.md` §12](DECISIONS.md) |

Per-domain rows for every fold are in `data/processed/training_domains.parquet` and
`data/processed/nn_baseline_domains.parquet`; the trained weights for all 38 runs are in
`data/processed/checkpoints/`, both the validation-selected and the budget model in each file.
