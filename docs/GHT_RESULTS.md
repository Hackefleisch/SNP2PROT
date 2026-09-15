# The GHT-SELEX arm — what it measured

**Written 2026-09-15**, from the run described in [`GHT_PLAN.md`](GHT_PLAN.md). The companion to
[`ML_RESULTS.md`](ML_RESULTS.md), which is the PBM arm's reading. Generated tables live in
[`reports/ght_results.md`](../reports/ght_results.md); this is the interpretation.

## The headline, in three sentences

**The DNA encoder reaches and passes a PWM fitted on the protein's own peaks — 0.929 auPRC
against 0.902, winning on 27 of 33 TFs — and every one of the 33 proteins stays above chance
(auROC 0.54-0.91) when held out entirely, which is rung 4 of the plan's ladder.** The
protein-blind control says how to read that: conditioning on the protein is worth **+0.233** auPRC
in Setting 1 and **+0.012** under a whole-neighbourhood holdout, so the architecture works and the
protein axis is used, but **what transfers to a protein with no relative left in training is
thin** — the same conclusion the PBM arm reached about the protein representation, from the other
direction. The strongest single result is rung 5: **the PBM arm's own protein tower, frozen,
matches the tower trained here and beats it by +0.041 auROC on `aliens`, the one negative set only
a protein can answer** — one protein space, two assays, two DNA vocabularies, no co-training.

**Two further results the plan did not expect.** The **per-residue protein tower**, which
`GHT_PLAN.md` §6 defaulted off on capacity grounds, is better on **10 of 10** held-out-protein
folds at **half** the parameters (§7), +0.029 auPRC on both `G1` and `G2` — and the capacity
reasoning was inverted: at 33 proteins the *pooled* tower is the expensive one.

**And one methodological finding that changes how any of these numbers may be quoted.** Deranging
the panel's protein vectors — same axis, same distribution, wrong protein — costs **nothing** in
Setting 1 (0.931 against 0.929) and **0.108 auPRC** under a held-out TF, where it falls below the
protein-blind control. A derangement is a bijection, so with every TF in training the model just
learns the permuted assignment: **Setting 1 cannot distinguish a protein representation from a
protein index.** Only Setting 2 can, and it does.

### What to claim, and what not to

| claim | evidence | verdict |
|---|---|---|
| we ingested their benchmark, their splits, their baselines | §1-2 | **yes**, rung 1-2 |
| parity with a per-TF PWM on held-out chromosomes | 0.929 vs 0.902, 27/33 TFs | **better than parity**, rung 3 |
| the DNA encoder is a learned PWM bank | median motif match 0.891 vs a 0.715 null, 33/33 | **yes**, §4.5 |
| generalises to a held-out TF where no per-TF method runs | 33/33 above chance; +0.040/+0.047 over motif transfer | **yes**, rung 4 |
| …and the margin is largest where identity transfer has nothing to copy | Spearman(identity, margin) −0.388 | **yes**, §5.1 |
| …but it is *not* mostly the protein | protein-blind reaches 0.564 of `G2`'s 0.576 | **state this out loud** |
| the protein is a *representation*, not an index | deranging it costs nothing under `C1` and 0.108 auPRC under `G1` | **yes, and only Setting 2 can show it** |
| one protein space serves both assays | frozen PBM tower matches and beats on `aliens` | **yes**, rung 5, §6 |
| the per-residue tower is worth turning on | better on 10 of 10 Setting 2 folds, at half the parameters | **yes**, §7 — one seed |
| anything about point mutations | none — this arm has no variants | **no** |

---

## What was built

| step | artifact | report |
|---|---|---|
| the benchmark | `data/raw/mex_archipelago/`, Zenodo `10.5281/zenodo.15555251` | `PROVENANCE.md` |
| the panel | `data/interim/ght/panel.parquet` — 33 TFs | [`ght_panel.md`](../reports/ght_panel.md) |
| the windows | `data/interim/ght/windows/` — 9,211,249 windows | [`ght_windows.md`](../reports/ght_windows.md) |
| the budgets | measured before the grid, not guessed | [`ght_preflight_C1-all.md`](../reports/ght_preflight_C1-all.md) |
| the baselines | PWM in Setting 1, motif transfer in Setting 2 | [`ght_results.md`](../reports/ght_results.md) |
| the grid | `data/processed/ght/ght_folds.parquet` | [`ght_results.md`](../reports/ght_results.md) |

---

## 1. The panel is 33, not 47, and the reason is worth a slide

`GHT_PLAN.md` §4 counted the usable panel by subtracting the 92 TFs labelled `C2H2 ZF` from the
benchmark's 139, giving 47, and §11 then warned in the same breath: *do not string-match family
labels; that produced a wrong panel count once already.* Running every construct through
`snp2prot.domains` — the admission path every other source in this project passed — gives **39
admitted, 33 after `D10`**.

The gap is not zinc fingers. It is:

| reason | count |
|---|---:|
| `repeat_array` — more than one Pfam hit of one family | 85 |
| `no_domain` — no hit in any family on `dbd_families.txt` | 13 |
| `c2h2_array_d10` — one `zf-C2H2` hit where the curators count 2-25 fingers | 6 |
| `mixed_families` — two DNA-binding families in one construct | 2 |

The 13 `no_domain` rejections are the interesting ones and 11 of them are **a coverage gap in our
own family list, not an absent domain**: they have Pfam hits, and none of those hits is on
`data/external/pfam/dbd_families.txt`, which was generated from what UniPROBE and CIS-BP curate.
Those two corpora are mouse- and yeast-heavy; the human Codebook panel reaches `FLYWCH`,
`Myb_DNA-bind_4/5`, `CGGBP1_N`, `ZNF704_C` and `DUF4772`, which they do not. The other two,
`FAM200B` and `ZNF518B`, have no Pfam hit anywhere above the gathering threshold — a genuine
absence rather than a gap.

One of the 13 is a latent defect rather than a gap. `MKX`'s only hit is `Homeobox_KN`, which
`snp2prot.domains.FAMILY_ALIASES` already maps to `Homeodomain` — but `call_domain` filters on the
family list **before** applying the alias, so an alias can only ever fire for a family that is
already admitted. Flagged rather than fixed: the list is shared with the PBM corpus and editing it
re-decides 1,338 domains.

**33 proteins is thin, and the plan already called 47 thin.** Every Setting 2 number below is
reported per fold as well as pooled.

---

## 2. The budgets were measured first, and the answer is the opposite of the PBM arm's

`GHT_PLAN.md` §8.1 asks for one fold run well past where it looks done, with the test metric
logged at every evaluation and used to set the budget. Run on `C1/all` for 20,000 steps:

| step | train loss | test auPRC |
|---:|---:|---:|
| 1,000 | 0.210 | 0.911 |
| 2,000 | 0.192 | 0.923 |
| 3,000 | 0.150 | 0.928 |
| **4,250** | — | **0.930 (peak)** |
| 8,000 | 0.123 | 0.926 |
| 14,250 | 0.086 | 0.926 |

**The PBM arm's lesson holds, with the sign reversed.** There, the training loss flattened at
~6,000 while held-out AUPR kept climbing for another 9,000 steps. Here the training loss keeps
falling from 0.19 to 0.09 long after held-out auPRC has peaked and begun to drift down. Two
opposite shapes, one conclusion: **neither curve predicts the other, and only the held-out one is
evidence.** The budget is 5,000 steps.

Wall clock, on one RTX 3070: ~41 ms/step, 1.7 s per evaluation pass, ~5 minutes per fold, 4.2 GB
peak. The whole grid — 11 folds × 3 seeds, plus the controls — is a few hours, which is what makes
the controls in §4 affordable at all.

---

## 3. What the negative sets actually ask

Three negative sets, three different questions, and they are **not comparable on auPRC**: the
positive fraction is 0.339 on `shades` and 0.091 on the 10:1-subsampled `random` and `aliens`, so
those are three different scales. auROC is the one number invariant to the imbalance, and it is
what every cross-set comparison below is read on.

| set | what it is | what it tests |
|---|---|---|
| `shades` | a flank 450-750 bp from the same summit, ~2:1 | local discrimination — can it find the site rather than the neighbourhood |
| `random` | GC-matched random genomic regions | **composition**: `shades` is not GC-matched and `random` is, so a model living off base composition scores well on the first and badly on the second |
| `aliens` | *other TFs' peaks*, GC-matched | **the protein**: the DNA is demonstrably bindable, so nothing but the protein can separate it from a positive |

`aliens` is the sharpest of the three, and [`ght_cobinding.md`](../reports/ght_cobinding.md)
says by how much. Measured over the test chromosomes:

- **49% of a TF's alien windows are a peak of another TF in our own 33-TF panel** — windows the
  model was shown in training as positives, for a different protein. A model that answers *does
  this window look bound* must get those wrong.
- **38% of a TF's own positives are also bound by another panel TF** (median; 8% to 93%). That is
  the size of the protein-blind shortcut on the positive side: co-bound promoters and enhancers
  where the answer is the same for several proteins at once.

The second number is why a Setting 1 auPRC means little on its own, and why the control in §4 is
not a strawman but a genuinely strong model.

---

## 4. Setting 1 — the DNA encoder against the PWM, and against itself blindfolded

Every TF is in training; eleven chromosomes are held out. The competitor here is a PWM — the
single best of that TF's own top-20 motifs, selected on the training chromosomes — and
ArChIPelago, an ensemble of such motifs feeding a random forest, cited under `D11`.

| | auPRC (`shades`) | auROC (`shades`) | auROC (`random`) | auROC (`aliens`) |
|---|---:|---:|---:|---:|
| two-tower | **0.929 ± 0.000** | 0.963 | 0.968 | 0.947 |
| PWM, best of that TF's own | 0.902 | 0.938 | 0.942 | 0.935 |
| protein-blind control | 0.697 | 0.818 | 0.806 | 0.661 |
| chance | 0.339 | 0.500 | 0.500 | 0.500 |

**Parity was the goal and the margin is better than parity**: +0.027 auPRC over the fitted PWM,
and the model wins on **27 of 33** TFs individually, median +0.021. That is a distributed win
rather than a lopsided one, which is what the per-TF table exists to check.

**No composition shortcut.** `GHT_PLAN.md` §12 named it as the realistic failure mode, since
`shades` is not GC-matched and `random` is. The model scores 0.963 auROC on `shades` and 0.968 on
GC-matched `random` — indistinguishable. It is not living off base composition.

**And the protein is doing heavy lifting here.** The blind control reaches 0.818 auROC on
`shades`, which is the size of the "does this window look bound" signal, and falls to **0.661 on
`aliens`** against the real model's 0.947. Asked to say *whose* peak it is looking at, a model
with no protein input cannot — which is exactly right, and is the clearest evidence in the whole
arm that the protein axis is being used.

### 4.5 The learned filters are the known motifs

A `Conv1d` weight over a one-hot encoding **is** a position weight matrix — `weight[c, b, j]` is
what base `b` at offset `j` contributes to channel `c`, which is a log-odds column up to an
additive constant that no ranking can see. So the §5 claim that this encoder is the learned
generalisation of the PWM baseline is checkable rather than rhetorical.

Measured on the `C1/all` checkpoint ([`ght_filters.md`](../reports/ght_filters.md)): every one of
the 33 panel TFs' selected motifs has a match in the filter bank, **median similarity 0.891**,
range 0.641–0.968. Against the null the same best-of-256 search returns on the *same filters with
their columns permuted* — same weights, same distribution, no order — median **0.715**, and
**33 of 33 TFs beat their own permuted best.**

The widths are all used, which is the §5 decision measured rather than argued: 8 TFs' best match
sits in the 8 bp bank, 13 in the 12 bp bank, 9 in 16 bp and 3 in 20 bp. A single-width bank would
have thrown away the last two groups.

---

## 5. Setting 2 — held-out TF, where no per-TF method can run

Both axes held out at once: a held-out protein, on held-out chromosomes. A PWM and a random
forest over PWM hits are fitted per TF and have nothing to compute here. What can be asked is the
field's own move — transfer the motif of the most identical characterised DNA-binding domain,
which is how CIS-BP infers a motif for an uncharacterised TF.

| regime | two-tower | 1NN motif | 5NN motif | chance |
|---|---:|---:|---:|---:|
| `G1` random TFs | **0.650 ± 0.077** | 0.610 | 0.599 | 0.339 |
| `G2` components at ≥ 0.5 identity | **0.576 ± 0.037** | 0.529 | 0.527 | 0.339 |

auPRC on `shades`, three seeds, macro-averaged over held-out proteins. The margins are **+0.040**
and **+0.047** against `1NN`, **+0.051** and **+0.050** against `5NN`.

**Every one of the 33 proteins is above chance when held out** — auROC 0.54 to 0.91 — which is
the fallback ladder's rung 4 reached. But the mean is not the interesting part.

### 5.1 The margin lives where identity transfer has nothing to copy

Joined per **(held-out protein, fold)**, because `G1` and `G2` leave a protein with different
neighbours and that difference is the experiment:

| nearest identity | pairs | from | two-tower | 1NN motif | delta |
|---|---:|---|---:|---:|---:|
| 0.00 – 0.35 | 49 | `G1`, `G2` | 0.571 | 0.508 | **+0.062** |
| 0.35 – 0.50 | 6 | `G1`, `G2` | 0.667 | 0.631 | **+0.036** |
| 0.50 – 0.70 | 9 | `G1` | 0.759 | 0.810 | **−0.051** |
| 0.70 – 1.00 | 2 | `G1` | 0.864 | 0.796 | **+0.068** |

Spearman against nearest-neighbour identity: **1NN +0.680**, **two-tower +0.559**, **delta
−0.388**, over 66 pairs. Read those three together:

1. Motif transfer works, and works better the closer the neighbour is. That is why it is the
   right bar, not a strawman.
2. **The model has not escaped the identity dependence** — its own score tracks identity too.
3. What changed is the *margin*, which grows as the neighbour gets more distant.

**And it goes negative where the neighbour is close.** 43 of 66 pairs are wins; the losses cluster
in the 0.50–0.70 band, which exists only under `G1` — `G2` removes whole identity components, so a
protein held out there cannot have a close relative left in training. The sharpest single case is
`SRY` under `G1/fold-4`, where `SOX15` at 60% identity stays in training: copying `SOX15`'s motif
scores 0.72 and the model scores 0.38. **When a near-twin is available, copying it beats learning
a general map.** That is the honest shape of the result, and it is a better slide than a uniform
win would have been.

The same statement in its sharpest form: the 10 (protein, fold) pairs where the held-out TF is the
**only member of its family in the panel** score 0.477 against motif transfer's 0.387 (+0.090);
the 56 with a relative score 0.639 against 0.602 (+0.037).

### 5.2 One protein the model cannot infer, and it is instructive

`LEF1` holds out at auROC 0.54 — barely above chance — yet in Setting 1, where `LEF1` is in
training, the filter it learns matches `LEF1`'s published motif at **r = 0.97**, the best match in
the whole bank (§4.5). So the DNA tower can represent that motif perfectly well. What it cannot do
is *derive* it from the `LEF1` domain it has never seen. That is the protein tower's limit stated
as cleanly as this data can state it.

### 5.3 The control that decides what to claim

`protein = constant` is the identical model with every TF handed the same vector — structurally
protein-blind. It is not a strawman: a median 38% of a TF's peaks are co-bound by another panel
TF (§3), so "does this window look bound" has a good answer.

Beside it runs a second control. `protein = shuffled` **deranges** the panel's vectors: the
protein axis is present, carries the same distribution, and is wrong. It separates *uses the
protein as a representation* from *uses it as an index* — which is the project's whole premise.

| holdout | `shades` auPRC | | | `aliens` auROC | | |
|---|---:|---:|---:|---:|---:|---:|
| | model | blind | **wrong** | model | blind | **wrong** |
| `C1` chromosomes | 0.929 | 0.697 | **0.931** | 0.947 | 0.661 | **0.947** |
| `G1` random TFs | 0.650 | 0.580 | **0.542** | 0.685 | 0.566 | **0.572** |
| `G2` whole neighbourhoods | 0.576 | 0.564 | **0.551** | 0.606 | 0.556 | **0.569** |

### 5.4 Setting 1 cannot tell a representation from an index

**Under `C1` a deranged protein axis scores 0.931 against the real model's 0.929**, and 0.947
against 0.947 on `aliens`. That is the control working, not failing: a derangement is a
*bijection*, every TF is in training, so the model simply learns the permuted assignment.

**Setting 1 is therefore blind to the difference between a protein representation and a protein
index.** No number from it can support the claim this project exists to make. That is the whole
reason Setting 2 exists, and it is worth saying on a slide.

**Under a held-out TF the ordering inverts and the claim becomes measurable.** `shuffled` falls
**below** `constant`, which falls below the real model: 0.542 < 0.580 < 0.650 on `G1`, 0.551 <
0.564 < 0.576 on `G2`. A wrong embedding is *worse than no embedding*, because the model
faithfully applies a map it learned for a different protein. Three orderings, one conclusion: the
embedding carries protein-specific information, and it is not a free per-TF parameter.

**Three readings, and the third is the one to take to the talk.**

1. **In Setting 1 the protein is doing heavy lifting**, and most visibly where it must: on
   `aliens` the blind model falls to 0.661 against 0.947.
2. **On `aliens` the protein earns its keep at every level of holdout**, including `G2`, and on
   every one of its five folds.
3. **On `shades`, transfer to a protein with no family neighbour left in training is weak.** `G2`
   gains **+0.012** auPRC and two of its five folds are negative. The +0.047 margin over `1NN` is
   real — and so is the fact that a model given no protein at all reaches 0.564 where the real
   model reaches 0.576.

So the architecture works and the protein axis is used. **What does not yet transfer is the
protein representation itself, once every relative has been removed.** That is the same conclusion
the PBM arm reached from the other direction — there the model returns the wild type's answer for
a single-residue variant ([`ML_RESULTS.md`](ML_RESULTS.md) §9.3, `T36`). Two arms, two assays, one
open problem, and the agreement between them is itself worth reporting.

**One asymmetry when reading the gaps**: the model ran at three seeds per fold and the control at
one, so a per-fold gap carries the control's single-run noise. `G2`'s within-cell spread on the
model is 0.003-0.012 auPRC, the same size as two of the five `G2` gaps.

---

## 6. Rung 5 — the PBM arm's own protein tower, frozen

The two arms' protein towers are **literally the same module**: `ProteinTower(1280, 256)`, four
tensors, same names, same shapes. That is not a coincidence — [`ML_PLAN.md`](ML_PLAN.md) §4.2
fixed the shared width `D = 256` across arms so that no comparison between them could be a
comparison of widths, and the same decision is what makes the weights portable.

So it can be lifted across. Below, the protein tower is loaded from
`A1_P3_all_seed20260819_lam30.pt` — a PBM run over 1,338 domains against 32,896 8-mers — and
**frozen**. The protein representation is fixed by a different assay on a different DNA
vocabulary; only the genomic DNA tower learns.

| holdout | metric | trained here | **frozen PBM tower** | protein-blind |
|---|---|---:|---:|---:|
| `C1` | `shades` auPRC | 0.929 | **0.931** | 0.697 |
| `C1` | `aliens` auROC | 0.947 | **0.949** | 0.661 |
| `G1` | `shades` auPRC | 0.650 | **0.649** | 0.580 |
| `G1` | `aliens` auROC | 0.685 | **0.712** | 0.566 |
| `G2` | `shades` auPRC | 0.576 | **0.563** | 0.564 |
| `G2` | `aliens` auROC | 0.606 | **0.646** | 0.556 |

**It matches, and on the hardest test it wins.** On `aliens` — the set only a protein-conditioned
model can answer — the frozen PBM tower beats the tower trained on this data by **+0.027** under a
random TF holdout and **+0.041** under a whole-neighbourhood holdout.

**The caveat, measured rather than waved away.** 8 of the 33 panel domains appear verbatim in the
PBM corpus (`ELF3`, `FLI1`, `FOSL2`, `GABPA`, `GCM1`, `LEF1`, `MYF6`, `NR1H4`), so for those the
PBM tower was fitted on that protein's own 8-mer behaviour. That is not a leak of GHT labels — the
PBM arm has never seen a genomic window — but it is not an unseen protein either. Split on it,
`aliens` auROC under a held-out TF:

| | trained here | frozen PBM tower |
|---|---:|---:|
| 25 domains the PBM tower has **never seen** | 0.671 | **0.699** |
| 8 it has | 0.570 | **0.622** |

The transfer holds on both halves, so the caveat does not explain it.

**Why frozen is better, and it is the capacity argument landing.** The frozen tower has **zero**
trainable protein parameters against 330,496 for the one trained here, and a Setting 2 fold has
about 22 training proteins. `GHT_PLAN.md` §6.2 predicted a capacity problem and located it in the
per-residue conv tower at ~3.4k parameters per protein; §2 above measured that even the *linear*
tower is ~10k per protein at this panel size. Freezing it removes the problem entirely — and the
representation that replaces it was fitted on 1,338 proteins rather than 22.

**This is rung 5 of the fallback ladder, and it is the strongest single claim the arm makes.** One
protein space, two assays, two DNA vocabularies, no co-training.

---

## 7. The per-residue protein tower, which the plan defaulted off — and should not have

[`GHT_PLAN.md`](GHT_PLAN.md) §6 built a convolution over ESM-2's per-residue vectors behind a
config flag and defaulted it **off on capacity grounds**: ~158k parameters against 47 proteins is
~3.4k each, where [`TRAINING.md`](TRAINING.md) §5 calls the PBM arm's 245 the central constraint.

Two things that reasoning did not anticipate. The panel is **33**, so the *pooled* tower is
already ~10k parameters per training protein (§2). And the conv tower is **smaller** — 160,672
against 330,496 — because the pointwise reduction to 64 dimensions happens before any convolution.

| holdout | metric | pooled vector | **conv over residues** |
|---|---|---:|---:|
| `C1` | `shades` auPRC | 0.929 | **0.929** |
| `C1` | `aliens` auROC | 0.947 | **0.946** |
| `G1` | `shades` auPRC | 0.650 | **0.679** |
| `G1` | `aliens` auROC | 0.685 | **0.714** |
| `G2` | `shades` auPRC | 0.576 | **0.605** |
| `G2` | `aliens` auROC | 0.606 | **0.622** |

**Better on 10 of 10 Setting 2 folds**, by +0.001 to +0.067 auPRC, never worse; identical under
`C1`. Ten out of ten with no fold going the other way is strong at two seeds against the pooled
tower's three — the seed counts are not equal, so the *size* of the margin is provisional even
though its sign is not.

**And it helps exactly where a pooled vector should hurt.** Mean pooling over a domain's residues
is a summary; a bank of kernels at widths 3, 7 and 15 asks whether particular local patterns are
present — the homeodomain's `WFQNRR`, the bZIP basic region before its leucine heptad, the Cys/His
spacing of `zf-C4`. Under `C1`, where every protein is in training, that distinction buys nothing.
Under a held-out protein it buys something on every fold.

`GHT_PLAN.md` §6.3 warned against letting `reports/unpooling_check.md` decide whether this tower
got built, and was right to: the pre-flight is contested, n = 86, 56 of them homeodomains. **The
tower is worth building on its own measurement, and this is that measurement.** The right default
is now `ght.model.protein.tower: residue`, subject to more seeds.

**It is also a live lead for the PBM arm.** `T36` asks for a protein tower that can resolve a
point mutation; the pooled representation is what destroys that signal (`ML_RESULTS.md` §9.3).
This is the first evidence in the project that an un-pooled tower helps a real task rather than
merely preserving a correlation.

---

## 8. What this does not show

**It says nothing about point mutations, which is the project's central claim.** There are no
variants in this arm at all. On the PBM arm, where there are 173, the model returns the wild
type's answer for a single-residue variant and cannot detect a lost interaction
([`ML_RESULTS.md`](ML_RESULTS.md) §9.3). That result stands; nothing here touches it, and nothing
here should be used to suggest otherwise.

**33 proteins is thin.** A Setting 2 fold holds out six or seven. The per-fold spread is the
number to argue with, not the mean.

**The construct and the protein disagree.** Many peak files come from a full-length construct
while the model is fed the isolated padded domain. `panel.parquet` carries `insert_composition`
so a reader can see which; nothing in the pipeline can fix it.

**The PWM baseline may be advantaged.** Of the 33 selected motifs, 17 are `AFS.*` — derived from
GHT-SELEX itself — and the release does not say whether motif discovery was restricted to the
training chromosomes. If it was not, the Setting 1 comparison is conservative rather than
flattering. Re-deriving the motifs would abandon `D11` and with it the comparability that makes
the whole benchmark worth using.

**A chromosome split is not a repeat split.** A repeat family present on both sides can be
memorised. That is the benchmark's own design (`D5`), taken unchanged so the numbers remain
comparable to theirs, and it is a property of every published number on this benchmark rather
than of ours alone.

**`random` and `aliens` auPRC is not comparable to the published numbers.** Those sets are
subsampled 10:1 here against the released 100:1, which moves chance auPRC from 0.0099 to 0.091.
auROC is unaffected and is what every cross-set comparison uses.

---

## 9. What the plan asked for and what was not built

| `GHT_PLAN.md` | state |
|---|---|
| §5 multi-width conv, global max over both strands | built, and §4.5 checks it recovers real motifs |
| §6 pooled protein tower as default | built; it is the headline |
| §6.2 per-residue conv tower behind a flag | built, run, and **it should be the default** — §7 |
| §6.3 BLOSUM62 severity feature | **not built, and not applicable**: it scores a substitution, and this arm has no variants |
| §6.4 `T36`'s non-linear tower (`protein.hidden > 0`) | **not run.** One config value; §7 makes it less urgent than the per-residue tower |
| §7 BCE primary, `mu` contrastive term optional | built; `mu = 0` is exactly plain BCE and is the headline |
| §7 null anchor unused | **absent**, not untrained — see `DECISIONS.md` §14 |
| §8.1 step and wall-clock pre-flight | run first, before anything was trained |
| §8.2 log the split, both weights, the code stamp | built, via the PBM arm's own `tracking` and `Trainer.save` |
| §9 PWM baseline, 1NN/5NN | built; the NN form is motif transfer, see `DECISIONS.md` §14.5 |
| §9 ArChIPelago | cited under `D11`, not rerun |
| §11 the 10 clean extras outside the benchmark | **not taken** — `TODO.md` `T40` |
| §12 `PROVENANCE.md` rows | written; three URLs are missing and are `TODO.md` `T39` |

Two things were built that the plan does not ask for, because the measurement said they were
needed: the **protein-blind and protein-shuffled controls** (§3 says why a Setting 1 number is
unreadable without them) and the **`aliens` negative set** (`DECISIONS.md` §14.2).

---

## Appendix — provenance of every number here

| number | source |
|---|---|
| panel counts, rejections | [`reports/ght_panel.md`](../reports/ght_panel.md) ← `scripts/build_ght_panel.py` |
| window counts, chromosome split | [`reports/ght_windows.md`](../reports/ght_windows.md) ← `scripts/build_ght_windows.py` |
| co-binding, `aliens` composition | [`reports/ght_cobinding.md`](../reports/ght_cobinding.md) ← `scripts/check_ght_cobinding.py` |
| step and wall-clock budgets | [`reports/ght_preflight_C1-all.md`](../reports/ght_preflight_C1-all.md) ← `scripts/preflight_ght.py` |
| baselines | `data/processed/ght/ght_baselines.parquet` ← `scripts/run_ght_baselines.py` |
| every model number | `data/processed/ght/ght_folds.parquet`, `ght_tfs.parquet` ← `scripts/run_ght_grid.py` |
| the tables below them | [`reports/ght_results.md`](../reports/ght_results.md) ← `scripts/make_ght_report.py` |
| the figures | `results/figures/ght_*.png` ← `scripts/make_ght_figures.py` |
| the talk page | `results/ght_story.html` ← `scripts/make_ght_story.py` |

Raw inputs: `PROVENANCE.md`, rows `mex_archipelago`, `genomes`, `codebook_ghtselex`,
`codebook_htselex`. Three of those carry **no source URL** — the download logs were empty and
rule 1 forbids guessing one; `TODO.md` `T39` is the ask.
