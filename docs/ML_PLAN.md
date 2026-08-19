# ML plan — contrastive TF–DNA binding, for the talk

Status: **plan of record, not yet executed.** Written 2026-08-19 from the owner's outline.
§1–§9 are the plan. §10 is Claude's review of it — kept separate so the plan and the criticism
never get confused for each other.

Review items move **out** of §10 and into the plan when the owner accepts them, leaving a dated
pointer behind so the reasoning stays findable. Accepted so far, both 2026-08-19:

| review item | landed in | what changed |
|---|---|---|
| §10.1 | **§8.1** | the nearest-neighbour lookup became the bar, and step 0 of the build order |
| §10.2 | **§4.1** | the DNA encoder switched from a language model to one-hot + revcomp pooling |
| §10.3 | **§2.3** | the SELEX k-mer format settled: 8-mers only, on the shared canonical space |
| §10.4 | **§2.3, §7** | SELEX data is exhaustive, so its negatives are `assayed_unbound`; reliability gets measured, not assumed |
| §10.5 | **§7** | the 45.8% cross-lab agreement goes on the transfer figure as a ceiling |
| §10.6 | **§6.1** | three protein-axis regimes specified: P1 homeodomain, P2 variant clusters, P3 fractions |
| §10.7 | **§5, §7** | DNA-side holdout dropped entirely; transfer sets collapse to easy/hard on the protein |
| §10.8 | **§1, §4.2, §4.3** | C2 reworded as inductive bias; ensemble carries spread not mean; all arms one width |
| §10.9 | *discarded* | shared-lab context; two consequences kept in `T31` |
| §10.10 | **§3.1** | dense in-memory training shape. Its per-residue recommendation was later overruled in favour of pooled |
| §10.11 | **§3, §3.1** | two-tower + multi-positive InfoNCE; batch shape follows; `null` anchor for all-negative rows |
| §10.12 | **§5.2** | UMAP: 3 highlighted-protein panels + a promiscuity panel + the NN-identity histogram |
| §10.13 | **§4.2, §8.2** | TransBind is a re-implemented benchmark baseline; ESM-DBP joins as fourth arm A4 |
| — | **§3.1** | protein embeddings pooled, not per-residue (owner, 2026-08-19), for cross-arm comparability |
| §10.14 | **§5.3** | positives-first metrics everywhere, computed **per protein**, macro-averaged |
| §10.15 | **§9** | SELEX frames and closes the talk; PBM carries the middle; §5 becomes desirable, not required |
| §10.16 | **§9.2**, `T28`, `T33` | second config file, MLflow behind a wrapper, seeds in config, deps in `dev` |

**All sixteen are resolved** — fifteen accepted in whole or part, §10.9 discarded, nothing left
open. Each §10 entry is now a dated pointer to where its content lives; the reasoning stays there
because several of the decisions only make sense with it.

Related: [TODO.md](../TODO.md) for open items this raises, [DECISIONS.md](DECISIONS.md) for what
is already settled, [RESULTS.md](RESULTS.md) for what the PBM dataset actually contains.

---

## 1. The goal

**Preliminary results, fast, for a talk.** Not a paper, not a finished method. The work is
scoped to the big cornerstones and to producing a story that holds together in front of an
audience. Where a choice is between "more rigorous" and "gets a defensible figure this month",
this project takes the second — but records what it gave up.

Two claims to establish:

- **C1 — Contrastive learning is a promising avenue for predicting TF–DNA binding, and it can
  infer the effect of point mutations that lie outside its own training corpus.**
- **C2 — A structure-derived representation is a better *inductive bias* for this task than
  sequence alone.** Phrased as bias, not information, and deliberately: the structures are
  *predicted from sequence*, so a structural embedding contains nothing the sequence did not
  already contain. Any gain is that the structure makes the relevant geometry easier to learn.
  A biologist audience accepts that framing readily and would rightly push back on the stronger
  one.

Everything below exists to support one of those two, or to make them presentable.

---

## 2. The second dataset — Codebook SELEX

### 2.1 Source

The SELEX data from **Jolma, Laverty et al. 2026, *Nature*, "An expanded codebook of human
transcription factor DNA-binding specificity"**, DOI
[`10.1038/s41586-026-10798-9`](https://doi.org/10.1038/s41586-026-10798-9) — verified against
Crossref 2026-08-19, 55 authors, Hughes and Kulakovskiy senior. Published 2026-08-05.

This is the peer-reviewed version of the bioRxiv preprint already in
[REFERENCES.md](REFERENCES.md) as `10.1101/2024.11.11.618478`. The PDF is **not yet
downloaded** and it blocks this phase — see §9.

Deposit: **<https://codebook.ccbr.utoronto.ca/v2/index_v2.php>**.

The paper's data is the required core. **More data from the same portal may be included**, but
the paper's SELEX must be in there. The portal was inspected on 2026-08-19 and lists ChIP-seq,
GHT-SELEX, HT-SELEX, PBM, PWMs and accessory ATAC/H3K9me3 data, distributed as `.tar.gz` and
`.xlsx`, with raw reads in SRA/GEO. Exact accessions to be read off the paper's data
availability statement, never guessed (rule 1).

### 2.2 What it is and is not

As far as we understand it, this holds **very similar data to the current PBM dataset but with
no proteoform data** — that is, it varies on the DNA axis and not on the protein axis. There are
no designed single-residue variants of one domain.

Collaborators have said they can turn the SELEX data into **exhaustive binding k-mers for each
variant**. If that lands, we can build a dataset closely parallel to the PBM one, at least on
the DNA axis.

### 2.3 What we ask the collaborators for — settled 2026-08-19

**8-mers. Only 8-mers.** That decides §10.3 and makes the transfer experiment (§7) run on a
single shared DNA encoder.

The exact spec to hand over:

| | |
|---|---|
| k | **8**, nothing else |
| space | the same **32,896** non-redundant 8-mers the PBM tables already use |
| canonicalisation | the **lexicographically smaller** of `(s, revcomp(s))` |
| score | a **continuous** value per 8-mer, not a binary call |
| coverage | **exhaustive — every one of the 32,896 8-mers scored**, positives *and* negatives |

**Why 8 and not longer.** Universal PBM is 8-mer-based and always will be, so an 8-mer SELEX
table makes the two assays *directly comparable* rather than merely joinable — the encoder,
the DNA space and the split definitions are all shared. A longer k would also reintroduce the
`dna_len` leak this project restructured itself to remove (`TODO.md` N5): length alone would
encode assay identity and with it the label prior.

**Why the canonicalisation matters.** Verified against the tables 2026-08-19: all **19 sources
carry one identical 8-mer set**, and every stored 8-mer is the lexicographically smaller member
of its reverse-complement pair (the 256 palindromes satisfy this trivially). If collaborators
canonicalise the other way, or not at all, the same physical k-mer arrives under a different
string and the join fails **silently** — the table validates, the row count looks plausible, and
half the k-mers simply never match. Cheap to state now, expensive to discover later.

**Why exhaustive.** Confirmed with the collaborators 2026-08-19: they can score **all** 8-mers,
not just the enriched ones. That is what makes the SELEX table structurally the same object as a
PBM table — one complete row per (domain, 8-mer) — and it is what lets §7 compare the two at all.

It also settles the `neg_provenance` question (rule 4). An exhaustively scored SELEX negative was
**measured and scored below threshold**, which is precisely what `assayed_unbound` means in
`schema.VALID_NEG_PROVENANCE` — not `selection_absent`, which is for a sequence that merely never
turned up. The two assays stay distinguishable through `assay` and `score_type`, so nothing is
pooled indistinguishably and the rule is satisfied.

**The residual caveat, kept deliberately.** A SELEX negative is still *derived from depletion in
a selection*, while a PBM negative is a direct spot measurement. With exhaustive enumeration that
gap narrows a lot, but it does not vanish: enrichment estimates get noisy for 8-mers that were
poorly covered in the input library. **SELEX may simply be less reliable on the negative side**,
and that is worth holding in mind rather than assuming away. §7 says how to check it empirically.

**Why a continuous score.** So we set our own cutoff from `configs/thresholds.yaml` rather than
inheriting theirs. Their binarization would be invisible in the stored table and would make the
SELEX labels non-self-describing, which the schema forbids (`threshold_pos` / `threshold_neg` on
every row).

**The cost, stated plainly.** SELEX can resolve motifs longer than 8 bp, and restricting to
8-mers discards that. This is a deliberate trade for comparability, not an oversight: PBM
truncates long motifs the same way, so both datasets are limited *identically*, which is what the
transfer experiment needs. Worth saying out loud if anyone asks why a SELEX dataset is being read
at PBM's resolution.

### 2.4 Consequence

SELEX gives us: a second, independent assay; a much larger human TF panel; and a held-out set
for the transfer experiment in §7. It does **not** give us protein-axis depth — that remains the
PBM dataset's unique contribution, and it is why the variant-holdout experiment (§6) exists on
the PBM side only.

---

## 3. The model — contrastive learning

A **simple contrastive learning approach**: two encoders, one for the protein and one for the
DNA site, trained so that a binding (domain, site) pair is close in a shared space and a
non-binding pair is far apart.

Deliberately simple. The point of the talk is that the *approach* is promising and that the
*embedding modality* matters, not that a particular architecture is the best one.

**Decided 2026-08-19: a two-tower model with multi-positive InfoNCE.** Protein tower, DNA tower,
scored by similarity in a **shared space**. The objective is not an implementation detail — it is
what makes §5.2's figures exist. A two-tower model puts proteins and DNA in *one* space, so the
UMAP is a direct read-out of the trained model; a cross-attention model with a BCE head
(TransBind's shape, §8.2) has no joint space to project and the figure could not be built from it.
That is also why TransBind stays a **benchmark baseline** and never becomes our architecture.

### 3.1 Training shape — the whole dataset fits in memory

**Accepted 2026-08-19.** The merged table is a **complete matrix**: 1,338 domains × 32,896 8-mers,
every cell present, and the row count is an exact invariant — 44,014,848 = 1,338 × 32,896
(`TODO.md` `N1`). It is *stored* as 44 M rows in a 974 MB parquet. It does not have to be
*trained* that way.

Measured budget, against the RTX 3070's 8 GB:

| tensor | shape | dtype | size |
|---|---|---|---:|
| labels | 1,338 × 32,896 | `int8` | **42.0 MiB** |
| DNA embeddings | 32,896 × 256 | `fp32` | 32.1 MiB |
| protein embeddings, pooled | 1,338 × 1280 | `fp32` | 6.5 MiB |
| logits, 64 domains × every 8-mer | 64 × 32,896 | `fp32` | 8.0 MiB |

Everything is resident at once, so there is **no dataloader, no row shuffling and no I/O in the
training loop** — which turns the experiment grid (5-fold × 2 regimes × 3 arms, plus §6.1's
PBM-only regimes) from days into an afternoon. That is the difference between having preliminary
results for the talk and not (§1).

#### Batch shape — settled with the objective, 2026-08-19

**Shuffling is on the protein axis:** each step draws a random batch of domains (64, say), never a
single one. **The 8-mer axis stays complete** — each domain in the batch is scored against all
32,896. With InfoNCE that is the exact form, not a shortcut: we hold the complete label matrix, so
approximating a negative set we actually know buys nothing, and the denominator being the full
8-mer space makes the **training objective identical to §5.3's per-protein ranking metric**.

Positives per domain are skewed — **median 49, mean 69, max 236, 5th percentile 1**
(`reports/label_health.md`) — so normalise the loss per domain, or a 236-positive domain will
outweigh a 1-positive one by construction.

**Fallback if it proves unstable:** *sampled softmax* — keep all of a domain's positives,
subsample the negatives (e.g. 4,096 of ~32,847) per step. Under full sweeps **every DNA embedding
takes a gradient at every step**, so there is no stochasticity at all on the DNA side; sampling
restores it at no cost worth counting. A knob, not a redesign.

#### All-negative rows are training signal, and need one mechanism to become it

A domain with no positives should push **every** 8-mer away from itself — real signal, and worth
keeping. But **plain multi-positive InfoNCE cannot express it**: the loss averages over the
positive set, so an empty positive set makes the average vacuous and the row contributes no
gradient at all. It is not that the row is useless; it is that the standard loss has nowhere to
put it.

**The mechanism: a learned `null` anchor on the DNA side**, present in every denominator. For a
domain with positives it is just one more negative. For a domain with none, **the null is the
target** — which drives every real 8-mer's similarity down and is exactly the intended behaviour,
under one uniform loss. (One-line variant: add `λ · logsumexp_k s(p, k)` for all-negative rows,
i.e. a fixed anchor at score 0.)

**On this dataset the mechanism applies to exactly 54 records, and 20 of them are spoken for.**
Every domain is scored against all 32,896 8-mers, so all-negative rows are precisely the 54 that
`reports/label_health.md` already flags — and they are not interchangeable:

| verdict | n | what it is | use |
|---|---:|---|---|
| `dead_variant` | 20 | variant that measurably lost binding, **with a control in its own series** | **evaluation only** — the C1 money rows (`T30`) |
| `no_evidence` | 34 | no positives **and no control to read that against** | unverifiable: may be a true non-binder, may be a failed experiment |

Training on the 20 would consume the sharpest C1 evidence we have. Training on the 34 risks
teaching the model that a perfectly good TF binds nothing, because without a control there is no
way to tell a real non-binder from an assay that failed. So the mechanism is worth building — it
will matter more on SELEX, where all-negative proteins may be commoner — but on PBM it should be
pointed at the 34 at most, and that is `T30`'s decision to make.

**Embeddings are a precomputed lookup, not part of the model.** 1,338 domains through the PLM and
32,896 8-mers through the DNA encoder are one-off inferences, cached to disk. The PLM never runs
inside the training loop.

**Protein embeddings are pooled to one vector per domain — decided 2026-08-19.** Mean-pooled over
the sequence, as TransBind does, giving a fixed 1280-d vector regardless of domain length.

The deciding reason is **cross-arm comparability**, and it outranks the alternative. A structural
embedder produces one vector per structure, so if the sequence arms were per-residue and the
structure arms were per-structure, A1/A4 and A2/A3 would feed *differently shaped towers* — and
the A1 → A2 comparison that carries claim **C2** would confound modality with architecture. One
vector per domain everywhere keeps the four arms genuinely interchangeable. It also makes the
TransBind baseline (§8.2) exactly controlled rather than nearly so, and it is the more universal
choice for anything added later.

**The cost is real and lands on C1, so measure it before trusting it.** Pooling dilutes a
single-residue change: a variant differs from its reference by 1-5 residues out of a median 77.
The worst case is that the domain vector moves by about one part in 77 — but that bound assumes
only the mutated residue's representation changes, and **ESM-2 is contextual**, so a substitution
shifts its neighbours' representations too. The true attenuation is an empirical question, not
1/77.

So it gets a **pre-flight check, before any training**: embed all 173 variants and their cluster
references, and compare variant-to-reference distance against the distance between unrelated
domains of the same family. If variants are not measurably separated from their own wild type in
the pooled embedding, **no §6.1 regime can produce a meaningful C1 number on that arm**, and it is
far better to know that in an afternoon than after the grid has run. The same check §4.1 specifies
on the DNA side, pointed at the axis that actually carries the claim.

If pooling does prove too lossy, the fallback that keeps every one of the reasons above is
**attention pooling** — still one fixed-width vector per domain, but learned rather than uniform.
Recorded in `TODO.md` `T32`.

**Two masks, not filters.** The gray band (`label == -1`, 1.7% of rows) is *kept in the matrix and
excluded from the loss* — it is absent evidence, not a negative. Records flagged by
`label_health.usable()` are masked the same way. Both are boolean masks over a resident matrix, so
neither costs anything.

---

## 4. The embeddings

One DNA representation, four protein representations — **1 × 4 = 4 embedding combinations**, all
four run through every experiment. The protein arms are named **A1 … A4** so they do not collide
with §6.1's `P1` / `P2` / `P3` split regimes.

**A1 and A4 are buildable today**; A2 and A3 wait on the structure ensembles (§4.3, §9). That is
the whole reason the build order starts with sequence embeddings — and with A4 added, the
unblocked phase carries a real comparison of its own rather than a single arm waiting for company.

### 4.1 DNA — one method

**Decided 2026-08-19: plain one-hot with reverse-complement mean pooling. No DNA language
model.** One encoder, one arm — this axis is not what the talk is arguing about. The argument
that got here is §10.2; the short form:

**8 × 4 one-hot into a small CNN.** One-hot is an *input parameterisation*, not an embedding.
The continuity lives in the encoder, which is Lipschitz, so 8-mers a single base apart land near
each other — and *how* near is **learned from the binding data** rather than imposed. That is the
right way round, because the similarity that matters here is not generic base chemistry but
*which substitutions preserve binding*, and that is position- and protein-dependent: a G→A change
can be neutral at one position of a motif and abolish binding at another. No fixed chemical prior
can express that.

**Reverse-complement mean pooling.** Measured 2026-08-19 against the tables: the 8-mer space is
**revcomp-collapsed**. There are 32,896 = (4⁸ + 4⁴)/2 distinct 8-mers, 256 of them palindromic,
and **zero** revcomp pairs with both strands stored. Each stored 8-mer is one arbitrary
representative, and the array measured both strands together, so the E-score is inherently
strand-symmetric. A naive encoder would have to learn that equivalence from data it never sees,
and would carry an arbitrary strand asymmetry the assay does not have. So: run the encoder over
both strands and **mean-pool**, making the embedding revcomp-symmetric by construction. Standard
practice in the DeepBind/Basset lineage, and it removes a real artefact rather than a
hypothetical one.

**What was rejected.** A genomic DNA language model: 8 bp is 1–3 BPE tokens for models trained on
512 bp–1 Mb of context, so the embedding would be dominated by tokenizer segmentation and two
8-mers one base apart can tokenise completely differently. 

A **learned lookup table** over the 32,896 8-mers was also rejected, but the argument for that
weakened on 2026-08-19 and is recorded honestly rather than left overstated. It was originally
*fatal*: no parameter sharing, so nothing transfers to an unseen 8-mer, which sinks a DNA-motif
holdout. **That regime is now dropped** (§5, §10.7) and splits are on proteins only — so every
8-mer appears in training and a lookup table would in fact train perfectly well. What is left is
weaker but still decisive enough: a table costs 32,896 x d parameters against a small CNN's few
thousand, it forecloses a DNA-side regime if one is ever wanted again, and a compositional encoder
is the better inductive bias for a space shared with the protein tower. One-hot + CNN stays; the
reason is now parsimony, not impossibility.

**Verify the continuity rather than assuming it.** After training, check whether 1-substitution
neighbours are close in the learned DNA embedding. 32,896 points is small enough to do
exhaustively in seconds, it tests the concern that motivated this choice, and it makes a slide:
*"the model learned which substitutions are equivalent, and here they are."*

**Variants deferred**, not dropped — chemical feature channels, alternative pooling functions,
DNAshape, and symmetry-by-augmentation are all one-encoder swaps against a fixed pipeline. They
are collected in `TODO.md` **T32**, to run after the first result exists.

### 4.2 Protein — three arms, from two modalities

**Revised 2026-08-19.** The third arm now carries the ensemble's **spread**, not its mean.

| # | modality | what it is |
|---|---|---|
| **A1** | sequence | one protein language model, the most standard available |
| **A2** | structure | embedding of the **single best** predicted structure |
| **A3** | structure | **A2 ⊕ the ensemble's spread**, projected back to the common width |
| **A4** | sequence | **ESM-DBP** — a PLM domain-adapted to DNA-binding proteins |

P1 was renamed A1 and so on, to stop the protein arms colliding with §6.1's `P1`/`P2`/`P3`
split regimes.

**Why spread rather than mean.** The ensemble mean is nearly redundant with A2 — both are "a
central conformation" — so a mean-vs-best comparison would confound two similar quantities and
answer nothing. Spread is the part of an ensemble that a single structure genuinely cannot carry:
*how much does this domain move, and where*. If a point mutation destabilises a domain, the
spread is where it shows and the mean is what erases it — which matters directly for **C1**, since
the variants are the hard case.

**A2 and A3 are a controlled pair.** A3 is A2 *plus* the spread, same embedder, same structures.
The only difference between them is the ensemble term, so the A2 → A3 delta isolates exactly the
question "does conformational spread carry signal beyond one structure?" That is a cleaner
experiment than the original best-vs-mean, and it is the one that speaks to **C2**.

**Spread alone is a diagnostic, not an arm.** A pure-spread representation may not identify the
protein at all — two unrelated domains can have similar flexibility profiles. Run it once as a
check (does a nearest-neighbour classifier on spread vectors recover domain identity?) rather than
as a headline arm; if it does retain identity, that is a finding, and if it does not, A3 is still
the right construction.

**A4 separates the modality from the pretraining corpus — accepted 2026-08-19.** A1 and A4 are the
*same modality* differing only in what the language model was trained on, so the A1 → A4 delta
answers a question the other arms cannot: when a protein-side gain appears, is it because the
representation is protein-*sequence*, or because the model was specialised to DNA-binding proteins?
Without A4 those two are confounded in A1.

ESM-DBP is ESM-2 fine-tuned on 170,264 non-redundant DNA-binding protein sequences from UniProtKB
(Zeng et al. 2024, *Nat Commun* 15:7838). Its width is **1280, the same as ESM-2**, so A1 and A4
also happen to enter the common projection on equal terms.

It is cheap in the way that matters here: sequence-only, so it needs nothing from a collaborator
and lands in build step 1 alongside A1 (§9).

**All four arms are one pooled vector per domain, and share one output dimensionality.** The first
half is what makes the arms comparable at all (§3.1); the second is what stops "structure wins"
from being "structure had more width".

**All four arms share one output dimensionality.** Each tower ends in a learned projection to the
same width `D`, so "structure wins" can never be "structure happened to be 1280-d against a 320-d
PLM". Raw embedding widths differ (a PLM is typically 1280, structural embedders vary), so the
projection layers differ in size — **report each tower's trainable parameter count alongside its
result**, since equal `D` controls the shared space but not the tower capacity feeding it.

### 4.3 Structures

Structures have to be built. **This task is handed to a colleague in the same lab**, who knows
the project; how the ensembles are generated is hers to decide or ours to discuss together. Only
two points are pinned from this side, both in `TODO.md` `T31`: a **fixed ensemble size**, and the
option of scoping the run to the **variant-bearing clusters only** (281 domains rather than 1,338)
— which would restrict arms A2/A3 to a subset of the protein axis, and so has to be decided
jointly rather than assumed.

The expectation is **10–50 structures per sequence** — an ensemble, not a single model.

**Three requirements the spread depends on:**

1. **A uniform ensemble size.** A standard deviation estimated from 10 samples is both noisier and
   differently biased than one from 50. If ensemble size varies per domain, the spread feature
   partly encodes *how many structures that domain happened to get*, which correlates with nothing
   biological and would leak straight into A3. Ask for a fixed `N` for every domain, or subsample
   every ensemble down to the smallest one before computing spread. Either is fine; silently mixed
   sizes are not.
2. **A common superposition.** Spread is only meaningful after the ensemble is structurally
   aligned; per-residue spread additionally needs a shared residue frame, which the domains already
   have via `mut_positions`' reference numbering.
3. **The per-structure embeddings, not just the structures.** Spread is computed by embedding each
   member and taking the elementwise standard deviation across the ensemble — so the embedder runs
   `N` times per domain, once per structure. Cheap, but it is `N` × 1,338 inferences, worth sizing
   before the run.

---

## 5. Experiments on SELEX

**Revised 2026-08-19: the talk is about the protein axis, and only the protein axis.** The
DNA-side holdout is dropped — here, on PBM, and in the transfer sets. That decides §10.7 by
removing the question rather than answering it, and it buys a sharper story for a 
biologist audience.

**The whole protein axis is in play at all times.** No subsetting: all 1,338 domains (SELEX: its
full panel), every experiment.

**5-fold cross-validation, in two split regimes**, each reported separately:

| regime | held out | tests |
|---|---|---|
| **S1 — random across proteins** | random *domains* | ceiling; the easy case |
| **S2 — similarity-based** | groups of similar domains | generalisation to genuinely new proteins |

**S1 splits on proteins, not on pairs.** Splitting (domain, 8-mer) pairs at random would put the
same domain on both sides and measure almost nothing. It is still the optimistic regime — variants
of one cluster can land either side — and that is its job.

**S2 groups by connected component, decided 2026-08-19.** The unit held out is not a CD-HIT
cluster but a **connected component** of the domain graph — an edge wherever two domains are
within `cluster.max_edits` (5) at `cluster.min_overlap` (0.6). That is what makes "completely new
proteins" true rather than nearly true: near-identical domains that the greedy algorithm assigned
to *different* clusters are held out together instead of one leaking into training (`T27`, settled
by this decision — `DECISIONS.md` §2).

Grouping is **evaluation-time only**. The stored clusters do not change; `wt_id`, `mut_positions`
and the cluster inventory are untouched.

**The chaining objection was tested and does not apply.** Connected components are single-linkage,
which this project rejected *for clustering* because it produced a 35-domain blob at 5 edits. That
measurement predates `T25`'s overlap floor. Re-measured 2026-08-19 with `min_overlap = 0.6`:

| | components | largest |
|---|---:|---:|
| connected components | **1,154** | **8 domains** |
| CD-HIT clusters (unchanged) | 1,165 | 8 domains |

Only **11 components merge more than one cluster**, the largest being 6 domains from 2 clusters.
No blob. The overlap floor is what killed the chaining, and the historical objection was measured
before it existed.

Chaining would in any case be the *safe* direction here: for **splitting**, over-grouping removes
more from training than strictly necessary but never leaks, whereas for **clustering** it makes
false claims about relatedness. That asymmetry is why this is right at evaluation time and wrong
in `build_clusters.py`.

Cost: **23 s**, 111,992 within-family alignments — and free in practice, since §8.1 builds the
same distances for the nearest-neighbour baseline.

### 5.1 Expected result

**Generally good performance on all three embedding combinations under the random split, and
degrading performance once the leaks are removed.**

The expected cause is a **lack of overlap between proteins** — the panel is broad and shallow,
so holding out a family removes anything the model could have transferred from.

**That degradation is only interpretable against the NN-lookup baseline (§8.1)**, which produces
the same shape. The claim is not that performance drops — it is *by how much, relative to
copying the nearest neighbour's profile*. The baseline goes on every figure as a horizontal
line.

### 5.2 Visualisation — settled 2026-08-19

A **UMAP of the shared protein–DNA space** the two-tower model produces (§3). Four panels:

| panel | what it shows |
|---|---|
| 1-3 | **all proteins as big dots**; one selected protein highlighted per panel, with **DNA in the background coloured binding / non-binding for that protein** |
| 4 | every DNA point coloured by **how many proteins bind it** |

Panels 1-3 are the owner's own second proposal and it **replaces** the first. Colouring DNA by
"its binding protein" is ill-defined the moment a site binds several, and selecting one protein
per panel removes the ambiguity rather than arbitrating it. It also reads instantly from the back
of a room.

Panel 4 turns that same promiscuity from a plotting nuisance into a result about the DNA landscape
— which is worth a slide of its own given §6.2 expects PBM's landscape to look more transient than
SELEX's.

**A UMAP of a contrastively-trained space shows clusters almost by construction. It illustrates;
it does not evidence.** So the claim *"proteins do not overlap, and that is why performance
degrades"* gets a **quantitative companion**: the distribution of nearest-neighbour sequence
identity between each held-out domain and its closest training domain, per split regime. That
histogram *explains* the degradation curve instead of gesturing at it, and it is nearly free — it
reads off the same 1,338 × 1,338 distance matrix §8.1 builds in step 0.

### 5.3 Metrics — positives first, everywhere

**Decided 2026-08-19, and it applies to §5, §6 and §7 alike.** The dataset is **466:1 negative to
positive** ([RESULTS.md](RESULTS.md)). Any metric that rewards correct negatives is measuring the
class balance rather than the model: calling every pair non-binding scores **99.8% accuracy** and
is worthless. Metrics must put the emphasis on the positives, at every point in the project.

| | metric | why |
|---|---|---|
| **primary** | **AUPR** / average precision | the standard answer at extreme imbalance; the baseline AUPR is the positive rate, so the number is honest about the prior |
| secondary | **precision@k**, **recall at fixed precision** | reads directly as "of the k 8-mers we called, how many bind" — the question a biologist actually has |
| secondary | **Spearman vs. the raw E-score** | uses the continuous score behind the label, so it does not throw away the gray band's information |
| **never headline** | AUROC, accuracy | both are inflated to near-uninformative at 466:1; two models far apart in usefulness differ in the third decimal |

AUROC may be *reported* alongside, since some readers expect it — but never as the number a claim
rests on. Note that TransBind (§8.2) reports AUPR next to AUROC for exactly this reason, so the
choice is also the comparable one.

**Computed per protein, decided 2026-08-19.** For each held-out domain, rank all 32,896 8-mers,
score that ranking, then **macro-average across domains** — every protein counts once, whatever
its positive count. Never pooled: the real question is *"does the model recover this protein's
motif"*, which is a retrieval problem per protein, and pooling hides exactly the per-protein
failures §5's S2 regime exists to expose.

Two consequences of going per-protein, both worth handling on the figure rather than in a
footnote:

**Report the distribution, not only the mean.** Positive counts run from 1 to 236 (median 49), and
a domain with a single positive has an essentially binary AUPR — did that one 8-mer rank high or
not. A violin or histogram of per-protein AUPR says far more than its average, and it is the
figure that shows *which* proteins fail rather than that some do.

**Per-protein AUPR is not comparable across split regimes on its own.** A random-baseline AUPR is
the domain's own positive rate, which varies from 1/32,896 to 236/32,896 — so regimes that hold
out different domain mixes (`P1` holds out homeodomain; `S2` holds out components) carry different
baselines, and a raw difference between them mixes that in. **The NN-lookup baseline solves it**:
computed on the *same* held-out domains, it absorbs the domain-mix effect, so the model-minus-
baseline gap is comparable across regimes even when the raw numbers are not. One more reason §8.1
belongs on every figure (§8.1) rather than in a table consulted once.

---

## 6. Experiments on the PBM dataset

Then we pull out the dataset we just built.

**The talk must start with SELEX.** The target audience is biologists, and the project requires
us to use SELEX — it is needed as a story element, not just as data.

On PBM, **repeat the exact same experiments as §5** — same 5-fold CV, same three regimes, same
three embedding combinations — and **additionally**:

- **two variant-holdout regimes**, §6.1. These carry claim **C1**; they are possible on PBM and
  impossible on SELEX, because only PBM has proteoform depth.

### 6.1 The three protein-axis variants — specified 2026-08-19

On top of **S1** and **S2** from §5, which run unchanged here, PBM adds three regimes that SELEX
cannot support. All counts measured from `data/interim/clusters/clusters.parquet`: **1,338
domains in 1,165 clusters — 1,057 singletons and 108 variant-bearing clusters holding 173
variants across 29 families.**

**The 1,057 singleton clusters stay in training throughout.** They are the backbone that keeps the
protein axis at full width in every regime; only the designated group is ever removed.

| | held out | test | train | train variants |
|---|---|---:|---:|---|
| **P1** | all `Homeodomain` | **428** domains | **910** | 89, across 28 families |
| **P2** | the 66 non-homeodomain variant-bearing clusters | **155** domains (89 variants) | **1,183** | 84, all homeodomain |
| **P3** | a fraction of all 173 variants | 173 / 86 / 43 | 1,165 / 1,252 / 1,295 | 0 / 87 / 130 |

#### P1 — hold out the whole homeodomain family

*Can an entire unseen family be predicted, having seen how single residues matter in other folds?*
Homeodomain is 428 of 1,338 domains and carries 84 of the 173 variants, so this is the one split
where the held-out family is large enough to be a real test and the training set is not dominated
by it. Training keeps 89 variants spread over 28 other families.

Unaffected by `T27` either way — clustering runs `same_family_only=True`, so every cluster (and
every connected component) lies inside one Pfam family; a family holdout removes both members of a
near-identical pair together.

#### P2 — hold out the non-homeodomain variant clusters

*Given the family is known and homeodomain's mutation structure has been seen, can mutation
effects be predicted in other folds?* Scoped 2026-08-19 to **only the variant-bearing clusters**,
not the whole families: all 1,057 singletons stay, so the model **has** seen these families —
just never a mutated member of them. That makes P2 a mutation-transfer test, not a family
holdout, and it is deliberately the gentler mirror of P1.

Together P1 and P2 bracket the transfer question from both ends: P1 learns from **variant breadth**
(89 variants, 28 families) to predict one thick family; P2 learns from **variant depth in one fold**
(84 variants, homeodomain) to predict thin ones.

#### P3 — drop all, half, or a quarter of the variants

*How much mutation supervision is actually needed?* Remove a fraction of the 173 variants from
training and test on the removed ones — a data-efficiency curve, with the **all** row as the
zero-shot extreme where training holds wild types only.

**The `all` row makes the §8.1 baseline the exact null hypothesis for C1.** With no variants in
training, the nearest training neighbour of any held-out variant *is its own wild type*, one edit
away — so NN-lookup predicts **precisely the wild-type profile**, which is the hypothesis *"the
mutation has no effect."* Beating it is claim C1, directly. The 20 `dead_variant` records
(`TODO.md` `T30`) are where the gap is widest: for a variant that measurably lost binding, the
wild-type profile is maximally wrong.

**Repeat `half` and `a quarter` over several seeded draws.** Which variants get dropped matters at
n = 86 and n = 43, so a single draw is a sample of one. Seed via `snp2prot.utils.seed` and report
the spread.

**Report stratified by family, always.** 84 of 173 variants are homeodomain, so a pooled variant
number is a homeodomain number wearing a general claim's clothes.

### 6.2 Expected result

The same UMAP, where we expect to see **more clustered proteins** (the PBM panel has dense
variant clusters that SELEX lacks) and a **more transient DNA landscape**.

---

## 7. Transfer — SELEX as a held-out test set for PBM

Finally, use the **SELEX data as a holdout test set for a model trained on PBM**, split by how
much of it the training set had already seen:

| set | overlap with PBM training data |
|---|---|
| **easy** | the protein was seen in training |
| **hard** | the protein was not |

**Revised 2026-08-19 from three sets to two, and the reason is structural rather than a
simplification.** §2.3 puts SELEX on the *same* 32,896 canonical 8-mer space as PBM, so **DNA
overlap is 100% by construction**. The original "hard" set — no overlap at all — is therefore
*empty and unreachable*, and the original "medium" is the only genuinely hard case there is. The
two sets above are what the design can actually contain, so they are named for what they test.
Do not present a DNA-novelty axis; there isn't one.

For **at least one** of our methods trained on PBM, show how well it transfers to a different
experiment. This is the finale.

**Scored with §5.3's metrics**, positives first — that matters more here than anywhere, because
the negatives are the part of the SELEX table we trust least (§2.3).

**Measure the negative reliability rather than assuming it.** The **easy** set — proteins *and*
8-mers present in both datasets — is a direct cross-assay agreement measurement, for free. It
answers "how much do PBM and SELEX agree on what does *not* bind?" before any model result is
read. If agreement on negatives is much worse than on positives, that is a property of the assay
pair and must be reported as one, not absorbed into the model's score.

**Put the ceiling on the figure.** `reports/overlap.md` measured it: for the 47 domains stored by
two *PBM* sources, the median pair agrees on **45.8%** of the 8-mers either called positive —
against 70-72% for replicates within one source. That is the noise floor between two labs running
the **same** assay. Across two different assays it will be worse. Drawing that line on the
transfer figure turns a mediocre-looking number into a result and pre-empts the obvious
objection.

---

## 8. Baselines

### 8.1 The nearest-neighbour lookup — the bar

**Accepted into the plan 2026-08-19.** `CLAUDE.md` already names it non-negotiable: *"The
nearest-neighbour baseline is the bar. Any model that does not beat it under
leave-one-cluster-out has learned nothing transferable."* It is implemented at
`snp2prot.baselines.nn_lookup`, currently a docstring stub.

**It is not a strawman — it is the incumbent method.** Transferring a motif to an
uncharacterised TF by DBD sequence identity is how CIS-BP does inference, and Weirauch et al.
2014 established the per-family identity thresholds for exactly that. That paper is
`weirauch2014`, one of this dataset's own 19 sources. So this is not a baseline we invented to
clear easily; it is the standard approach in the field, and it is the first question a
biologist audience asks — *"isn't this just copying the most similar protein's motif?"*

**Selecting the neighbour.** Percent identity over the aligned domain, from
`snp2prot.align.edit_profile`: `identity = 1 - n_edits / aligned_len`. The `overlap` guard is
mandatory, not optional — free terminal gaps are precisely what let two unrelated domains look
three edits apart, which produced 31 bad cluster memberships before it was caught (`T25`,
`DECISIONS.md`). Use the already-tuned `cluster.min_overlap` from `configs/thresholds.yaml`.
This matters most under leave-one-family-out, where the nearest neighbour is cross-family by
construction and so is the degenerate case by construction.

**What it predicts.** For a held-out domain `i`, take `j = argmin(D[i, train_idx])` and copy
domain `j`'s **E-score** profile — not its binary labels. The E-scores give a ranking over all
32,896 8-mers, which is the same output shape the contrastive model produces, so AUPR /
precision@k / Spearman are directly comparable between the two without any special-casing.

Primary form is **k = 1**: the pure lookup table, and the thing to beat. A distance-weighted
average over the top-k neighbours is one extra line and a slightly stronger bar; run it as a
secondary if time allows.

**Cost — measured 2026-08-19, not estimated.** `edit_profile` runs at **0.17 ms** per pair on
this corpus (median domain 77 aa). All-vs-all over 1,338 domains is **894,453 pairs ≈ 2.6
minutes single-core**, done once and cached as a 1,338 × 1,338 float matrix. Given that matrix,
`fit_predict` is roughly twenty lines.

The matrix is **shared infrastructure, not baseline overhead**. Three consumers need exactly these
distances: this baseline; **§5's S2 grouping**, which builds its connected components from them;
and §10.12's nearest-neighbour-identity histogram, which explains the degradation curve. Build it
once, in step 0.

**Why it cannot be replaced by inspecting the model's own curve afterwards.** The NN signature is
a *shape* — strong on random splits, degrading on holdouts. That shape is not diagnostic, because
every imperfectly-generalising model has it; a model that learned real biophysics still does
worse on unseen families than on seen ones, since that is what a generalisation gap is. What
separates the two cases is the *level*, and a level cannot be read off a single curve.

Concretely: leave-one-family-out AUPR comes back at 0.22. If the baseline scores 0.21, the model
is a lookup table. If it scores 0.05, that is the talk. Same number, opposite conclusions, and
nothing in the model's own curve distinguishes them. The slide cannot be written without the
baseline — so "build the model, then check whether it looks like NN" needs the very thing it was
trying to avoid building.

**Run it before any model is trained.** This is the stronger argument for doing it first. It is
the cheapest available test of whether the splits do what they are believed to do: if NN-lookup
scores near-ceiling under S2, the split groups leak. That was `T27` — near-identical domains in
different clusters — now handled by connected-component grouping (§5), and the baseline is the
cheapest check that the fix worked. Finding that out
in an afternoon beats finding it out after 45 training runs on splits that were never valid.

**Where it appears.** As a horizontal line on **every** figure in §5, §6 and §7. Not a table
entry consulted once.

### 8.2 Published comparators

Optional addition, if available: include other published methods as a baseline. The one named is:

- **TransBind** — Basnet & Cheng 2026, *NAR Genomics and Bioinformatics* 8(2):lqag047, DOI
  [`10.1093/nargab/lqag047`](https://doi.org/10.1093/nargab/lqag047). Protein-aware deep
  learning integrating DNA sequence with ESM-DBP protein embeddings via cross-attention;
  evaluated on 690 ChIP-seq experiments, 161 TFs, and supports label-zero-shot prediction for
  unseen TFs. Code at <https://github.com/jianlin-cheng/TransBind>, archived at
  `10.5281/zenodo.19462292`.
  **Filed 2026-08-19** as `docs/papers/basnet2026_transbind.pdf`.
  **Role, decided 2026-08-19: a benchmark baseline only** — never our architecture, and it
  produces no UMAP (§3: a cross-attention BCE head has no joint space to project).
  ⚠️ Its released weights operate on 101 bp genomic regions and **cannot** be run on our 8-mers,
  so the baseline is a **re-implementation of its architecture** — a PLM protein embedding plus
  cross-attention over the DNA encoder — trained on *our* data. Label it as a re-implementation
  on every figure, never as "TransBind's reported numbers".

**A4 makes this baseline close to controlled.** TransBind's protein encoder *is* ESM-DBP, so
running it against arm A4 holds the protein representation fixed and varies the architecture —
cross-attention + BCE against two-tower + InfoNCE. That is a far more interpretable comparison
than A1 vs TransBind, which would confound architecture with pretraining corpus.

**The comparison is now clean.** TransBind mean-pools ESM-DBP into one 1280-d vector and §3.1
does the same (decided 2026-08-19), so A4 and the re-implementation share both the protein encoder
*and* its pooling. The only thing varying between them is the architecture — cross-attention + BCE
against two-tower + InfoNCE — which is exactly what the baseline should isolate.

## 9. Order of work

**Strategically, start with the PBM ML part, sequence embeddings only.**

Reasons, as given by the owner:

1. It has **no collaborator dependency** — the PBM dataset exists today and the sequence
   embeddings need nothing from anybody.
2. **A lot of collaborators are on vacation.** Both the SELEX k-mer conversion and the
   structure ensembles are blocked on people who are not available.
3. Everything built for it is **easily transferable to SELEX** once that data lands.
4. It is **extendable by simply switching the embedding** to include structure.

So the build order and the talk order are deliberately **opposite**:

```
BUILD ORDER                        TALK ORDER
0. NN-lookup baseline              1. SELEX as framing  — the assay, why
1. PBM + sequence arms (A1, A4)          biologists care, what the panel covers
2. + structure                     2. PBM as the workhorse — the regimes,
3. SELEX, when it lands               and the variant story that only it has
                                   3. SELEX again as the transfer target
                                      — the finale, PBM-trained
```

**Step 0 comes before any model is trained** (§8.1). It sets the bar every later number is read
against, and it validates the splits before the training runs are spent on them.

**SELEX opens and closes, PBM does the work — settled 2026-08-19 (§10.15).** The requirement to
use SELEX is met by the framing and the transfer finale; it does not require SELEX to carry the
experimental middle. That is what lets the build order start with PBM and never build anything
twice, and it moved SELEX's own §5 experiments from *required* to *desirable*.

**What cutting them would cost, stated so the trade is visible:** §6.2 expects PBM's proteins to
look more clustered and its DNA landscape more transient *than SELEX's*. That comparison needs
§5 actually run on SELEX. If time forces a cut, §5 is the first thing to go — PBM carries the same
regimes plus the variant ones — but the SELEX-vs-PBM UMAP contrast goes with it.

### 9.2 Configuration, tracking and reproducibility — decided 2026-08-19

**A second config file, not an extension of `thresholds.yaml`.** `CLAUDE.md` already anticipated
this and gave the reason: changing a threshold invalidates the dataset, every report and every
provenance row; changing a learning rate does not. Two things with different blast radii do not
belong in one file. New file, read by a module mirroring `snp2prot.thresholds` — that pattern
already works and is about thirty lines.

**Seeds live in the experiment config**, not in code and not in call sites. Every fold assignment,
every `P3` variant draw and every UMAP is seeded from there and the seed is recorded with the run.
`snp2prot.utils.seed` stays as the mechanism.

**Experiment tracking: MLflow, local backend, behind a thin wrapper.**

The grid is real — per arm on PBM: 5 folds `S1` + 5 folds `S2` + `P1` + `P2` + 7 `P3` runs
(3 fractions, with `half` and `a quarter` repeated over seeded draws) ≈ **19 runs per arm, ~76
across four arms**, and past 110 if §5's SELEX experiments run too, before the TransBind baseline. That is well past what a
notebook and a filename convention can keep straight.

| tool | verdict |
|---|---|
| **MLflow** | **adopt.** Local-first — a plain `mlruns/` directory, no server, no cloud — which matches a project whose data is git-ignored and whose provenance is deliberately self-contained. Logs params, metrics and artifacts, which is exactly the need. |
| Weights & Biases | equivalent in function and nicer to look at, but cloud-hosted. Fine if the owner prefers it — the wrapper below makes it a one-line swap. |
| Hydra | **not yet.** Its value is config composition and multirun sweeps. The grid here is *enumerated*, not searched, and a loop over configs covers it. Revisit if configs start nesting. |
| Optuna | **not yet.** It solves hyperparameter *search*. Nothing here is being searched — the arms and regimes are fixed. Add it only if tuning becomes a real bottleneck. |

**Put a thin `snp2prot.tracking` wrapper in front of whichever is chosen.** The owner's brief was
"MLflow, W&B, whatever", and a wrapper is what makes that "whatever" cheap later instead of a
rewrite.

**Log the split, not just the hyperparameters.** This is the project's own standard applied to
modelling: rule 5 makes every row carry the cutoffs actually applied, so a table describes its own
binarization. A run is the same — a metric is meaningless without knowing exactly which domains
were held out. **Record the held-out domain list (or its hash) as a run artifact**, alongside the
regime name, the fold index and the seed. Given §5's connected-component grouping and §6.1's three
regimes, "which domains were in test" is not reconstructible from the regime name alone.

**Dependencies all go in `dev` for now**, not a separate `ml` extra. They will grow as
implementation proceeds; splitting out a runtime subset is a question for publication or
deployment, and answering it now would be guessing.

**`results/figures/` and `results/checkpoints/` already exist** in `snp2prot.config`
(`FIGURE_DIR`, `CHECKPOINT_DIR`) and are fine to use — and free to be reworked if the tracking
layout suggests better ones. Nothing builds a path by hand either way.

### 9.1 Blocking, right now

- `docs/papers/jolma2026_expanded-codebook.pdf` is **not downloaded**. It blocks the SELEX
  binarization rule, which has to come from its methods section and not be guessed
  (`docs/papers/README.md`, known gaps).
- The SELEX k-mer conversion is with collaborators.
- The structure ensembles are with a colleague.

None of these block step 1 of the build order, which is the point of doing it first.

---

---

# 10. Review — risks and proposed changes

**Claude's, not the owner's. Nothing here is decided or acted on.** Items worth turning into
tracked decisions are cross-referenced to [TODO.md](../TODO.md).

Ordered by how much damage each does if ignored.

## 10.1 The NN-lookup baseline is missing — **accepted 2026-08-19, now §8.1**

Raised as a review item, argued, and **accepted into the plan**. The method, the neighbour
selection rule, the measured cost and the reasoning are now **§8.1**, and it is **step 0 of the
build order** in §9.

Kept here as a pointer rather than deleted, because the reasoning is why it is in the plan: the
degradation curve the whole story rests on is *the same shape* a nearest-neighbour lookup
produces, so without the baseline's *level* the headline result cannot be distinguished from a
lookup table — and the field's incumbent method for this task (CIS-BP motif transfer by DBD
identity, Weirauch 2014) **is** that lookup.

## 10.2 A DNA language model is the wrong tool for 8-mers — **accepted 2026-08-19, now §4.1**

Raised, agreed, and **the method was switched**. The DNA LM is out; the encoder is plain one-hot
with reverse-complement mean pooling, specified in **§4.1**. Deferred variants are `TODO.md`
**T32**.

Kept as a pointer because the reasoning is load-bearing. Two findings came out of the discussion
and are recorded in §4.1: that a one-hot input still yields a *continuous learned* embedding (the
continuity belongs to the encoder, not the parameterisation), and — measured against the tables —
that the **8-mer space is reverse-complement collapsed**, 32,896 = (4⁸ + 4⁴)/2 with 256
palindromes and zero both-strand pairs. The second was not in the original review and turned out
to matter more than the chemical-similarity question that prompted it.

The one thing §10.2 argued for that did **not** survive: it proposed running the DNA LM as a
second comparison arm. The owner wants a single DNA method for now, so there is no second arm.

## 10.3 Both datasets must speak one DNA vocabulary — **accepted 2026-08-19, now §2.3**

Raised and **settled immediately: only 8-mers.** The full spec to hand to collaborators is
**§2.3**; `TODO.md` `T28` is now the task of sending it rather than the question of what to ask.

One thing the review did not have and §2.3 now does: the **canonicalisation rule**. Checking the
tables showed all 19 sources share one identical 32,896-entry 8-mer set, and each entry is the
**lexicographically smaller** member of its reverse-complement pair. That is the detail that would
have failed silently — a differently-canonicalised table joins cleanly on nothing.

## 10.4 SELEX negatives are weaker than PBM negatives — **resolved 2026-08-19**

**Largely dissolved by the data being exhaustive.** The review assumed "exhaustive binding
k-mers" meant a positives-only list; the collaborators can in fact score **every** 8-mer,
positives and negatives. That makes a SELEX negative *measured and scored below threshold* —
`assayed_unbound` under `schema.VALID_NEG_PROVENANCE`, not `selection_absent` — so rule 4 is
satisfied and the finale is not built on a label-semantics mismatch. Recorded in **§2.3**.

**What survives, and is kept on purpose.** A SELEX negative is still derived from depletion in a
selection rather than from a direct spot measurement, and enrichment is noisy for 8-mers poorly
covered in the input library. So **SELEX may be less reliable on the negative side.** §7 now says
to *measure* that instead of assuming either way: the transfer experiment's **easy** set is a
free cross-assay agreement check, and the 45.8% cross-lab ceiling from `reports/overlap.md` goes
on the figure. That also absorbs §10.5.

The review's proposed mitigation — score the finale on positives-only ranking — is superseded by
the broader metric decision in **§5.3**, which applies positives-first metrics everywhere rather
than only here.

## 10.5 Set the ceiling before showing the transfer numbers — **accepted 2026-08-19, now §7**

Folded into §7 alongside §10.4: the measured **45.8%** median cross-lab agreement from
`reports/overlap.md` is drawn on the transfer figure as an explicit ceiling.

## 10.6 The variant experiment — **specified 2026-08-19, now §6.1**

The owner supplied the two regimes the review asked for: **V1**, leave one variant-bearing family
out, and **V2**, keep every family and drop a fraction of the variants. Both are in **§6.1** with
measured feasibility counts.

**Two things resolved that the review had wrong or missing.**

- **`T27` does not threaten P1.** Clustering is family-restricted, so the near-identical
  cross-cluster pairs are all *within* one family and a family-level holdout removes them
  together. The review filed `T27` against "the split regime the degradation story rests on"
  without distinguishing cluster-level from family-level holdout. It bit the similarity regime
  (§5) only — and that is now settled there by connected-component grouping.
- **V2 at "drop all variants" turns §8.1 into the exact null hypothesis for C1.** With no variants
  in training, NN-lookup's nearest neighbour for every held-out variant is its own wild type, so
  the baseline predicts the wild-type profile — i.e. *"the mutation has no effect."* Beating it is
  C1. Neither the review nor the original plan had noticed this; it is the cleanest result
  available here and it costs nothing extra.

**What still stands from the review.** The axis is thin and lopsided: **84 of 173 variants are
homeodomain**, so §6.1 requires results stratified by family. `Homeodomain` is the only
near-balanced V1 split (84 out, 89 retained); every other family trains mostly on homeodomain.
And `T30` — whether the 20 `dead_variant` records become a dedicated C1 evaluation set — is still
open, and §6.1 now shows exactly why they matter.

## 10.7 "Left-out DNA motifs" is not yet a definition — **dissolved 2026-08-19**

**The DNA-side holdout was dropped rather than defined.** The talk rests entirely on the protein
axis (§5), so there is no DNA-motif regime to specify, on PBM or on SELEX. The review's proposal —
cluster the 8-mer space and hold out whole clusters — is not needed and was not adopted.

One consequence the review had not connected: because §2.3 puts both datasets on the same 32,896
8-mer space, **DNA overlap between PBM and SELEX is 100%**, so §7's three-way easy/medium/hard
split was already unreachable. It is now two sets, keyed on the protein alone. Recorded in §7.

If a DNA-side regime is ever wanted again, the grouping problem the review raised is still real:
random 8-mer holdout is not leak-free, because a held-out 8-mer typically differs by one base from
something in training.

## 10.8 C2 as phrased is not establishable — **accepted in full 2026-08-19**

All three parts adopted.

- **The claim is now phrased as inductive bias** (§1, C2): *a structure-derived representation is
  a better inductive bias for this task than sequence alone*. The owner confirms this was the
  intended meaning. Predicted structures carry no information beyond the sequence they were
  predicted from, so a gain can only be that the geometry is easier to learn — and stated that
  way it survives contact with a biologist audience.
- **The ensemble is represented by its spread, not its mean** (§4.2, arm A3). The review proposed
  spread as an *addition*; the owner prefers it as the *replacement*, which is better: the
  ensemble mean is nearly redundant with the single best structure, so best-vs-mean would have
  confounded two similar quantities. A2 and A3 are now a controlled pair differing only in the
  spread term.
- **All arms share one output dimensionality** (§4.2), so "structure wins" cannot be "structure
  had more width". Equal `D` controls the shared space but not tower capacity, so parameter counts
  get reported too.

**Two things §4.2/§4.3 add that the review did not have.** Spread alone may not identify the
protein at all, so pure-spread is specified as a diagnostic rather than an arm. And spread needs a
**uniform ensemble size** — a standard deviation from 10 samples is noisier and differently biased
than one from 50, so mixed sizes would let A3 encode how many structures a domain happened to get.
That requirement is now in `T31`.

The review's warning that point mutants are the hard case for structure still stands and is not
resolved by any of this: AF2-family predictors return near-identical backbones for single-residue
variants. The spread formulation is the best available answer — if a mutation destabilises a
domain, the spread is where it shows — but whether it is enough is exactly what A2 → A3 measures.

## 10.9 Pin the structure vocabulary with the colleague — **discarded 2026-08-19**

Dropped by the owner. The colleague is in the same lab and shares the project's context, so the
review's written-spec-and-vocabulary framing was solving a problem that does not exist here; those
decisions are hers or are made jointly.

Two points from it survive in `T31` because they are consequences of decisions made *on this side*
rather than requests about her method: a **fixed ensemble size** (arm A3 represents an ensemble by
its spread, so mixed `N` would leak ensemble size into the feature), and the reminder that scoping
to the **variant-bearing clusters** would cut the run to 281 domains — with the caveat that A2/A3
could then no longer run the full-axis regimes unless every arm is restricted to the same subset.

## 10.10 An engineering note that changes the schedule — **accepted 2026-08-19, now §3.1**

Adopted; the owner had the same shape in mind. The dense-matrix formulation, the measured memory
budget and the precomputed-embedding lookup are **§3.1**.

**A recommendation this item made was subsequently overruled, and correctly.** The review argued
protein embeddings must stay **per-residue**, because pooling attenuates a single substitution and
C1 depends on that signal. The owner chose **pooled** on 2026-08-19, for a reason the review had
not weighed: a structural embedder yields one vector per structure, so per-residue sequence arms
against per-structure structure arms would have made A1 → A2 confound modality with architecture —
and that comparison is claim **C2**. Cross-arm comparability outranks per-arm sensitivity here.

The review also overstated the cost. "One part in 77" is a worst case assuming only the mutated
residue's representation changes; ESM-2 is **contextual**, so a substitution shifts its neighbours
too, and the true attenuation is empirical. §3.1 turns it into a **pre-flight measurement** rather
than an assumption, and names attention pooling as the fallback that would keep every one of the
owner's reasons intact.

## 10.11 The contrastive objective and the UMAP figures are coupled — **accepted 2026-08-19, now §3**

Accepted as written: **two-tower with multi-positive InfoNCE**, giving a shared protein–DNA space.
The owner confirms this was the intent and the reason the UMAP was chosen — a cross-attention BCE
model has no joint space, so there would be no figure. TransBind therefore stays a benchmark
baseline (§8.2) rather than a candidate architecture.

Batch shape followed from it and is settled in **§3.1**: shuffle the protein axis, keep the 8-mer
axis complete, since with InfoNCE the full denominator is exact *and* makes the training objective
identical to §5.3's per-protein ranking metric.

**One correction to the review, from the owner.** §3.1 previously said an all-negative row was
unusable. It is not — pushing every 8-mer away from a protein is real signal. What is true, and
narrower, is that *plain* multi-positive InfoNCE cannot express it, because averaging over an empty
positive set yields no gradient. §3.1 now specifies a **learned `null` anchor** as the mechanism,
and records the catch: on PBM this applies to exactly 54 records, of which 20 are the
`dead_variant` C1 evaluation set and 34 are unverifiable without a control.

## 10.12 The two UMAP proposals — **accepted in full 2026-08-19, now §5.2**

All of it adopted: the owner's second proposal replaces the first, plus the fourth panel colouring
DNA by how many proteins bind it, plus the nearest-neighbour-identity histogram as the
*quantitative* companion to a figure that would otherwise show clusters by construction.

## 10.13 TransBind as a baseline — **resolved in full 2026-08-19**

TransBind is a **benchmark baseline only**, re-implemented (its released weights cannot run on
8-mers) and labelled as such — §8.2.

**ESM-DBP accepted as a fourth protein arm, A4** (§4.2). It is the same modality as A1 with a
different pretraining corpus, so it separates *"the protein-side gain is the modality"* from
*"the gain is the pretraining domain"* — confounded in A1 alone. Sequence-only, so it costs nothing
in schedule and lands in build step 1.

A consequence worth having: since TransBind's own protein encoder is ESM-DBP, **A4 vs the
re-implementation holds the protein representation fixed and varies the architecture**, which is
a much more interpretable baseline comparison than A1 would have given. The residual difference the
review flagged here — TransBind mean-pools — **disappeared** when the owner moved our own
embeddings to pooled on 2026-08-19 (§3.1), so A4 and the baseline now share encoder *and* pooling,
and only the architecture varies.

## 10.14 Metric choice — **accepted in full 2026-08-19, now §5.3**

Positives-first metrics everywhere, and — the part left open until now — **computed per protein**,
macro-averaged across held-out domains, never pooled.

§5.3 adds two consequences the review had not drawn out: the per-protein **distribution** matters
more than its mean when positive counts run 1-236, and per-protein AUPR is **not comparable across
split regimes** on its own, because each regime holds out a different domain mix with a different
random baseline. The NN-lookup baseline fixes the second by being evaluated on the same held-out
domains.

## 10.15 A story-order alternative — **accepted 2026-08-19, now §9**

Adopted. SELEX frames the talk and closes it as the transfer target; PBM carries the experimental
middle. The owner's reason is speed: it satisfies the "the project requires SELEX" constraint
without requiring SELEX to be experimentally complete first, so nothing is built twice and the
critical path stays on the data that already exists.

§9 records the consequence the review did not: this demotes §5's SELEX experiments from required
to desirable, and if they are cut the §5-vs-§6 UMAP contrast in §6.2 is what is lost.

## 10.16 Smaller points — **resolved 2026-08-19**

| point | outcome |
|---|---|
| Codebook publishes PBM with E-scores | **check it, unhurried** — `TODO.md` `T33`. Acquisition stays closed; this is a look, not a reopening |
| collaborator-derived k-mers need a provenance decision | folded into `T28`, the task that collects them |
| modelling dependencies | all in `dev` for now; a runtime subset is a publication/deployment question |
| a second config file | **yes** — §9.2 |
| existing `results/` paths | fine to use, free to rework — §9.2 |
| seeding | moved into the experiment config and the tracker — §9.2 |

The one thing the owner added beyond the review: **experiment tracking**. §9.2 recommends MLflow
on a local backend behind a thin wrapper, and argues against Hydra and Optuna *for now* — the grid
is enumerated rather than searched, so neither earns its complexity yet. It also adds the
requirement the review did not have: **log the held-out domain list**, because with
connected-component grouping and three PBM regimes, "which domains were in test" cannot be
reconstructed from a regime name.
