# Decisions and resolved items

The record of what was **decided**, what was **fixed**, and what was **excluded** — the
counterpart to [`TODO.md`](../TODO.md), which holds only what is still open. An item leaves
`TODO.md` by arriving here.

This is a log, not a specification. How the dataset is built is
[`docs/METHODS.md`](METHODS.md); what `dbd_seq` is and what the admission policy excludes is
[`docs/DOMAIN_POLICY.md`](DOMAIN_POLICY.md); the owner's brief is
[`docs/TFDNA_MERGE_BRIEF.md`](TFDNA_MERGE_BRIEF.md).

Entries before 2026-08-14 carry their original `#N` from the archived open-items list
([`reports/archive/OPEN_ITEMS_2026-08-14.md`](../reports/archive/OPEN_ITEMS_2026-08-14.md)) so
older commit messages and code comments still resolve.

---

## 1. Scope

### 2026-08-18 — `D1`, `T18`, `T22`: the protein table reconnected, and its limits labelled
**`D1` decided by the owner: 25-34% full-length coverage is acceptable.** The domain sequence
is the representation this dataset is for; a full-length protein is a bonus where it resolves
and never a filter. Domain-level and construct-level work covers all 1,334 domains.

**`T18` — the protein table lost 111 domains to a substring test.** `build_protein_table.py`
located a stored domain inside its construct with `dbd in construct`, which fails by
construction for any canonical sequence that borrowed flank from a reference protein
(`snp2prot.canonical`). Those are now placed by alignment, with the span read off the
alignment so an indel cannot shift the window. The same fallback was added to
`proteins.locate_in_protein`, which had the same defect one level up: an engineered variant
appears verbatim in no wild-type protein, so **every variant** was failing to map onto its
full-length sequence.

| | before | after |
|---|---:|---:|
| domains with a protein record | 1,223 of 1,334 (92%) | **1,303 (98%)** |
| located by substring / by alignment | 1,223 / 0 | 1,223 / **80** |
| domains with a full-length sequence | 439 | 453 |
| domains located inside that sequence | 363 | **445** |

The 31 that remain unlocated are constructs whose stored domain matches nothing their source
publishes: 27 in CIS-BP, plus 2 `PNAS08` AP2, 1 `SCI09` GCM and 1 `SHO18A` bZIP.

**`T22` — `BAR15A:SIX6` pointed at the wrong protein, and it is not alone.** The stored
accession was `Q6P051`, a 305 aa TrEMBL entry titled *"SIX6 protein (Fragment)"* from a cDNA
clone; reviewed `SIX6_HUMAN` is `O95475` at 246 aa. The stored `SIX6_REF` domain is a
substring of both, at offset 182 in the TrEMBL entry and **123** in the reviewed one — exactly
the 59-residue shift by which our numbering disagreed with the literature. Corrected through
`ACCESSION_OVERRIDES` in `build_protein_table.py`, a corrections map that carries the evidence
for each entry and holds this one alone.

The audit `T22` asked for found the general case: **41 of the accessions the deposits publish
are unreviewed (TrEMBL) and 9 are flagged fragments.** Rather than chase them, the protein
table now records `uniprot_reviewed` and `uniprot_fragment` per row, so a full-length sequence
carries its own provenance and nothing downstream has to assume Swiss-Prot.

### 2026-08-18 — `D3`: the reference is the medoid, and it is not the seed
**Decided by the owner.** A cluster now carries two named domains instead of one. The **seed**
is what greedy assignment compared candidates against — longest first, and it decides
*membership*. The **reference** is the medoid, the member with the smallest total distance to
the others, and it decides the *coordinate frame* `mut_positions` is expressed in.

**Why one field could not do both.** Inside a variant series every domain is the same padded
length, so the seed fell out of the alphabetical order of the amino-acid string. In
`C:BAR15A:HOXD13` that elected `HOXD13_S316C`: the wild type was stored as a one-mutation
variant, all seven disease variants as two-mutation ones, and **position 50 — S316C's own
substitution — was written into all seven siblings' `mut_positions`**. Any per-position
analysis of this corpus would have found a recurrent hot spot that is one variant's own
change. `C:BAR15A:ARX` had the same pattern at position 14. Across a paralogue merge it looked
different and was the same defect: BAR15A's human HOXB7 series framed against `Cell08`'s mouse
copy, every disease variant reading three edits with two of them mouse-human differences.

**Why the medoid and not the seed order.** The medoid cannot seed the clustering — it is only
defined once the members are known, and the members are only known once something has seeded.
Resolving that circularity means an iterative k-medoids refinement, which is a different
algorithm, changes cluster *membership* rather than only its description, and would reopen a
settled decision whose blast radius is the splits. Longest-first also suits the data: with
free terminal gaps a clipped copy matches a longer seed, which is the direction the padding
artefact runs.

**Effect.** 26 of 108 multi-member clusters take a new reference; 113 domains change their
`n_mut_from_wt` and `mut_positions`; the total edits reported across cluster members falls
from 480 to **398**, so 17% of the recorded protein-axis variation was an artefact of the
frame. Cluster membership, counts, label health and the merge resolution are all unchanged:
1,161 clusters, 108 multi-member, 173 variants, 20 `dead_variant` / 34 `no_evidence`, 48
records dropped at merge.

**One invariant weakens.** Members are within `max_edits` of the seed by construction, but the
triangle inequality allows `2 × max_edits` from the reference. Measured: worst case 5 of 5,
zero clusters over. `build_clusters.py` now checks and reports rather than assuming.

A reference is a frame, not a claim about ancestry — for a cluster of genuine paralogues there
is no wild type and the medoid is simply the most central member. Implementation:
`clusters.medoid` / `clusters.references`, `seed` and `reference` columns in the inventory,
4 new tests.

### 2026-08-18 — `D2`: `MAR17A:Esrrb` is excluded, and unresolved residues are now a rejection
**Decided by the owner.** The domain carried two `X` residues in an 89 aa zf-C4 window — the
depositor's uncertainty, not a property of the protein. A sequence model would embed the
ambiguity as a token and a structure predictor cannot place the residue at all.

Implemented as a **general rule rather than a named exclusion**: `domains.call_domain` rejects
any padded window containing a character outside the 20 standard residues
(`unresolved_residue`, now in `DOMAIN_POLICY.md`), and `schema.validate` refuses such a
`dbd_seq` outright, mirroring the non-ACGT check on `dna_seq`. `MAR17A:Esrrb` was the only
domain in the corpus affected; any future one is caught at parse time with a reason.

Corpus after the single-source rebuild: 1,334 domains, 1,382 records, 45,462,272 rows. The
full rebuild later the same day settled at **1,338 domains, 1,386 records, 45,593,856 rows** in
1,165 clusters (`docs/METHODS.md` §9.3) — the difference is the 13 stale sources catching up
with the current parsers, not this decision. Variant counts are untouched: `Esrrb` was a
singleton.

**The rebuild also exposed a stale table, which is worth recording.** `MAR17A` on disk
predated a parser change: re-running it dropped `Esrrb` *and added* `MAR17A:Mlx` (HLH, 78 aa),
which the current parser admits and the committed table did not. That is exactly what `TODO.md`
`N6` warns about, and `MAR17A` now joins the list of sources whose stored form is newer than
the rest. A full rebuild would move other numbers too.

### 2026-08-18 — `T25`: free terminal gaps need an overlap floor, and 31 memberships fail it
**Decided by the owner and applied.** `cluster.min_overlap: 0.6` is new in
`configs/thresholds.yaml`: a domain joins a cluster only if the alignment matched residues over
at least that fraction of the shorter sequence.

**What was wrong.** With terminal gaps weighted zero, an arbitrarily long prefix of one
sequence and suffix of the other can be discarded for free, and `align.edit_profile` then
counts edits only over what remains. That is no longer a distance between two domains — it is
the edit distance of the best-matching window after free trimming, and it has no triangle
inequality. Two unrelated Myb domains of 72 and 60 aa align on six residues, score +14 against
−22 for the honest alignment, and report **3 edits where a charged alignment reports 59**.

**How much of the corpus it touched.** Every one of the 202 cluster memberships was measured:
**31 of them (15%), across 18 clusters, were formed on 3-10% overlap** — Myb 13, AP2 8, HLH 7,
Homeodomain 3. Genuine memberships align at ≥90%. Nothing falls between 10% and 90%, so the
floor is not a tuned number; any value in that band rejects the same 31.

**What it cost, and what it bought.**

| | before | after |
|---|---:|---:|
| clusters | 1,133 | **1,162** |
| clusters holding a variant | 122 | **108** |
| variant domains | 202 | **173** |
| variants that are homeodomain | 87 of 202 | 84 of 173 |

The correction falls hardest on the families the breadth claim leant on — **Myb 16 → 4 and
AP2 9 → 1** — so "the protein axis is no longer the bottleneck" was overstated by more than
the totals suggest. `mut_positions` and `n_mut_from_wt` for the 31 were positions in a frame
the variant had never been aligned to; they are now correct or gone.

It also closes two demonstrated leaks: an HLH domain and a Myb domain each sat in one cluster
on a spurious match while a genuine ≤5-edit relative sat in another, which under
leave-one-cluster-out puts near-identical sequences on both sides of a split.

**Why the floor is safe.** It is measured against the *shorter* sequence, so the case free
terminal gaps exist for — the same domain clipped to less padding, ROG18A's 83 aa bare against
105 aa padded — is wholly contained and scores 1.0. The subsumption pass that folds
short-window duplicates together now demands 1.0 exactly, which is what "one is the other with
different padding" means.

Implementation: `overlap` on `align.EditProfile`, `min_overlap` in `snp2prot.clusters.cluster`,
5 new tests. Label health and the merge resolution were rebuilt on the new clusters and are
unchanged: 20 `dead_variant` / 34 `no_evidence`, 48 records dropped at merge.

### 2026-08-18 — `T15` and `D4`: one record per domain, the one with more positives
**Decided by the owner.** 47 domains are stored by more than one source, 95 records in all,
and their labels do not agree — the median pair on 46% of the 8-mers either called positive,
two pairs on almost nothing (`C:LIN14B:NAP` 12 against 126 sharing none; `C:Cell08:Tlx2` 20
against 14 sharing three). Left alone they hand a sequence model identical input with two
different labels, inside one cluster, where no split can separate them.

**The record with more positives wins — unless one candidate belongs to a variant series, in
which case the series wins.** 48 of the 95 records are dropped at merge, 1,579,008 rows
(3.5% of the corpus). The rule is applied at *merge*, never at parse: `data/interim/` keeps
every measurement, so `reports/overlap.md` still measures the noise floor from the full
evidence and lists every resolution.

Why more positives: a PBM fails by missing binding, not by inventing it — a weak array
compresses its E-score distribution and calls fewer positives (`METHODS.md` §5.1) — so the
deeper measurement is the more informative one. Intersecting the two instead would inherit the
worse array's sensitivity everywhere and would leave `C:LIN14B:NAP` with no positives at all.

**The series exception was added after measuring what the plain rule does.** It flips 6
clusters, all of them wild types: `BAR15A`'s `ARX_REF` has 188 positives against `Cell08`'s
206, so more-positives would take `Cell08`'s copy and leave ARX's five `BAR15A` variants to be
compared against a wild type from another lab on another array — across a 46% cross-source
noise floor that dwarfs any single-residue effect. Where one candidate's source also supplies
the other domains of that cluster, that source wins.

Two properties made the rule safe to apply corpus-wide, both checked and both worth re-checking
whenever a source is added: **no duplicate pair disagrees about whether the protein binds at
all** (all 95 records carry positives, so the rule never discards a measured non-binding), and
**there are no ties** (1.8x margin at the median, 10x at the worst).

This also settles `T15`'s general question — the corpus stores **one row per sequence**, not
one per (sequence, protein) — and with it the species contradiction under `C:Cell08:Tlx2`:
`Cell08` wins on 20 positives against 14, so the domain is stored as *Mus musculus*.

Implementation: [`src/snp2prot/merge.py`](../src/snp2prot/merge.py), 6 tests, and the
resolution table in [`reports/overlap.md`](../reports/overlap.md).

### 2026-08-18 — `T20`: the cutoff stays at `E >= 0.45`
**Decided by the owner**, on the sweep in [`reports/threshold_review.md`](../reports/threshold_review.md).
`configs/thresholds.yaml` is unchanged: `positive: 0.45`, `negative: 0.35`,
`per_experiment: true`.

The sweep found 0.45 is the **maximum over every absolute cutoff tried** — median cross-source
Jaccard 0.520, against 0.509 at 0.40 and 0.500 at 0.47 — so loosening it adds 8-mers the two
labs disagree about more often, not less. A rank-matched rule scores 8-12 points higher, and
that was rejected on three grounds:

1. Its advantage is partly mechanical. Forcing equal counts removes the 1.6x positive-count
   asymmetry from the denominator, so what it measures is how much of the cross-lab
   disagreement is *sensitivity* rather than *ranking*. That is a diagnosis, not a rule.
2. "Every protein binds exactly N 8-mers" is not a fact about biology. Specificity breadth
   genuinely differs between families.
3. It would break `neg_provenance`, and with it rule 4. Negatives would become "outside this
   experiment's top N" — a ranking artifact — rather than measured non-binding. On the 18
   `BAR15A` variants that lost binding it would invent ~100 positives each, destroying the
   signal the corpus exists to carry.

**What the E-score does and does not guarantee, since this is what the decision turns on.** It
is a rank statistic — for each 8-mer, the Mann-Whitney AUC of the probes carrying it against
the rest, within the brightest half of the array, shifted to [-0.5, +0.5]. Every monotone
change of laser power, protein concentration, antibody or scanner cancels, which is why one
cutoff can span 19 sources and two decades of arrays. `E >= 0.45` says "at most about three of
this 8-mer's ~32 probes are less than convincing". What it does not equalise is **sensitivity**:
a noisy array shuffles ranks, every shuffle costs AUC, and the distribution compresses. Two
labs on one protein differ 1.6x in positive count (`Hoxa2`: 165 against 17). That asymmetry is
a measured property of the corpus, recorded in `docs/METHODS.md`, and it caps the agreement any
model can be expected to reach — the same fact as the 45.8% noise floor.

### 2026-08-18 — `T21`: label health is a side table, not a threshold change
**Decided by the owner.** `Cell09` and `LIU18B` are not re-binarized, not excluded, and not
left silently indistinguishable from measured non-binding. Instead
[`data/interim/label_health/`](../data/interim/label_health/) records, per
`(dbd_seq, source_dataset)` record, whether it has positive evidence — and filtering on it is a
**training** decision made at featurization, exactly as cluster-size restriction is (§ above,
2026-08-14). The dataset stores evidence; what a model is fed is chosen later.

Three verdicts, of 1,383 records:

| verdict | records | rows | meaning |
|---|---:|---:|---|
| `ok` | 1,329 | | at least one 8-mer at or above the cutoff |
| `dead_variant` | 20 | 657,920 | no positives, but another record of the same cluster **and source** has them — its own series is the control, so this is measured non-binding |
| `no_evidence` | 34 | 1,118,464 | no positives and no such control: indistinguishable from an experiment too weak to see |

`no_evidence` is 2.5% of the corpus and is what `label_health.usable(df)` drops; 18 of the 20
`dead_variant` records are `BAR15A` variants that lost binding, and keeping them is the point
of drawing the distinction at all.

Two things the build turned up on the way:

- **13 of the 54 silent records are not weak — their replicates disagreed.** They reach 0.45
  on their stored score, which is the replicate mean, while the label is `LABEL_GRAY` because
  the two experiments disagreed. `n_at_cutoff` sits beside `n_pos` in the table for this
  reason.
- **`frac_silent`, not `median_max_escore`, is the source-level statistic.** A source that
  deliberately contains dead variants drags its own median down; `Cell09` stands out on both
  (71% silent, median 0.428) but `BAR15A` at 20% silent would look bad on the median alone.

Report: [`reports/label_health.md`](../reports/label_health.md), rebuilt by
`scripts/build_label_health.py` (12 s). No label, threshold or parsed table was touched.

### 2026-08-18 — Kock et al. 2024 excluded; PBM acquisition is closed
**Decided by the owner**, after the source was acquired, screened and rescored. Full record in
[`reports/kock2024_excluded.md`](../reports/kock2024_excluded.md). The corpus stays on its two
source families — UniPROBE (18 accessions) and CIS-BP / Weirauch 2014 — and nothing further is
sought.

What was on offer: **67 new domains**, 87 of 93 clones admitted, variants 202 → ~269, and 13
singleton clusters becoming variant-bearing. Real, and the reason this was not an easy call.

What decided it, in one sentence: **the deposit publishes no E-scores, and neither of the two
ways in gives labels comparable with the rest of the corpus.**

- Its own statistic, `affinityEstimate`, has no transferable scale — the value at our
  `E = 0.45` boundary runs 10.82-12.14 across four proteins while the `E = 0.35` boundary runs
  10.08-11.16, so one protein's positive threshold is another's negative. `affinityQ` is a
  detection statistic, calling 6-13% of 8-mers positive against our 0.3%. Entering on those
  terms meant a second `score_type` and two new cutoffs (the 2026-08-14 decision below, now
  moot).
- Recomputing E-scores from the raw `.gpr` scans was built and validated against `BAR15A`
  rescored from Barrera's own scans: the **ranking** reproduces (Spearman 0.78-0.91 where a
  protein binds; rank-matched positives at Jaccard 0.57-0.77, against 0.70-0.72 for replicates
  within one source) but the **values** do not — `E >= 0.45` selects 0-32 8-mers where the
  stored values select 129-204, and the count-matching cutoff is 0.312-0.343. That is one
  recalibrated number and a mixed-provenance corpus, 19 sources on published E-scores and one
  on ours.

Neither price was worth 67 domains in a dataset whose whole value is that every row is
comparable. Removed with the decision: `T11`, `T19`, `T24`, `src/snp2prot/rawpbm.py`,
`scripts/validate_escores.py`, `tests/test_rawpbm.py`, `data/external/pbm_design/` and its
`PROVENANCE.md` rows. Kept, because they are defects in our own data that the screening
exposed: `T22` (`BAR15A:SIX6` carries a TrEMBL fragment accession, 59 residues off) and `T23`
(padded windows have never been checked against the reference proteome). `T20` is untouched
and now stands on its own.

### 2026-08-25 — `T33`: the Codebook PBM panel is excluded, on story grounds

**Checked, and it is good data — which is exactly the problem.** Jolma, Laverty et al. 2026
*Nature* ([`10.1038/s41586-026-10798-9`](https://doi.org/10.1038/s41586-026-10798-9)) does run
PBMs, and runs them on the same platform we already parse: *"we analysed proteins on two
different **universal PBM arrays (HK and ME)**"* — the two array designs `weirauch2014` uses,
standard protocol, so the 8-mer E-scores on the portal would need no new binarization rule.

Three findings decided it.

**1. It would break the SELEX transfer experiment.** The study assayed **393 proteins across up
to five assays**, drawn from one pool:

| assay | TFs assayed | motifs identified |
|---|---:|---:|
| GHT-SELEX | 392 | 179 |
| HT-SELEX | 392 | 184 |
| ChIP–seq | 373 | 179 |
| SMiLE-seq | 299 | 84 |
| **PBM** | **173** | **63** |

So **at least 172 of the 173 PBM proteins were also SELEX-assayed.** [`ML_PLAN.md`](ML_PLAN.md)
§7 splits the SELEX test set into *easy* (protein seen in training) and *hard* (not seen).
Training on the Codebook PBM would move nearly the whole SELEX panel into *easy*, empty the hard
set, and reduce the finale to measuring the same protein on two platforms — which is the
transfer claim itself.

**2. It adds no protein-axis depth, which is the axis that matters.** The paper contains **zero**
occurrences of "mutant", "substitution", "point mutation", "missense", "alanine" or "wild type",
and no engineered constructs. Its 41 mentions of SNPs and 21 of alleles are all **DNA-side** —
allele-specific binding at genomic SNPs in ChIP–seq and GHT-SELEX peaks. The only within-protein
multiplicity is construct architecture (full-length versus DBD) and protein source, which our
policy treats as replicates of one domain, not a variant series. The only within-family
multiplicity is *"the handful of paralogues analysed"* — SP140/SP140L, DACH1/DACH2,
CAMTA1/CAMTA2, ZXDA/ZXDB/ZXDC — natural relatives with no wild-type control, and far enough
apart to land in separate clusters anyway.

So the yield would be **singleton proteins**: 63 motifs, of which 39 are novel Codebook TFs and
the rest controls we likely already hold. [`ML_RESULTS.md`](ML_RESULTS.md) §3.2 is the reason
that settles it — the one fold where the model beat the baseline is the one with 155 same-family
training domains, so depth is what predicts success and 39 more singletons buy none of it.

**3. The 8-mers there are PBM 8-mers, not SELEX-derived.** The paper never uses the terms
"8-mer" or "E-score"; those come from the PBM deposit. The SELEX→k-mer conversion the
collaborators have described is a separate capability and remains `T28`'s ask.

**What is not excluded.** The paper stays the source of record for the SELEX phase
([`ML_PLAN.md`](ML_PLAN.md) §2), and its PBM panel remains a legitimate *external validation*
set if one is ever wanted — held out entirely, never trained on. That use has none of the
problems above.

### 2026-08-14 — The dataset is PBM only
**Decided by the owner.** B1H and SNP-SELEX are out; the corpus grows by extending PBM
coverage instead. Then merge, then modelling. Tier 4 held-out sets are still wanted.

Consequences, all of them simplifications:
- every stored row is an 8 bp site, so the assay-identity leak through `dna_len` cannot occur
  while this holds (`TODO.md` `N5`);
- **padded-20 bp vs. common-core stops being a question** (was `#3`) — there is nothing of a
  different length to reconcile;
- **the SNP-SELEX OBS cutoff stops being a question** (was `#5`) — it was the last blocker on
  a phase that no longer exists;
- `b1h:` and `snp_selex:` in `configs/thresholds.yaml` become dead config (`TODO.md` `T8`).

### 2026-08-14 — upbm Q-values may enter as their own `score_type`
**Superseded 2026-08-18: Kock was excluded, so no second `score_type` exists.** Kept because
the reasoning is what the exclusion was weighed against.

If the Kock deposit turns out to publish only upbm affinity/contrast/specificity Q-values and
no probe-level or E-score data, it is still admitted — with a new `score_type` and its own
cutoffs, rather than being skipped to protect the single-E-score scale.

**The cost is real and was accepted knowingly.** Every stored row is currently `pbm_escore` on
one scale with one pair of cutoffs, and the validator makes any value outside [-0.5, 0.5] a
hard error precisely because reading the wrong column once produced a 4:1 ratio (§5, `#36`).
A second score type means a second cutoff pair in `configs/thresholds.yaml`, and by this
project's own convention that invalidates the dataset, every report and every provenance row.
Weighed against 122 alleles in 30 designed series — the single largest addition of protein-axis
depth identified — and the depth won. Tracked as `TODO.md` `T11`.

### 2026-08-14 — Homeodomain may go to ~58% of domains; depth beats balance
**Superseded 2026-08-18.** The projection was made before CIS-BP landed 868 domains; measured
against the corpus as it now stands, admitting Kock would have moved homeodomain from 29.7%
to 33.3%. Moot either way — the source was excluded.

Admitting Kock takes homeodomain from 47.9% to roughly 58% of domains and from 67% to ~84% of
point variants. Accepted as a straight consequence of the earlier call that family imbalance
is fixed by parsing more rather than down-sampling (`#10`), and because multi-member clusters
roughly double — 28 to ~58 — which is the corpus's scarcest property.

It does make the transfer question harder to answer, not easier: whether single-residue
sensitivity generalises to another fold is what `TODO.md` `N2` is about, and this deepens the
fold the corpus is already richest in. That was the trade, made with the numbers in view.

### 2026-08-14 — Clusters are formed by sequence distance, not only by construct lineage
Two natural paralogues one substitution apart previously sat in separate singleton clusters,
because `wt_id` came from construct naming — a gene folder and its insert names. They now
merge. Worth roughly a dozen new multi-member clusters, at no acquisition cost, and it
subsumes the distance-0 case that `TODO.md` `T4` describes.

`wt_id` thereby stops meaning "one reference and the variants engineered from it" and starts
meaning "domains within *k* edits of each other". Two parameters follow that this decision did
not settle, and neither has a defensible default — the threshold, and which member of a merged
paralogue pair is the reference that `mut_positions` is expressed against. Open as `D3`.

### 2026-08-14 — Construct architecture is recorded beside the table, not in it
Flank length, affinity tag and expression system are perfectly confounded with source study,
so they go in the protein-side companion table keyed by construct — not as columns in the
22-column row schema. Keeps the schema frozen and rewrites no parser, while still letting a
modeller condition on the covariate or hold it out.

Noted at the time: the confound is already partly live. `dbd_seq` is the padded envelope
clipped where the construct ends, so flank length is readable off `dbd_seq` length today —
ROG18A's bare domains run 83-85 aa against 95-105 aa for the padded ones.

**Closed 2026-08-17 as recorded rather than as fully built.** `data/interim/proteins/` carries
`construct_seq` with `dbd_start` / `dbd_end` marking the padded envelope inside it, so the
flank actually present on each side is exact and derivable per construct — the part of the
covariate that varies within a source and is therefore the part worth having. Affinity tag and
expression system are *not* recorded and are not planned: both are constant within a source
and would be read off `source_dataset` anyway, so a column would restate the source name with
extra steps. If a source ever varies its tag across constructs, this reopens.

### 2026-08-14 — Cluster-size restriction is a training decision, not a dataset one
Was `#2`: whether to admit only DBDs with ≥5 variants. **The dataset stores everything; the
modeller filters.** Only 9 of 425 clusters hold 5 or more domains, so the rule would have cut
the protein axis to almost nothing — but that is beside the point, which is that a dataset
should not bake in a training-set choice. What remains is the convenience of filtering by
cluster size cheaply: `TODO.md` `T3`.

### 2026-08-17 — `T5b` done: one organism convention, enumerated not inferred
The convention is the **UniProt-style binomial** — `Genus species`, no strain, no
abbreviation, no hybrid marker — and anything that is not an organism is stored as the empty
string rather than as a word that looks like one. Corpus-wide: **137 distinct strings down to
132**, and 21 constructs now hold an empty organism where seven different placeholders used to
sit.

Handled in `snp2prot.metadata.organism`, wired into all three parser families
(`uniprobe_panels`, `bar15a`, `weirauch2014`):

| was | is | why |
|---|---|---|
| `C. elegans` (NAR11, 4) | `Caenorhabditis elegans` | abbreviated genus; every other deposit writes it out |
| `Acyrtosiphon pisum` (1) | `Acyrthosiphon pisum` | transposed `h`; the corpus carried both spellings for one organism |
| `Malus x domestica` (2) | `Malus domestica` | hybrid marker; UniProt and NCBI both index it without |
| `Chimera` (ROG18A, 12) | *empty* | engineered chimeras have no source organism |
| `N/A` (LIU18B, 2) | *empty* | reconstructed ancestors, likewise |
| `PBM CONSTRUCTS` (3) | *empty* | a stray heading from CIS-BP Table S6 |
| `None Available` | *empty* | UniPROBE's own missing-value string |

**The corrections are enumerated, not inferred from a pattern.** An alias table with a reason
per entry can be reviewed; a rule that rewrites any string of a given shape will eventually
rewrite a name that was correct. Every entry was found by surveying the distinct values
actually present.

Three names are deliberately left alone, and the tests pin each: `Sarsia sp. Long Island
Sound` is a genuinely unnamed species rather than a formatting variant, and trimming it to a
binomial would assert a species nobody has assigned; `Acanthamoeba polyphaga mimivirus` is a
virus, correctly named in three words; `Physcomitrella patens` is *Physcomitrium patens* under
current taxonomy, which is a revision of the name rather than a disagreement about how to
write it, and following it would put this module in the business of tracking taxonomy.

**Not fixed, because it is not a naming problem:** the domain shared by `Cell08` and
`weirauch2014` is byte-identical and stored as *Mus musculus* by one and *Homo sapiens* by the
other. One of them is wrong about the protein, and no normalisation can decide which.

### 2026-08-17 — `T3` done: a cluster side table, not a `cluster_size` column
The open design question was a stored `cluster_size` column against a side table.
**Side table.** A column would answer size queries by Parquet predicate pushdown, but it
denormalises a per-cluster fact onto 45 million rows and requires touching `schema.py` and
every parser's output; the schema stays frozen and no parser was rewritten.

`data/interim/clusters/clusters.parquet` — one row per cluster: family, representative, domain
count, variant count, maximum edit distance, contributing sources, construct and row counts.
Written by `scripts/build_clusters.py`, which is the only pass that already knows the
assignment, and read through `snp2prot.clusters.load` / `ids_with_at_least` / `select`.

Two things settled along the way:

- **"Cluster size" means distinct canonical domains** — not rows, not constructs. One domain
  assayed by two laboratories counts once. Asserted in a test, because counting constructs
  instead would inflate 47 clusters.
- **Row counts are accumulated during the rewrite pass, not in a second scan.** The inventory
  needs a row count per cluster and the rows are already in memory; a second pass over 45
  million rows would cost more than everything else in the script.

Incidental fix: three places discovered the corpus with an ad-hoc `"/proteins/" not in path`
string test, which would silently not have covered the new `clusters/` directory.
`config.source_tables()` is now the single discovery helper and `config.NON_SOURCE_INTERIM`
the single list of companion tables.

### 2026-08-17 — Full rebuild, and two reproducibility defects it exposed
The owner waived the no-long-running-commands rule once for a complete rerun. Everything
reproduced: 45,495,168 rows, 1,335 domains, 1,133 clusters, 122 multi-member, largest 8, and
the audit sweep clean apart from the two known excluded accessions. Timings are recorded in
`METHODS.md` §10.2 — the full build is **13 minutes**, of which `build_dataset.py --all` is
**7 min 28 s**, not the 4-5 min previously claimed in `TODO.md`.

Two defects surfaced that only a fresh build could show:

- **`reports/clusters.md` was not stable under re-run.** It stated "after dropping *N*
  superseded short-window copies", where *N* is 17 on a fresh build and 0 on a repeat: it
  described the state of `data/interim/` when the script ran, not the dataset. A committed
  report whose diff depends on how many times a step was repeated cannot serve as evidence, so
  the line now states the dataset property and the run-dependent count goes to stdout only.
- **A single-source rebuild silently detaches that source from its clusters.**
  `build_dataset.py --source X` writes `wt_id` in lineage form while the rest of the corpus
  carries the cluster form. The table validates, the row count is right, and the domain simply
  leaves its cluster. `build_clusters.py` must follow any single-source rebuild; it is
  idempotent, and the cluster table was verified byte-identical across two consecutive runs.

The clusters report also gained the multi-domain family breakdown and the fifteen largest
clusters, which the inventory made free to produce.

### 2026-08-17 — Tier 4 dropped; `T7` discarded
**Decided by the owner.** The bHLH dimer held-out sets are no longer wanted. The two problems
`T7` was opened to resolve go with it, unresolved and now moot: a heterodimer is two chains
forming one binding unit and so fails admission condition 1, and MAX's substitutions are "in
and around" the DBD so condition 3 needed checking per variant.

Consequences: step 4 of the plan is gone, and the `tier4:` block in `configs/thresholds.yaml`
became dead config and was removed with the others (below). `data/testsets/` and
`config.TESTSET_DIR` are kept as empty scaffolding — any future held-out set belongs there
whatever it measures, and the directory costs nothing.

### 2026-08-17 — `T8` done: the dead threshold blocks are gone
`b1h:`, `snp_selex:` and `tier4:` are deleted outright from `configs/thresholds.yaml` — no
tombstone comments, no test asserting their absence. They configured phases 3, 4 and Tier 4,
all dropped; they carried null cutoffs and `TODO` markers for work that will not happen, which
reads as unfinished configuration rather than as deleted scope.

Safe because `thresholds.for_assay` looks blocks up by key and no caller asked for any of the
three. This file is where the removal is recorded; git history is where the content is. The
config file is the owner's control surface and should show only what is live.

### 2026-08-17 — `T5` done: unresolved organisms recovered from UniProt, or left empty
All 17 constructs carrying the literal `$species` were in `Cell09`, all with a Swiss-Prot
accession. `snp2prot.metadata.organism` reads the organism from the `OS=` field of the
canonical FASTA already cached under `data/external/uniprot/`: **13 resolve to *Caenorhabditis
elegans*, 4 do not** — their accessions are cached as zero-byte, UniProt's negative cache for
an entry it no longer serves — and those are stored as empty. Metadata only; no label, cluster
or boundary moved.

Two choices inside this worth stating:

- **Nothing here fetches.** The lookup is a local read of an existing cache, so parsing stays
  offline and a cold cache degrades to an empty organism rather than to a network call inside
  a parser.
- **The study's scope is not used as evidence.** Every affected construct is from one
  *C. elegans* panel, so filling the remaining 4 from the paper would be right today — and
  would be a rule that invents an organism the next time a mixed-species deposit has the same
  template defect. Rule 1 covers accessions and URLs; the same reasoning applies here.

Rebuild cost: `Cell09` alone (4 s), then `build_clusters.py` corpus-wide (79 s) because a
single-source rebuild resets that source's `wt_id` to its lineage form and desynchronises it
from the cluster ids everything else carries. **That is a general trap, not a one-off** — any
single-source rebuild needs the cluster pass after it. Verified identical afterwards: 1,335
domains, 1,133 clusters, 122 multi-member, largest 8.

What this does **not** fix is the other half of `T5`: `species` still mixes `C. elegans` with
full binomials across sources, and `reports/overlap.md` found one domain stored as *Mus
musculus* by `Cell08` and *Homo sapiens* by `weirauch2014`. Normalisation was not asked for
and was not done.

### 2026-08-17 — PBM acquisition closes after Kock; `T2` and `T2b` dropped
**Decided by the owner.** Extending PBM coverage ends with `T11` (Kock et al. 2024). Two
tasks go with it:

> Closed 2026-08-18: Kock was itself excluded (§1), so acquisition ended with the sources
> already parsed.

- **`T2` — survey other PBM deposits.** Individual GEO / ArrayExpress submissions and paper
  supplements, always subject to the same gate as everything else: does the deposit publish
  the assayed construct sequence? Dropped as not worth the yield now that CIS-BP has landed
  868 domains and the corpus spans 56 families.
- **`T2b` — resolve UniProt accessions for the CIS-BP domains.** Table S6 publishes gene name
  and species but no accession, so 884 domains have no full-length sequence. Resolving them
  from gene plus species across 124 organisms would not be clean, and a wrong mapping is worse
  than a missing one — the reasoning that left `PP15` unparsed. Dropped; the consequence is
  that full-length coverage stays at about a quarter of the corpus, which is what `D1` is
  about. Domain-level and construct-level representations are unaffected and complete.

### Earlier — B1H and C2H2 arrays dropped on structural grounds
`#1`, `#4`, `#19`, `#20`. Persikov's B1H varies a different subunit than the one that binds,
and a C2H2 zinc-finger array is 2-6 separate ~23-residue folds on flexible linkers, each with
its own Zn(2+) — not one continuous unit, so it fails admission condition 2. ~8,000 domains
excluded. Reversible via `domain.allow_repeat_arrays`, and the weak-negatives worry about B1H
(`#1`) became moot when the phase went.

---

## 2. What `dbd_seq` is, and how clusters are formed

### `#18` — `dbd_seq` is the Pfam envelope padded by 10 residues each side
Not the bare envelope, not what was on the array. Pfam's models start at the structural core;
for homeodomains that clips the N-terminal arm, which reaches into the DNA minor groove and
makes real base contacts. Three BAR15A variants (KLF11 R402Q, VSX1 G160D, SNAI2 D119E) mutate
residues 2-9 N-terminal of the Pfam start — trimmed to the bare envelope each collapses onto
its own wild type while carrying a different label, invisible to any model. 10 covers the
observed maximum of 9 with margin, and is set from the biology rather than fitted to those
three cases.

### `#44` — Cluster distance is alignment-based, with free terminal gaps
Distance was Hamming, which excluded every variant carrying an indel and every variant whose
padding had been clipped differently. It is now a pairwise alignment (BLOSUM62, affine gaps).
**Terminal gaps cost nothing**, because `dbd_seq` is the padded envelope *clipped where the
construct ends*: in ROG18A every bare domain is 83-85 aa while the padded ones run 95-105 aa,
purely from how much padding fitted. Charging for that reported 8 phantom indels between
FoxJ3 and its own chimera. Internal indels still count.

### `#45` — `mut_positions` is in the reference's coordinate frame
So every member of a cluster shares one coordinate system and maps onto one predicted
structure. A variant with a deletion may therefore name a position past its own length. The
validator bounds positions by the cluster's reference and rejects a cluster carrying variants
but no reference row. One entry per edited residue, so `len(mut_positions) == n_mut_from_wt`
always holds.

### `#47` — BAR15A refuses to rebase an allele across an indel
Its construct-frame positions are shifted onto the padded domain by a constant offset, which
an indel invalidates. No current allele triggers this; a future one is rejected with a reason
rather than silently mis-positioned.

### `#30` — ROG18A's engineered chimeras cluster with their parents
Six substitutions is close enough. `FoxN3_J3` differs from `FoxN3` by 6 substitutions in
97 aa — real protein-axis depth, which the corpus is short of. Implementing it exposed a much
larger defect, `#39` below.

---

### 2026-08-19 — `T27`: splits group by connected component, clustering is unchanged

**The problem.** CD-HIT greedy clustering is order-dependent: a domain joins the *first* seed
within `max_edits`, not the best. Two paralogues 1-5 edits apart can therefore land in different
clusters — 13 such pairs were measured across Homeodomain, Myb, AP2 and HLH, including
`BAR15A:VAX2`/`Cell08:Vax1` at 1 edit. Holding out one cluster left its near-twin in training,
which is leakage in exactly the regime that claims to test generalisation to new proteins.

**The decision.** Splits hold out **connected components** of the domain graph — an edge wherever
two domains are within `cluster.max_edits` (5) at `cluster.min_overlap` (0.6) — not CD-HIT
clusters. **Evaluation-time only: the stored clusters are untouched**, so `wt_id`,
`mut_positions` and the cluster inventory keep their meaning. Specified in
[`ML_PLAN.md`](ML_PLAN.md) §5 as the `S2` regime.

**Why this is not the single-linkage that was already rejected.** `snp2prot.clusters` rejects
single-linkage because it produced a 35-domain blob at 5 edits. **That measurement predates
`T25`'s overlap floor.** Re-measured 2026-08-19 with `min_overlap = 0.6`, over all 1,338 domains
(111,992 within-family alignments, 23 s):

| | groups | largest |
|---|---:|---:|
| connected components | 1,154 | 8 domains |
| CD-HIT clusters | 1,165 | 8 domains |

Only **11 components merge more than one cluster**, the largest 6 domains from 2 clusters. The
chaining the earlier rejection was based on does not occur once a membership must also cover 60%
of the shorter sequence.

**And chaining would be the safe direction here anyway.** For *splitting*, over-grouping removes
more from training than strictly necessary but never leaks. For *clustering*, it makes false
claims about which domains are variants of one another. That asymmetry is why connected components
are right at evaluation time and would still be wrong in `build_clusters.py`.

**Cost: none.** The same 1,338 x 1,338 distance matrix is built anyway for the nearest-neighbour
baseline (`ML_PLAN.md` §8.1, step 0 of the build order).

**Superseded in one parameter on 2026-08-19 — see `D5` below.** The grouping stands; the
threshold that defines an edge moved from `<= 5 edits` to `>= 0.5 identity`, because the
baseline measured how little the 5-edit version changed.

### 2026-08-19 — `D5`: the `S2` edge is 50% identity, not 5 edits

**Raised by the nearest-neighbour baseline, which [`ML_PLAN.md`](ML_PLAN.md) §8.1 put in step 0
to raise exactly this**: *"if NN-lookup scores near-ceiling under `S2`, the split groups leak"*.

**What the measurement said** ([`reports/nn_baseline.md`](../reports/nn_baseline.md)). Grouping
by connected component at 5 edits did the job `T27` specified — the share of held-out domains
with a >= 90% identical neighbour still in training fell from **24% under `S1` to 7%** — and it
cost the lookup **0.018 AUPR**, 0.769 to 0.751. The two facts are consistent because removing
near-twins was never where the performance came from:

| nearest training neighbour | share of `S2` holdouts | lookup AUPR |
|---|---:|---:|
| >= 0.9 identity | 7% | 0.947 |
| 0.7 - 0.9 | 42% | 0.927 |
| 0.5 - 0.7 | 26% | 0.778 |
| < 0.5 | 25% | 0.359 |

Copying a **70-90% identical** neighbour already scores 0.927. 5 edits on a 77 aa domain is
about 94% identity, so the old floor removed only the top row of that table.

**The decision, by the owner: group at `>= 0.5` identity.** The regime exists to be hard —
*"this should actually be something difficult for an ML algorithm to solve"* — and 0.5 is well
below the 60-70% band Weirauch et al. 2014 set for transferring a motif between TFs by DBD
identity. The parameter is `splits.s2_min_identity` in `configs/experiment.yaml`, not in
`thresholds.yaml`: the grouping is evaluation-time only and changes no stored table.

**It chains, and that was measured rather than assumed:**

| edge rule | components | largest | largest component's median internal identity |
|---|---:|---:|---|
| `<= 5` edits | 1,154 | 8 | — |
| `>= 0.9` identity | 1,103 | 8 | — |
| `>= 0.7` identity | 746 | 36 | — |
| **`>= 0.5` identity** | **396** | **273** | **0.377**, with 12% of its pairs at the floor |
| `>= 0.4` identity | 242 | 406 | — |

The 273-domain component is 64% of the homeodomain family, and it is a **chain** — most of its
members are not within 50% of each other, only of something that is. Two consequences, both
accepted deliberately:

1. **Chaining is the safe direction for a split.** Over-grouping removes more from training than
   strictly necessary but can never leak. The cost is statistical power, not validity — the same
   asymmetry that makes single linkage wrong in `build_clusters.py` and right here.
2. **One of the five `S2` folds is now close to a homeodomain holdout**, so `S2` and `P1`
   partly overlap and fold-to-fold spread will be wide. Per-fold numbers are reported and the
   spread is the honest reading; a pooled `S2` mean now hides more than it did.

### 2026-08-19 — `D6`: C1 is evaluated on the variants where the wild-type copy fails

**The problem, measured.** `P3/all` holds out every variant and leaves every wild type in
training, so the nearest-neighbour baseline copies each variant's own wild type — which *is* the
hypothesis "the mutation has no effect" ([`ML_PLAN.md`](ML_PLAN.md) §6.1). It scores **mean AUPR
0.928 and median 1.000**. For most of the 173 variants the wild-type profile simply is the right
answer, because most single substitutions do not measurably change which 8-mers a domain binds.
A model could be perfect and gain two points.

**The decision, by the owner: report C1 on the part where the baseline goes wrong.** A held-out
variant enters the **C1 evaluation set** when its wild-type copy scores below
`c1_set.max_baseline_aupr` (0.7), or when it scores nothing at all because the variant has no
positive 8-mer — a `dead_variant`, whose AUPR is undefined rather than zero.

**29 of 173 variants**: 18 that lost binding entirely and 11 whose profile changed enough that
their own wild type does not predict it, across 5 families (Homeodomain 17, Forkhead 5, zf-C4 4,
HLH 2, THAP 1). Built by `snp2prot.evaluation.c1`, written to
`data/processed/c1_variants.parquet`, fixed before any model runs so every arm is scored on the
same variants.

The cut at 0.7 sits in the sparse part of the distribution — 8 variants fall below 0.6 and 14
below 0.8 — so it is not perched on a cliff. `no_evidence` records are excluded as everywhere:
with no control there is no telling a true non-binder from a failed assay (`T21`).

**The caveat is part of the decision.** The set is *defined* by the baseline scoring badly on it,
so the baseline scores badly on it by construction and "the model beats the baseline here" is
close to vacuous. What makes it mean something is the model's **absolute** number on the set,
reported next to its number on the complement — a model that had merely learned to distrust wild
types everywhere would gain here and lose there, which is a shifted prior and not C1. And with
5 families carrying the whole set, any C1 claim from it is a claim about those folds.

**This settles half of `T30`.** The 20 `dead_variant` records are 18 of these 29 (two are cluster
references rather than variants, so no `P3` fold holds them out). Whether the 34 `no_evidence`
records ever become training signal through the `null` anchor is still open.

## 3. Composition decisions

### 2026-08-17 — `T4` closed: the label noise floor is measured, and it is not small
Forty-seven domains are stored by more than one source — 48 duplicate constructs, 49 source
pairs of one protein against one 32,896-8-mer set. The task had two jobs and both are now
discharged.

**Job 1, a shared split group, fell out of clustering.** These pairs are distance 0, so CD-HIT
puts them in one `wt_id` by construction and no split can separate them. Nothing further to do.

**Job 2, label agreement, is measured** by `reports.overlap_report`, regenerated into
[`reports/overlap.md`](../reports/overlap.md) by `scripts/make_overlap_report.py`. Two passes
over the interim tables, because only ~3M of 45M rows are replicated.

The number that matters is not the one that looks reassuring:

| measure | value |
|---|---|
| pooled hard-label agreement | 99.98% (300 of 1,572,985 calls differ) |
| **median positive-call Jaccard** | **0.458** |
| pooled positives agreed on | 2,863 of 5,816 (49.2%) |
| median E-score rank correlation | 0.608 |

Positives are under 0.5% of an 8-mer table, so agreeing on the negatives holds pooled
agreement near 100% however badly two labs agree about binding. Against the within-source
figures already in `METHODS.md` §6 — 70.4% and 72.2% positive-set Jaccard for byte-identical
distinct genes and for technical replicates — cross-laboratory agreement is some 25 points
worse. **A model that reproduces held-out positives much past this is reproducing a source.**

Nine pairs are flagged below Jaccard 0.20 or rho 0.30, reported and not repaired. One is not
noise on any reading: `C:LIN14B:NAP` — *Arabidopsis* ANAC092 / O49255, an identical 144 aa
stored domain in both `LIN14B` and `weirauch2014` — shares **zero** positive calls between the
two deposits, at rho 0.068. Twelve positives against 126, disjoint. Open as `D4`.

A second, smaller inconsistency surfaced with it: the domain shared by `Cell08` and
`weirauch2014` under `C:Cell08:Tlx2` is byte-identical but carries `species` as *Mus musculus*
in one source and *Homo sapiens* in the other. Metadata only, no label affected, folded into
`T5`.

### `#9` — The 318:1 negative:positive ratio stays in the stored table
Controlled when sampling training batches, not by discarding rows. A PBM low-E-score negative
is real evidence of non-binding and throwing it away destroys information. `pbm.positive`
stays at 0.45.

### `#10` — Family imbalance is fixed by parsing more panels, not by down-sampling
Homeodomain has fallen from 79.5% to 47.9% of domains this way.

### `#27`, `#33` — Family-share and family-size worries, closed
The brief's 40% homeodomain flag was written because "C2H2 will otherwise swamp everything";
C2H2 is now 4 domains (0.8%), so the threat it guarded against does not exist. Homeodomain is
also *more* internally diverse than the minority families — median pairwise identity 34%, only
1.1% of pairs ≥90% identical, against Forkhead's 20.3%. Nothing is being swamped. Which
families are large enough to hold out is now a note, `TODO.md` `N4`, not a decision.

### `#24`, `#32` — Protein-axis depth accepted as-is
81 point variants in 28 clusters. `#32` framed this as the project's central gap — three
families with zero point variants and 54 of 66 variants homeodomain. The `#39`/`#46` fixes
moved Forkhead from 6 variants to 18 and HLH from 0 to 3, so the premise no longer holds as
stated. **Owner: acceptable while point variants are not the majority of the data**; subsets
can be built later, and the question resurfaces on its own when modelling starts. Carried
forward as `TODO.md` `N2`.

---

## 4. Individual constructs

### `#7` — Keep all 120 BAR15A variant alleles, not the paper's 117
Trust what is on disk. The three extra are HOXD13 `I297V`, `N298S`, `Q325K`, whose files are
named `*_8mers_11111111.txt` — exactly what a naive glob drops, which is the likely origin of
117. HOXD13 survives the policy with all 8 domains.

### `#8`, `#41` — `PROP1_R112Q` kept as deposited
Its name says residue 112 Arg→Gln; the deposited sequence gives Arg→**Met** at that same
position. Position and wild-type residue are both right; only the replacement disagrees.
PROP1's other allele `R99Q` matches its name exactly, and the Barrera supplement lists R99Q
but **not** R112Q. Both Arg→Gln and Arg→Met are single-nucleotide changes from an Arg codon,
so nothing available settles which the clone really was.

It matters more than most: construct position 57 maps to **homeodomain position 43, the first
residue of helix 3 — the recognition helix** that reads bases in the major groove.

Kept, because the allele name is stored nowhere in the schema. `dbd_seq` is the sequence and
`mut_positions` (53, within the padded domain) is derived by diffing against the reference,
never parsed from the name. The stored claim — "this domain, with M at position 53, binds
these 8-mers" — is true whichever name is right. **The one exposure is a future join on allele
name `R112Q` that assumes Gln; anyone doing that must read this entry first.**

`PITX2_T114P` remains excluded: its insert shows no substitution at all. Two further anomalies
vanished when their genes were rejected as mixed-family.

---

## 5. Defects found and fixed

Each of these was caught by a distribution looking wrong, **never by a failing validator**.
That is why per-source sanity summaries exist, and why `scripts/audit_sources.py` now sweeps
the invariants they violated. The live checklist is `TODO.md` `N3`.

### `#34` — The HMM set covered 9 families; UniPROBE spans 86 domain labels
Constructs from bZIP, HMG_box, GATA, T-box, Zn_clus, IRF and others were rejected as
`no_domain` because the HMM was missing, not because nothing was there. Fixed by scanning the
full Pfam-A library (30,134 families), with the admission policy evaluated over a DNA-binding
whitelist (`data/external/pfam/dbd_families.txt`, 41 families at the time, 65 after `T15`
below) so a DBD beside an unrelated
domain is still admitted. **340 → 468 domains, 7 → 31 families.**

### `#35` — UniPROBE labels predate several Pfam renames
`Homeobox`→`Homeodomain`, `Fork_head`→`Forkhead`, `MADS`→`SRF-TF`, `E2F_TDP`→`WHD_E2F_TDP`,
`Zn2Cys6`→`Zn_clus`, `BRIGHT`→`ARID`, `NR`→`zf-C4`, `AP-2`→`TF_AP-2`, `PBX`→`PBC`,
`RFX`→`RFX_DNA_binding`. Each resolved by Pfam accession lookup, never by string guessing.
Missing `Homeobox` alone would have excluded the largest family.

### `#36`, `#37` — The E-score column was read by position and was wrong for three accessions
`usecols=[0, 2]` holds for BAR15A/EMBO10/PNAS13/Cell08 but not for GR09 and SCI09 (headerless,
median intensity at index 2, E-score at index 3) or RAD13A (`8-mer / 8-mer / Median / Z-score`
— no E-score at all). The corpus briefly showed a 4:1 negative:positive ratio with `raw_score`
up to 776,106. The column is now identified by header name where present, else by being
bounded to [-0.5, 0.5] **and** taking negative values (p/q-values cannot). Files with no such
column, or several, are rejected with a reason. A validator guard makes the class of bug
impossible to miss: any `pbm_escore` outside [-0.5, 0.5] is a hard schema error.

### `#39`, `#40` — One gene folder is not one protein
The costliest defect so far. UniPROBE archives encode the construct at path level 2 and detail
pages can carry several insert sequences, but the parser grouped by gene and took the first
insert — **merging distinct engineered proteins as replicates**. ROG18A's `FoxJ3_N3/`, six
different chimeras, became one protein; NAR11's `HLH-1/`, a point-mutant series L13R/L13T/L13V,
became one; LIN14B's `ANAC092_DBD` and `ANAC092_FL` were averaged together; GD13's two CLAMP
constructs likewise. Fixed by keying inserts on their full construct name and splitting
experiments against those names. **ROG18A 7 → 15 domains, NAR11 1 → 4**, and it produced the
corpus's first two non-homeodomain variant clusters.

### `#42` — A bare `Insert sequence` label collapsed 88 GR09 genes into one cluster
Keying inserts by construct name (needed for `#39`) assumed the label always carries one.
GR09's reads simply "Insert sequence", so the key came out empty and every GR09 gene shared
it, producing a single `GR09:` cluster of 45 distinct Zn_clus domains. Caught by the
cluster-size distribution; the validator passed. Fixed by falling back to the gene name for
bare and plasmid-named labels, with a regression test.

### `#46` — Evidence the alignment places gaps correctly
After `#44`, ROG18A's 15 domains form 3 clusters (FoxJ3 with 7 variants, FoxN3 with 4, FoxN2
with 1) instead of mostly singletons. The edit positions are compositionally exact: `_6aa`
gives 6 consecutive positions, `_loop` gives 8, and `_6aa+loop` gives precisely their union.

### `#51` — PNAS08's archive is flat
Files sit at the top level with no gene folder, so the gene came out as a filename and nothing
matched. Fixed by deriving the gene from the filename stem when there is no directory.
**PNAS08 goes from 0 to 2 domains (65,792 rows).**

---

## 6. Sources and constructs excluded

### `#58` — `Path10`: tables truncated to the enriched end
Its files carry a valid E-score column but only 341-1,391 rows each, cut at E ≥ 0.25, against
the 32,896 of a full universal PBM. It therefore holds **no non-binding evidence at all**, and
applying our cutoffs would label its enriched 8-mers between 0.25 and 0.35 as non-binding —
the exact inverse of the truth. 18 *Plasmodium* AP2 proteins excluded; deregistered, reason in
`PROVENANCE.md`.

### `#59` — `GD09`: excluded twice over
Its tables hold 32,895 of the 32,896 non-redundant 8-mers (the palindrome `CCTTAAGG` is
missing), and the fully crossed design means a protein measured against a different 8-mer set
cannot join the others. Independently, its single protein NSY-7 — deposited twice, as `NSY-7`
and `Nsy-7` — returns no Pfam hit anywhere in its 294-residue construct. The completeness
error now distinguishes "truncated to the enriched end" (Path10, ~4% of the set) from
"near-complete but *n* short" (GD09, one short).

### `#52` — `PP15` must stay unmatched
Its archive holds both *Arabidopsis thaliana* and *A. lyrata* experiments, but its detail pages
carry only the *thaliana* sequence. A prefix match would have assigned the thaliana sequence to
the lyrata ortholog — a mis-attribution of exactly the kind the policy exists to prevent. Left
unparsed deliberately.

### `#50` — 22 archive entries are protein complexes
`Myc_Max`, `Kay_Jra`, `Da_Twi`, ten `HLH-2_*` heterodimers (Cell09), nine `CSL/NOTCH/MAML`
complexes (PO10). Two or three chains form one binding unit, so none can satisfy admission
condition 1. UniPROBE publishes no sequence for them, so they were already excluded — but by
accident. Now an explicit rejection category in [`DOMAIN_POLICY.md`](DOMAIN_POLICY.md).
This was the reasoning that made Tier 4 (bHLH dimers) doubtful; Tier 4 was dropped outright on 2026-08-17 (§1).

### `#53` — 10 genes have a detail page with no sequence
Cell08 4, MAR17A 5, GR09 1 — Hoxa3, Nkx3-1, Six6, Arid5b, Cebpa, E4f1, Xbp1 among them.
UniPROBE simply does not publish those. Reported per source by the parser.

### `#38` — 12 of 30 registered accessions contribute nothing
Eleven rejected entirely by the domain policy (C2H2 arrays, multi-domain, or families outside
the DBD whitelist); the rest for the data reasons above. All reported per source, none silent.
Attrition is tabulated in [`METHODS.md`](METHODS.md) §9.1.

---

## 7. Build performance

### 2026-08-14 — The build was profiled, then made ~6x faster without changing a byte
Was `TODO.md` `T6`. `build_dataset.py --all` took ~25 minutes. Measured on real source
builds, two things dominated and neither was doing useful work.

**The Pfam-A load was 100% fixed cost.** Scanning 3 sequences took 19.5 s; scanning 100 took
19.4 s. The library is 2.2 GB of *text*, reopened once per accession, so a 30-source build
spent about ten minutes re-parsing one unchanging file to do ~0.1 s of scanning.

**Validation re-checked constant columns 32,896 times per construct.** `dbd_seq` has 489
distinct values in the corpus but every check ran per row, and
`df.apply(_positions_in_range, axis=1)` built a pandas Series for each of 16.6 M rows to test
a handful of integers — 50 s of a 139 s corpus validation on its own.

Four changes, in descending order of payoff:

| # | change | effect |
|---|---|---|
| 1 | **Press the library** (`scripts/press_pfam.py`) | 19.5 s → 3.5 s per source |
| 2 | **Cache the loaded library per process** (`domains._library`) | → 0.5 s once, then ~0.6 s per source |
| 3 | **Validator axis checks on distinct values**, weighted back to row counts | validate 162.7 s → 27.0 s over the corpus |
| 4 | **Vectorized `add_pair_ids` payload** | 2.8x; it runs twice per build, in `coerce` and again in `validate` |

**Nothing about the dataset changed, and that was checked rather than assumed:**
- all 18 sources revalidated — **every report string, every stat and every `pair_id`
  identical** to the previous implementation;
- six sources fully rebuilt and compared against the committed Parquet — **byte-identical**;
- the pressed library checked against the unpressed one by moving the `.h3*` files aside:
  identical calls, 169 hits over 150 constructs;
- 8 regression tests added. They pin *behaviour*, not speed — specifically that the
  row-weighted counts in every validator message still equal what the per-row implementation
  produced, since that weighting is the one thing a future edit could silently drop.

**Pressing is a per-machine setup step**, not a repo artifact: the `.h3*` files are 2.6 GB of
git-ignored derived data beside the library, and re-downloading Pfam-A means pressing again.
`domains._library` says so out loud when it is handed a large unpressed library rather than
silently running six times slower.

**Parallelism was deliberately left out** and is `TODO.md` `T9`. The serial fixes were worth
more at a fraction of the risk: parallelism means restructuring `build_dataset.py` and
`reconcile_replicates`, and process pools make failures harder to attribute — which matters
for a build whose per-source rejection reasons are load-bearing. Threads were measured at
**1.0x** (the CSV parser holds the GIL); processes at **7.6x**.

---

## 8. Questions closed without action

| was | resolution |
|---|---|
| `#6` Are UniPROBE 8-mer files licence-gated? | No. Direct download via `filehandler.php`; no click-through, cookie or MTA. |
| `#11`, `#15`, `#16` Competing `dbd_seq` conventions; EMBO10/PNAS13 lacking a DBD field; the 719 aa `Jumeau` outlier | All dissolved by the domain policy. Every domain is `pfam_hmmer_padded`; `dbd_seq` length no longer correlates with source. Jumeau is 107 aa. |
| `#12` Noyes 2008 accession | Not in UniPROBE — it is B1H, not PBM. |
| `#13` Which panels supply forkhead / ETS / bZIP | Identified and all parsed: EMBO10, PNAS13, ROG18A, MAR17A, SCI09, GR09, SHO18A. |
| `#14` How should identical-domain pairs be split? | Reframed, not closed. The 11 same-source Cell08 pairs it described are no longer byte-identical — the 10 aa padding separated them (Evx1/Evx2 now differ at 12 of 77 positions, Lhx2/Lhx9 at 3). The real case is 17 **cross-source** pairs: `TODO.md` `T4`. |
| `#17` Panels add breadth, not depth | Confirmed and accepted — see `#24`. |
| `#21` Cost of the domain policy | Measured: 11.81M → 9.15M rows at the time. |
| `#28` Parse `GB11`? | Moot — its detail pages carry no sequence at all, so it is unparseable regardless. HLH reached 35 domains by other means. |
| `#31` SHO18A detail pages 500ing | Worked around via `detailsDef.php` with the site's own ids, verified against SHO18A/Apt. No ids guessed. |
| `#48` Did construct-splitting drop data anywhere? | No — verified clean across all 30 accessions, and the check is now part of `audit_sources.py`. |
| `#49` Does SCI09 lose anything to its rejected combined files? | No. All 104 genes keep at least one readable per-replicate file; the 105 rejected 20-column files are redundant. |

---

## 9. Source candidates screened and rejected

### 2026-08-14 — a literature sweep for PBM sources with designed protein variation
A separate research session surveyed the literature for PBM datasets carrying reference plus
mutant alleles of the same DBD. Screening its results against the admission policy and
against what is already parsed:

| candidate | verdict |
|---|---|
| **Kock et al. 2024, *Nat Commun* 15:3110** — 30 HD allelic series, 122 alleles | **acquired, screened, excluded 2026-08-18** — no E-scores published and no comparable scale; [`reports/kock2024_excluded.md`](../reports/kock2024_excluded.md). |
| Rogers et al. 2019, *Mol Cell* 74:245 | **already parsed as `ROG18A`** — 15 domains in 3 clusters, including the `N3(J3-6aa)` swap the sweep highlighted. |
| Liu et al. 2018, *eLife* 7:e34594 | **already parsed as `LIU18B`**, but only 2 of its constructs. The lead is the Dryad deposit, not the paper: `T14`. |
| Ibrahim et al. 2013, *Genome Res* 23:2091 — HOXD13 wt/Q325R/Q325K | **mostly held already.** `Q325K` is in BAR15A, which carries HOXD13 as 8 domains / 7 variants. Only `Q325R` would be new — one domain. |
| **Siggers et al. 2014, *Mol Cell* 55:640** — Msn2/Msn4/Com2/Usv1/Rgm1 | **rejected.** All are C2H2 zinc-finger arrays: separate ~23-residue folds on flexible linkers, failing admission condition 2. These are the exact grounds on which Phase 3 was dropped and ~8,000 domains excluded (§1). Reversible only via `domain.allow_repeat_arrays`. |
| "Mine Berger 2008 for pairs at ≤3 substitutions" | **measured, and smaller than proposed.** Over all 489 domains: 31 pairs at ≤3 edits, but 23 are cross-source and largely re-find the 17 `T4` duplicates plus their variant halos. Genuinely new same-source paralogues at ≤3 edits: 8 pairs. Now handled generally by `T12`. |
| Chu 2012; Noyes 2008; Aditham 2021; gcPBM paralogue panels | **out by assay**, as the sweep itself noted — B1H, microfluidic affinity, or an array design that does not yield the 32,896 8-mers. |

**Identifiers from that sweep are unverified** and are recorded in `TODO.md` as such. Rule 1
applies: none of them may be turned into a download URL by pattern, and each must be confirmed
against the source before use.

---

## 10. The CIS-BP gate, answered

### 2026-08-14 — CIS-BP publishes the assayed construct sequence
`TODO.md` `T1` rested on one question: does CIS-BP publish the sequence that was physically on
the array, or only its own DBD annotation of the full-length protein? Attributability — the
first admission condition — depends on the answer, and it decided whether CIS-BP was the
largest available expansion or a dead end.

**It publishes the construct.** Weirauch et al. 2014's own methods: *"All inserts were
sequence verified in full. Insert sequences and other information are available in Table
S6."* The project site describes Table S6 as *"DBD clone source material. This spreadsheet
provides information on the clone source material, the experimental construct sequences, and
the clone source contributors."*

Verified at every link rather than taken on the paper's word:

| link | contents |
|---|---|
| `TabS6_DBD_clone_information.xlsx`, sheet *Experimental constructs* | 1,032 rows with `Plasmid ID`, `Insert AA`, `#Flanking AAs`, `Backbone`, `Tag location` |
| GEO `GSE53348` | 2,064 samples, titled `pTH####_<HK\|ME>_8mer_<n>` — the plasmid ID is the join key |
| `GSM1291226` (`pTH1294_HK_8mer_593`) | `ID_REF / VALUE / E-Score / Z-Score`, **32,896 rows** |

The construct architecture is stated per row, so the confound `T13` exists to record is
published rather than inferred: 671 constructs carry 50 flanking residues, 96 carry 15, and
**265 carry none at all** — bare DBDs from an Agilent oligo pool, described in the methods as
*"These constructs did not contain any AAs flanking the DBD."*

**Measured yield through our own admission policy**, by scanning all 1,032 inserts and calling
each with `call_domain`: **731 admitted (71%)** — 240 rejected `no_domain`, 49 `repeat_array`,
12 `mixed_families` — giving 722 distinct sequences of which **677 are new**. That takes the
protein axis from 489 domains to roughly 1,166, across 106 species instead of 24, and fills
the families the corpus is thinnest in (Myb_DNA-binding 58 against our 2, AP2 39 against 3,
GATA 27 against 2, WRKY 22 against 1).

It also brings the depth the corpus lacks, because the paper's first selection strategy
deliberately populated nine DBD identity bins up to **90-99.99%**: within the new domains
alone there are 61 pairs within 3 edits and 103 within 5, led by Myb_DNA-binding — a family
with no variant depth here today.

**Nothing about the existing configuration has to change**, which is what separates this from
Kock: the same E-score statistic on the same scale, with the paper quoting `E > 0.45` as its
own significance cutoff — identical to `pbm.positive` — and HK/ME dual arrays per plasmid,
which `per_experiment: true` and `reconcile_replicates` already handle.

Two limits recorded honestly. The 240 `no_domain` rejections come from a stricter bar than
CIS-BP's own (gathering thresholds over full Pfam-A, against their 81 models at Eval < 0.01)
and should be spot-checked rather than assumed correct. And this is Weirauch 2014's own 1,032
constructs, **not** the ~2,294 TFs with PBM data the database aggregates: for aggregated
entries the construct sequence belongs to the contributing study, and much of that is UniPROBE
already held.

### 2026-08-14 — The whitelist is built from every source database's vocabulary, not one
`#34` fixed the HMM library but left the *whitelist* derived from a single database:
`build_dbd_family_list.py` read UniPROBE's `Domain` labels, so it could only ever name
families UniPROBE happened to publish. Screening CIS-BP exposed the same failure a second
time — of 240 constructs rejected as `no_domain`, **only 33 had no Pfam hit at all**; the
other 207 had a perfectly good hit belonging to a family we had no word for. Plant TCP, Dof,
NAC and SBP, Doublesex DM, the Rel homology domain of NF-kB.

The script now reads **both** UniPROBE's `Domain` labels and CIS-BP's `Pfam ID` column, each
being that database's own curation of what binds DNA, and each resolved against Pfam-A `NAME`
fields rather than trusted verbatim. **41 -> 65 families**, 24 added, none removed. The output
records which database vouched for each family.

Six CIS-BP labels predate Pfam renames and could not be resolved by string match. Rather than
guess, each was resolved **by the sequence**: scan the constructs CIS-BP gave that label to
against full Pfam-A and read off which family actually hits them. The method reproduces the
three existing hand mappings (`Homeobox`, `Fork_head`, `E2F_TDP`) exactly, which is why it was
trusted for the rest, and every result is confirmed by the target's own Pfam description —
`DUF573` -> `GeBP-like_DBD` ("DBD domain"), `RHD` -> `RHD_DNA_bind` and not `RHD_dimer`,
`EIN3` -> `EIN3_DNA-bd` and not `EIN3_N`.

Effect, measured: CIS-BP admission **731 -> 896**, and on sources already parsed **477 -> 490**
constructs with nothing lost — `LIN14B` alone goes from 1 admitted construct to 12, because
NAC is a family the corpus previously could not see. A rebuild follows (`TODO.md` `T16`).

Pinned by `tests/test_domains.py::test_whitelist_covers_the_families_recovered_from_cisbp`,
because this failure mode is silent: a missing family looks exactly like an absent domain.

### 2026-08-14 — Corpus rebuilt on the widened whitelist
`T16`. `build_dataset.py --all` re-run after `T15`; same 18 of 30 accessions contribute, all
validator-clean, and the twelve that yield nothing do so for their existing recorded reasons.

| | before | after |
|---|---:|---:|
| rows | 16,645,376 | **17,040,128** |
| distinct domains | 489 | **501** |
| clusters | 425 | **437** |
| Pfam families | 30 | **32** |
| binding / non-binding / excluded | 51,387 / 16,323,689 / 270,300 | 51,642 / 16,713,577 / 274,909 |
| admitted constructs | 568 (72.9%) | 581 (74.6%) |
| rejected `no_domain` | 55 | 42 |

Every added domain traces to the whitelist: `LIN14B` 1 → 12 constructs (11 `NAM`) and `SCI09`
68 → 69 (`MH1`). Only those two sources' reports changed, which is the cleanest available
evidence that nothing else moved.

The row-count invariant still holds exactly — **17,040,128 = 518 × 32,896**, where 518 is 501
distinct domains plus the 17 measured by two studies each. Every `dna_len` is 8. Cluster
structure is otherwise untouched: still 81 variants in 28 clusters, and the new NAC domains are
all singletons, so this bought breadth rather than depth.

The protein table was rebuilt with it: 501 domains, canonical UniProt for 462 (92%), domain
located in the full-length sequence for 354 (71%).

### 2026-08-14 — `weirauch2014` parsed; the corpus nearly triples
`T1` closed. 29,080,064 rows from 884 domains, validator-clean with no warnings, joined to
Table S6 on plasmid ID. The corpus goes from 501 domains to **1,364**, 32 families to **56**,
24 organisms to **138**, and 17,040,128 rows to **46,120,192**.

**What it bought, and what it did not.** Families holding ten or more domains went from 8 to
**27**, which turns leave-one-family-out from a handful of proteins into a real test and
largely settles the concern recorded as `#33`. Cross-source duplicate domains went from 17 to
38, so the label-agreement check of `T4` now has twice the evidence and, for the first time,
compares two different laboratories rather than two deposits from one. But **protein-axis depth
is untouched**: all 81 variants still sit in 28 clusters, and 1,293 of 1,321 clusters hold a
single domain. This was breadth, exactly as predicted before parsing.

**Two costs, both recorded rather than absorbed.** Full-length UniProt coverage fell from 71%
to 26%, because Table S6 publishes gene and species but no accession (`T2b`). And the source
carries three construct architectures — 671 with 50 flanking residues, 96 with 15, **265 with
none** — so its stored domains are systematically shorter than a UniPROBE construct of the
same protein. That is `T13`'s confound, now present in the data rather than anticipated.

Three pieces of structure came out of the work:

- **`parsers/_pbm.py`.** Everything that depends on the assay rather than the distributor —
  E-score identification, completeness, binarization, replicate combining, frame assembly —
  now has one home. `_uniprobe.py` keeps only UniPROBE's file layout and re-exports the rest,
  and five sources were re-parsed to confirm the extraction is byte-identical.
- **A latent bug, found by the extraction.** A header naming the E-score column used to
  return early, skipping the row-count check entirely, so a *headed* but truncated table would
  have been admitted silently; only headerless sources like `Path10` were ever caught. The
  check is now uniform. It matters here: GEO's tables are headed.
- **`build_protein_table.py` is no longer UniPROBE-shaped.** It derives its sources from the
  parser registry, which was right, but assumed every source publishes HTML detail pages —
  registering this parser made it raise `FileNotFoundError`. Construct lookup is now per
  source.

---

## 11. The canonical domain sequence

### 2026-08-17 — `dbd_seq` becomes construct-independent
Until now `dbd_seq` was the Pfam envelope plus 10 residues, **clipped wherever the assayed
construct happened to end**. So the same domain from two labs could be stored as two different
strings — VENTX is 76 aa from BAR15A and 57 aa from CIS-BP's zero-flank construct. That is an
artefact of cloning, and it reaches the model as if it were biology.

**Decided: the stored sequence is a project-internal canonical form.** Not "canonical" in
UniProt's sense — a definition local to this dataset, whose defining property is *stability*:
the same domain in gives the same string out, however truncated the input was. Length still
varies with the biology, because an insertion or deletion genuinely is a different sequence.

The window is **the Pfam envelope ± 10 residues**, and the padding stays: it is part of the
canonical definition rather than an addition to it, and now that it is taken from a reference
it is always fully available instead of surviving by luck. VSX1 G160D mutates 6 residues
before the Pfam start and currently survives only because BAR15A's construct happens to be
long enough.

Why the stored residues, not a fixed-width alignment: a match-state-only representation
**collapses VSX1 G160D onto its own reference** — identical sequence, different label — which
is the exact failure the padding was introduced to prevent (§2, `#18`).

### 2026-08-17 — A reference is consulted only where the construct falls short
Measured over 1,477 admitted constructs: **1,120 (76%) already cover envelope ± 10**, and for
those the window extracted from the construct is provably the string the full-length protein
would give. 37 canonical windows already come out byte-identical across sources
(`BAR15A:ARX_REF` = `Cell08:Arx`, 77 aa).

The remaining 357 need an external reference, and every one of them carries an identifier — 68
a UniProt accession, 289 a CIS-BP Gene ID. Resolution through UniProt cross-references was
verified working for Ensembl, FlyBase and Araport. Roughly 121 of the CIS-BP identifiers are
standard cross-references; the rest are assembly-specific and will partly fail.

**Expect to lose on the order of 10%.** That is accepted: a construct that cannot be placed is
discarded rather than force-fitted, which is the `PP15` rule (§6) — mis-attribution is worse
than absence.

The reference **extends**, it never replaces. The construct's own residues are kept, carrying
whatever mutations were engineered into them; the reference supplies only flanking residues the
construct is missing. So a shorter sequence caused by a real deletion or a true protein
terminus stays shorter — that is the thing that was assayed.

### 2026-08-17 — Placement is accepted at ≤5 edits outside the envelope
100% identity is impossible: the corpus is full of deliberate point mutants. Measured on the
462 constructs where a reference is already held, edits between construct and reference are
0 at the median and:

| ceiling | accepted | with ≥99% of the construct aligned |
|---|---:|---:|
| ≤1 | 94.6% | 92.4% |
| ≤2 | 96.3% | 94.2% |
| **≤5** | **97.6%** | **95.2%** |
| ≤10 | 97.8% | 95.2% |

Five sits at the knee — ten buys 0.2% more. The 11 it rejects are genuinely wrong, not
marginal: `Hoxc11` at 111 edits over 57% coverage, `Etv4` at 36% coverage.

**Edits are counted outside the Pfam envelope only.** Inside it, mutations are the subject of
the dataset and are kept unconditionally; outside is where residues are borrowed, so that is
what must match. A flat ceiling would reject ROG18A's chimeras, which differ from their parent
by 6-8 substitutions while being perfectly well placed. Mismatches and gaps are counted
separately: `Mlx` and `Rfx3` show 0-1 mismatches with 25-54 gaps, which is an alternative
isoform rather than a wrong protein.

### 2026-08-17 — Constructs differing outside the canonical window are discarded, not merged
Of 54 stored domains produced by more than one construct, 24 come from byte-identical
constructs and 26 differ only in how much flank each lab cloned — both benign. **Four carry a
genuine substitution outside the stored window** and were being silently reconciled as
replicates of one protein:

| substitutions outside | stored | constructs |
|---:|---:|---|
| 4 | 76 aa | `Cell09:HLH-25` / `Cell09:HLH-27` — two distinct *C. elegans* genes |
| 2 | 70 aa | `weirauch2014:pTH8163` / `pTH9718` |
| 1 | 106 aa | `MAR17A:Foxc1` / `BAR15A:FOXC1_REF` / `weirauch2014:pTH2673` |
| 1 | 77 aa | `Cell08:Msx3` / `MAR17A:Msx3` |

**Owner's decision: bin them.** A difference the model cannot see, attached to measurements
that may differ because of it, is a confound — the same reasoning as admission condition 3,
applied symmetrically. Until now condition 3 was enforced only in `bar15a.py`; in the panel
and CIS-BP parsers it was vacuous, because `mut_positions` is computed *from* the stored
sequences and so nothing could fall outside by construction. This makes it a real check
everywhere.

### 2026-08-17 — Clusters are built on canonical domains, never on reference proteins
Tempting, since references are construct-independent and a cluster representative need not be
a dataset row. But DNA-binding domains are conserved while the rest of the protein diverges.
Among pairs whose domains are within 5 edits, median domain identity is **96%** while median
full-length identity is **68%**, and the extremes are stark: **Irx3 and Irx4 have 94% identical
domains and 32% identical proteins.** Any protein-level clustering separates them, after which
a model trained on Irx3 predicts Irx4 for free — precisely the leakage clusters exist to
prevent. Also Hoxa10/Hoxd10 at 45%, Vax1/VAX2 at 46%.

The reference's job is to make the domain canonical and to establish protein identity for
deduplication. It is not the clustering unit.

### 2026-08-17 — CD-HIT greedy incremental is the clustering algorithm
Not invented here: the published greedy incremental algorithm (Li & Godzik 2006; MMseqs2
`--cluster-mode 2`). Sort by decreasing length, longest becomes the representative, and each
remaining sequence is compared **only to representatives**. Three properties earn it the job:

* the longest member becomes the representative, so the least-clipped form is the reference;
* comparison is never transitive, which kills chaining — single-linkage produced 35-domain
  blobs at 5 edits, against 8 at 1 edit;
* every member is within *k* of the representative, which is the invariant `mut_positions`
  needs (`#45`).

Implemented on `snp2prot.align` rather than by installing MMseqs2: the aligner already has the
right semantics, including free terminal gaps, and a full pairwise scan of the corpus takes
16 seconds. Adding a C++ dependency to a 16-second job is not warranted; the algorithm is
reused, which is the part that matters.

### 2026-08-17 — Built, and the corpus rebuilt on it
`snp2prot.canonical`, `snp2prot.references` and `snp2prot.clusters`, wired into all three
parsers, with `scripts/build_clusters.py` assigning `wt_id` corpus-wide.

| | before | after |
|---|---:|---:|
| rows | 46,120,192 | **45,495,168** |
| distinct domains | 1,364 | **1,335** |
| clusters | 1,321 | **1,133** |
| clusters holding >1 domain | 28 | **122** |
| variants | 81 | **202** |

Protein-axis depth roughly tripled, and stopped being a homeodomain story: 87 of 202 variants
are homeodomain against 54 of 81 before, with Myb 16, bHLH 15, forkhead 13, AP2 9, zf-C4 7,
SAND 5, RFX 5, TCP 5. That is the corpus finally able to ask whether sensitivity to a single
residue transfers between folds.

**Engineered constructs keep a short window** (owner, 2026-08-17). ROG18A's FoxJ3/FoxN3
chimeras are synthetic and cloned as bare domains: no natural protein exists to extend them
from, and their termini are not truncations. Discarding them cost the corpus its only
non-homeodomain variant depth, which was the wrong trade. A short entry that duplicates a
fully-canonical copy of the same domain is dropped at clustering instead, so no two rows
describe one domain at two lengths.

**Clustering threshold: 5 edits** (owner). 1,133 clusters, largest 8 — against the 35-domain
blob single-linkage produced at the same threshold, which is CD-HIT's non-transitivity earning
its place.

Four defects were found and fixed while doing this, each recorded because each was silent:

* **Stale interim output.** `build_dataset.py` left a source's previous Parquet on disk when
  that source stopped yielding, and it was still being read as part of the corpus. Two
  accessions were affected. The build now deletes it.
* **Cluster ids collided.** Naming a cluster after its representative's lineage `wt_id` is
  readable but not unique — `ROG18A:FoxJ3` names six different representatives, so six
  clusters collapsed into one, undoing the clustering for exactly the chimeras that had just
  been rescued. Collisions now take a suffix.
* **The validator assumed clusters are source-local.** They are not: a `Cell08` variant
  routinely belongs to a cluster represented by a `BAR15A` domain, which a per-source
  validator cannot see. It reported every such cluster as a missing reference. Now a warning,
  with the bounds check done corpus-wide — where it passes for all 202 variants.
* **An O(rows) alignment.** `build_clusters.py` computed an edit profile per row rather than
  per distinct domain: 45 million alignments instead of 1,335.

**Known gap:** the protein-side table covers 1,224 of 1,335 domains, because it locates a
domain by substring within its construct and a reference-extended canonical sequence is no
longer a substring of it. Tracked as `T18`.

---

## 12. The modelling audit — 2026-08-27

A defect-hunting review of everything under `src/snp2prot/` that touches the model, the
baseline or the metrics. Eighteen findings; all closed. **Every number in
[`reports/training.md`](../reports/training.md),
[`reports/nn_baseline.md`](../reports/nn_baseline.md) and
[`docs/ML_RESULTS.md`](ML_RESULTS.md) is superseded by the changes below and none of them have
been rerun yet.** Read those three against this section until they are.

### 12.1 The comparison was not like for like

**`D7` — the nearest-neighbour baseline now copies binary calls, not E-scores.** It copied the
neighbour's continuous profile, and the model trains on thresholded labels, so the two were
never compared on equal information. Measured, same neighbour and same overlap guard:

| fold | E-score copy | binary copy | of the baseline's score |
|---|---:|---:|---:|
| `S1/fold-0` | 0.7863 | 0.4721 | 40% was the ordering |
| `S2/fold-0` | 0.3912 | 0.1452 | 63% |
| `P1/Homeodomain` | 0.0244 | 0.0057 | 76% |
| `P3/all` | 0.9278 | 0.6603 | 29% |

Between 29% and 76% of the old baseline came from the ordering *within* the copied profile —
information the model is never given. **The owner's reading settled it**: a universal-PBM
E-score is a rank-enrichment statistic against background, read by the field at a cutoff and
stored here at 0.45 / 0.35, not a graded affinity. Its ordering is not a quantity the assay
reports, which is why the dataset stores a binary label at all. The continuous form was
therefore not a harder bar but a different question, and it is gone rather than kept as a
second column. **The whole PBM corpus is binary: training, baseline and evaluation alike.**

Consequence to expect on the rerun: the model lost to the baseline on 33 of 38 runs against the
continuous form. Against the matched bar it should win on most. That is a fairer comparison,
not a softer one.

**`D8` — Spearman is deleted as a model metric**, for the same reason: it scored a prediction
against the raw E-score ordering. `metrics.spearman` survives as a utility for
`scripts/check_pooling.py`, which correlates embedding displacement against edit count.

**The C1 training-pool exclusion is removed.** It was applied to the model and never to the
baseline, so in 18 of 19 folds the baseline could copy from 10-29 domains the model had been
denied — 358 domain-slots across the grid. It also bought nothing: C1 is scored in exactly one
place, `run_grid`'s `P3/all` section, and under `P3/all` every variant is held out by
construction. **If a C1 number is ever wanted from another regime the mask has to come back.**

**The C1 cut moved 0.7 → 0.25, and the set 29 → 41 variants.** Not because any variant became
harder: under `P3/all` the lookup makes the *same* prediction either way — "this variant binds
what its wild type binds" — and copying the continuous profile simply earned extra credit for
*ordering* the wild type's sites, which the binary copy does not. The same prediction now scores
median 0.776 where it scored 1.000, so a cut of 0.7 had drifted onto the median and was selecting
83 of 173 variants. The criterion that chose 0.7 is unchanged — put the cut where the
distribution is sparse — and on the new scale that is 0.25, the single-variant bin between
0.20 and 0.25. The 11 variants the old cut selected now score 0.000–0.330, every one still in the
tail; the two orderings agree at rank correlation 0.614. The set is larger than 29 because binary
scoring exposes variants whose positive *set* differs from the wild type's, which the ordering
credit used to partly rescue — and those are exactly the variants C1 is about.

### 12.2 A metric for the domains AUPR cannot reach

**`D9` — `metrics.suppression`, and the null anchor is not a threshold.** 20 held-out records
have no positive 8-mer — every one a `dead_variant` — so AUPR, AUROC and R@P0.5 are undefined
for them, and they are the sharpest evidence the corpus holds for **C1**.

Counting "predicted positives" needs a decision threshold. The null anchor was designed as one
and **does not work as one**: measured on a trained checkpoint it sits at −24.3 where the
negatives' median is −26.4 and the positives' median is +9.6, with a median of **15,159 of
32,460 8-mers scoring above it**. The reason is structural — for the ~1,318 domains that have
positives it is one more negative being pushed down, and the only force raising it comes from
the 20 that have none. It loses 1,318 to 20. Its first job is real and it stays for that: a
domain with no positives gets a target, so it contributes a gradient instead of nothing.

So the metric compares the variant with its own wild type down the **protein axis**, needing no
threshold and no continuous value: *of the sites the wild type binds, the fraction the model
ranks lower in the variant*. Ranks, not scores — comparing scores would hand 1.0 to a model
that merely scores the variant lower everywhere. **1.0** it saw the mutation abolish binding ·
**0.5** the sites moved at random · **0.0** it predicts the variant exactly like the wild type,
which is what the lookup does by construction. Scored only where the wild type stayed in
training, and never on `no_evidence` (`T21`).

First measurement, `A1` on `P3/all`: **model 0.4921, baseline 0.0000, chance 0.5000.** The
model is at chance on the C1 claim.

### 12.3 The ceilings were not ceilings

**`D10` — `ridge_probe` and `kernel_probe` deleted; `rank_ceiling` kept and reframed.** The
ridge probe was documented as *"strictly upper-bounds any model whose protein side is a linear
map"* and was **violated on 5 of the 8 rows of its own report, by up to 68%**. It minimises
squared error on E-scores and is scored by AUPR, two objectives that demonstrably disagree — on
`S2/fold-0` the MSE-optimal ridge scored 0.4324 where the AUPR-optimal one scored 0.4736.

The general rule, which is what decided the section:

> Every probe of this kind is a **lower** bound on its model class — one estimator's score, when
> the class can always contain a better member. A lower bound licenses a conclusion only when
> the number comes out **high**: it can rule a component *out* as the constraint, never *in*.

`ridge_probe` and `kernel_probe` came out low in §4's use, so they licensed nothing.
`rank_ceiling` comes out high (0.9167 at `D = 256` against 0.275 mean on `S2`), so "the shared
space is not the binding constraint" holds. The conclusion the deleted probes were used for —
*would a stronger or non-linear tower help* — is now answered by an ablation instead (`T36`).

The two-tower model at `protein_hidden = 0` **is** a linear protein map fitted with a ranking
loss, and on `P1` it scored 0.0588 where the ridge scored 0.0286. The objective was the lever,
not the tower.

### 12.4 Training

**`D11` — the step budget is 15,000, not 6,000, and early stopping is off.** The old value
rested on *"the training loss keeps falling to ~0.01, so everything past that is fitting the
training proteins"*, which is exactly backwards. Measured to 15,000 steps with stopping
disabled, comparing windows of three evaluations *within* a single run:

| fold | mean test AUPR @5.5–6.5k | @13.5–15k | gain |
|---|---:|---:|---:|
| `P3/all` | 0.8563 | 0.8798 | +0.024 |
| `S2/fold-3` | 0.2275 | 0.2457 | +0.018 |
| `S1/fold-2` | 0.6584 | 0.6863 | +0.028 |

All three still climbing at 15,000. The training loss does reach ~0.005 — the training proteins
are memorised by about step 6,000 — but held-out AUPR keeps improving for another 9,000 steps.
**Memorising the training set and generalising to unseen proteins are not coupled here.** Cosine
LR decay was measured as the cheaper alternative and is slightly *worse* at equal budget (0.8438
against 0.8458 on `P3/all`), so the budget is the lever, not the schedule.

`patience: 8` fired on 5 of 38 runs and all five were `S2` — `A1 S2/fold-4` stopped at step
2,750 with its best at 750, while `S2/fold-3` was still gaining at 12,000. It is now `0`.
**Both models are kept per run** — the validation-selected checkpoint and the model at the
budget — and scored side by side, so the selection loss is a reported number per fold rather
than an assumption.

**`D12` — the protein tower LayerNorms its input.** ESM-2's pooled vectors have norms 4.83–9.85
and ESM-DBP's 0.76–1.24, so with a bias on the projection `‖b‖/‖Wx‖` was **0.13 for `A1` and
1.16 for `A4`**: for one arm the bias was a correction, for the other it outweighed the signal
and set the output direction. The cost was not a handicap but an erasure — `A4`'s embeddings are
the better separated (mean pairwise cosine 0.752 against `A1`'s 0.874), yet the untrained tower
emitted **0.8930 for both arms, identical to four decimals**. With the LayerNorm they come
through at 0.891 and 0.759. L2-normalising the input instead is a trap: it sets `‖Wx‖ ≈ 0.26`
beside `‖b‖ = 0.251` and reproduces the same collapse for *both* arms.

This did not rescue `A4` on `P1`, which was the hypothesis: it moved 0.0035 → 0.0044 against a
random-ranking 95th percentile of 0.0036, with AUROC 0.551 and no meaningful signal. `A4` is
better at in-family discrimination and worse out of family, consistent with its embeddings being
more family-clustered (within-Homeodomain cosine 0.963 against 0.731 across, versus `A1`'s 0.933
and 0.874).

**`D13` — the temperature clamp moved off the forward pass, and weight decay off the
non-weights.** `logit_scale` and `null` were being decayed with the convolution weights; neither
is a weight. And `score()` clamped in the forward pass, which makes the gradient exactly 0 above
the ceiling: with decay the parameter oscillated across the boundary, without it froze at
4.60882 to five decimals with no force acting on it at all. `clamp_temperature()` now pins the
value after each step, CLIP's own arrangement, and the gradient stays live.

**None of this changes a reported number, and that is the point.** `score` is `scale × cosine`
and every metric ranks *within* one domain, so no positive scalar can reorder anything —
verified to 8 decimals, and end to end 0.8558 clamped at 100 against 0.8554 running free to 148.
What the temperature does set is *what the model learns from*: at the converged scale a median
of **2 of 32,460 negatives carry half the negative-side gradient**, against ~1,500 at the
initial temperature. That concentration is not evidence being discarded — a negative's gradient
share is proportional to how wrongly close it still sits, so the rest are satisfied constraints
— and the owner's decision is to keep the temperature low deliberately, for a space where
non-binding is pushed far away.

### 12.5 Things that were true and unrecorded

* **`recall_at_precision` returned `nan` for a ranking that never reached the target**, and
  `macro_average` dropped those domains. The published `P1` baseline recall was **0.130 over 24
  of 412 positive-bearing domains**; over all of them it is **0.008**. A failure is not an
  absence: at depth 1 the precision is 1.0 whenever the top 8-mer binds, so never reaching 0.5
  means the top hit is wrong and nothing recovers. It now scores 0, and `macro_average` reports
  `n_scored_<metric>` for every metric so no figure can again be a mean over an unstated subset.
* **Nothing showed chance level.** A per-protein AUPR is anchored to the protein's own positive
  rate, which ranges 0.0015–0.0034 across the folds. `metrics.chance_aupr` and
  `metrics.random_baseline` (200 sampled draws, giving a 95th percentile) are now reported with
  every fold, and a result at or below that percentile is flagged **at chance**. On the old grid
  that caught exactly one row: `A4` on `P1`, at 0.0035 against a null of 0.0035.
* **`nn_lookup` ties.** The documented rule — ties break to the lowest training row — was untrue
  in **51% of the 358 tied cases** because `np.argpartition` promises no order among equal
  values. Worth up to 0.008 AUPR and liable to change with the numpy version. Now a stable sort.
* **`OVERSHOOT` landed in the same commit as the report it fixed.** `S2/fold-1` carved 38% of
  its pool to validation where every sibling carved 15%; the guard fixed it and the grid was
  never rerun, so two rows of `reports/training.md` came from a `validation_split` that no
  longer exists. Runs now record `code.commit` and `code.dirty`, and the report says when its
  rows were not all produced by the same tree.
* **Trained weights were thrown away.** A grid was 38 numbers and no models. Both checkpoints
  per run now go to `data/processed/checkpoints/` and to MLflow.
* **`model.seed` split out of `splits.seed`.** Which domains are held out and how the towers are
  initialised are unrelated choices; one number was doing both. `run_grid.py --seeds N` varies
  only the second, default 1. Nothing has measured the spread yet — `T37`.
* **The alignment guarantee was documented where it was not enforced.** `KmerMatrix`'s docstring
  claimed the row order was "checked on load"; it was not, and `Trainer` indexed the matrix and
  the embeddings together on the strength of a convention neither asserted.
  `corpus.require_aligned` now enforces it in the library and in every script — comparing
  sequences, since a reordered corpus passes both a length check and a set check.
* **`snp2prot.training` had no tests.** It is the module that decides which domains a model
  sees and which checkpoint is kept. It now has 15, and the suite went 262 → 320.
