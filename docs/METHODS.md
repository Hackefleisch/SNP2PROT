# Methods

Construction of a paired transcription-factor DNA-binding domain × DNA-site dataset with
binary binding labels.

This document records what was done and why. It is written to be sufficient for an
independent reimplementation. Parameters quoted here are the values in
[`configs/thresholds.yaml`](../configs/thresholds.yaml); statistics are from the most recent
complete build (see §10).

---

## Glossary

Terms are used throughout in the following senses.

| term | meaning |
|---|---|
| **transcription factor (TF)** | A protein that binds specific DNA sequences and regulates transcription of nearby genes. |
| **domain** | A part of a protein that folds into a stable shape on its own and carries one function. Proteins are built from domains the way a program is built from functions: they recur across otherwise unrelated proteins, and one can often be studied in isolation. A transcription factor typically has a DNA-binding domain plus separate domains for other roles. Throughout this document "domain" means one stored `dbd_seq` — one protein entry in the dataset. |
| **DNA-binding domain (DBD)** | The compact, independently folding part of a TF that contacts DNA. The rest of the protein does other work and is not represented here. Stored as `dbd_seq`. |
| **construct** (or *insert*) | The protein fragment actually cloned and put on the array. It is usually a DBD plus some flanking residues, and is not the full-length protein. Domains are cut from the construct, so the construct bounds what can be stored. |
| **family** | A group of DBDs sharing a fold and a recognition mechanism — homeodomain, forkhead, bZIP, C2H2 zinc finger. Assigned here from Pfam, and stored as `dbd_family`. Families are the coarsest unit at which generalization can be tested. |
| **Pfam / HMM** | Pfam is a library of protein-family models; each family is a profile hidden Markov model that scores how well a stretch of sequence matches that family. Used here to locate domains and assign families. |
| **envelope** | The span of a construct that a Pfam model matches, in residue coordinates. `dbd_seq` is this envelope plus padding (§3.3). |
| **reference** (or *wild type*) | The unmodified version of a domain, against which its variants are compared. `n_mut_from_wt = 0`. |
| **variant** | A version of a reference domain carrying one or more sequence changes — a point substitution, or an engineered chimera. |
| **cluster** | A reference domain together with all of its variants, keyed by `wt_id`. The unit at which train/test partitions must be made: splitting two variants of one protein across the partition does not test generalization. A cluster of size 1 is a reference with no variants. |
| **paralogue** | Two genes descended from a duplication within one organism. Their DBDs are sometimes byte-identical while the rest of the protein has diverged (§6). |
| **site** (or *8-mer*) | The DNA sequence whose binding is measured — here always 8 base pairs. Stored as `dna_seq`. |
| **PBM** | Protein-binding microarray. A single experiment scores one protein against every possible DNA sequence on the array at once, giving a complete profile rather than a handful of measurements. |
| **E-score** | The PBM enrichment statistic: a rank-based measure bounded to [−0.5, 0.5], comparable across experiments, high for bound sequences. Thresholded to produce the binary label. |
| **replicate** | A repeat measurement of the same protein, possibly on a different array design. Reconciled rather than averaged (§6). |

## 1. Design objective and its consequences

The dataset is intended to support models that are sensitive to variation on the protein and
the DNA axis simultaneously. Two requirements follow, and both constrain the construction.

**The protein-side unit must be attributable.** A measurement is only usable if the stored
protein sequence is the thing that produced it. This is stricter than "the protein was a
transcription factor", and it excludes several classes of experiment that are otherwise
perfectly sound (§3.2).

**The DNA axis must be uniform across sources.** Binding-site length is confounded with assay
identity, and assay identity is confounded with the positive rate. Pooling sources of
differing site length therefore leaks the label prior. The dataset is restricted to universal
protein-binding microarray (PBM) experiments, in which every protein is scored against the
same complete set of 32,896 non-redundant 8-mers, giving a fully crossed design with no
length variation at all.

## 2. Data source and acquisition

All primary data were obtained from UniPROBE (Hume et al. 2015, *NAR* 43:D117;
<https://thebrain.bwh.harvard.edu/uniprobe/>). The download index lists 36 accessions. For
each, the contiguous 8-mer archive `<ACC>_contig8mers.zip` was retrieved, together with one
HTML detail page per gene, which carries the assayed clone insert sequence, the Pfam domain
label, the UniProt accession and the source organism.

Every retrieved file is recorded in [`PROVENANCE.md`](../PROVENANCE.md) with its URL,
accession, UTC retrieval time, byte size and SHA-256. Detail pages are recorded per accession
through a `_manifest.tsv` carrying the same fields per file. Raw files are never modified;
all parsing writes to a separate tree.

Three classes of source were not usable and are recorded as such rather than substituted:

- two accessions advertise an archive the server does not hold;
- three publish detail pages carrying no sequence, so the assayed construct is unknown;
- one publishes 8-mer tables truncated to the enriched end (341-1,391 rows cut at E ≥ 0.25,
  against 32,896 for a complete experiment), which carry no non-binding evidence and would
  invert the labels of the protein's strongest sites if thresholded.

The resulting attrition is quantified in §9.1.

Files are additionally required to contain the complete non-redundant 8-mer set, since the
design is fully crossed and a partial table cannot contribute negatives.

Supporting reference data: the complete Pfam-A HMM library (release of 2026-01-22, 30,134
families) and canonical UniProt sequences, both cached under `data/external/` and recorded in
`PROVENANCE.md`.

## 3. Protein-side annotation

### 3.1 Domain identification

Each **assayed construct** — not the full-length protein — was scanned against the complete
Pfam-A library using HMMER3 through pyhmmer 0.12.1, with each family's curated gathering
threshold as the inclusion criterion.

Scanning the construct rather than the protein is deliberate. UniPROBE's own domain label
describes the full-length protein, and for constructs carrying only part of a protein the two
disagree: constructs annotated `Homeobox, POU` were found to contain both domains within the
assayed insert, which the protein-level label does not establish.

Because a whole-library scan also returns domains unrelated to DNA binding, the admission
policy is evaluated over a restricted set of **65 DNA-binding families**. A family absent from
that set is indistinguishable from an absent domain, so its derivation matters: it is taken
from the source databases' own curation of what binds DNA — UniPROBE's `Domain` labels and
CIS-BP's per-construct `Pfam ID` — each resolved against Pfam-A `NAME` fields rather than
trusted verbatim.

Deriving it from a single database proved insufficient twice. The first set covered only the
nine families present in the earliest sources; the second covered 41, all of them families
UniPROBE happens to publish, which silently discarded constructs carrying plant TCP, Dof, NAC
and SBP domains, Doublesex `DM`, and the Rel homology domain of NF-κB. Adding CIS-BP's
vocabulary took the set to 65 and recovered 13 UniPROBE constructs that had been rejected as
`no_domain`, eleven of them NAC.

Both vocabularies are stale in places. Labels predating a Pfam rename are resolved by
accession lookup where the accession is published (`Homeobox` → `Homeodomain`, `Fork_head` →
`Forkhead`, `MADS` → `SRF-TF`, `E2F_TDP` → `WHD_E2F_TDP`, `Zn2Cys6` → `Zn_clus`, `BRIGHT` →
`ARID`, `AP-2` → `TF_AP-2`, `PBX` → `PBC`, `RFX` → `RFX_DNA_binding`), and otherwise **by the
sequence**: the constructs carrying a given label are scanned against the full library and the
family that actually hits them is read off (`CXC` → `TCR`, `DUF260` → `LOB`, `DUF573` →
`GeBP-like_DBD`, `EIN3` → `EIN3_DNA-bd`, `RHD` → `RHD_DNA_bind`, `zf-Dof` → `Zn_ribbon_Dof`).
That procedure reproduces the accession-derived mappings exactly where both are available,
which is the basis for trusting it where only one is.

Pfam hits outside this set are recorded but do not affect admission: a DNA-binding domain
adjacent to an unrelated domain remains solely responsible for the DNA interaction.

### 3.2 Admission policy

A construct enters the dataset only if its stored domain satisfies three conditions.

1. **Sole responsibility.** The stored subunit alone produced the measured interaction.
2. **One continuous region.** The subunit is a single contiguous stretch of sequence.
3. **Containment of variation.** Every sequence change that alters binding lies inside the
   stored region.

The conditions follow from the intended use. A model that predicts a structure, locates the
centre of the DNA-contacting residues and embeds the residues around it requires a subunit
that folds as a unit, whose spatial arrangement is determined, and within which the
experimental contrast is visible.

Four rejection categories follow:

| category | condition | description |
|---|---|---|
| mixed families | 1, 2 | Pfam hits from two DNA-binding families in one construct (paired + homeodomain; POU-specific + homeodomain). Both domains contact DNA; the data do not attribute a given 8-mer's score to either, and their relative geometry is undetermined in the absence of DNA. |
| repeated arrays | 2 | Multiple hits of one family, characteristically C2H2 zinc-finger arrays of two to six ~23-residue folds joined by flexible linkers, each requiring its own Zn²⁺. Not one continuous unit. |
| protein complexes | 1 | Two or three distinct polypeptide chains forming one binding unit (bHLH heterodimers; Fos/Jun; CSL–NOTCH–MAML ternary complexes). No single sequence is responsible. |
| no domain | — | No Pfam hit above the gathering threshold anywhere in the construct. |

Condition 3 is enforced per variant: a substitution falling outside the stored window is
rejected, because the variant would otherwise be sequence-identical to its reference while
carrying a different label.

Rejections are reported per source and never repaired.

### 3.3 Domain boundaries

The stored sequence, `dbd_seq`, is the Pfam envelope **extended by 10 residues on each
side**, clipped where the construct ends.

The padding is required by condition 3. Pfam models begin at the structural core; for
homeodomains this excludes the N-terminal arm, which inserts into the DNA minor groove and
makes base contacts. Empirically, four of 91 point substitutions in one panel fall outside the
bare envelope, three of them 2–9 residues N-terminal of the Pfam start. Trimmed to the bare
envelope, those variants become identical to their own references. A 10-residue margin covers
all observed cases and is set from the position of the N-terminal arm rather than fitted to
them. The condition is re-checked for every variant at parse time.

A companion protein table stores each domain at four levels — bare Pfam envelope, padded
domain, assayed construct and full-length UniProt sequence, with offsets — so that
representation can be varied downstream without re-parsing.

## 4. Assignment of sequences to measurements

UniPROBE archives encode the assayed construct as a path element, and a single gene directory
may contain several distinct proteins: engineered chimeras, point-mutant series, or separate
domain and full-length constructs. Detail pages correspondingly carry one insert sequence per
construct.

Experiments are therefore matched to constructs by name, longest match first, rather than
grouped by gene. Constructs are keyed by the full name the page gives. Where a page labels its
sequence by plasmid or gives no name, the gene name is used.

Two situations are resolved conservatively. Where a gene's experiments match no construct on
its own page, the data are not assigned. Where an archive contains experiments for orthologues
of two species but the page carries only one species' sequence, no assignment is made, since a
name-based match would attach one orthologue's sequence to the other's measurements.

## 5. Label assignment

PBM E-scores are rank-based enrichment statistics bounded to [−0.5, 0.5] by construction.
Labels are assigned per 8-mer: E ≥ 0.45 → 1 (binding), E ≤ 0.35 → 0 (non-binding), and
intermediate values → −1, retained in the table and excluded from training.

Thresholds are applied **per experiment**, never to pooled scores, because array designs
differ in intensity scale. The applied cutoffs are written into every row, so a table is
self-describing.

Negatives from a PBM are assayed-and-unbound rather than merely unobserved, and are recorded
as such in a mandatory `neg_provenance` field, so that they remain distinguishable from
weaker forms of negative evidence should other assay types be added.

The E-score column is identified by content, not position: column order, column count and the
presence of a header row all vary between accessions, and one accession publishes no E-score
at all. The column is taken from a header name where one is given, and otherwise as the unique
numeric column bounded to [−0.5, 0.5] that also takes negative values — a criterion that
distinguishes an E-score from intensities, z-scores and p-values. Files in which no such
column exists, or in which several do, are rejected with a reason.

## 6. Replicate reconciliation

Where an accession provides several experiments for one construct, each is binarized
independently and the labels are then combined: agreement is retained, disagreement is
assigned to the excluded band. Scores are never averaged before thresholding, since replicates
may sit on different array designs whose intensity scales differ by an order of magnitude.

Distinct genes whose stored domains are byte-identical are treated the same way. Such cases
arise from recent gene duplication, where the domain is invariant while the flanking sequence
has diverged. Their agreement is comparable to that of technical replicates of a single gene
(median positive-set Jaccard 70.4% against 72.2%), supporting the treatment.

### 6.1 Cross-source replicates and the label noise floor

Reconciliation above operates within a source. Forty-seven domains are stored by more than one
distributor, giving 49 pairs in which two laboratories assayed the same padded domain against
the same 32,896 8-mers. These are not reconciled — each source's rows are retained as
deposited — but they are measured, because the rate at which two independent measurements of
one protein disagree is the noise floor against which any model's accuracy has to be read.

Pooled hard-label agreement is 99.98% (300 of 1,572,985 calls). That figure is uninformative:
positives are under 0.5% of an 8-mer table, so agreeing on the negatives dominates it. On the
calls the dataset exists to make, agreement is far weaker — the median pair agrees on 45.8% of
the 8-mers either of them called positive, the pooled figure is 49.2%, and the median Spearman
correlation between the two laboratories' E-scores is 0.61.

The comparison with §6 is the point: within a source, byte-identical domains and technical
replicates agree at 70.4% and 72.2% positive-set Jaccard. Across sources the same quantity is
45.8%. Agreement falls by some 25 points when the two measurements come from different
laboratories rather than from one, so a substantial share of the disagreement tracks
laboratory and array design rather than the protein — and a model that reproduces held-out
positives much beyond this level is reproducing a source, not a binding preference.

Nine pairs fall below a positive-set Jaccard of 0.20 or a rank correlation of 0.30 and are
flagged rather than corrected; one — the *Arabidopsis* NAC domain ANAC092, stored by both
Lindemose et al. 2014 and CIS-BP — shares no positive call at all between the two deposits
(Jaccard 0.000, rho 0.07) despite an identical stored domain. Per-pair figures are regenerated
into `reports/overlap.md`.

## 7. Cluster construction

A cluster (`wt_id`) is a reference domain together with its variants, and is the unit at which
train/test partitions must be made: a partition that separates two variants of one protein
does not test generalization.

Distance between a variant and its reference is computed from a global pairwise alignment
(BLOSUM62, affine gaps), so variants may carry insertions and deletions. **Terminal gaps are
free.** This is required rather than convenient: `dbd_seq` is the padded envelope clipped by
the construct, so two constructs of the same domain differ at the termini by however much
padding fitted. In one panel the bare domains are uniformly 83–85 residues while the padded
ones range 95–105; charging for those terminal gaps would report indels between a protein and
its own chimera. Internal indels are counted.

`n_mut_from_wt` is the number of edited residues and `mut_positions` lists their positions in
the **reference's** coordinate frame, one entry per edited residue, so that all members of a
cluster share one coordinate system.

## 8. Schema and validation

Every parser emits an identical 22-column row schema, and every table passes a validator
before it is written. The validator enforces, among other checks: unambiguous nucleotides and
a maximum site length; a mandatory negative-provenance field on every negative; agreement
between the recorded mutation count and positions; mutation positions bounded by the cluster's
reference; deterministic row identity and its uniqueness; and E-scores within their
theoretical bounds.

Row identity is a hash of the domain sequence, site sequence, assay and stringency, computed
by the schema rather than supplied by parsers.

A separate audit sweeps all accessions and checks properties that a per-table validator cannot
see: that every experiment file in an archive is claimed by exactly one construct, that every
construct is emitted or rejected with a reason, that no two constructs share a key, that every
domain is measured against the complete 8-mer set exactly once, that the 8-mer set is identical
across sources, and that no cluster is a size outlier.

## 9. Dataset composition

### 9.1 Acquisition and attrition — UniPROBE

The UniPROBE download index lists 36 accessions. Two advertise a contiguous 8-mer archive the
server does not hold, and one further accession provides no detail pages, leaving **33
archives retrieved**. Of those, three publish detail pages carrying no sequence at all, so the
assayed construct is unknown and no measurement can be attributed to a protein.

The remaining accessions supply **779 constructs carrying a sequence**. Each was scanned
against Pfam-A and put through the admission policy of §3.2:

| outcome | constructs | share |
|---|---:|---:|
| admitted | 581 | 74.6% |
| rejected — repeated array | 106 | 13.6% |
| rejected — mixed families | 50 | 6.4% |
| rejected — no domain above threshold | 42 | 5.4% |

Admission of a construct is necessary but not sufficient. A further 20 admitted constructs
belong to accessions whose measurements could not be used: one publishes 8-mer tables
truncated to the enriched end, and the rest provide no experiment file that resolves to a
complete, readable 8-mer table. Twelve registered accessions therefore contribute nothing,
each reported with its reason rather than passed over.

Of the 561 admitted constructs in usable accessions, **501 distinct domain sequences** remain;
the difference is constructs that resolve to an identical stored sequence, chiefly paralogues
whose domains are invariant and separate array designs of one protein. These are reconciled
rather than deduplicated (§6).

### 9.2 A second distributor: CIS-BP / Weirauch 2014

UniPROBE is not the only repository of universal-PBM data. Weirauch et al. (2014) cloned and
assayed 1,032 DNA-binding domains across 124 organisms, depositing the measurements in GEO
(`GSE53348`, 2,064 samples) and the **assayed insert sequences** in a supplementary
spreadsheet (Table S6). That spreadsheet is what makes the source admissible: attributability
requires the sequence physically on the array, and each GEO sample is titled
`<plasmid>_<HK|ME>_8mer_<n>`, so the plasmid identifier joins measurement to sequence. The
join is exhaustive — 1,032 plasmids on each side, none unmatched.

Its measurements need no special treatment: the same Wilcoxon E-score on the same
[-0.5, 0.5] scale, and the paper adopts `E > 0.45` as its own significance threshold, which is
the cutoff already in use here. The two array designs (HK, ME) carry different probe sequences
and are treated as replicates under §6.

Two properties do differ and are recorded rather than absorbed. Constructs come in **three
architectures** — 671 with 50 endogenous flanking residues, 96 with 15, and 265 with none —
so a zero-flank construct yields a shorter stored domain than a UniPROBE construct of the same
protein, because the padding of §3.3 is clipped by the construct. And the positive rate is
lower (0.16% against 0.31%), which is expected: this panel reaches far more organisms and more
weakly-binding factors.

Of its 1,032 constructs, 896 satisfy the admission policy, giving 884 distinct domain
sequences of which 21 were already held from UniPROBE.

### 9.3 The dataset as it stands

Build of 2026-08-14, from 19 contributing sources.

| property | value |
|---|---|
| rows | 46,120,192 |
| distinct DNA-binding domains | 1,364 |
| clusters | 1,321 |
| Pfam families | 56 |
| source organisms | 138 |
| distinct 8-mers | 32,896 (identical in every source) |
| binding / non-binding / excluded | 97,320 / 45,238,991 / 783,881 |
| negative:positive ratio | 465:1 |
| `dbd_seq` length | 30-377 aa (median 76) |

The length range is set by the domains themselves, not by an error: the shortest are AT-hooks,
a ~10-residue motif that the padding of §3.3 extends to about 30, and the longest is a plant
GRAS domain of ~350 residues. A 30-residue AT-hook is a motif rather than a fold, which is
worth knowing before it reaches a structure predictor.

Composition by family:

| family | domains | family | domains |
|---|---:|---|---:|
| Homeodomain | 440 | GATA | 29 |
| bHLH (`HLH`) | 103 | ETS | 26 |
| bZIP | 79 | HMG box | 26 |
| Zn2Cys6 (`Zn_clus`) | 74 | NAC (`NAM`) | 26 |
| Forkhead | 72 | WRKY | 23 |
| zf-C4 | 62 | TCP | 21 |
| Myb | 60 | Dof (`Zn_ribbon_Dof`) | 18 |
| AP2 | 44 | ARID | 17 |

**27 families hold ten or more domains**, against eight before this source, with a further 29
below that. Leave-one-family-out is a real test on 27 of them rather than on a handful.

Protein-side depth remains concentrated. Each row below is one cluster *size*: how many
clusters contain exactly that many domains, and how those domains divide into references and
variants.

| domains per cluster | clusters of this size | domains in them | references | variants |
|---:|---:|---:|---:|---:|
| 1 | 1,293 | 1,293 | 1,293 | 0 |
| 2 | 8 | 16 | 8 | 8 |
| 3 | 6 | 18 | 6 | 12 |
| 4 | 5 | 20 | 5 | 15 |
| 5 | 4 | 20 | 4 | 16 |
| 6 | 2 | 12 | 2 | 10 |
| 7 | 1 | 7 | 1 | 6 |
| 8 | 2 | 16 | 2 | 14 |
| | **1,321** | **1,402** | **1,321** | **81** |

The columns are related exactly. *Domains in them* is simply the size times the number of
clusters — the six clusters of size three account for 18 domains. Every cluster contains
**exactly one reference** (verified: all 1,321), so *references* equals the cluster count and
*variants* is the remainder:

> variants = domains − clusters

Hence a cluster of size 1 is a reference with nothing to compare it to and contributes no
variants, while the two clusters of size 8 contribute 7 variants each.

**98% of clusters hold a single domain.** The 81 variants sit in the 28 clusters of size two
or more, and half of them in the eight clusters of size five or more. The largest are
`BAR15A:HOXD13` and `ROG18A:FoxJ3` at eight domains each, then `BAR15A:FOXC1` at seven.

Domains covered (1,402) exceeds distinct domains (1,364) because 38 domains appear in more
than one cluster: the same sequence assayed by two studies enters each study's cluster
independently, which is what makes the cross-source label agreement in `overlap.md`
measurable. Adding CIS-BP raised that from 17 to 38, so the agreement check now has more than
twice the evidence — and for the first time it compares two different laboratories rather than
two deposits from the same one.

Variants are predominantly homeodomain (54), then forkhead (18), zf-C4 (5), bHLH (3) and
paired domain (1). Two of the three families carrying variants outside the homeodomain arose
only once distinct constructs sharing a gene directory were separated (§4): a bHLH
point-mutant series and a set of forkhead chimeras.

The companion protein table carries all 1,364 domains, each located inside the construct it
was cut from. Full-length coverage is now the exception rather than the rule: Table S6
publishes no UniProt accession, so a canonical sequence resolves for 462 domains (34%) and the
domain is located within it for 354 (26%). Domain-level and construct-level representations
are available for every one; full-protein for about a third.

## 10. Reproducibility

The dataset is rebuilt from the raw archives by a single command, and every parameter that
affects its content lives in one configuration file. Reports summarizing binarization, cluster
inventory and validation are regenerated per source and version-controlled, so a change in the
dataset appears as a diff.

All statistics in §9 are from the build of 2026-08-14 and reflect the procedures described
above; no method change is outstanding. Consistency between the tables, the protein table and
the reports is asserted after each build: every corpus domain appears in the protein table with
offsets that index the stored sequence, every source carries the complete 8-mer set exactly
once, and every source has a full set of current reports.

## References

Sources, with DOIs verified against Crossref, are listed in
[`REFERENCES.md`](REFERENCES.md). The accession survey, with per-accession citations and
family composition, is in [`UNIPROBE_ACCESSIONS.md`](UNIPROBE_ACCESSIONS.md). The admission
policy is stated in full, with its exclusion audit, in
[`DOMAIN_POLICY.md`](DOMAIN_POLICY.md).
