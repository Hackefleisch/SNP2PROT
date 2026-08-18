# Methods

Construction of a paired transcription-factor DNA-binding domain × DNA-site dataset with
binary binding labels.

This document records what was done and why, in enough detail for an independent
reimplementation. Parameters quoted here are the values in
[`configs/thresholds.yaml`](../configs/thresholds.yaml); all statistics are from the complete
rebuild of **2026-08-17** (§10.2). Where a choice had alternatives, the alternative and the
measurement that decided between them are given — the full record is
[`DECISIONS.md`](DECISIONS.md), and this document cites its sections rather than repeating
them.

It is written to be readable without a background in structural biology or in machine
learning. The glossary defines every domain-specific term used; §1 states what the dataset is
for, which is what makes the rest of the choices follow.

---

## Glossary

| term | meaning |
|---|---|
| **transcription factor (TF)** | A protein that binds specific DNA sequences and regulates transcription of nearby genes. |
| **domain** | A part of a protein that folds into a stable shape on its own and carries one function. Proteins are built from domains the way a program is built from functions: they recur across otherwise unrelated proteins, and one can often be studied in isolation. Throughout this document "domain" means one stored `dbd_seq` — one protein entry in the dataset. |
| **DNA-binding domain (DBD)** | The compact, independently folding part of a TF that contacts DNA. The rest of the protein does other work and is not represented here. Stored as `dbd_seq`. |
| **construct** (or *insert*) | The protein fragment actually cloned and put on the array. Usually a DBD plus some flanking residues; not the full-length protein. Constructs bound what can be observed, and different laboratories choose different amounts of flank. |
| **family** | A group of DBDs sharing a fold and a recognition mechanism — homeodomain, forkhead, bZIP, C2H2 zinc finger. Assigned from Pfam and stored as `dbd_family`. The coarsest unit at which generalization can be tested. |
| **Pfam / HMM** | Pfam is a library of protein-family models; each family is a profile hidden Markov model scoring how well a stretch of sequence matches that family. Used to locate domains and assign families. |
| **envelope** | The span of a construct that a Pfam model matches, in residue coordinates. |
| **canonical domain** | The stored sequence: the Pfam envelope ± 10 residues, made independent of how much flank the construct happened to carry (§3.4). Project-internal — unrelated to UniProt's sense of "canonical". |
| **reference** (or *wild type*) | The unmodified version of a domain, against which its variants are compared. `n_mut_from_wt = 0`. Also used in §3.4 for the full-length protein a construct is extended from. |
| **variant** | A version of a reference domain carrying one or more sequence changes — a point substitution, or an engineered chimera. |
| **cluster** | A reference domain together with the domains within a fixed distance of it, keyed by `wt_id`. The unit at which train/test partitions must be made (§7). |
| **paralogue** | Two genes descended from a duplication within one organism. Their DBDs are sometimes byte-identical while the rest of the protein has diverged. |
| **site** (or *8-mer*) | The DNA sequence whose binding is measured — here always 8 base pairs. Stored as `dna_seq`. |
| **PBM** | Protein-binding microarray. One experiment scores one protein against every possible DNA sequence on the array at once, giving a complete profile rather than a handful of measurements. |
| **E-score** | The PBM enrichment statistic: a rank-based measure bounded to [−0.5, 0.5], comparable across experiments, high for bound sequences. Thresholded to produce the binary label. |
| **replicate** | A repeat measurement of the same protein, possibly on a different array design. Reconciled rather than averaged (§6). |
| **distributor** | A repository that publishes another laboratory's measurements. Two are used here, UniPROBE and CIS-BP; "source" (`source_dataset`) means one accession or deposit within a distributor. |

---

## 1. Design objective and its consequences

The dataset supports models that are sensitive to variation on the protein axis and the DNA
axis **simultaneously**: many single-residue variants of the same DBD, each scored against
many single-base variants of a site. Three requirements follow, and each one excludes
otherwise sound data.

**The protein-side unit must be attributable.** A measurement is usable only if the stored
protein sequence is the thing that produced it. This is much stricter than "the protein was a
transcription factor", and it is enforced as an explicit three-condition policy (§3.2).

**The DNA axis must be uniform across sources.** Binding-site length is confounded with assay
identity, and assay identity is confounded with the positive rate, so pooling assays of
differing site length leaks the label prior — length alone would let a model infer which assay,
and therefore which prior, a row came from. The dataset is restricted to **universal PBM**
experiments, in which every protein is scored against the same complete set of 32,896
non-redundant 8-mers. The design is fully crossed and there is no length variation at all.

**The protein axis must be dense.** If every DBD in the set is distant from every other, a
model can memorize a per-protein profile and score well without learning anything transferable
— a lookup table with extra steps. Density is why clusters exist (§7) and why the evaluation
protocol holds out whole clusters rather than random rows.

### 1.1 Scope decisions and what they exclude

Three bodies of data were considered and dropped. They are recorded because their absence is a
design choice, not an oversight ([`DECISIONS.md`](DECISIONS.md) §1).

| excluded | reason |
|---|---|
| bacterial one-hybrid (B1H), e.g. Persikov et al. | The assay varies a different subunit than the one that binds; the varied sequence is not the sequence responsible for the measurement (condition 1). ~8,000 domains. |
| C2H2 zinc-finger arrays | Two to six separate ~23-residue folds on flexible linkers, each with its own Zn²⁺ — not one continuous region (condition 2). |
| SNP-SELEX | 19 bp sites. Admitting it would reintroduce the length/assay confound above. |
| bHLH heterodimer held-out sets | Two chains forming one binding unit (condition 1). Dropped 2026-08-17. |

The consequence is a **PBM-only corpus**: every stored row is an 8 bp site, so the length
confound cannot occur while that holds. It returns the moment any non-PBM assay is admitted.

---

## 2. Data sources and acquisition

Two distributors of universal-PBM data are used.

**UniPROBE** (Hume et al. 2015, *NAR* 43:D117). The download index lists 36 accessions. For
each, the contiguous 8-mer archive `<ACC>_contig8mers.zip` was retrieved together with one
HTML detail page per gene, which carries the assayed clone insert sequence, the Pfam domain
label, the UniProt accession and the source organism.

**CIS-BP / Weirauch et al. 2014**. 1,032 DBDs cloned and assayed across 124 organisms, with
measurements in GEO (`GSE53348`, 2,064 samples) and the **assayed insert sequences** in
supplementary Table S6. The insert sequences are what make the source admissible:
attributability requires the sequence physically on the array, not a database's annotation of
the full-length protein. Each GEO sample is titled `<plasmid>_<HK|ME>_8mer_<n>`, so the plasmid
identifier joins measurement to sequence; the join is exhaustive, 1,032 plasmids on each side
with none unmatched ([`DECISIONS.md`](DECISIONS.md) §10).

Every retrieved file is recorded in [`PROVENANCE.md`](../PROVENANCE.md) with its URL,
accession, UTC retrieval time, byte size and SHA-256; detail pages through a per-accession
`_manifest.tsv` carrying the same fields. Raw files are never modified — all parsing writes to
a separate tree — and no identifier, URL or filename is ever constructed by pattern. A download
that fails is reported, never substituted with a similar dataset.

Supporting reference data, cached under `data/external/`: the complete **Pfam-A HMM library**
(30,134 families, HMMER3/f 3.3 format) and canonical **UniProt** sequences fetched on demand
and cached one FASTA per accession. The UniProt cache has a `PROVENANCE.md` row; **the Pfam-A
library does not yet**, which is a gap in the provenance record rather than in the data — the
whole admission policy rests on that library, so its download URL, release and checksum are
owed (`TODO.md` T10).

Three classes of source were unusable and are recorded as such rather than replaced:

- two accessions advertise an archive the server does not hold;
- three publish detail pages carrying no sequence, so the assayed construct is unknown;
- one publishes 8-mer tables truncated to the enriched end (341–1,391 rows cut at E ≥ 0.25,
  against 32,896 for a complete experiment). These carry no non-binding evidence and would
  invert the labels of that protein's strongest sites if thresholded.

A file is additionally required to contain the complete non-redundant 8-mer set. The design is
fully crossed, so a partial table cannot contribute negatives; one accession was excluded for
being a single 8-mer short of 32,896. Attrition is quantified in §9.1.

---

## 3. Protein-side annotation

### 3.1 Domain identification

Each **assayed construct** — not the full-length protein — is scanned against the complete
Pfam-A library using HMMER3 through pyhmmer 0.12.1, with each family's curated **gathering
threshold** as the inclusion criterion.

Scanning the construct rather than the protein is deliberate. A distributor's domain label
describes the full-length protein, and for constructs carrying only part of a protein the two
disagree: constructs annotated `Homeobox, POU` were found to contain both domains inside the
assayed insert, which the protein-level label does not establish.

A whole-library scan also returns domains unrelated to DNA binding, so the admission policy is
evaluated over a restricted set of **65 DNA-binding families**. A family absent from that set
is indistinguishable from an absent domain, so its derivation matters: it is taken from the
source databases' own curation of what binds DNA — UniPROBE's `Domain` labels and CIS-BP's
per-construct `Pfam ID` — each resolved against Pfam-A `NAME` fields rather than trusted
verbatim (`scripts/build_dbd_family_list.py`).

Deriving it from a single database proved insufficient twice. The first set covered only the
nine families in the earliest sources. The second covered 41, all of them families UniPROBE
happens to publish, which silently discarded constructs carrying plant TCP, Dof, NAC and SBP
domains, Doublesex `DM`, and the Rel homology domain of NF-κB — of 240 CIS-BP constructs
rejected as `no_domain`, only 33 had no Pfam hit at all; the other 207 had a good hit belonging
to a family the whitelist had no word for. Adding CIS-BP's vocabulary took the set to 65.

Both vocabularies are stale in places. Labels predating a Pfam rename are resolved by
accession lookup where the accession is published (`Homeobox` → `Homeodomain`, `Fork_head` →
`Forkhead`, `MADS` → `SRF-TF`, `E2F_TDP` → `WHD_E2F_TDP`, `Zn2Cys6` → `Zn_clus`, `BRIGHT` →
`ARID`, `AP-2` → `TF_AP-2`, `PBX` → `PBC`, `RFX` → `RFX_DNA_binding`), and otherwise **by the
sequence**: the constructs carrying a given label are scanned against the full library and the
family that actually hits them is read off (`CXC` → `TCR`, `DUF260` → `LOB`, `DUF573` →
`GeBP-like_DBD`, `EIN3` → `EIN3_DNA-bd`, `RHD` → `RHD_DNA_bind`, `zf-Dof` →
`Zn_ribbon_Dof`). That procedure reproduces the accession-derived mappings exactly where both
are available, which is the basis for trusting it where only one is.

Five Pfam model names denoting one biological family are merged so that `dbd_family` is a
family and not a model id (`bZIP_1`/`bZIP_2` → `bZIP`, `zf-H2C2_2` → `zf-C2H2`, `Homeobox_KN` →
`Homeodomain`, `SOXp` → `HMG_box`). The list is deliberately minimal: `TF_AP-2` (mammalian
TFAP2) and `AP2` (plant AP2/ERF) merely look alike by name and are **not** merged.

Pfam hits outside the DNA-binding set are recorded but do not affect admission: a DNA-binding
domain adjacent to an unrelated domain remains solely responsible for the DNA interaction.

### 3.2 Admission policy

A construct enters the dataset only if its stored domain satisfies three conditions.

1. **Sole responsibility.** The stored subunit alone produced the measured interaction.
2. **One continuous region.** The subunit is a single contiguous stretch of sequence.
3. **Containment of variation.** Every sequence change that alters binding lies inside the
   stored region.

The conditions follow from the intended use. A model that predicts a structure, locates the
centre of the DNA-contacting residues and embeds the residues around it requires a subunit that
folds as a unit, whose spatial arrangement is determined, and within which the experimental
contrast is visible. Condition 3 has a second, sharper justification: a variant whose
substitution falls outside the stored window is *sequence-identical to its own reference while
carrying a different label*, which is an unlearnable contradiction rather than a hard example.

Four rejection categories follow:

| category | condition | description |
|---|---|---|
| mixed families | 1, 2 | Pfam hits from two DNA-binding families in one construct (paired + homeodomain; POU-specific + homeodomain). Both contact DNA; the data do not attribute a given 8-mer's score to either, and their relative geometry is undetermined without DNA. |
| repeated arrays | 2 | Multiple hits of one family, characteristically C2H2 zinc-finger arrays. Not one continuous unit. |
| protein complexes | 1 | Two or three distinct chains forming one binding unit (bHLH heterodimers; Fos/Jun; CSL–NOTCH–MAML ternary complexes). No single sequence is responsible. 22 archive entries. |
| no domain | — | No Pfam hit above the gathering threshold anywhere in the construct. |

One refinement prevents a false rejection. Two Pfam models of the same region — `bZIP_1` and
`bZIP_2` both matching residues 18–74 — are one domain seen twice, not an array. Hits that
overlap are collapsed before the policy is applied, and only genuinely separate hits count as
repeats.

Rejections are reported per source and never repaired. Both strict rules are reversible through
`domain.allow_repeat_arrays` and `domain.allow_mixed_families`, which readmit the C2H2 arrays
and the PAX/POU constructs respectively; both change what the dataset *is*, which is why they
are configuration rather than code.

### 3.3 Domain boundaries: the padded envelope

The stored sequence is the Pfam envelope **extended by 10 residues on each side**
(`domain.padding_aa`).

The padding is required by condition 3. Pfam models begin at the structural core; for
homeodomains this excludes the N-terminal arm, which inserts into the DNA minor groove and
makes base contacts. Empirically, four of 91 point substitutions in one panel fall outside the
bare envelope, three of them 2–9 residues N-terminal of the Pfam start. Trimmed to the bare
envelope those variants become identical to their own references. Ten residues covers every
observed case with margin, and is set from the position of the N-terminal arm rather than
fitted to the four. The condition is re-checked for every variant at parse time.

This is also why the stored form is the **residues in sequence order** rather than a
fixed-width profile alignment: a match-state-only representation discards insert columns, which
collapses the same variants onto their own references — the exact failure the padding prevents.

### 3.4 The canonical domain: removing the construct from the sequence

The padded envelope alone is not enough, because it is **clipped wherever the construct
happens to end**. The same domain from two laboratories was therefore stored as two different
strings — VENTX is 76 residues from one deposit and 57 from a zero-flank construct in another.
That difference is cloning, not biology, and it reaches a model as though it were biology.

The stored sequence is consequently a **project-internal canonical form**
(`snp2prot.canonical`), whose defining property is *stability*: the same domain in gives the
same string out, however truncated the input was. Length still varies with the biology, because
an insertion or deletion genuinely is a different sequence.

Three routes produce it, tried in order.

1. **From the construct.** If the construct already spans envelope ± 10, the window is a
   slice of it and no external data is consulted. This is provably the string the full-length
   protein would have given. Measured when the canonical form was introduced, **1,120 of 1,477
   admitted constructs (76%)** took this route.
2. **Extended from a reference protein.** The remainder are extended using the
   full-length protein — 357 constructs on that measurement. The reference **only ever
   extends; it never replaces** — every position
   the construct covers is taken from the construct, carrying whatever mutations were
   engineered into it, and the reference supplies only the flank the construct is missing. A
   short sequence caused by a real deletion or a true protein terminus therefore stays short,
   because that is what was assayed.
3. **Left short, and flagged.** Engineered chimeras have no natural protein to extend from —
   one panel's FoxJ3/FoxN3 hybrids are synthetic and were cloned as bare domains. Their
   termini are not truncations; nothing was cut off. Discarding them would have removed the
   corpus's only non-homeodomain variant depth, so they are kept, and the resulting length
   duplication is resolved at clustering instead (§7.1).

**References are resolved by cross-reference, never by gene name** (`snp2prot.references`).
Each construct needing one carries an identifier, but not the same kind: UniPROBE publishes a
UniProt accession, CIS-BP a `Gene ID` from whichever database that organism's genome project
used — Ensembl, FlyBase, Araport/TAIR, RefSeq, SGD, WormBase, dictyBase, NCBI GeneID, and a
tail of assembly-specific identifiers that resolve nowhere. Searching UniProt for a gene symbol
plus an organism is tempting and wrong: `AT1G28420` resolves by cross-reference to a
1,705-residue protein and by gene-name search to a 317-residue one. A construct that cannot be
resolved is discarded rather than approximated.

A placement is accepted only if the construct aligns into the reference at **≥95% coverage**
with **≤5 edits outside the Pfam envelope**. Both numbers are measured, not chosen
([`DECISIONS.md`](DECISIONS.md) §11):

| edit ceiling | constructs accepted | with ≥99% aligned |
|---|---:|---:|
| ≤1 | 94.6% | 92.4% |
| ≤2 | 96.3% | 94.2% |
| **≤5** | **97.6%** | **95.2%** |
| ≤10 | 97.8% | 95.2% |

Five sits at the knee — ten buys 0.2% more — and what it rejects is not marginal (`Hoxc11` at
111 edits over 57% coverage). **Edits are counted outside the envelope only**: inside it,
mutations are the subject matter of the dataset and are kept unconditionally, whereas outside is
where residues are borrowed, so that is what must match. A flat ceiling would reject the
engineered chimeras, which differ from their parent by 6–8 substitutions while being perfectly
well placed. Mismatches and gaps are counted separately, because 0–1 mismatches with 25–54 gaps
is an alternative isoform rather than a wrong protein. The 0.95 coverage bound is likewise set
at an observed gap in the distribution: rejections sit at 0.987/0.986/0.986 and then at 0.842
and below, with nothing between.

**Constructs differing outside the canonical window are discarded, not merged.** Of 54 stored
domains produced by more than one construct, 24 came from byte-identical constructs and 26
differed only in cloned flank — both benign — but four carried a genuine substitution outside
the stored window and were being silently reconciled as replicates of one protein. A difference
the model cannot see, attached to measurements that may differ because of it, is a confound;
this is condition 3 applied symmetrically.

A **companion protein table** (`data/interim/proteins/`) stores each domain at four levels —
bare Pfam envelope, padded canonical domain, assayed construct, and full-length UniProt
sequence, with offsets between them — so representation can be varied downstream without
re-parsing. Construct architecture (how much flank each construct carried) is recorded here
rather than as a row column, because flank length, affinity tag and expression system are
perfectly confounded with the source study: as row columns a model could learn "which
laboratory made this" instead of "what does this residue do".

---

## 4. Assignment of sequences to measurements

UniPROBE archives encode the assayed construct as a path element, and **a single gene directory
may contain several distinct proteins**: engineered chimeras, point-mutant series, or separate
domain and full-length constructs. Detail pages correspondingly carry one insert sequence per
construct. This was the costliest defect found in the project — treating a gene directory as
one protein silently averaged different proteins as replicates of each other.

Experiments are therefore matched to constructs **by name, longest match first**, rather than
grouped by gene, and constructs are keyed by the full name the page gives. Where a page labels
its sequence by plasmid or gives no name, the gene name is used; an empty key would fuse every
such gene in an accession into one cluster.

Two situations are resolved conservatively. Where a gene's experiments match no construct on
its own page, the data are not assigned. Where an archive contains experiments for orthologues
of two species but the page carries only one species' sequence, no assignment is made, since a
name-based match would attach one orthologue's sequence to the other's measurements.

**Organism names** are normalised to one convention across deposits
(`snp2prot.metadata.organism`): the UniProt-style binomial — capitalised genus, lower-case
epithet, no strain, no abbreviation, no hybrid marker. Anything that is not an organism is
stored as the empty string rather than as a word that looks like one. Two distinct problems are
handled.

*Missing organisms.* A UniPROBE detail page can carry the unsubstituted template placeholder
`$species`, a defect in the deposit. Where the construct carries an accession whose canonical
FASTA is already cached, the organism is read from that record's `OS=` field; where it does
not, the field is left empty. Nothing here fetches, so parsing stays offline and a cold cache
degrades to an empty organism rather than to a wrong one. The study's scope is deliberately not
used as evidence: all 17 affected constructs came from one *C. elegans* panel, so guessing
would have been right in this case and would have been a rule that invents an organism the next
time a mixed-species deposit has the same defect. 13 of the 17 resolved. Seven further strings
that occupy the field without naming an organism are treated the same way — `N/A` for
reconstructed ancestors, `Chimera` for engineered constructs with no source organism, a stray
spreadsheet heading, and UniPROBE's own `None Available`.

*Inconsistent spellings.* Deposits disagree about how to write a name that is not itself in
dispute. Corrections are **enumerated, not inferred**: an alias table with a reason per entry
can be reviewed, whereas a rule that rewrites any string matching a pattern will eventually
rewrite a name that was correct. Three entries, each found by surveying the values actually
present — an abbreviated genus (`C. elegans`), a transposed letter in *Acyrthosiphon pisum*,
and a hybrid marker (`Malus x domestica`). Three values are deliberately **not** normalised:
`Sarsia sp. Long Island Sound` is a genuinely unnamed species rather than a formatting variant;
`Physcomitrella patens` has been *revised* to *Physcomitrium patens*, which is taxonomy rather
than orthography; and `Acanthamoeba polyphaga mimivirus` is a virus, correctly named in three
words.

---

## 5. Label assignment

PBM E-scores are rank-based enrichment statistics bounded to [−0.5, 0.5] by construction.
Labels are assigned per 8-mer:

| E-score | label | meaning |
|---|---|---|
| ≥ 0.45 (`pbm.positive`) | 1 | binding |
| ≤ 0.35 (`pbm.negative`) | 0 | non-binding |
| between | −1 | no call; retained in the table, excluded from training |

Thresholds are applied **per experiment**, never to pooled scores, because array designs differ
in intensity scale and a pooled cutoff silently reweights arrays. The applied cutoffs are
written into every row (`threshold_pos`, `threshold_neg`), so a stored table is self-describing
and no cutoff literal appears in parser code.

Negatives from a PBM are **assayed-and-unbound**, not merely unobserved, and are recorded as
such in a mandatory `neg_provenance` field on every negative row. This is enforced by the
validator. It exists so that PBM negatives stay distinguishable from weaker forms of negative
evidence — a sequence merely absent from a selection experiment is not evidence of non-binding
— should another assay type ever be admitted.

### 5.1 What a fixed cutoff does and does not guarantee

The E-score is a rank statistic: for each 8-mer, the Mann-Whitney statistic of the probes
carrying it against those that do not, computed within the brightest half of the array and
shifted to [−0.5, +0.5]. Because only ranks enter, every monotone change of laser power,
protein concentration, antibody labelling or scanner cancels. That invariance is what lets one
cutoff span 19 sources and two decades of arrays, and `E ≥ 0.45` states, concretely, that at
most about three of an 8-mer's ~32 probes are less than convincing.

What the statistic does **not** equalise is sensitivity. A noisy or weak array shuffles the
probe ranking slightly, each shuffle costs the statistic, and the whole distribution compresses
toward zero. Measured on the 47 domains held by two sources, the median pair differs **1.6×**
in how many positives it calls at the same cutoff, and the tail is worse — `Hoxa2` is 165
against 17. Rank-matched labelling removes that asymmetry and raises cross-source agreement by
8–12 points, but it asserts that every protein binds a fixed number of 8-mers, and it would
turn a negative into "outside this experiment's top N" rather than measured non-binding. The
fixed cutoff was kept for that reason (`docs/DECISIONS.md`, 2026-08-18); the asymmetry is
recorded here as a property of the data rather than removed by definition.

**At the extreme, an experiment sees nothing at all.** 54 of 1,383 domain records have no
positive 8-mer. They are of two kinds, and the difference matters because the second kind is
the signal this dataset exists to carry:

| kind | records | what it is |
|---|---:|---|
| `dead_variant` | 20 | no positives, but another record of the same cluster **and source** has them — its own series is the control, so this is measured non-binding |
| `no_evidence` | 34 | no positives and no such control: indistinguishable from an experiment too weak to see anything |

18 of the 20 `dead_variant` records are `BAR15A` variants that lost binding. Of the 34
`no_evidence` records, 12 are `Cell09`, whose per-domain maximum E-score has a median of 0.428
— the source's typical domain never reaches the cutoff at its single best 8-mer. A further 13
of the 54 are silent only because their replicates disagreed: they reach 0.45 on the stored
replicate mean while the label is a no-call.

This is recorded per record in `data/interim/label_health/`, alongside label counts, the best
E-score, and how many 8-mers reach the cutoff regardless of label. **No label depends on it.**
Excluding unsupported records is a training decision, applied at featurization via
`label_health.usable`, on the same principle as cluster-size restriction.

**The E-score column is identified by content, not position.** Column order, column count and
the presence of a header row all vary between accessions, and one accession publishes no
E-score at all. The column is taken from a header name where one is given, and otherwise as the
unique numeric column bounded to [−0.5, 0.5] that also takes negative values — a criterion that
distinguishes an E-score from intensities, z-scores and p-values. Files in which no such column
exists, or in which several do, are rejected with a reason. Reading the wrong column by
position once produced a plausible-looking but wrong positive rate, which is why the check is
now structural.

---

## 6. Replicate reconciliation

Where an accession provides several experiments for one construct, each is binarized
independently and the labels are then combined: **agreement is retained, disagreement is
assigned to the no-call band**. Scores are never averaged before thresholding, since replicates
may sit on different array designs whose intensity scales differ by an order of magnitude.

Distinct genes whose stored domains are byte-identical are treated the same way. Such cases
arise from recent gene duplication, where the domain is invariant while the flanking sequence
has diverged. Their agreement is comparable to that of technical replicates of a single gene
(median positive-set Jaccard 70.4% against 72.2%), which is what supports the treatment.

### 6.1 Cross-source replicates and the label noise floor

Reconciliation above operates within a source. **47 domains are stored by more than one
source**, giving 49 pairs in which two laboratories assayed the same canonical domain against
the same 32,896 8-mers. These are *not* reconciled — each source's rows are retained as
deposited — but they are measured, because the rate at which two independent measurements of
one protein disagree is the noise floor against which any model's accuracy has to be read.
Regenerated into [`reports/overlap.md`](../reports/overlap.md).

| measure | value |
|---|---|
| pooled hard-label agreement | 99.98% (300 of 1,572,985 calls differ) |
| **median positive-call Jaccard** | **0.458** |
| pooled positives agreed on | 2,863 of 5,816 (49.2%) |
| median E-score rank correlation | 0.608 |

**The first row is not the informative one.** Positives are under 0.5% of an 8-mer table, so
agreeing on the negatives pins pooled agreement near 100% however badly two laboratories agree
about binding. On the calls the dataset exists to make, the median pair agrees on 45.8% of the
8-mers either of them called positive.

The comparison with §6 is the point: within a source, byte-identical domains and technical
replicates agree at 70.4% and 72.2% positive-set Jaccard; across sources the same quantity is
45.8%. Agreement falls by some 25 points when the two measurements come from different
laboratories rather than from one, so a substantial share of the disagreement tracks laboratory
and array design rather than the protein — and **a model that reproduces held-out positives
much beyond this level is reproducing a source, not a binding preference.**

Nine pairs fall below a positive-set Jaccard of 0.20 or a rank correlation of 0.30 and are
flagged rather than corrected. One is outside anything noise explains: the *Arabidopsis* NAC
domain ANAC092, stored by both Lindemose et al. 2014 and CIS-BP with a byte-identical
144-residue domain, **shares no positive call at all** between the two deposits (12 positives
against 126, disjoint; rank correlation 0.068). One of the two measurements is probably
attached to the wrong construct; both rows currently sit in one cluster, so no split separates
them. Open for the project owner, not resolved here.

---

## 7. Cluster construction

A cluster (`wt_id`) is the unit at which train/test partitions must be made: a partition that
separates two variants of one protein does not test generalization, because the model can copy
an answer rather than derive it. A cluster is **not** a claim that its members are
interchangeable — every member remains its own training row, and a cluster never collapses
them.

Clusters are built on **canonical domains, never on reference proteins**. Protein-level
clustering is tempting because references are construct-independent, but DBDs are conserved
while the rest of the protein diverges: among pairs whose domains are within 5 edits, median
domain identity is 96% while median full-length identity is 68%. Irx3 and Irx4 have 94%
identical domains and 32% identical proteins; any protein-level clustering separates them,
after which a model trained on Irx3 predicts Irx4 for free — precisely the leakage clusters
exist to prevent.

### 7.1 Subsumption

Route 3 of §3.4 leaves a small number of domains stored at a short window while another source
holds the same domain at full width. Two domains that are **zero edits apart under free
terminal gaps are the same domain**, so all but the longest are dropped before clustering. In
the 2026-08-17 build this removed **17 constructs (559,232 rows)** and is the reason no two
stored rows describe one domain at two lengths.

### 7.2 Distance

Distance comes from a global pairwise alignment (BLOSUM62, affine gaps), so variants may carry
insertions and deletions. **Terminal gaps are free.** This is required rather than convenient:
even after canonicalisation a short-window entry ends where its construct did, and charging for
those terminal gaps would report indels between a protein and its own chimera. Internal indels
are counted.

`n_mut_from_wt` is the number of edited residues and `mut_positions` lists their positions in
the **reference's** coordinate frame — one entry per edited residue, so `len(mut_positions) ==
n_mut_from_wt` always holds and every member of a cluster shares one coordinate system that
maps onto one predicted structure. A variant carrying a deletion can therefore name a position
past its own length.

### 7.3 Algorithm and threshold

**CD-HIT greedy incremental clustering** (Li & Godzik 2006; equivalently MMseqs2
`--cluster-mode 2`), reimplemented on the project's aligner rather than by adding a C++
dependency to a job that takes 16 seconds. Sequences are sorted by decreasing length, the
longest becomes a cluster representative, and each remaining sequence is compared **only to
representatives**, joining the first it is close enough to or founding a new cluster.

Three properties earn it the job:

* the longest member becomes the representative, so the least-clipped form is the reference;
* **comparison is never transitive**, which kills chaining. Single-linkage at the same
  threshold produced a 35-domain blob — A near B near C, with A and C unrelated, is not a
  cluster;
* every member is within *k* of the representative, which is exactly the invariant
  `mut_positions` needs.

Its known weakness is order dependence: a sequence joins the first representative it matches
rather than its best. Ties are broken deterministically by length then sequence, so a rebuild
reproduces the same clusters. Only sequences of the same Pfam family are compared, which is
both correct and a large saving.

**The threshold is 5 edits** (`cluster.max_edits`, owner's decision). It gives 1,133 clusters
with a largest of 8 domains.

### 7.4 Cluster identity

A cluster is named after its representative's own lineage name, so it stays recognisable —
`C:BAR15A:HOXD13`, with the `C:` prefix marking it cluster-derived rather than
construct-derived. That name is **not unique**: one gene directory held six engineered chimeras
that now sit in six different clusters, all of which would otherwise carry the id
`ROG18A:FoxJ3` and collapse back into one, undoing the clustering for precisely the constructs
it matters most for. Collisions take a numeric suffix; the 2026-08-17 build has 8 such ids.

### 7.5 The cluster inventory

Selecting on cluster size — "train only on DBDs with at least five variants" — is a
training-time choice, not a dataset one: the dataset stores everything and the modeller
filters. What the dataset owes the modeller is the ability to make that selection cheaply.
Answering it from the row tables means grouping 45 million rows and counting distinct
sequences, roughly a minute of work to learn 1,133 numbers.

A per-cluster side table (`data/interim/clusters/clusters.parquet`, one row per cluster:
family, representative, domain count, variant count, maximum edit distance, contributing
sources, construct and row counts) is therefore written by the same pass that assigns the
cluster ids, and read through `snp2prot.clusters.load` / `ids_with_at_least` / `select`.

It is a side table rather than a `cluster_size` column on every row: a column would be a
Parquet predicate pushdown, but it would denormalise a per-cluster fact onto 45 million rows
and require changing the schema and every parser's output. **"Cluster size" means distinct
canonical domains** throughout — not rows and not constructs, so one domain assayed by two
laboratories counts once.

---

## 8. Schema and validation

Every parser emits an identical **22-column row schema**, and every table passes a validator
before it is written. The validator enforces, among other checks:

- unambiguous nucleotides (`ACGT` only) and a maximum site length. Sites are **never padded**
  in the stored table — padding is a modelling decision, and rejecting non-ACGT characters is
  what keeps it out;
- a mandatory `neg_provenance` on every negative;
- agreement between `n_mut_from_wt` and the length of `mut_positions`;
- mutation positions bounded by the cluster's reference;
- E-scores within their theoretical bounds — reading the wrong column once produced values
  outside [−0.5, 0.5], so this is a hard error;
- deterministic row identity and its uniqueness. `pair_id` is a hash of domain sequence, site,
  assay and stringency, computed by the schema rather than supplied by parsers. `pair_id` and
  `dna_len` are derived, never hand-set.

Clusters are corpus-wide, so a per-source validator cannot see a cluster whose reference row
lives in another source; that check is a warning locally and is performed corpus-wide by the
clustering pass.

A separate **audit** (`scripts/audit_sources.py`) sweeps all accessions and checks properties a
per-table validator cannot see: that every experiment file in an archive is claimed by exactly
one construct, that every construct is emitted or rejected with a reason, that no two
constructs share a key, that every domain is measured against the complete 8-mer set exactly
once, that the 8-mer set is identical across sources, and that no cluster is a size outlier.

Every defect the audit exists to catch was originally found by a distribution looking wrong,
never by an exception or a failing test. The recurring ones — no header row; a column count
differing from the previous source; the E-score in a different column or absent; one gene
directory holding several engineered proteins; a bare or plasmid-named insert label; a flat
archive with no gene directories; a table truncated to the enriched end or *n* rows short of
32,896 — are the checklist any new source is screened against.

---

## 9. Dataset composition

### 9.1 Acquisition and attrition — UniPROBE

The download index lists 36 accessions. Two advertise a contiguous 8-mer archive the server
does not hold and one provides no detail pages, leaving **33 archives retrieved**. Of those,
three publish detail pages carrying no sequence at all, so the assayed construct is unknown and
no measurement can be attributed to a protein.

The remaining accessions supply **779 constructs carrying a sequence**. Each was scanned
against Pfam-A and put through the admission policy of §3.2:

| outcome | constructs | share |
|---|---:|---:|
| admitted | 581 | 74.6% |
| rejected — repeated array | 106 | 13.6% |
| rejected — mixed families | 50 | 6.4% |
| rejected — no domain above threshold | 42 | 5.4% |

Admission is necessary but not sufficient. A further 20 admitted constructs belong to
accessions whose measurements could not be used: one publishes truncated 8-mer tables, and the
rest provide no experiment file that resolves to a complete, readable table. **Twelve
registered accessions therefore contribute nothing**, each reported with its reason rather than
passed over. Restricted to accessions that do contribute, the policy admits 561 of 734
constructs (76.4%): 97 repeated arrays, 49 mixed families, 27 with no domain.

### 9.2 A second distributor: CIS-BP / Weirauch 2014

Its measurements need no special treatment: the same Wilcoxon E-score on the same [−0.5, 0.5]
scale, and the paper adopts `E > 0.45` as its own significance threshold, which is the cutoff
already in use here. The two array designs (HK, ME) carry different probe sequences and are
treated as replicates under §6.

Of its **1,032 constructs, 896 satisfy the admission policy (86.8%)** — 65 repeated arrays, 57
with no domain, 14 mixed families — giving 868 stored domains.

Two properties differ from UniPROBE and are recorded rather than absorbed. Constructs come in
**three architectures** — 671 with 50 endogenous flanking residues, 96 with 15, and 265 with
none — which is precisely the variation §3.4 exists to remove from the stored sequence. And the
positive rate is lower (0.16% against 0.31%), which is expected: this panel reaches far more
organisms and more weakly-binding factors.

One limit is recorded honestly: this is Weirauch 2014's own 1,032 constructs, **not** the
~2,294 TFs with PBM data the CIS-BP database aggregates. For aggregated entries the construct
sequence belongs to the contributing study, and much of that is UniPROBE, already held.

### 9.3 The dataset as it stands

Build of **2026-08-17**, from 19 contributing sources.

| property | value |
|---|---|
| rows | 45,495,168 |
| constructs | 1,383 |
| distinct DNA-binding domains | 1,335 |
| clusters | 1,133 |
| Pfam families | 56 |
| distinct 8-mers | 32,896 (identical in every source) |
| binding / non-binding / no call | 95,840 / 44,626,532 / 772,796 |
| negative:positive ratio | 466:1 (positive rate 0.214% of calls made) |
| `dbd_seq` length | 30–378 aa (median 77) |
| source organisms | 131, plus 21 constructs with none recorded |

The high negative:positive ratio is a property of the assay and is retained deliberately.
Universal PBM scores every protein against every 8-mer, and a transcription factor binds a
small minority of them; discarding negatives to balance the classes would discard the assay's
main contribution, which is dense, positive evidence of *non*-binding.

The length range is set by the domains themselves, not by an error: the shortest are AT-hooks,
a ~10-residue motif that padding extends to about 30, and the longest is a plant GRAS domain.
A 30-residue AT-hook is a motif rather than a fold, which is worth knowing before it reaches a
structure predictor. The 21 constructs without an organism are reconstructed ancestors,
engineered chimeras, and accessions UniProt no longer serves — none of which has one.

**The row count is an exact invariant.** 45,495,168 = 1,383 × 32,896: every construct
contributes exactly one full 8-mer table. If a rebuild's row count is not a clean multiple of
32,896, rows were dropped or duplicated. It is the fastest single sanity check available.

Composition by family — 27 families hold ten or more domains, with a further 29 below that, so
leave-one-family-out is a real test rather than a handful of proteins:

| family | domains | family | domains |
|---|---:|---|---:|
| Homeodomain | 427 | GATA | 29 |
| bHLH (`HLH`) | 100 | ETS | 26 |
| bZIP | 79 | NAC (`NAM`) | 26 |
| Zn2Cys6 (`Zn_clus`) | 74 | HMG box | 26 |
| Forkhead | 69 | WRKY | 23 |
| Myb | 59 | TCP | 20 |
| zf-C4 | 58 | ARID | 17 |
| AP2 | 41 | Dof (`Zn_ribbon_Dof`) | 17 |

Protein-side depth. Each row is one cluster *size*: how many clusters hold exactly that many
domains, and how those domains divide into references and variants.

| domains per cluster | clusters | domains in them | references | variants |
|---:|---:|---:|---:|---:|
| 1 | 1,011 | 1,011 | 1,011 | 0 |
| 2 | 88 | 176 | 88 | 88 |
| 3 | 13 | 39 | 13 | 26 |
| 4 | 7 | 28 | 7 | 21 |
| 5 | 8 | 40 | 8 | 32 |
| 6 | 2 | 12 | 2 | 10 |
| 7 | 3 | 21 | 3 | 18 |
| 8 | 1 | 8 | 1 | 7 |
| | **1,133** | **1,335** | **1,133** | **202** |

The columns are related exactly. Every cluster contains **exactly one reference**, so
*references* equals the cluster count and

> variants = domains − clusters

A cluster of size 1 is a reference with nothing to compare it to. **122 clusters hold more than
one domain**, 34 hold three or more and 14 hold five or more; 75 clusters draw on more than one
source. The largest is `C:BAR15A:HOXD13` at 8 domains, then `C:BAR15A:CRX`, `C:BAR15A:FOXC1`
and `C:PNAS08:PF14_0633` at 7.

Variants by family — Homeodomain 87, Myb 16, bHLH 15, forkhead 13, AP2 9, zf-C4 7,
RFX 5, TCP 5. The point is not the total but the spread: with variants in homeodomain only, the
project's central question — *does sensitivity to single-residue change transfer across folds?*
— is unanswerable. Eight families now carry variants. It is still thin.

The companion protein table covers **1,224 of 1,335 domains**, each located inside the
construct it was cut from. Full-length coverage is the exception rather than the rule:
Table S6 publishes no UniProt accession, so a canonical full-length sequence resolves for 439
domains (33%) and the domain is located within it for 333 (25%). Domain-level and
construct-level representations are available for every domain; full-protein for about a third.

---

## 10. Reproducibility

### 10.1 Determinism

Every parameter that affects dataset content lives in one configuration file,
`configs/thresholds.yaml`. Reports summarizing binarization, cluster inventory, validation and
cross-source overlap are regenerated and version-controlled, so a change in the dataset appears
as a diff.

Three properties were verified rather than assumed on the 2026-08-17 rebuild:

- **The clustering pass is idempotent.** It strips the `C:` marker to recover the lineage name
  underneath, so running it twice produces a byte-identical cluster table. Verified by
  comparison.
- **Committed reports do not depend on how many times a step was run.** `reports/clusters.md`
  previously stated "after dropping *N* superseded copies", where *N* describes the state of
  `data/interim/` when the script ran — 17 on a fresh build, 0 on a repeat. It now states the
  dataset property instead.
- **A single-source rebuild must be followed by the cluster pass.** `build_dataset.py --source
  X` writes `wt_id` in its lineage form while every other source on disk carries the cluster
  form, so rebuilding one source silently detaches it from its cluster. The table still
  validates and the row count is still right. Always re-run `build_clusters.py` afterwards.

### 10.2 The build, in order

Measured on the reference machine (24 threads, 62 GB RAM) for the 2026-08-17 build.

| # | command | time | produces |
|---|---|---:|---|
| 0 | `python scripts/press_pfam.py` | once per machine | binary Pfam index; without it every source pays ~20 s re-reading 2.2 GB of text |
| 1 | `python scripts/build_dataset.py --all` | **7 min 28 s** | `data/interim/<source>/`, 19 of 31 registered sources |
| 2 | `python scripts/build_clusters.py` | **1 min 27 s** | `wt_id`, `n_mut_from_wt`, `mut_positions` corpus-wide; cluster inventory; `reports/clusters.md` |
| 3 | `python scripts/build_protein_table.py` | 4 s | `data/interim/proteins/` (warm UniProt cache) |
| 4 | `python scripts/make_reports.py --source <S>` ×19 | 1 min 33 s | per-source binarization, cluster inventory, validation |
| 5 | `python scripts/make_overlap_report.py` | 6 s | `reports/overlap.md` |
| 5b | `python scripts/build_label_health.py` | 12 s | `data/interim/label_health/`, `reports/label_health.md` |
| 6 | `python scripts/audit_sources.py` | 2 min 22 s | cross-source invariant sweep |

Steps 1–6 total about **13 minutes**. Step 1 is dominated by CIS-BP, which is 29 of the 45.5
million rows; within it, one source (2.3 GB of text) is the critical path, so parallelising
across sources alone cannot go below what that one source costs. The build is currently serial;
processes give ~7.6× on this workload where threads give nothing, because the CSV parser holds
the interpreter lock.

Software: Python 3.11.13, pandas 3.0.5, pyarrow 25.0.1, pyhmmer 0.12.1 (HMMER3), Biopython
1.88. Test suite: `python -m pytest -q`, 169 tests, under a second.

### 10.3 Consistency assertions after a build

Every corpus domain appears in the protein table with offsets that index the stored sequence;
every source carries the complete 8-mer set exactly once; every source has a full set of
current reports; the row count is a clean multiple of 32,896.

---

## 11. Known limitations

Stated because a user of the dataset needs them, not because they are open work items.

1. **The label noise floor is high.** Two laboratories measuring the same domain against the
   same 8-mers agree on 45.8% of their positive calls (§6.1). This bounds what any evaluation
   on this data can mean.
2. **One replicate pair is irreconcilable.** ANAC092 (§6.1) shares no positive call between its
   two deposits, and both rows sit in one cluster.
3. **Protein-axis depth is thin and homeodomain-heavy.** 202 variant domains across 122
   clusters; 1,011 of 1,133 clusters hold a single domain, and 87 of the 202 variants are
   homeodomain.
4. **Full-length protein sequence is available for about a third of domains** (§9.3), almost
   entirely because one source publishes no UniProt accession. Domain-level and
   construct-level work is unaffected; full-protein embeddings are not.
5. **One domain is two different organisms depending on the source.** Organism *spelling* is
   normalised (§4), but the domain under `C:Cell08:Tlx2` is byte-identical in two sources and
   stored as *Mus musculus* by one and *Homo sapiens* by the other. A homeodomain identical
   across mouse and human is entirely possible, so these may be two correct records — but their
   labels also disagree (Jaccard 0.097), and one deposit holding the wrong protein would
   explain both facts at once. Unresolved.
6. **One domain carries two unresolved `X` residues** (an 89-residue zf-C4 domain, 1 of 1,335).
   Harmless to a sequence model; a genuine problem for structure prediction.
7. **The canonicalisation route is not recorded per row.** Whether a stored domain came from
   its construct, from a reference extension, or was left short is known at parse time but is
   not a stored column; `dbd_source` records the annotation *method*, not the route.
8. **Cluster membership is order-dependent by construction.** Greedy clustering assigns a
   sequence to the first representative it matches, not the best. Deterministic, but not the
   unique optimal partition.

---

## References

Sources, with DOIs verified against Crossref, are listed in
[`REFERENCES.md`](REFERENCES.md). The accession survey, with per-accession citations and family
composition, is in [`UNIPROBE_ACCESSIONS.md`](UNIPROBE_ACCESSIONS.md). The admission policy is
stated in full, with its exclusion audit, in [`DOMAIN_POLICY.md`](DOMAIN_POLICY.md). Every
decision cited here, with the measurement behind it, is in [`DECISIONS.md`](DECISIONS.md).
