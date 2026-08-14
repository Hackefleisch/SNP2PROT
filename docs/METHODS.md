# Methods

Construction of a paired transcription-factor DNA-binding domain × DNA-site dataset with
binary binding labels.

This document records what was done and why. It is written to be sufficient for an
independent reimplementation. Parameters quoted here are the values in
[`configs/thresholds.yaml`](../configs/thresholds.yaml); statistics are from the most recent
complete build (see §10).

---

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
policy is evaluated over a restricted set of 41 DNA-binding families. That set was derived
from the domain labels UniPROBE itself publishes, each resolved against Pfam-A `NAME` fields
by accession lookup; several UniPROBE labels predate Pfam renames (`Homeobox` → `Homeodomain`,
`Fork_head` → `Forkhead`, `MADS` → `SRF-TF`, `E2F_TDP` → `WHD_E2F_TDP`, `Zn2Cys6` → `Zn_clus`,
`BRIGHT` → `ARID`, `AP-2` → `TF_AP-2`, `PBX` → `PBC`, `RFX` → `RFX_DNA_binding`). Pfam hits
outside this set are recorded but do not affect admission: a DNA-binding domain adjacent to an
unrelated domain remains solely responsible for the DNA interaction.

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

From 30 accessions retrieved, 779 constructs carried a sequence. 479 distinct domains were
admitted; the remainder were either rejected under §3.2 or resolved to a domain sequence
identical to another construct's.

| property | value |
|---|---|
| rows | 16,316,416 |
| distinct DNA-binding domains | 479 |
| clusters | 424 |
| Pfam families | 31 |
| source organisms | 22 |
| distinct 8-mers | 32,896 (identical in every source) |
| binding / non-binding / excluded | 50,655 / 15,998,460 / 267,301 |
| negative:positive ratio | 316:1 |
| `dbd_seq` length | 43–216 aa (median 77) |

Composition is dominated by homeodomains (234 domains), followed by forkhead (54), bHLH (35),
zf-C4 (28), ETS (22), HMG box (21) and Zn₂Cys₆ (17), with a further 24 families represented by
ten or fewer domains each.

Protein-side depth is concentrated: 72 variants across 26 clusters, with 381 clusters
containing a single domain. Variants are predominantly homeodomain (54), with forkhead (9),
zf-C4 (5), bHLH (3) and paired-domain (1) representation.

## 10. Reproducibility

The dataset is rebuilt from the raw archives by a single command, and all parameters that
affect its content live in one configuration file. Reports summarizing binarization, cluster
inventory and validation are regenerated per source and version-controlled, so that a change
in the dataset is visible as a diff.

Statistics in §9 are from the build of 2026-08-13. Two subsequent method changes — the
alignment-based cluster distance of §7 and the recovery of one accession whose archive has a
flat layout — are implemented but not reflected in those figures; they increase variant and
domain counts and do not alter any other procedure.

## References

Sources, with DOIs verified against Crossref, are listed in
[`REFERENCES.md`](REFERENCES.md). The accession survey, with per-accession citations and
family composition, is in [`UNIPROBE_ACCESSIONS.md`](UNIPROBE_ACCESSIONS.md). The admission
policy is stated in full, with its exclusion audit, in
[`DOMAIN_POLICY.md`](DOMAIN_POLICY.md).
