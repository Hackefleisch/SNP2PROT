# Kock et al. 2024 — screened, rescored, and excluded

**Decision, 2026-08-18, by the owner: the source does not enter the corpus.** The dataset stays
on its two source families, UniPROBE (18 accessions) and CIS-BP / Weirauch 2014, and PBM
acquisition is closed.

This report is the record of why. It replaces the two working reports produced during the
investigation and is written so the question does not have to be reopened from scratch: the
identifiers are verified, the yield is measured, and the one thing that failed is stated
precisely enough to recognise if it ever changes.

Everything downloaded for the investigation has been removed — the raw scans, the probe design
under `data/external/pbm_design/`, and the rescoring code. Nothing in `data/` or
`configs/thresholds.yaml` was touched at any point.

---

## 1. The source

Kock KH, Kimes PK, Gisselbrecht SS, … Bulyk ML. *DNA binding analysis of rare variants in
homeodomains reveals homeodomain specificity-determining residues.* **Nat Commun 15:3110
(2024)**, `10.1038/s41467-024-47396-0`, PMID 38600112, PMC11006913.

| what | verified value |
|---|---|
| PBM data | GEO `GSE233827` — 122 samples, platform `GPL33450`, AMADID #030236 (the ROG18A design) |
| processed 8-mer table | Harvard Dataverse `10.7910/DVN/FDQHCF`, CC0 — `SupplementaryData4.txt`, 581,717,822 B |
| construct sequences | Supplementary Data 3, `41467_2024_47396_MOESM6_ESM.xls`, 84,480 B |
| `10.5281/zenodo.10460649` | **not data** — the `upbm` R package (Kimes & Song) |

Supplementary Data 4 is 131 `(set, tf, allele)` groups × exactly 32,896 8-mers:
94 records under `This study` (93 clones, 27 reference alleles and 66 single-substitution
variants across 27 genes) and 37 under `Barrera et al., Science (2016)`.

## 2. What it would have added

Screened through this project's own admission policy — full Pfam-A at gathering thresholds,
`domains.call_domain`, `canonical.canonicalise`, then clustered against the live 5-edit
threshold:

| | |
|---|---|
| clones in Supplementary Data 3 | 93 |
| admitted | **87** (6 rejected: the whole PAX4 series, `PAX` + `Homeodomain` in one construct) |
| already held byte-identical | 20 |
| **genuinely new domains** | **67** |
| joining an existing cluster | 84 of 87; the 3 NKX2-6 domains would found one new cluster |
| clusters gaining domains | 25, of which **13 are singletons today** and would become variant-bearing |
| variants | 202 → ~269 (+33%) |
| homeodomain share | 29.7% → 33.3% |
| rows | +3,059,328 (+6.7%) |

Every one of the 66 variants is exactly one substitution from its own reference clone,
verified by alignment. The DNA axis is exact: 32,896 8-mers per allele, matching our stored
`dna_seq` strings directly, no reverse-complement remapping.

**This was a real gain and it is what makes the exclusion a cost, not a free win.** 13
singleton clusters becoming variant-bearing is the corpus's scarcest property.

## 3. Why it does not enter

### 3.1 The published scores are on a scale no fixed cutoff can cross

The deposit contains **no E-scores at any level**. GEO holds `.gpr` scanner files; the
per-sample `VALUE` and Supplementary Data 4 both carry `affinityEstimate`, a per-allele log
intensity, with `affinityQ` / `contrastQ` / `specificityQ` alongside.

`affinityEstimate` ranks agree well with our E-scores on shared proteins (Spearman 0.83–0.93),
but it has no transferable scale. Measured on four alleles held by both:

| allele | `affinityEstimate` at our `E = 0.45` | at our `E = 0.35` |
|---|---|---|
| SIX6-REF | 10.82 | 10.36 |
| SIX6-T165A | 10.93 | 10.38 |
| HOXD13-REF | 11.61 | 11.01 |
| PITX2-REF | 12.14 | 10.65 |

The positive boundary spans 10.82–12.14 and the negative boundary 10.08–11.16 — **one
protein's positive threshold is another's negative threshold**. Across all 130 alleles the
per-allele maximum runs 8.42 (`PROP1-F88S`, a variant that abolishes binding) to 15.70
(`ISX-REF`). `affinityQ` is a detection statistic, not a stringency one: at the paper's own
`Q < 10⁻⁶` it calls 6–13% of 8-mers positive against our 0.3%.

Admitting on those terms needed a second `score_type` and two new cutoffs, and would have
given this source a positive rate 20–40× every other source's — the label-prior leak the
project exists to avoid, wearing a different hat.

### 3.2 Recomputing E-scores from raw scans works, but not well enough

The alternative was to compute the corpus's own statistic from the raw scans. It was built and
validated: **32 arrays, 14 `BAR15A` alleles** (the whole ARX and HOXD13 series) rescored from
Barrera's own scans, which sit in the same GEO series, and compared with the UniPROBE values
already in `data/interim/BAR15A/`.

| | live alleles (11) | dead variants (3) |
|---|---|---|
| Spearman, 8-mers stored above 0.3 | 0.778 – 0.905 (median 0.851) | n/a |
| rank-matched positive set, Jaccard | 0.565 – 0.766 (median 0.639) | n/a |
| highest recomputed E-score | 0.448 – 0.493 | 0.325 – 0.427 |
| 8-mers at `E >= 0.45` | **0 – 32** | 0 |
| stored positives at `E >= 0.45` | **129 – 204** | 0, 0, 2 |

The ranking reproduces — 0.64 median rank-matched agreement sits between the corpus's
cross-source floor (0.458) and its within-source replicates (0.70–0.72, `reports/overlap.md`).
**The values do not.** The cutoff on the recomputed scale that selects as many 8-mers as the
stored source called positive is 0.312–0.343 across all 11 live alleles, median **0.321**.

So the route cost one recalibrated cutoff, near 0.32–0.35, derived rather than chosen. Better
than §3.1 by a wide margin, and it kept the property a rank rule destroys — the three dead ARX
variants top out at 0.33–0.43 and yield 0–2 positives. **Not good enough to admit a source
on:** it would have made the corpus mixed-provenance, 19 sources on published E-scores and one
on ours, with its own cutoff, and the residual disagreement with UniPROBE is not small.

Three findings from that work, in case it is ever revisited:

- **Keep saturated spots.** 274 spots on one ARX array are saturated and **260 of them are in
  the top 1,000 by signal**; masking them cost the top 8-mer 12 of its 18 probes. The E-score
  is a rank statistic, so a ceiling-clipped spot still carries the fact that matters.
- **Normalise Cy3 as a sequence-model residual, not a ratio.** Raw division drops agreement
  from 0.489 to 0.123: probe composition drives both dsDNA yield and homeodomain binding, so
  dividing by observed Cy3 removes real signal. Regressing Cy3 on mono- and dinucleotide
  composition and dividing by the unexplained part gives 0.595 → 0.621.
- **Join the design to a scan by `(Column, Row)`, never by probe ID.** The ID sets match
  exactly and name different probes; joining by ID gives rho = -0.003 and looks clean.

Two limits bound the method and neither is fixable from this deposit. **More replicates do not
restore the range** — five `HOXD13_REF` arrays give 11 positives against the stored 144 — most
likely because Barrera's masliner multi-power scan stitching cannot be reproduced from one
scan per array. And **array QC has no working gate**: upper-tail width separates a good slide
from a bad one (3.5–4.1 against 1.1–1.7 across eight chambers of each), but a variant that has
lost binding looks identical to a bad array by that measure, so filtering on it would delete
exactly the negative evidence the corpus most needs.

### 3.3 A third of the deposit is our own data

37 of the 131 records are labelled `Barrera et al., Science (2016)`: nine genes we already
hold as `BAR15A`, reprocessed through `upbm`. Not a second experiment — 1,217,152 rows of a
source agreeing with itself, which would have corrupted `reports/overlap.md`, the one report
that measures cross-source agreement. Any future attempt must parse `This study` only.

## 4. What the investigation left behind, and what it took away

**Kept, because they are defects in our own data that this work exposed** — both are open
tasks in `TODO.md`:

- **`T22`** — `BAR15A:SIX6` carries `Q6P051`, a 305 aa TrEMBL *"SIX6 protein (Fragment)"*,
  not reviewed `SIX6_HUMAN` `O95475` (246 aa). Our `full_seq` frame runs 59 residues off the
  literature's numbering because of it. Every other gene checked exactly: 17 variants across
  CRX, HESX1, MSX2, NKX2-5, PITX2, PROP1 and VENTX, all delta 0.
- **`T23`** — a stored `HOXD13` domain differs from `P35453` at one residue inside the ±10
  padding. Probably clone reality, never audited corpus-wide.

**Removed with this decision:** `T11` (acquire and parse Kock), `T19` (E-scores from raw
arrays), `T24` (parse-time gotchas held for `T11`), `src/snp2prot/rawpbm.py`,
`scripts/validate_escores.py`, `tests/test_rawpbm.py`, `data/external/pbm_design/` and its
`PROVENANCE.md` rows.

**`T20` is unaffected and now stands alone.** The cutoff question — absolute `E >= 0.45`
against a rank-matched rule — was raised by `reports/threshold_review.md` on the corpus we
already hold, not by this source. Nothing about excluding Kock settles it.

## 5. What would have to be true to revisit

1. The lab, or anyone, publishing E-scores for these arrays — the whole problem is that the
   deposit stops at `affinityEstimate`.
2. Multi-power scans per array appearing, which is the one input that could close the range
   gap in §3.2.
3. `T20` resolving in favour of a rank-matched rule at a fixed stringency, which would make a
   recomputed source's compressed scale irrelevant — the labels would be set by rate, not by
   value.

Short of those, the answer stands: 67 domains and 13 newly variant-bearing clusters are not
worth a second scale in a dataset whose entire value is that every row is comparable.
