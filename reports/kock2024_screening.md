# Kock et al. 2024 — acquisition screening (`T11`, exploratory pass)

Generated 2026-08-18. **Nothing has been downloaded into `data/raw/` and no code has been
written.** This is the screening `T11` step 1 asks for: what is actually in the deposit, what
it would cost to admit, and every issue that would arise. The decision is the owner's.

All work happened in a scratch directory. Every identifier below was resolved and the file
behind it opened; nothing is repeated on the strength of the earlier research session.

---

## 1. Identifiers, verified — and one correction

| what | verified value | how |
|---|---|---|
| paper | Kock KH, Kimes PK, Gisselbrecht SS, … Bulyk ML. *DNA binding analysis of rare variants in homeodomains reveals homeodomain specificity-determining residues.* **Nat Commun 15:3110 (2024)** | Crossref + Europe PMC |
| DOI | `10.1038/s41467-024-47396-0` | Crossref |
| PMID / PMCID | `38600112` / `PMC11006913` | Europe PMC |
| PBM data | **GEO `GSE233827`** — 122 samples, platform `GPL33450` "Universal PBM 8x60k (Bulyk lab design)", AMADID #030236 | the paper's Data availability, then the GEO record |
| processed 8-mer data | Harvard Dataverse **`10.7910/DVN/FDQHCF`**, CC0 — `SupplementaryData4.txt` (581,717,822 B, md5 `e507f327a4af871715697a46beaf23f4`) and `SourceData.xlsx` (389,135,137 B) | Dataverse API |
| construct sequences | **Supplementary Data 3** = `41467_2024_47396_MOESM6_ESM.xls`, 84,480 B | Europe PMC / BioStudies `S-EPMC11006913` |

**Correction to `TODO.md` `T11`:** the reported Zenodo deposit `10.5281/zenodo.10460649` is
**not data**. DataCite resolves it to *"upbm (R package)"* by Kimes & Song — the analysis
software, `IsSupplementTo github.com/pkimes/upbm`. The paper's Data availability names only
GEO `GSE233827` and the Dataverse DOI. **`GSE233827` is missing from `T11` entirely**, and it
is the more important of the two.

`T11` was right about `Nat Commun 15:3110`, right about the Dataverse DOI, right about the
array design (AMADID #030236, the ROG18A design), and right that constructs are DBD + 15 aa
flanks published in Supplementary Data 3.

## 2. What the deposit actually contains

**Supplementary Data 4** (= Dataverse `SupplementaryData4.txt` = `MOESM7`, same byte count) is
one tab-separated table, 4,309,377 lines:

```
set  tf  allele  seq  affinityEstimate  affinityQ  contrastDifference  contrastQ  contrastResidual  specificityQ
```

- 131 `(set, tf, allele)` groups × **exactly 32,896 8-mers each**, no exceptions;
- **`set` has two values** — `This study` (94 records) and `Barrera et al., Science (2016)`
  (37 records). §5 is about the second one;
- `contrastQ` / `specificityQ` are `NA` for all 1,217,152 Barrera-set rows.

**Supplementary Data 3** holds **93 clone sequences** for the new alleles only — 27 reference
alleles and 66 single-substitution variants across 27 genes, with a mutagenesis-primer sheet.
No sequences are published for the Barrera set. (94 SD4 records against 93 clones because
`NKX2-6-REF` was assayed in two series and appears as `NKX2-6` and `NKX2-6series2`.)

**GEO `GSE233827`** deposits `.gpr` scanner files only — four per sample, each covering eight
chambers of one slide. The GSM-level `VALUE` is defined as *"affiniteEstimate, an 8-mer-level
affinity score"*: the same statistic as SD4, not a second one.

> **There are no E-scores anywhere in this deposit**, at any level.

## 3. Protein axis — yield through the admission policy

All 93 clones were scanned against full Pfam-A at gathering thresholds and put through
`domains.call_domain` and `canonical.canonicalise` exactly as a parser would.

| | count |
|---|---|
| clones in Supplementary Data 3 | 93 |
| **admitted** | **87** |
| rejected `mixed_families` | 6 — the whole PAX4 series |
| canonical domains byte-identical to a domain we already hold | 20 |
| **genuinely new domains** | **67** |

- **PAX4 is out.** Its 234 aa clone carries `PAX:1-121` *and* `Homeodomain:164-219` — the
  PAX + homeodomain rejection already named in `CLAUDE.md`. This removes the one gene `T11`
  step 2 called "absent" from the corpus, and it takes 5 variants with it (F177V, F177Y,
  R200K, R200Q, R227W).
- **Every one of the 66 variants is exactly one substitution from its own REF clone**,
  same length, no indels. Verified by alignment, not assumed.
- Flanks are 15/15 in 43 of 87 clones but range 8–23. Three HOXB9 clones carry only 8 aa of
  C-flank, so the ±10 window cannot be filled from the construct (`origin=construct_short`);
  they need the UniProt reference, which is what `canonical` already does.
- The Pfam envelope can start one residue apart in a REF and its own variant (HOXB9-REF at 6,
  HOXB9-R186Q at 7). Harmless because `dbd_seq` is canonical, but it is the reason it must be.

## 4. Cluster impact — mostly depth, almost no breadth

Each canonical domain was compared to all 340 Homeodomain cluster representatives with the
project's own aligner and the live 5-edit threshold.

- **84 of 87 join an existing cluster**; the 3 NKX2-6 domains are 10–11 edits from anything and
  found **one new cluster**.
- **25 existing clusters gain domains**, 64 genuinely new ones between them.
- **13 of those 25 are singletons today** — `Cell08:Alx3`, `Hoxa2`, `Hoxa4`, `Hoxa7`, `Hoxb6`,
  `Hoxb9`, `Hoxc9`, `Hoxc10`, `Alx4`, `Cart1`, `Bapx1`, `Prrx1`, `weirauch2014:SHOX` — and
  would become variant-bearing.

| | now | after |
|---|---|---|
| domains | 1,335 | ~1,402 |
| variants | 202 | ~269 (+33%) |
| clusters holding a variant | 122 | ~136 |
| homeodomain share of domains | 29.7% | **33.3%** |
| rows | 45,495,168 | +3,059,328 (+6.7%) |

**The homeodomain-imbalance worry recorded in `docs/DECISIONS.md` is stale.** It projected
47.9% → ~58%; that was before CIS-BP landed 868 domains. The real cost today is 3.6
percentage points.

The 13 SIX6 alleles all sit 4–5 edits from `C:weirauch2014:six3`, i.e. right at the threshold.
They join at `max_edits: 5` and would found their own cluster at 4 — worth knowing before the
inventory is read as stable.

## 5. The Barrera set is our own data, reprocessed

37 of the 131 SD4 records are labelled `Barrera et al., Science (2016)`: ARX, HESX1, HOXD13,
ISX, NKX2-5, PBX4, PITX2, PROP1, SIX6 — 9 genes we already hold as `BAR15A`. These are not a
second experiment. They are Barrera's raw arrays pushed through `upbm`.

**Admitting them would put two rows on every `(domain, 8-mer)` pair that look like two sources
and are one experiment.** `reports/overlap.md` exists to measure cross-source agreement and
the 45.8% noise floor comes out of it; 1,217,152 rows of a source agreeing with itself would
corrupt exactly that number. **Parse `This study` only.**

Kept *outside* the training table, though, that set is unusually valuable: it is the same raw
data scored two ways, which measures pipeline noise directly (§6).

## 6. The blocker — no cutoff reproduces our label convention

Four alleles are held by both us and Kock with a byte-identical `dbd_seq`, so E-score and
`affinityEstimate` can be compared 8-mer by 8-mer. **All 32,896 8-mers matched our stored
`dna_seq` directly, with no reverse-complement remapping** — the DNA axis is exact, which
closes `T11` step 3.

Rank agreement is high — Spearman **+0.83 to +0.93** — so the two statistics are measuring the
same thing. What is missing is a scale.

**A global cutoff on `affinityEstimate` cannot work.** The value sitting at our `E = 0.45`
boundary, across four proteins:

| allele | at `E = 0.45` | at `E = 0.35` |
|---|---|---|
| SIX6-REF (Kock's own array) | 10.82 | 10.36 |
| SIX6-T165A (Kock's own array) | 10.93 | 10.38 |
| HOXD13-REF | 11.61 | 11.01 |
| PITX2-REF | 12.14 | 10.65 |

The positive boundary spans 10.82–12.14 and the negative boundary spans 10.08–11.16 across the
Barrera-set and own-array versions: **the ranges overlap, so one protein's positive threshold
is another's negative threshold.** Over all 130 alleles the per-allele maximum runs 8.42
(`PROP1-F88S`) to 15.70 (`ISX-REF`) — a dead variant's entire distribution lies below a strong
binder's floor. That is real biology, and it is why the E-score is rank-based to begin with.

**`affinityQ` is a detection statistic, not a stringency one.** At the paper's own `Q < 10⁻⁶`
it calls 1,939–4,439 8-mers positive per allele (6–13%), against ~100 (0.3%) under `E > 0.45`.
Our positives are a strict subset — containment 1.000 in every case tested — so the two rules
are nested, not contradictory. But a `label` column mixing 0.3% and 10% positive rates hands a
model the source identity for free, which is the `dna_len` leak wearing a different hat.

**A rank-matched rule is the only one that lands in the right place.** Taking the top *N* by
`affinityEstimate`, with *N* the count our E-score rule gives:

| allele | vs Kock's reprocessing of Barrera | vs Kock's own array |
|---|---|---|
| SIX6-REF | 0.70 | 0.52 |
| SIX6-T165A | 0.79 | 0.68 |
| HOXD13-REF | 0.51 | 0.59 |
| PITX2-REF | 0.78 | 0.75 |

(Jaccard against our BAR15A positives.) Measured against the corpus's own noise floor —
70–72% within-source replicates, 45.8% median cross-source — **this is as close as two PBM
sources ever get here.** The left column is the same raw data scored twice, so 0.51–0.79 is
the price of the pipeline alone.

The catch: a pure quantile rule would manufacture ~100 positives out of background for
`PROP1-F88S`, a variant that abolishes binding — destroying the exact signal the dataset
exists to carry. So an admission rule needs **two parameters**: a per-allele significance gate
(`affinityQ`) *and* a stringency calibration (rank or rate). Neither has a published
convention, and both land in `configs/thresholds.yaml` with the blast radius that implies
(`docs/DECISIONS.md`, 2026-08-14).

### Computing real E-scores ourselves — nothing is blocked

The first version of this report called that route blocked, then licence-gated. **Neither is
true.** Every input is public and now in hand:

| input | status |
|---|---|
| raw scanner output for every Kock array | GEO `GSE233827`, four `.gpr.gz` per sample (~3 MB), one array chamber per file |
| raw scanner output for the *Barrera* arrays | in that same series, dated `151124`-`151202` — the scans behind the `BAR15A` E-scores we already hold |
| probe design for AMADID #030236 | GEO `GPL34105` supplementary, 1.2 MB, **downloaded** to `data/external/pbm_design/` |
| the E-score algorithm | Berger et al. 2006; ~80 lines, prototyped, ~1 s per array |

The one real trap is the join key, and it cost a wrong conclusion in the previous version of
this section. A `.gpr` names its spots `dBr_11450_v_Jan07`; the design names them
`dBr_11450_Jan07`. **All 41,944 de Bruijn IDs match as a set, and they are not the same
probes** — 35 of 62,976 grid positions agree. Joining by ID yields a clean-looking table of
noise (rho = -0.003 against UniPROBE). Joining by `(Column, Row)` is correct, and the proof is
biological rather than statistical:

| array | check | result |
|---|---|---|
| CREB1 (`GSM8024474`, 2023) | probes carrying the CRE `TGACGTCA` | **18.6x** the array median (0.73x when joined by ID) |
| CREB1 | top 8-mers by recomputed E-score | `TGACGTCA`, `ATGACGTA`, `TGACGTAA` — the canonical CRE |
| HOXD13-I297V (`GSM7437471`, 2015, a Barrera-era scan) | top 8-mers | `CCATAAAA`, `CCCATAAA`, `CTCATAAA`, `CCAATAAA` |

One design covers both eras: Kock's arrays and Barrera's are the same v14.

**What is left is preprocessing, not access.** The prototype — background-subtract, drop
saturated spots, score the top half of the ranked list — reaches Spearman **0.49** against
UniPROBE's published E-scores for the same array, motif clearly right, dynamic range
compressed (max 0.457 against 0.495). The gap is Cy3 normalization against the
double-stranding control array (public, one per sample in the same series), spatial
detrending, spot masking, and combining replicate arrays. Details and the measurements are in
[`data/external/pbm_design/README.md`](../data/external/pbm_design/README.md).

**The validation loop needs nothing we do not have.** Recompute `BAR15A` from its raw scans,
compare against `data/interim/BAR15A/`, and point the code at Kock only once it reproduces
them. If it does, **Kock enters on the existing scale under the existing cutoffs** — no new
`score_type`, no new thresholds, nothing invalidated. That is `TODO.md` `T19`.

## 7. Everything else that would bite

1. **`HOXC9-K195R` is byte-identical to `Cell08`'s mouse `Hoxc9` wild type.** One `dbd_seq`,
   two identities, two labels from two assays — `T15` with a variant on one side of it. 9 of
   the 20 exact duplicates are human-Kock against mouse-`Cell08`, so `T15`'s species conflict
   goes from one case to ten.
2. **Our SIX6 full-length is the wrong protein.** `BAR15A:SIX6_*` carries `Q6P051`, a 305 aa
   TrEMBL *"SIX6 protein (Fragment)"*, not reviewed `SIX6_HUMAN` `O95475` (246 aa). Kock's
   variant numbering lands 59 residues off ours because of it, and Kock is right. Every other
   gene checks out exactly: 17 variants across CRX, HESX1, MSX2, NKX2-5, PITX2, PROP1 and
   VENTX all land on the residue the label names, delta 0. This is a defect in **our** table
   (`full_seq` only, not `dbd_seq`), found by this screening.
3. **One Kock clone disagrees with UniProt outside the domain.** `HOXD13-REF` has `T` where
   `P35453` has `D`, at protein position 261 — 6 residues N-terminal of the Pfam start and so
   *inside* our ±10 padding. Our `BAR15A` domain carries the same `T`, so both labs used one
   clone source; but the stored `dbd_seq` differs from the reference proteome at one residue.
4. **SD4 is already replicate-collapsed** — one value per `(allele, 8-mer)`. `per_experiment`
   binarization and `reconcile_replicates` have nothing to work on, and no within-source
   replicate agreement can be computed for this source. Replicate counts survive only as
   metadata (1 for 4 alleles, 2 for 52, 3 for 27, 4 for 9, 6 for 1). Four alleles rest on a
   single array.
5. **Supplementary Data 3 is legacy `.xls`**, which needs `xlrd`; the project reads `.xlsx`
   through `openpyxl` and has no `.xls` reader. A new dependency for one 84 KB file.
6. **`T11`'s cross-lab replicate hope was aimed at the wrong genes.** Barrera's VENTX variants
   are `E101K` and `R143C`, not an Arg→Gln, and Kock's PROP1 variant is `R99H`, not `R99Q`.
   The real cross-lab replicates are **`SIX6-H141N` and `SIX6-T165A`** — assayed by both labs,
   byte-identical domains, and the corpus's first cross-lab replicate of a *variant*. Two
   more are same-position-different-substitution: CRX `R41` (Barrera `R41Q`/`R41W`, Kock
   `R41L`) and PROP1 `R99` (Barrera `R99Q`, Kock `R99H`).
7. **`PROP1 R112Q` stays unresolved** (`docs/DECISIONS.md` §4) — Kock does not assay it.
8. Licensing is clean: Dataverse CC0, GEO public, no click-through. 581 MB into `data/raw/`,
   which is git-ignored.

## 8. The options

| | what | cost | gets |
|---|---|---|---|
| **A** | Parse `This study` only, new `score_type`, positives = `affinityQ` significant **and** top-*k* per allele | two new cutoffs in `thresholds.yaml`; label semantics differ from every other row; invalidates dataset + reports + provenance by this project's own convention | 67 domains, +33% variants, 13 singletons become variant-bearing |
| **B** | As A, but calibrate *k* so this source's positive **rate** matches the corpus PBM rate | same, minus one free parameter, plus an explicit "rate, not threshold" claim in METHODS | same |
| **C** | **Compute E-scores from the raw `.gpr`s** | a probe-level module (~150 lines with the normalization steps), validated against `BAR15A` before use. No new data access — the design is downloaded and everything else is public | the same data **on the existing scale, under the existing cutoffs** — no config change, nothing invalidated |
| **D** | Drop | none | nothing; PBM acquisition closes empty |

**Recommendation: C.** It was ranked third in the first version of this report on the strength
of a blocker that does not exist. Everything it needs is public and downloaded, the algorithm
is ~80 lines and prototyped, and — decisively — the *validation* is free, because `GSE233827`
deposits the raw arrays behind `BAR15A`, whose UniPROBE E-scores we already hold. The code
either reproduces them or it does not, and we find that out before a single new row is
written. C also costs nothing in config churn: A and B both put a second label convention into
a corpus whose entire argument rests on having one.

Fall back to B only if the normalization work fails to reproduce UniPROBE.

Three things are worth doing regardless of the decision: record the
same-raw-data-two-pipelines agreement (§6) in `docs/METHODS.md` as an independent read on the
noise floor, fix the `SIX6` full-length accession (§7.2), and read
[`threshold_review.md`](threshold_review.md) — the cutoff question this screening raised turns
out to have a defect of its own behind it, in `Cell09`.
