# Brief: build a merged two-axis (TF-protein × DNA) binary binding dataset

> Source of record, pasted verbatim by the project owner on 2026-08-13. Where this brief's
> file paths differ from the repo, the repo wins — see the mapping table in `CLAUDE.md`.

---

## 0. What we are building and why

I am building a training set for an ML model that predicts whether a transcription factor
binds a DNA site, and that must be **sensitive to variation on both axes simultaneously**:

- **Protein axis** — DNA-binding domain (DBD) sequence, ideally many single-residue variants
  of the same WT domain, forming dense clusters.
- **DNA axis** — short binding-site sequences, ideally many single-base variants.

Constraints, all deliberate:

- **Label is binary.** Bind / does not bind. I set the thresholds; expose them as config,
  never hard-code them in the middle of a parser.
- **Negatives are wanted, in bulk.** A 50–100:1 negative:positive ratio is fine and expected.
- **DNA ≤ 20 bp.** Trim longer sequences; never pad a short one and pretend it has flanks.
- **DBD only** is fine for the protein side; full-length is not needed.
- **Related proteins.** The point is generalization across protein sequence space. If the
  proteins in the set are all distant from each other, the model becomes a lookup table and
  the whole exercise is worthless. Prioritize sources with point mutants and close paralogs.

Do not optimize for a model yet. Phase 0–6 is acquisition, parsing, and merging only.

---

## 1. Ground rules for you (Claude Code)

1. **Never invent an accession, URL, or file name.** Every one below came from literature and
   some are unverified. If a download 404s or the file layout differs from what I describe,
   **stop, report exactly what you found, and ask** — do not silently substitute a different
   file or a "similar" dataset.
2. **Log everything.** Maintain `PROVENANCE.md` with one row per downloaded file: source URL,
   accession, download timestamp, byte size, sha256, and a one-line description of what it
   contains. This is going into a publication; provenance is not optional.
3. **Raw stays raw.** `data/raw/<source>/` is append-only and never edited. All parsing writes
   to `data/interim/<source>/`. Merged output goes to `data/processed/`.
4. **One parser per source**, in `src/parsers/<source>.py`, each exposing
   `parse() -> pd.DataFrame` conforming to the unified schema in §3. No cross-source logic
   inside a parser.
5. **Checkpoint often.** Parquet, not CSV, for anything over ~10⁵ rows.
6. Some sources need **manual download** (license click-through, no direct link). When you hit
   one, write the exact steps into `data/raw/<source>/HOWTO.md` and move on to the next source
   rather than blocking.
7. If a dataset turns out to be unusable under our constraints, say so plainly and record why
   in `PROVENANCE.md`. A short honest catalogue beats a padded one.

---

## 2. Sources, in priority order

Tier 1 and 2 are the backbone. Tier 3 adds DNA-side depth. Tier 4 is held-out test data —
**do not merge Tier 4 into the training table.**

### Tier 1 — PBM 8-mer E-scores (best fit: real negatives, 8 bp DNA, point-mutant clusters)

- **UniPROBE** — `https://thebrain.bwh.harvard.edu/uniprobe/`
  - Barrera et al. 2016, accession **`BAR15A`** — 41 reference + 117 variant TF alleles,
    mostly single missense changes in the DBD. This is the single most valuable source we have:
    matched WT/mutant pairs over an identical 8-mer space.
  - Homeodomain panels: Berger et al. 2008 (~168 mouse HDs) and Noyes et al. 2008 (Drosophila).
    Look up their UniPROBE accessions on the site rather than guessing.
  - Also pull any forkhead / ETS / bZIP family panels you find — we want family clusters.
  - Each experiment yields all **32,896 non-redundant 8-mers** with an E-score and a z-score.
  - Note: probe-sequence downloads require accepting an academic-use license, and the Barrera
    universal-array *designs* need an MTA. The **8-mer score files** are the ones we need;
    confirm whether they are behind the same click-through and write a `HOWTO.md` if so.

### Tier 2 — C2H2 zinc-finger bacterial one-hybrid (densest protein point-mutation axis)

- **Persikov et al. 2015, NAR 43:1965** — `https://zf.princeton.edu/b1h/`
  - Randomized ZF libraries screened against **all 64 possible 3 bp targets**, with the variable
    finger in the F2 or F3 position of a three-finger protein whose flanking fingers have fixed
    known specificities.
  - Selections were run at **2 mM and 10 mM 3-AT** (low / high stringency) — two natural
    binarization levels. Keep the stringency as a column; do not collapse it.
  - The site is only 3 bp of variable target, but sits in a defined 9 bp three-finger context.
    Reconstruct and store the **full 9 bp site**, with the variable triplet's offset recorded.
- **Najafabadi et al. 2015, Nat Biotechnol** —
  `https://hugheslab.ccbr.utoronto.ca/supplementary-data/C2H2_B1H/`
  - `OLS.info` (zinc-finger sequences, ~4.6 MB) and `B1H.s_scores` (score matrix, ~12 MB),
    covering ~8,138 distinct natural C2H2-ZF domains.

⚠️ Negatives from both of these are **selection-absence**, not assayed-and-unbound. Tag them
distinctly (§3, `neg_provenance`) — I may down-weight or drop them later.

### Tier 3 — DNA-side depth (many single-base variants, WT proteins only)

- **SNP-SELEX / GVAT** — `https://renlab.sdsc.edu/GVATdb/`, raw reads ENA **`PRJEB9797`**,
  code `github.com/ren-lab/snp-selex`.
  - 270 TFs × 95,886 noncoding variants; oligos are 40 bp — **trim to ±9 around the variant
    for a 19 bp window**. Keep OBS and PBS as raw scores.
- **Codebook** (optional, for extra motif-centered windows): PWMs on Zenodo
  `10.5281/zenodo.15667805`; intervals GEO `GSE280248` (ChIP-seq), `GSE278858` (GHT-SELEX);
  raw SRA `PRJEB78913`, `PRJEB76622`, `PRJEB61115`. Genomic fragments are far longer than 20 bp,
  so this only enters via motif-centered extraction — treat as lower priority.

### Tier 4 — held-out quantitative test sets (DO NOT TRAIN ON THESE)

- **BET-seq** (Le et al. 2018, PNAS): Figshare article **`5728467`**, code
  `github.com/FordyceLab/BET-seq`. 2 proteins (Pho4, Cbf1) × ~10⁶ sequences of the form
  `NNNNN-CACGTG-NNNNN` = **16 bp**, with real ΔG. Perfect length fit, tiny protein axis.
- **STAMMP** (Aditham et al. 2021, Cell Systems): ~210 Pho4 variants × 9 oligos, >1,800 Kd.
  Supplementary tables + `fordycelab.com/data`; code `github.com/FordyceLab/ProcessingPack-STAMMP`.
- **MAX / k-STAMMP** (Hastings et al. 2025, Nat Commun 16:636): 240 single-aa substitutions in
  and around the MAX DBD × 7 motif variants. Source Data + Supplementary Data 1–3; OSF `osf.io/jmz8t`.

### Supporting metadata (needed by every parser)

- **CIS-BP** — `https://cisbp.ccbr.utoronto.ca/` for TF → DBD sequence and family assignment.
- **Pfam/InterPro + HMMER** for DBD boundary extraction where CIS-BP has no entry. Record the
  method used per protein in a `dbd_source` column — hand-curated and HMM-derived boundaries
  are not interchangeable.
- **UniProt** for canonical sequences and to resolve variant notation (e.g. `p.Arg233Cys`) into
  an actual mutated DBD string. Verify every variant maps to the residue the paper claims;
  off-by-one errors from isoform numbering are the most likely silent bug in this whole project.

---

## 3. Unified schema

Every parser emits exactly these columns:

| column | type | notes |
|---|---|---|
| `pair_id` | str | deterministic hash of (dbd_seq, dna_seq, assay, stringency) |
| `dbd_seq` | str | amino acids, DBD only, uppercase |
| `dbd_family` | str | from CIS-BP/Pfam (e.g. `Homeodomain`, `C2H2 ZF`, `bHLH`) |
| `dbd_source` | str | how the boundaries were determined |
| `wt_id` | str | cluster key: the WT/reference domain this variant belongs to |
| `n_mut_from_wt` | int | Hamming distance to `wt_id` sequence; 0 for the WT itself |
| `mut_positions` | str | comma-separated, 1-based within the DBD; empty for WT |
| `protein_id` | str | UniProt/gene identifier where known |
| `species` | str | |
| `dna_seq` | str | observed bases only, ≤20 bp, no padding characters |
| `dna_len` | int | |
| `dna_context` | str | `full` / `core_only` — was real flanking sequence observed? |
| `label` | int | 1 bind, 0 non-bind, `-1` for the discarded gray band |
| `raw_score` | float | E-score, s-score, OBS, ΔG… whatever the source gives |
| `score_type` | str | e.g. `pbm_escore`, `b1h_sscore`, `snpselex_obs`, `dG_kcal` |
| `threshold_pos` / `threshold_neg` | float | the cutoffs actually applied to this row |
| `assay` | str | `PBM`, `B1H`, `SNP-SELEX`, `MITOMI`, … |
| `stringency` | str | e.g. `2mM_3AT`, `10mM_3AT`, or empty |
| `neg_provenance` | str | `assayed_unbound` / `selection_absent` / `synthetic` — **required for label 0** |
| `source_dataset` | str | `BAR15A`, `persikov2015`, … |
| `source_file` | str | relative path under `data/raw/` |

Two hard rules:

- **`neg_provenance` is mandatory on every negative.** PBM low-E-score negatives are genuine
  evidence of non-binding; a sequence merely absent from a selection is not. If I later find
  these were pooled indistinguishably, the dataset is dead.
- **Do not pad `dna_seq` in the stored table.** Padding is a modeling decision (§5), not a
  parsing one. Store what was measured.

---

## 4. Binarization (config-driven, in `config/thresholds.yaml`)

Defaults to start with — all overridable:

- **PBM**: E-score ≥ **0.45** → 1; ≤ **0.35** → 0; in between → `-1` (excluded, but kept in the
  table so I can revisit). Binarize **per experiment**, never on pooled values — HK and ME array
  designs shift E-scores relative to each other.
- **B1H**: use recovery at the two 3-AT stringencies. Start with "recovered at 10 mM" → 1,
  "not recovered at 2 mM" → 0, and the middle → `-1`. Sanity-check this against the paper's own
  treatment before committing.
- **SNP-SELEX**: threshold on OBS for the bound call; keep PBS separately as the allele-difference
  signal. Report what fraction of oligos survive.
- **Tier 4**: convert Kd/ΔG to binary only at evaluation time, with the cutoff as a swept parameter.

Emit a `reports/binarization_summary.md` per source: n_pos, n_neg, n_gray, pos rate, and the
distribution of `raw_score` with the cutoffs drawn on it.

---

## 5. Merge and the two failure modes I care about

**5a. Length leaks assay identity.** PBM gives 8 bp, B1H 9 bp, SNP-SELEX 19 bp — and each source
has a different positive rate. A model can infer the assay from `dna_len` alone and shortcut the
task. Mitigation, implemented at featurization time and **not** in the stored table:

- centre every site in a fixed 20 bp window, pad with a masked neutral token, and pass an explicit
  `dna_context` flag so the model knows padding is absent-evidence rather than observed base;
- **and** produce a second, stricter variant of the dataset restricted to a common 8–10 bp core
  across all sources. Build both. I want to compare.

**5b. The lookup-table risk.** Report these before any model work:

- **Cluster inventory** (`reports/cluster_inventory.md`): per `wt_id`, the number of variants,
  the distribution of `n_mut_from_wt`, and per family, the number of clusters. Flag any family
  that contributes >40% of rows — C2H2 will otherwise swamp everything.
- **Split protocol**, implemented in `src/splits.py`:
  - *leave-one-variant-out within cluster* — tests DNA-side and fine protein-side generalization;
  - *leave-one-cluster-out* and *leave-one-family-out* — tests real protein-space generalization.
  Both must be produced; report metrics on both, always separately.
- **Nearest-neighbour baseline** (`src/baselines/nn_lookup.py`): for a held-out DBD, copy the
  binary profile of the most sequence-similar training DBD (% identity over the aligned DBD).
  This is the number any future model must beat under leave-one-cluster-out. If it doesn't, we
  have a lookup table with extra steps, and I want to know that early and cheaply.

Also emit `reports/overlap.md`: how many (dbd_seq, dna_seq) pairs appear in more than one source,
and whether their labels agree. Cross-source disagreement rate is a direct read on label noise.

---

## 6. Phases

- **Phase 0** — repo scaffold, `PROVENANCE.md`, `config/thresholds.yaml`, schema dataclass +
  a validator that every parser output must pass. No downloads yet.
- **Phase 1** — UniPROBE / Barrera `BAR15A`. Get one source end-to-end through the validator
  before touching a second. Report cluster inventory for it alone.
- **Phase 2** — remaining UniPROBE family panels.
- **Phase 3** — Persikov B1H + Najafabadi C2H2.
- **Phase 4** — SNP-SELEX, trimmed.
- **Phase 5** — merge, overlap report, splits, NN baseline.
- **Phase 6** — Tier 4 test sets, kept in a physically separate directory.

Stop and check in with me at the end of each phase with the summary reports. Don't run ahead to
modeling — Phase 6 is the end of this brief.

---

## 7. Open questions I have not decided (flag, don't resolve)

- Whether to keep the B1H data at all, given its weak negatives.
- Whether to restrict the final set to DBDs with ≥5 variants (sharper clusters, much smaller set).
- Whether the padded-20 bp or common-core variant becomes primary.

If you hit evidence that bears on any of these while parsing, surface it rather than deciding.
