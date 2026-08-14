# SNP2PROT — working notes for Claude

## What this project is

A **training set** for a model that predicts whether a transcription factor binds a DNA
site, built so that it varies on **both axes at once**:

- **protein axis** — DNA-binding domain (DBD) sequence, with many single-residue variants of
  the same WT domain forming dense clusters;
- **DNA axis** — binding sites ≤20 bp, with many single-base variants.

The label is **binary**. Negatives are wanted in bulk (50–100:1 neg:pos is expected, not a
class-imbalance bug to fix). Full-length protein is not needed; DBD only.

The full specification is [docs/TFDNA_MERGE_BRIEF.md](docs/TFDNA_MERGE_BRIEF.md) — the
owner's brief, stored verbatim. **Read it before touching anything under `src/snp2prot/`.**
Where it and this file disagree about *paths*, this file wins (see the mapping below); where
they disagree about *intent*, the brief wins.

**Current phase: 2 complete. 18 UniPROBE PBM sources parsed and validator-clean —
16,645,376 rows, 489 domains in 425 clusters, 30 Pfam families, 24 organisms, every protein
scored against the same 32,896 8-mers. Phase 3 is DROPPED (see below); Phase 4 is next.**

**Read [docs/METHODS.md](docs/METHODS.md) first** — it records how the dataset was built and
why, in enough detail to reimplement. [docs/DOMAIN_POLICY.md](docs/DOMAIN_POLICY.md) states
the admission rules; [reports/OPEN_ITEMS.md](reports/OPEN_ITEMS.md) tracks decisions.

> ### `dbd_seq` is the PADDED Pfam domain
> Not the bare Pfam envelope, and not the sequence that was on the array. It is the Pfam
> envelope plus 10 residues each side. Three BAR15A variants mutate residues 2-9 residues
> N-terminal of the Pfam start, and without the padding each collapses onto its own wild type
> with a different label. **Read [docs/DOMAIN_POLICY.md](docs/DOMAIN_POLICY.md) before using
> the column.** Settings live in the `domain:` block of `configs/thresholds.yaml`.

## The admission policy — three conditions

A construct enters the dataset only if the stored `dbd_seq` satisfies all three. They follow
from what the data is for: predict a structure, find the centre of the DNA-contacting
residues, embed the residues around it.

1. **Sole responsibility** — that subunit alone produced the measured interaction.
2. **One continuous region** — one contiguous stretch, not fragments scattered through the
   protein.
3. **The variation lies inside it** — every change that alters binding falls within the
   stored region, or a sequence model sees identical inputs with different labels.

Enforced in `snp2prot.domains`, which rejects mixed-family constructs (PAX+homeodomain,
POU+homeodomain), C2H2 zinc-finger arrays, constructs with no Pfam hit, and any variant whose
mutation escapes the padded window. Rejections are reported, never silently repaired.

## The failure this project exists to avoid

If the DBDs in the set are all distant from each other, a model can memorize per-protein
profiles and score well without learning anything — a lookup table with extra steps. Two
concrete consequences that shape the code:

1. **Never pool sources without checking `dna_len`.** PBM gives 8 bp, B1H 9 bp, SNP-SELEX
   19 bp, and each has a different positive rate, so length alone leaks assay identity and
   with it the label prior. Mitigated at *featurization* time, never in the stored table.
2. **The nearest-neighbour baseline is the bar.** `snp2prot.baselines.nn_lookup` copies the
   binding profile of the most similar training DBD. Any model that does not beat it under
   leave-one-cluster-out has learned nothing transferable.

## Layout

```
docs/TFDNA_MERGE_BRIEF.md   the owner's spec, verbatim — source of record
docs/REFERENCES.md          dataset -> paper -> Crossref-verified DOI, plus reading order
docs/UNIPROBE_ACCESSIONS.md all 36 UniPROBE accessions, citations, family survey
docs/DOMAIN_POLICY.md       what dbd_seq is, the padding, and what gets excluded  <- READ THIS
docs/METHODS.md             publication-quality account of how the dataset was built
docs/papers/                paper PDFs (git-ignored); README.md there is the manifest
docs/papers_inbox/          the owner drops papers here; Claude identifies and files them
PROVENANCE.md               one row per raw file: URL, accession, timestamp, size, sha256
configs/thresholds.yaml     the ONLY config file; all binarization cutoffs live here
reports/                    generated markdown, committed; the phase-boundary deliverables

src/snp2prot/
  schema.py       unified 22-column row schema + validate(); the gate every parser passes
  thresholds.py   reads configs/thresholds.yaml
  config.py       every path in the project; nothing builds a path by hand
  domains.py      Pfam/HMMER annotation + the three-condition admission policy
  proteins.py     the protein-side companion table (bare/padded domain, construct, full-length)
  parsers/        one module per source, each exposing parse() -> pd.DataFrame
    _uniprobe.py  machinery shared by every UniPROBE accession (formats differ per accession)
  metadata/       CIS-BP / Pfam+HMMER / UniProt lookups shared by all parsers
  reports.py      binarization summary, cluster inventory, overlap report
  splits.py       leave-one-variant / -cluster / -family out
  baselines/      nn_lookup — the number a model must beat
  evaluation/     metrics for scoring the baseline per split regime; Phase 5
  data/           featurization (padded20 vs common_core); Phase 6+

data/raw/<source>/       append-only, never edited          (git-ignored)
data/interim/<source>/   per-source parsed Parquet          (git-ignored)
data/interim/proteins/   protein table: domain at 4 levels + UniProt full-length
data/external/pfam/      Pfam HMMs for boundary annotation
data/external/uniprot/   cached canonical sequences
data/processed/          merged training table              (git-ignored)
data/testsets/           Tier 4 held-out sets, kept physically apart  (git-ignored)
```

### Path mapping vs. the brief

The brief was written against a generic layout. These are the same things under this repo's
names — **do not create the left-hand paths**, that would fork the structure in two.

| brief says | this repo uses |
|---|---|
| `config/thresholds.yaml` | `configs/thresholds.yaml` |
| `src/parsers/<source>.py` | `src/snp2prot/parsers/<source>.py` |
| `src/splits.py` | `src/snp2prot/splits.py` |
| `src/baselines/nn_lookup.py` | `src/snp2prot/baselines/nn_lookup.py` |
| Tier 4 "separate directory" | `data/testsets/` |

## Rules that are not negotiable

1. **Never invent an accession, URL, or file name — or a DOI.** Several in the brief are
   unverified, and one plausible-looking DOI for Persikov 2015 is a different paper. If
   a download 404s or a file's layout differs from the description, **stop and report exactly
   what was found** — do not substitute a "similar" dataset, and do not guess a URL pattern.
2. **Raw stays raw.** `data/raw/` is append-only. Parsers read from it and write to
   `data/interim/`. Nothing edits a raw file, ever.
3. **Every raw file gets a `PROVENANCE.md` row** — URL, accession, UTC timestamp, bytes,
   sha256, one-line description. Use `scripts/record_provenance.py`. This is publication
   evidence, not a convenience log.
4. **`neg_provenance` is mandatory on every `label == 0` row.** A PBM low-E-score negative is
   real evidence of non-binding; a sequence merely absent from a B1H selection is not. If
   these are ever pooled indistinguishably, the dataset is dead. The validator enforces it.
5. **No cutoff literals in parser code.** Read `snp2prot.thresholds`; write the applied values
   into each row's `threshold_pos` / `threshold_neg` so a table is self-describing.
6. **Never pad `dna_seq` in the stored table.** Padding is a modeling decision. The validator
   rejects non-ACGT characters for this reason.
7. **One parser per source, no cross-source logic inside it.** Merging, dedup and
   reconciliation live in the merge step, not in a parser.
8. **Blocked ≠ stuck.** A source behind a licence click-through or MTA gets a
   `data/raw/<source>/HOWTO.md` with exact manual steps, a row in the PROVENANCE manual queue,
   and then you move to the next source. Do not burn effort brute-forcing a download.
9. **Stop at each phase boundary** and hand back the reports. Do not run ahead to modeling.
10. **Never run a command that takes more than a few minutes.** Hand over the exact command
    instead, with how long it takes and what to check in the output. `scripts/build_dataset.py
    --all` is ~25 min and is always the owner's to run; single small sources and test subsets
    are fine to execute here.
11. **Never `git commit` unsolicited.** The owner reviews work as an uncommitted diff;
    committing removes the review surface. Finish, run tests and lint, leave the tree dirty.

## Conventions

- **One config file: `configs/thresholds.yaml`.** There is no `default.yaml` and no per-run
  config, because at this stage a "run" *is* a dataset build and the thresholds are the only
  thing that varies. When modeling arrives (post-Phase 6) its hyperparameters get their own
  file — deliberately not merged into this one, because changing a threshold invalidates the
  dataset, every report and every provenance row, while changing a learning rate does not.
  Two things with different blast radii do not belong in one file.
- **Paper PDFs go in `docs/papers/`**, git-ignored, named `<firstauthor><year>_<slug>.pdf`
  (`_supp`, `_supp-<what>`, `_fig<N>` for the rest). The manifest is `docs/papers/README.md`;
  `scripts/check_papers.py` reports what is missing, unlisted, or waiting in the inbox.
- **The owner adds papers by dropping them in `docs/papers_inbox/`** under whatever name they
  downloaded with. Filing them is Claude's job, and it is done by **reading each file's
  embedded metadata and first page — never by trusting the download name.** Elsevier `mmc`
  numbers have twice turned out to be main articles rather than supplements, and the two 2008
  Cell homeodomain papers have PIIs that run opposite to their DOI order. After filing, add a
  manifest row and leave the inbox empty.
- **Parquet, not CSV**, for anything over ~10⁵ rows. Checkpoint often.
- Parser output goes through `schema.coerce(df)`, then `schema.validate(df, source)`, and the
  report's `raise_if_failed()` before any write to `data/interim/`.
- `pair_id` and `dna_len` are **derived** — `coerce()` computes them. Never hand-set them.
- **Cluster distance is alignment-based, so variants may carry insertions and deletions.**
  `snp2prot.align.edit_profile` aligns a variant to its reference (BLOSUM62, affine gaps) and
  counts every edited residue. **Terminal gaps are free**: `dbd_seq` is the padded envelope
  *clipped* where the construct ends, so two versions of one domain differ at the termini for
  reasons that have nothing to do with the protein — in ROG18A every bare domain is 83-85 aa
  while the padded ones run 95-105 aa. Internal indels still count. Settings live in the
  `cluster:` block of `configs/thresholds.yaml`.
- `mut_positions` is 1-based **in the reference's frame**, not the variant's, so every member
  of a cluster shares one coordinate system and maps onto one predicted structure. A variant
  with a deletion can therefore name a position past its own length. One entry is emitted per
  edited residue, so `len(mut_positions) == n_mut_from_wt` always holds. Isoform-offset
  off-by-ones remain the most likely silent bug — every variant must land on the residue the
  paper claims.
- Type hints on public functions; `from __future__ import annotations` at module top.
- Line length 100; `ruff check` + `ruff format` clean before a commit.

## Commands

```bash
uv venv --python 3.11 .venv && uv pip install --python .venv -e ".[dev]"
.venv/bin/python -m pytest -q
.venv/bin/python scripts/build_dataset.py --all            # parse -> validate -> interim
.venv/bin/python scripts/make_reports.py --source BAR15A    # regenerate reports/
.venv/bin/python scripts/audit_sources.py                  # invariant sweep, ~4 min
.venv/bin/python scripts/build_protein_table.py            # protein-side companion table
.venv/bin/python -m ruff check . && .venv/bin/python -m ruff format .
.venv/bin/python scripts/record_provenance.py data/raw/<source>/<file> --url ... --desc ...
```

Note: system Python 3.10 has no `ensurepip`, so `python -m venv` produces a venv without pip.
Use `uv` (already installed at `~/.local/bin/uv`).

## Phases

| phase | scope | status |
|---|---|---|
| 0 | scaffold, `PROVENANCE.md`, thresholds config, schema + validator | **done** |
| 1 | UniPROBE / Barrera `BAR15A` end-to-end, cluster inventory for it alone | **done** |
| 2 | remaining UniPROBE family panels | **done** — Cell08, EMBO10, PNAS13, then SCI09, GR09, MAR17A, SHO18A, ROG18A to rebuild breadth after the policy. Survey in `docs/UNIPROBE_ACCESSIONS.md`; `GB11` (27 bHLH) is the best remaining candidate. |
| 3 | ~~Persikov B1H + Najafabadi C2H2~~ | **DROPPED** — C2H2 arrays fail condition 2; Persikov varies a different subunit than the one that binds. ~8,000 domains excluded. |
| 4 | SNP-SELEX, trimmed to a 19 bp window | next — screen every TF through the domain policy; many of its 270 are C2H2 and will be rejected |
| 5 | merge, overlap report, splits, NN baseline | |
| 6 | Tier 4 test sets, in `data/testsets/` | end of brief — bHLH dimers; MAX substitutions are "in and around" the DBD, so condition 3 needs checking per variant |

## Deferred to the owner — flag, do not resolve

Tracked in [reports/OPEN_ITEMS.md](reports/OPEN_ITEMS.md). Currently: whether to keep B1H at
all given its weak negatives; whether to require ≥5 variants per DBD; whether padded-20 bp or
common-core becomes the primary dataset. If parsing turns up evidence bearing on any of these,
add it to that file rather than acting on it.
