# SNP2PROT

A merged, binary-labelled dataset of transcription-factor **DNA-binding domain × DNA site**
interactions, built to vary on both axes simultaneously: dense clusters of single-residue DBD
variants against short (≤20 bp) binding sites with single-base variants.

Specification: [docs/TFDNA_MERGE_BRIEF.md](docs/TFDNA_MERGE_BRIEF.md).
Agent working notes and conventions: [CLAUDE.md](CLAUDE.md).

## Status

**Phase 2 complete.** 18 UniPROBE PBM sources parsed and validator-clean:

| property | value |
|---|---|
| rows | 17,040,128 |
| DNA-binding domains | 501 in 437 clusters |
| Pfam families / organisms | 32 / 24 |
| DNA sites | all 32,896 non-redundant 8-mers, identical in every source |
| labels | 51,642 binding / 16,713,577 non-binding / 274,909 excluded |
| negative:positive | 324:1 |
| `dbd_seq` length | 43-216 aa (median 77) |

Largest families: Homeodomain 234, Forkhead 54, bHLH 35, zf-C4 28, ETS 23, HMG box 21,
Zn2Cys6 17, bZIP 13, NAC 11, T-box 11, then 22 further families of ten or fewer.

Protein-side depth is concentrated in 28 clusters holding 81 variants; **409 of 437 clusters
(94%) hold a single domain**.

**How the dataset was built, in full: [docs/METHODS.md](docs/METHODS.md).**
**What is open right now: [TODO.md](TODO.md).** What has been settled and why:
[docs/DECISIONS.md](docs/DECISIONS.md).

> **`dbd_seq` is the Pfam domain padded by 10 residues each side** — not the bare envelope.
> See [docs/DOMAIN_POLICY.md](docs/DOMAIN_POLICY.md) for the admission policy and what it
> excludes.

## Setup

```bash
uv venv --python 3.11 .venv
uv pip install --python .venv -e ".[dev]"
.venv/bin/python -m pytest -q
```

## Scripts

All are run from the repository root with the project's interpreter, `.venv/bin/python`.

| script | what it does | runtime |
|---|---|---|
| [`build_dataset.py`](scripts/build_dataset.py) | Parse a source, validate it, write `data/interim/<source>/` | seconds per source; **~4-5 min for `--all`**, provided the Pfam library is pressed |
| [`audit_sources.py`](scripts/audit_sources.py) | Sweep every source for the invariants past defects violated | ~4 min |
| [`make_reports.py`](scripts/make_reports.py) | Regenerate the committed markdown reports for one source | seconds |
| [`build_protein_table.py`](scripts/build_protein_table.py) | Build the protein-side table; fetches UniProt sequences, cached | ~2 min |
| [`record_provenance.py`](scripts/record_provenance.py) | Emit a `PROVENANCE.md` row for a raw file (size, sha256, timestamp) | instant |
| [`press_pfam.py`](scripts/press_pfam.py) | **Run once after downloading Pfam-A.** Presses the HMM library; without it every source pays 19.5 s to re-read 2.2 GB of text | ~25 s |
| [`build_dbd_family_list.py`](scripts/build_dbd_family_list.py) | Regenerate the DNA-binding Pfam family whitelist | seconds |
| [`check_papers.py`](scripts/check_papers.py) | Report which source papers are present, missing or unfiled | instant |

### Building the dataset

```bash
# once per machine — without it every source re-reads 2.2 GB of text HMM
.venv/bin/python scripts/press_pfam.py

# one source, to check a change
.venv/bin/python scripts/build_dataset.py --source BAR15A

# everything (~4-5 min; Cell08 and SCI09 dominate the reads)
.venv/bin/python scripts/build_dataset.py --all 2>&1 | tee /tmp/build.log
```

Expect a `PASS <source>: N rows` line per source and a final `N/30 sources built`. Sources
yielding nothing print `SKIPPED` with a reason — an accession whose every construct is
rejected by the domain policy is a normal outcome, not a failure.

### Checking it

```bash
.venv/bin/python scripts/audit_sources.py            # all sources
.venv/bin/python scripts/audit_sources.py --source GR09
```

Exits non-zero if any check fails. Warnings are expected (redundant files, genes UniPROBE
publishes no sequence for, protein complexes); any `ERROR` line needs attention before the
build is trusted. The audit reads archives, detail pages and the parsed tables, and does not
re-scan Pfam, so it is cheap to run before a rebuild.

### Regenerating reports

```bash
for a in $(ls data/interim | grep -v proteins); do
  .venv/bin/python scripts/make_reports.py --source "$a"
done
.venv/bin/python scripts/build_protein_table.py
```

Reports are committed, so a change in the dataset shows up as a diff.

### Adding a paper

Drop the PDF in [docs/papers_inbox/](docs/papers_inbox/) under whatever name it downloaded
with; it gets identified from its contents, renamed and filed. `check_papers.py` reports what
is still missing.

## Layout

| Path | Purpose |
| --- | --- |
| [docs/](docs/) | The owner's brief (verbatim) and [REFERENCES.md](docs/REFERENCES.md), the verified DOI list |
| [docs/papers_inbox/](docs/papers_inbox/) | Drop new papers here; they get identified, renamed and filed into `docs/papers/` |
| [PROVENANCE.md](PROVENANCE.md) | One row per raw file: URL, accession, timestamp, size, sha256 |
| [configs/thresholds.yaml](configs/thresholds.yaml) | The project's only config file — every binarization cutoff |
| [src/snp2prot/schema.py](src/snp2prot/schema.py) | The 22-column row schema and its validator |
| [src/snp2prot/parsers/](src/snp2prot/parsers/) | One module per source dataset |
| [reports/](reports/) | Generated markdown deliverables, committed to git |
| `data/` | Raw / interim / processed / testsets — all git-ignored |
| [scripts/](scripts/) | Command-line entry points; all logic lives in `src/` |
| [tests/](tests/) | pytest — `.venv/bin/python -m pytest -q` |

## The dataset schema

Every parser emits the same 22 columns; see [src/snp2prot/schema.py](src/snp2prot/schema.py)
for the authoritative definition. Two invariants matter more than the rest:

- **`neg_provenance` is mandatory on every negative.** An assayed-and-unbound negative (PBM
  low E-score) is real evidence; a sequence merely absent from a selection is not. They must
  stay distinguishable downstream.
- **`dna_seq` is never padded.** The table stores what was measured. Padding to a fixed window
  is a featurization choice, made in `snp2prot.data`.
