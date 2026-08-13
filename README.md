# SNP2PROT

A merged, binary-labelled dataset of transcription-factor **DNA-binding domain × DNA site**
interactions, built to vary on both axes simultaneously: dense clusters of single-residue DBD
variants against short (≤20 bp) binding sites with single-base variants.

Specification: [docs/TFDNA_MERGE_BRIEF.md](docs/TFDNA_MERGE_BRIEF.md).
Agent working notes and conventions: [CLAUDE.md](CLAUDE.md).

## Status

**Phase 2 complete.** Four UniPROBE PBM sources are parsed and validator-clean, 11.8M rows:

| source | citation | rows | clusters | DBDs | contributes |
|---|---|---:|---:|---:|---|
| `BAR15A` | Barrera 2016 | 5,263,360 | 41 | 160 | matched WT/point-mutant pairs |
| `Cell08` | Berger 2008 | 5,164,672 | 157 | 157 | mouse homeodomains |
| `EMBO10` | Wei 2010 | 723,712 | 22 | 22 | ETS family |
| `PNAS13` | Nakagawa 2013 | 657,920 | 20 | 20 | forkhead, across four species |

See [reports/](reports/) for per-source binarization summaries, cluster inventories and
validation output, and [reports/OPEN_ITEMS.md](reports/OPEN_ITEMS.md) for decisions waiting
on the project owner.

## Setup

```bash
uv venv --python 3.11 .venv
uv pip install --python .venv -e ".[dev]"
.venv/bin/python -m pytest -q
```

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
| [tests/](tests/) | pytest |

## The dataset schema

Every parser emits the same 22 columns; see [src/snp2prot/schema.py](src/snp2prot/schema.py)
for the authoritative definition. Two invariants matter more than the rest:

- **`neg_provenance` is mandatory on every negative.** An assayed-and-unbound negative (PBM
  low E-score) is real evidence; a sequence merely absent from a selection is not. They must
  stay distinguishable downstream.
- **`dna_seq` is never padded.** The table stores what was measured. Padding to a fixed window
  is a featurization choice, made in `snp2prot.data`.
