# SNP2PROT

A merged, binary-labelled dataset of transcription-factor **DNA-binding domain × DNA site**
interactions, built to vary on both axes simultaneously: dense clusters of single-residue DBD
variants against short (≤20 bp) binding sites with single-base variants.

Specification: [docs/TFDNA_MERGE_BRIEF.md](docs/TFDNA_MERGE_BRIEF.md).
Agent working notes and conventions: [CLAUDE.md](CLAUDE.md).

## Status

**Phase 2 complete, rebuilt under the domain policy.** Four UniPROBE PBM sources, parsed and
validator-clean:

| source | citation | rows | clusters | domains | `dbd_seq` length | families |
|---|---|---:|---:|---:|---|---|
| `BAR15A` | Barrera 2016 | 2,960,640 | 24 | 90 | 72–146 | Homeodomain, Forkhead, PAX, zf-C4 |
| `Cell08` | Berger 2008 | 4,835,712 | 147 | 147 | 58–91 | Homeodomain |
| `EMBO10` | Wei 2010 | 723,712 | 22 | 22 | 94–105 | ETS |
| `PNAS13` | Nakagawa 2013 | 625,024 | 19 | 19 | 95–110 | Forkhead |
| | **total** | **9,145,088** | **212** | **278** | | |

> **`dbd_seq` is the Pfam domain padded by 10 residues each side** — not the bare envelope.
> See [docs/DOMAIN_POLICY.md](docs/DOMAIN_POLICY.md), which also lists what is excluded and
> why Phase 3 was dropped.

A companion [protein table](src/snp2prot/proteins.py) carries each domain at four levels —
bare Pfam envelope, padded domain, assayed construct and full-length UniProt sequence — so
embeddings can be compared across representations from one build.

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
