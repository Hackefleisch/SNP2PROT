# SNP2PROT

A merged, binary-labelled dataset of transcription-factor **DNA-binding domain × DNA site**
interactions, built to vary on both axes simultaneously: dense clusters of single-residue DBD
variants against short (≤20 bp) binding sites with single-base variants.

Specification: [docs/TFDNA_MERGE_BRIEF.md](docs/TFDNA_MERGE_BRIEF.md).
Agent working notes and conventions: [CLAUDE.md](CLAUDE.md).

## Status

**Phase 2 complete.** Nine UniPROBE PBM sources, parsed and validator-clean under the
domain policy:

| source | citation | rows | domains | families |
|---|---|---:|---:|---|
| `Cell08` | Berger 2008 | 4,835,712 | 147 | Homeodomain |
| `BAR15A` | Barrera 2016 | 2,960,640 | 90 | Homeodomain, Forkhead, PAX, zf-C4 |
| `MAR17A` | Mariani 2017 | 789,504 | 24 | Forkhead, HLH, Homeodomain, zf-C4 |
| `SCI09` | Badis 2009 | 756,608 | 23 | ETS, Forkhead, HLH, Homeodomain |
| `EMBO10` | Wei 2010 | 723,712 | 22 | ETS |
| `PNAS13` | Nakagawa 2013 | 625,024 | 19 | Forkhead |
| `SHO18A` | Shokri 2019 | 559,232 | 17 | ETS, Forkhead, HLH, Homeodomain |
| `GR09` | Zhu 2009 | 230,272 | 7 | HLH, Homeodomain, zf-C2H2 |
| `ROG18A` | Rogers 2019 | 230,272 | 7 | Forkhead (incl. 4 engineered chimeras) |
| | **total** | **11,710,976** | **340** | 290 clusters |

> **`dbd_seq` is the Pfam domain padded by 10 residues each side** — not the bare envelope.
> See [docs/DOMAIN_POLICY.md](docs/DOMAIN_POLICY.md), which also lists what is excluded and
> why Phase 3 was dropped.

Every construct is admitted only if its stored subunit is solely responsible for the measured
interaction, forms one continuous region, and contains every change that alters binding.
340 of 638 constructs pass; the rest are C2H2 arrays or multi-domain and are reported, not
silently repaired.

A companion protein table ([data/interim/proteins/](src/snp2prot/proteins.py), 361 domains)
carries each domain at four levels — bare Pfam envelope, padded domain, assayed construct and
full-length UniProt sequence — so embeddings can be compared across representations.

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
