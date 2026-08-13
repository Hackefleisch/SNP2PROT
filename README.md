# SNP2PROT

> TODO: one-paragraph description — what question this project answers, and for whom.

## Status

Skeleton. Structure is in place; contents are filled in incrementally.

## Setup

```bash
python -m venv .venv && source .venv/bin/activate
pip install -e ".[dev]"
```

## Layout

| Path | Purpose |
| --- | --- |
| `src/snp2prot/` | Importable package: data, models, training, evaluation |
| `configs/` | Experiment configs (one file per run) |
| `scripts/` | Thin CLI entry points that call into `src/` |
| `notebooks/` | Exploration only — nothing load-bearing lives here |
| `data/` | Local data, git-ignored (`raw` is read-only by convention) |
| `results/` | Figures, tables, checkpoints, logs — git-ignored |
| `tests/` | pytest suite |

## Usage

```bash
# TODO
python scripts/train.py --config configs/default.yaml
```
