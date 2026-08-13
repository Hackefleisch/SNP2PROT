# SNP2PROT — working notes for Claude

Machine learning project. **Currently a skeleton** — most modules raise `NotImplementedError`
and configs are placeholders. Filled in step by step; don't treat TODOs as bugs to fix
unprompted.

## What this project is

TODO: fill in — the scientific question, the inputs, the prediction target, the success criterion.

## Layout

```
src/snp2prot/       importable package (the only place real logic lives)
  config.py         project paths + Config dataclass; single source of truth for paths
  data/             loading, preprocessing, splits
  models/           architectures; all satisfy models/base.py::BaseModel
  training/         training loops
  evaluation/       metrics
  utils/            seeding and shared helpers
configs/            one YAML per experiment; copy default.yaml, never edit it in place
scripts/            thin CLIs — argparse + a call into src/, nothing more
notebooks/          exploration only; nothing importable
data/               git-ignored (raw/ interim/ processed/ external/)
results/            git-ignored (figures/ tables/ checkpoints/ logs/)
tests/              pytest
```

## Conventions

- **Logic goes in `src/`.** Scripts and notebooks call into it; they don't define it.
- **Paths come from `snp2prot.config`.** No hardcoded or relative paths in modules.
- **`data/raw/` is read-only.** Every derived file is produced by a committed script.
- **Seed everything** via `utils.seed.set_seed` at the top of any entry point.
- **Configs are data, not code.** A run should be reproducible from its YAML plus the commit hash.
- Type hints on public functions; `from __future__ import annotations` at the top of modules.
- Line length 100, ruff for lint + import order.

## Commands

```bash
pip install -e ".[dev]"     # setup
pytest                      # tests
ruff check . && ruff format .
python scripts/train.py --config configs/default.yaml
```

## Not yet decided

Fill these in as they're settled — they change how code should be written:

- Framework (PyTorch / sklearn / JAX) — dependency list in `pyproject.toml` is empty on purpose.
- Experiment tracking (W&B / MLflow / plain CSV logs).
- Data source and how it's fetched.
- Headline evaluation metric and validation strategy (random split vs. grouped/held-out).
