"""The GHT-SELEX arm. `docs/GHT_PLAN.md` is the plan of record.

Deliberately a separate package from the PBM modelling path, and deliberately writing to
`data/interim/ght/` and `data/processed/ght/` rather than into the merged table. The two arms
share the **protein tower** and nothing else: a PBM row is one 8-mer scored against a complete,
shared vocabulary, while a GHT row is one 301 bp genomic window that exists for exactly one TF.
`GHT_PLAN.md` §7 is why the loss cannot be shared either.

Nothing here merges into `data/processed/training.parquet`.
"""
