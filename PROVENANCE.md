# Provenance

One row per file under `data/raw/`. Append-only: never edit or delete a row, and never edit
the file it describes. This table is publication evidence, not a convenience log.

Regenerate the hash/size columns for a file with:

```bash
python scripts/record_provenance.py data/raw/<source>/<file>
```

## Downloaded files

| source_dataset | file (rel. to `data/raw/`) | source URL | accession | downloaded (UTC) | bytes | sha256 | contents |
|---|---|---|---|---|---|---|---|
| _none yet_ | | | | | | | Phase 0 is scaffold-only; no downloads attempted. |

## Sources attempted and rejected

Record anything that turned out to be unusable under the project's constraints, with the
reason. A short honest catalogue beats a padded one.

| source | date | outcome | reason |
|---|---|---|---|
| _none yet_ | | | |

## Manual-download queue

Sources behind a licence click-through, MTA, or an interactive portal, which the agent
cannot fetch unattended. Each has a `data/raw/<source>/HOWTO.md` with exact steps.

| source | blocker | HOWTO | status |
|---|---|---|---|
| _none yet_ | | | |
