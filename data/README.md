# Data

Contents are git-ignored. Document here how each dataset is obtained so the tree is reproducible.

| Dir | Meaning |
| --- | --- |
| `raw/` | As downloaded. Treat as read-only — never write here. |
| `interim/` | Intermediate transforms. Cheap to regenerate. |
| `processed/` | Model-ready data. Produced by a script under `scripts/`. |
| `external/` | Third-party reference data (annotations, mappings). |

## Sources

- TODO: name, URL/accession, version, license, download command.
