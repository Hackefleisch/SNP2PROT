"""One parser per source dataset.

Contract — every module here exposes:

    def parse() -> pd.DataFrame

returning a frame that passes `snp2prot.schema.validate()`. Rules:

- A parser reads only from `data/raw/<source>/` and never writes there.
- A parser knows about exactly one source. No cross-source logic, no merging, no
  deduplication against other datasets — that all lives in `snp2prot.merge`.
- Cutoffs come from `snp2prot.thresholds`, never from a literal in the parser body.
- Finish with `schema.coerce(df)` then `schema.add_pair_ids(df)`.

Registry below is filled in as each source lands, so `scripts/build_dataset.py` can
iterate without importing modules that do not exist yet.
"""

from __future__ import annotations

from collections.abc import Callable

import pandas as pd

#: source_dataset -> zero-argument parse function. Populated in Phases 1-4.
REGISTRY: dict[str, Callable[[], pd.DataFrame]] = {}
