"""One parser per source dataset.

Contract — every module here exposes:

    def parse() -> pd.DataFrame

returning a frame that passes `snp2prot.schema.validate()`. Rules:

- A parser reads only from `data/raw/<source>/` and never writes there.
- A parser knows about exactly one source. No cross-source logic, no merging, no
  deduplication against other datasets — that all lives in `snp2prot.merge`.
- Cutoffs come from `snp2prot.thresholds`, never from a literal in the parser body.
- Finish with `schema.coerce(df)`, which fills the derived columns (`pair_id`, `dna_len`).

Registry below is filled in as each source lands, so `scripts/build_dataset.py` can
iterate without importing modules that do not exist yet.
"""

from __future__ import annotations

from collections.abc import Callable

import pandas as pd

from snp2prot.parsers import bar15a, uniprobe_panels, weirauch2014

#: source_dataset -> zero-argument parse function. Populated in Phases 1-4.
REGISTRY: dict[str, Callable[[], pd.DataFrame]] = {
    bar15a.SOURCE: bar15a.parse,
    weirauch2014.SOURCE: weirauch2014.parse,
    **{acc: uniprobe_panels.make_parser(acc) for acc in uniprobe_panels.PANELS},
}
