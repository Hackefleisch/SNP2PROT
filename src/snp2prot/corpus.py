"""One row per domain of the merged training table, with everything a split needs.

The merged table is 44 million rows and answers "which domains are there?" expensively.
Every question the modelling code asks about the protein axis — which family, which cluster,
which of them are variants, which carry positive evidence — is answered by 1,338 rows that
already exist in the companion tables. This module joins them into that single view and
nothing else; it computes no new fact.

Three tables feed it, each still the authority on its own column:

- `label_health` — one row per `(dbd_seq, source_dataset)` record, with the label counts and
  the `ok` / `dead_variant` / `no_evidence` verdict (`T21`);
- `merge.resolve` — which record of a duplicated domain survives into the merged table, so
  the domain set here is exactly the merged table's 1,338 and the `source_dataset` on a row
  is the source that actually supplied its measurements (`T15`/`D4`);
- the cluster inventory — `wt_id`, the cluster's reference domain and its size, from
  `scripts/build_clusters.py`.

**`is_variant` is derived from the cluster reference, not from `n_mut_from_wt`.** Both say the
same thing, and the reference is a 1,165-row lookup while `n_mut_from_wt` lives on 44 million
rows. A domain is a variant exactly when it is not its own cluster's reference — which is the
definition `docs/DECISIONS.md` `D3` settled when the reference became the medoid.
"""

from __future__ import annotations

from collections.abc import Sequence

import numpy as np
import pandas as pd

from snp2prot import clusters, label_health, merge

#: The domain view. `source_dataset` is the surviving record's source, `n_pos` its positive
#: count, `verdict` its label-health class, `cluster_size` the number of domains sharing its
#: `wt_id` — the column `docs/DECISIONS.md` (2026-08-14) makes a training-time selector.
DOMAIN_COLUMNS = (
    "dbd_seq",
    "dbd_family",
    "wt_id",
    "source_dataset",
    "protein_id",
    "gene",
    "n_pos",
    "verdict",
    "is_variant",
    "cluster_size",
)


def domains(records: pd.DataFrame | None = None) -> pd.DataFrame:
    """One row per domain in the merged table, sorted by `dbd_seq`.

    Sorted so that the row order is a stable identity: the distance matrix and the 8-mer
    matrix are both indexed by it, and neither stores anything that would let a mismatched
    rebuild be detected other than the domain list itself.
    """
    records = label_health.load() if records is None else records
    kept = merge.resolve(records)
    kept = kept[kept.keep]

    inventory = clusters.load()
    reference = dict(zip(inventory.wt_id, inventory.reference, strict=True))
    size = dict(zip(inventory.wt_id, inventory.n_domains, strict=True))

    out = kept.loc[:, ["dbd_seq", "dbd_family", "wt_id", "source_dataset", "protein_id", "gene"]]
    out = out.assign(
        n_pos=kept["n_pos"].to_numpy(),
        verdict=kept["verdict"].to_numpy(),
        is_variant=[reference.get(w) != s for w, s in zip(kept.wt_id, kept.dbd_seq, strict=True)],
        cluster_size=[size.get(w, 1) for w in kept.wt_id],
    )
    return out.loc[:, list(DOMAIN_COLUMNS)].sort_values("dbd_seq").reset_index(drop=True)


def trainable(domain_table: pd.DataFrame, keep_dead: bool = True) -> pd.Series:
    """Boolean mask of the domains a model is allowed to train on.

    The row-level `label_health.usable` filter, lifted to the domain axis: `no_evidence`
    records are dropped because a silent record with no control cannot be told apart from a
    failed assay, and `dead_variant` records are kept because a variant that measurably lost
    binding is the most informative negative in the corpus (`T21`, `ML_PLAN.md` §3.1).
    """
    drop = {label_health.NO_EVIDENCE}
    if not keep_dead:
        drop.add(label_health.DEAD_VARIANT)
    return ~domain_table["verdict"].isin(drop)


def reference_rows(domain_table: pd.DataFrame) -> np.ndarray:
    """For each domain, the row of its cluster's reference — or `-1` if it is its own.

    What `snp2prot.evaluation.metrics.suppression` needs: a dead variant is scored by comparing
    it with the wild type it came from, which is a *row* of the same tables rather than anything
    stored on the domain itself.
    """
    inventory = clusters.load()
    reference = dict(zip(inventory.wt_id, inventory.reference, strict=True))
    row_of = {seq: i for i, seq in enumerate(domain_table.dbd_seq)}
    out = np.full(len(domain_table), -1, dtype=np.int64)
    for i, (wt_id, seq) in enumerate(zip(domain_table.wt_id, domain_table.dbd_seq, strict=True)):
        found = reference.get(wt_id)
        if found is not None and found != seq:
            out[i] = row_of.get(found, -1)
    return out


def require_aligned(reference: Sequence, reference_name: str = "the corpus", **tables) -> None:
    """Raise unless every named table is in exactly the same domain order as `reference`.

    **The row order is the only thing joining these artifacts.** `domains()` sorts by `dbd_seq`
    and the 8-mer matrix, the distance matrix and each arm's embeddings are all built against
    that order; nothing else links a row of one to a row of another. A rebuild that changes the
    domain set — or merely its order — silently repairs into a table where every protein carries
    another protein's measurements, and every downstream number is quietly wrong rather than
    obviously broken.

    **Order, not membership.** A length check and a set check both pass on a corpus that was
    merely reordered, which is the realistic failure: the same domains, sorted after a rebuild
    that added and removed a few. Only comparing the sequences catches it.

    This lives here rather than in the loaders because a `KmerMatrix` cannot check itself against
    the corpus without importing this module, and a data container should not depend on the
    tables that describe it. So the check sits where two artifacts meet: `snp2prot.training`
    before it indexes them together, and each script after it loads them.
    """
    expected = [str(s) for s in reference]
    for name, table in tables.items():
        found = [str(s) for s in table]
        if found == expected:
            continue
        detail = f"{len(found)} domains against {len(expected)}"
        for row, (a, b) in enumerate(zip(found, expected, strict=False)):
            if a != b:
                detail = (
                    f"first difference at row {row}: "
                    f"{a[:24]!r} where {reference_name} has {b[:24]!r}"
                )
                break
        raise ValueError(
            f"the {name} table is not aligned with {reference_name} — {detail}.\n"
            "It was built for a different domain set; rebuild it before using these together."
        )
