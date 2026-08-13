"""DNA-binding domain annotation, and the admission policy for what may enter the dataset.

**`dbd_seq` is the Pfam envelope padded by `domain.padding_aa` residues on each side, not
the bare envelope.** See `configs/thresholds.yaml` and `docs/DOMAIN_POLICY.md` for why.

A construct is admitted only if the stored subunit satisfies all three conditions:

1. **Sole responsibility** — it alone produced the measured interaction. A construct whose
   Pfam hits span two families (PAX+homeodomain, POU+homeodomain) fails: two domains bind,
   and nothing in the data says which produced a given 8-mer's score.
2. **One continuous region** — not fragments scattered through the protein. A C2H2 array is
   2-6 separate ~23-residue folds on flexible linkers, each needing its own Zn(2+), so it
   fails even though every finger belongs to one family.
3. **The variation lies inside it** — a mutation outside the stored region makes a variant
   sequence-identical to its wild type while carrying a different label. Padding exists to
   satisfy this condition, and `audit_variant_positions` checks it directly.

Rejections are reported, never silently repaired.
"""

from __future__ import annotations

import functools
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import pyhmmer
from pyhmmer.easel import Alphabet, DigitalSequenceBlock, TextSequence

from snp2prot import thresholds
from snp2prot.config import PROJECT_ROOT


class Rejection(str):
    """Why a construct was not admitted. A plain string so it lands in reports unchanged."""


MIXED_FAMILIES = Rejection("mixed_families")
REPEAT_ARRAY = Rejection("repeat_array")
NO_DOMAIN = Rejection("no_domain")


@dataclass(frozen=True)
class DomainHit:
    family: str
    start: int  # 1-based, inclusive, in the scanned sequence
    end: int


@dataclass(frozen=True)
class DomainCall:
    """The admission decision for one construct."""

    sequence: str | None  # the padded domain, or None if rejected
    family: str
    start: int
    end: int
    hits: tuple[DomainHit, ...]
    rejection: Rejection | None = None

    @property
    def ok(self) -> bool:
        return self.rejection is None


def _name(x: Any) -> str:
    return x.decode() if isinstance(x, bytes) else str(x)


@functools.lru_cache(maxsize=1)
def load_hmms(hmm_dir: str | None = None) -> tuple:
    """Load every Pfam HMM once. Cached: parsing them per call dominates runtime."""
    cfg = thresholds.load()["domain"]
    d = Path(hmm_dir or cfg["hmm_dir"])
    if not d.is_absolute():
        d = PROJECT_ROOT / d
    files = sorted(d.glob("*.hmm"))
    if not files:
        raise FileNotFoundError(f"no Pfam HMMs under {d}")
    out = []
    for f in files:
        with pyhmmer.plan7.HMMFile(f) as fh:
            out.extend(list(fh))
    return tuple(out)


def scan(sequences: dict[str, str], hmm_dir: str | None = None) -> dict[str, list[DomainHit]]:
    """Pfam-scan many sequences at once, using each family's gathering threshold."""
    abc = Alphabet.amino()
    keys = list(sequences)
    block = DigitalSequenceBlock(
        abc,
        [TextSequence(name=k.encode(), sequence=sequences[k]).digitize(abc) for k in keys],
    )
    found: dict[str, list[DomainHit]] = {k: [] for k in keys}
    for hmm in load_hmms(hmm_dir):
        pipeline = pyhmmer.plan7.Pipeline(abc, bit_cutoffs="gathering")
        for hit in pipeline.search_hmm(hmm, block):
            for dom in hit.domains.included:
                found[_name(hit.name)].append(
                    DomainHit(_name(hmm.name), dom.alignment.target_from, dom.alignment.target_to)
                )
    return {k: sorted(v, key=lambda h: h.start) for k, v in found.items()}


def call_domain(
    sequence: str, hits: list[DomainHit], config: dict[str, Any] | None = None
) -> DomainCall:
    """Apply the admission policy to one construct and return the padded domain."""
    cfg = config or thresholds.load()["domain"]
    pad = int(cfg["padding_aa"])
    families = {h.family for h in hits}

    if not hits:
        if not cfg.get("allow_no_domain", False):
            return DomainCall(None, "", 0, 0, (), NO_DOMAIN)
        return DomainCall(sequence, "", 1, len(sequence), ())
    if len(families) > 1 and not cfg.get("allow_mixed_families", False):
        return DomainCall(None, ", ".join(sorted(families)), 0, 0, tuple(hits), MIXED_FAMILIES)
    if len(hits) > 1 and not cfg.get("allow_repeat_arrays", False):
        return DomainCall(None, next(iter(families)), 0, 0, tuple(hits), REPEAT_ARRAY)

    lo = max(1, min(h.start for h in hits) - pad)
    hi = min(len(sequence), max(h.end for h in hits) + pad)
    return DomainCall(sequence[lo - 1 : hi], next(iter(families)), lo, hi, tuple(hits))


def audit_variant_positions(call: DomainCall, mut_positions: list[int]) -> list[int]:
    """Condition 3: return the mutation positions that fall outside the stored region.

    A non-empty result means the variant would be sequence-identical to its wild type once
    trimmed, which is exactly the failure the padding exists to prevent.
    """
    if not call.ok:
        return list(mut_positions)
    return [p for p in mut_positions if not (call.start <= p <= call.end)]


def rebase_positions(call: DomainCall, mut_positions: list[int]) -> list[int]:
    """Re-express construct positions as 1-based positions within the stored `dbd_seq`."""
    return [p - call.start + 1 for p in mut_positions]
