"""The 301 bp genomic windows: read the benchmark's BEDs, cut them out of hg38, store tokens.

`GHT_PLAN.md` D2-D5. Positives are MEX peaks pooled across replicates, 301 bp centred on
summits; the negative sets and the chromosome split are the benchmark's, unchanged, so the
numbers are comparable to theirs by construction.

## What is stored

`uint8` base tokens, `A C G T` -> `0 1 2 3` and anything else `4`, one row per window. Not
one-hot: one-hot at `float32` is 4.8 kB per window against 301 B, which is the difference
between a corpus that lives on the GPU and one that does not. The encoder one-hots a batch.

## Four negative sets, and only one of them is stored in full

| set | what it is | stored |
|---|---|---|
| `positives` | MEX peaks, 301 bp on the summit | all |
| `shades` | a flank at [450, 750] bp either side of each summit, nominally 2:1 | all |
| `random` | GC-matched random genomic regions, 100:1 | **subsampled to `NEGATIVE_RATIO`:1** |
| `aliens` | *other TFs' peaks*, GC-matched, 100:1 | **subsampled to `NEGATIVE_RATIO`:1** |

The subsample is seeded and its realised ratio is stored on every window file, because an auPRC
is anchored to the positive fraction and a changed ratio changes it: at 100:1 chance auPRC is
0.0099, at 10:1 it is 0.091. **So auPRC on `random` here is not comparable to the published
100:1 numbers, and the cross-set comparison is done on auROC**, which is invariant to the
imbalance. That is what `random` is for in the first place (`GHT_PLAN.md` §12): `shades` is not
GC-matched and `random` is, so a model living off base composition scores well on the first and
badly on the second, and the gap between two auROCs is the diagnostic.

**`aliens` is kept although `GHT_PLAN.md` §3 does not single it out, and it is the most
informative negative set this arm has.** For a model fitted per TF — every baseline in §9 —
alien windows are ordinary negatives. For a *protein-conditioned* model they are not: the DNA is
demonstrably bindable, because it is some other TF's peak, so the only thing that can separate
it from a positive is the protein. A per-TF method cannot be asked this question at all.

## Absence of evidence

Every negative here is **0 by construction, not by measurement** (`GHT_PLAN.md` §3). No
experiment says the protein failed to bind. That is categorically different from a PBM negative,
where a low E-score is a measurement against that specific 8-mer, and it is why nothing in this
module writes into `data/processed/training.parquet` and why `neg_provenance` would have to be
distinct if they ever shared a table.
"""

from __future__ import annotations

import gzip
import zipfile
import zlib
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path

import numpy as np
import pandas as pd

from snp2prot.ght import config

#: The window width the benchmark centres on every summit. Checked, never assumed.
WINDOW = 301

#: Negative sets, in the order the plan ranks them: `shades` primary, `random` secondary, and
#: `aliens` the protein-conditioned control described above.
NEGATIVE_SETS = ("shades", "random", "aliens")
POSITIVE_SET = "positives"
ALL_SETS = (POSITIVE_SET, *NEGATIVE_SETS)

SPLITS = ("Train", "Test")

#: How many negatives per positive are kept from the 100:1 sets. 10 rather than 100 because
#: 40 million windows is 12 GB of tokens and hours of scoring for a control, while 10:1 leaves
#: ample power for the auROC the control is actually read on. Recorded on every file.
NEGATIVE_RATIO = 10
SUBSAMPLE_SEED = 20260915

#: `A C G T` -> `0 1 2 3`, everything else (`N`, IUPAC codes, soft-masked lower case handled by
#: upper-casing first) -> 4. A window containing any 4 is dropped and counted.
_TABLE = np.full(256, 4, dtype=np.uint8)
for _i, _b in enumerate(b"ACGT"):
    _TABLE[_b] = _i
    _TABLE[_b + 32] = _i  # soft-masked lower case


@dataclass(frozen=True)
class TFWindows:
    """One TF's windows, every split and every negative set in one file."""

    tf: str
    tokens: np.ndarray  # (n, 301) uint8 in 0..3
    label: np.ndarray  # (n,) int8
    split: np.ndarray  # (n,) 'Train' | 'Test'
    negset: np.ndarray  # (n,) '' for positives, else the set name
    chrom: np.ndarray  # (n,)
    start: np.ndarray  # (n,) int64, 0-based BED start
    summit: np.ndarray  # (n,) int64
    ratio: float  # the realised negatives-per-positive of the subsampled sets

    def rows(self, split: str | None = None, negset: str | None = None) -> np.ndarray:
        """Positions of the positives plus one negative set, for one split."""
        keep = np.ones(len(self.label), dtype=bool)
        if split is not None:
            keep &= self.split == split
        if negset is not None:
            keep &= (self.label == 1) | (self.negset == negset)
        return np.flatnonzero(keep)

    def save(self, path: Path | None = None) -> Path:
        p = path or config.window_file(self.tf)
        p.parent.mkdir(parents=True, exist_ok=True)
        np.savez(
            p,
            tf=np.asarray(self.tf),
            tokens=self.tokens,
            label=self.label,
            split=self.split,
            negset=self.negset,
            chrom=self.chrom,
            start=self.start,
            summit=self.summit,
            ratio=np.asarray(self.ratio),
        )
        return p

    @classmethod
    def load(cls, tf: str, path: Path | None = None) -> TFWindows:
        p = path or config.window_file(tf)
        if not p.exists():
            raise FileNotFoundError(f"no windows for {tf} at {p}\nrun scripts/build_ght_windows.py")
        with np.load(p, allow_pickle=False) as z:
            return cls(
                str(z["tf"]),
                z["tokens"],
                z["label"],
                z["split"],
                z["negset"],
                z["chrom"],
                z["start"],
                z["summit"],
                float(z["ratio"]),
            )


#: What `aliens.bed` writes in the summit column. Those windows are other TFs' peaks, so the
#: summit belongs to a different experiment and the release does not carry it forward.
NO_SUMMIT = "."


def read_bed(archive: zipfile.ZipFile, tf: str, split: str, name: str) -> pd.DataFrame:
    """One `chrom start end summit` BED out of `datasets.zip`, as a frame.

    `summit` is `-1` where the file writes `.`, which `aliens.bed` does throughout: an alien
    window is some other TF's peak and its summit is not this TF's to report. Stored as a
    sentinel rather than dropped, so a reader can tell "no summit" from "summit at 0".
    """
    raw = archive.read(f"GHTS/{split}/{tf}/{name}.bed").decode()
    rows = [line.split("\t") for line in raw.splitlines() if line]
    return pd.DataFrame(
        {
            "chrom": [r[0] for r in rows],
            "start": np.array([int(r[1]) for r in rows], dtype=np.int64),
            "end": np.array([int(r[2]) for r in rows], dtype=np.int64),
            "summit": np.array(
                [-1 if r[3] == NO_SUMMIT else int(r[3]) for r in rows], dtype=np.int64
            ),
        }
    )


def plan_windows(
    archive: zipfile.ZipFile,
    tfs: list[str],
    ratio: int = NEGATIVE_RATIO,
    seed: int = SUBSAMPLE_SEED,
) -> pd.DataFrame:
    """Every window to cut, with its TF, split, set and label — before any sequence is read.

    Built first and in full so that the single pass over hg38 can be driven by a chromosome-major
    order: the genome is a 3.1 GB gzip stream with no index, so it is read once and each
    chromosome is served while it is resident.

    The 100:1 sets are subsampled **per (TF, split)** with a seed derived from the TF name, so a
    rebuild of one TF reproduces its own draw without depending on the others.
    """
    frames = []
    for tf in tfs:
        for split in SPLITS:
            positives = read_bed(archive, tf, split, POSITIVE_SET)
            n_pos = len(positives)
            positives = positives.assign(label=np.int8(1), negset="", tf=tf, split=split)
            frames.append(positives)
            for name in NEGATIVE_SETS:
                bed = read_bed(archive, tf, split, name)
                if name != "shades" and len(bed) > ratio * n_pos:
                    # crc32, not hash(): str hashing is salted per process, so `hash(tf)` would
                    # give a different draw on every run and the file could not be reproduced.
                    rng = np.random.default_rng(
                        [seed, zlib.crc32(tf.encode()), SPLITS.index(split), ALL_SETS.index(name)]
                    )
                    keep = rng.choice(len(bed), ratio * n_pos, replace=False)
                    bed = bed.iloc[np.sort(keep)].reset_index(drop=True)
                frames.append(bed.assign(label=np.int8(0), negset=name, tf=tf, split=split))
    return pd.concat(frames, ignore_index=True)


def cut(plan: pd.DataFrame, genome: Path | None = None, progress=None) -> dict[str, np.ndarray]:
    """One pass over hg38, returning `(n, 301)` tokens in `plan` row order.

    `hg38.fa.gz` is a plain gzip, not BGZF, so there is no random access to buy — the file is
    decompressed once, one chromosome held at a time, and every window on that chromosome cut
    before the next begins. Peak memory is the longest chromosome, ~249 MB.

    A window that runs off the end of a chromosome, or that carries a base other than `ACGT`, is
    left as all-`4` and the caller drops it: neither is expected (the benchmark excludes `N` and
    blacklist regions) and silently repairing one would hide a coordinate-system error.
    """
    path = genome or config.HG38
    width = int((plan.end - plan.start).mode().iloc[0])
    tokens = np.full((len(plan), width), 4, dtype=np.uint8)
    # Positions, not index labels: `tokens` is a fresh array and the two coincide only while the
    # plan carries a RangeIndex, which is not a property worth relying on.
    position = np.arange(len(plan))
    by_chrom: dict[str, np.ndarray] = {
        str(chrom): position[block] for chrom, block in plan.groupby("chrom").indices.items()
    }
    starts = plan.start.to_numpy()
    ends = plan.end.to_numpy()
    wanted = set(by_chrom)

    current: str | None = None
    buffer: list[bytes] = []

    def flush() -> None:
        if current is None or current not in wanted:
            return
        sequence = np.frombuffer(b"".join(buffer), dtype=np.uint8)
        bases = _TABLE[sequence]
        rows = by_chrom[current]
        lo, hi = starts[rows], ends[rows]
        fits = (lo >= 0) & (hi <= len(bases)) & (hi - lo == width)
        index = lo[fits][:, None] + np.arange(width)[None, :]
        tokens[rows[fits]] = bases[index]
        if progress is not None:
            progress(current, int(fits.sum()), int((~fits).sum()))

    with gzip.open(path, "rb") as handle:
        for line in handle:
            if line.startswith(b">"):
                flush()
                current = line[1:].split()[0].decode()
                buffer = []
                continue
            if current in wanted:
                buffer.append(line.rstrip())
    flush()
    return tokens


def build(
    tfs: list[str],
    archive_path: Path | None = None,
    genome: Path | None = None,
    ratio: int = NEGATIVE_RATIO,
    progress=None,
) -> tuple[list[TFWindows], pd.DataFrame]:
    """`(per-TF windows, a summary frame)`. Windows carrying a non-`ACGT` base are dropped."""
    archive = zipfile.ZipFile(archive_path or config.DATASETS_ZIP)
    plan = plan_windows(archive, tfs, ratio=ratio)
    tokens = cut(plan, genome=genome, progress=progress)

    bad = (tokens == 4).any(axis=1)
    plan = plan.reset_index(drop=True).assign(_bad=bad)
    position = np.arange(len(plan))
    out: list[TFWindows] = []
    rows = []
    for tf, block_rows in plan.groupby("tf", sort=False).indices.items():
        block = plan.iloc[block_rows]
        keep = position[block_rows][~block._bad.to_numpy()]
        sub = plan.iloc[keep]
        n_pos = int((sub.label == 1).sum())
        n_sub = int(sub.negset.isin({"random", "aliens"}).sum())
        out.append(
            TFWindows(
                tf=str(tf),
                tokens=tokens[keep],
                label=sub.label.to_numpy(dtype=np.int8),
                split=sub["split"].to_numpy(dtype=np.str_),
                negset=sub.negset.to_numpy(dtype=np.str_),
                chrom=sub.chrom.to_numpy(dtype=np.str_),
                start=sub.start.to_numpy(dtype=np.int64),
                summit=sub.summit.to_numpy(dtype=np.int64),
                ratio=float(n_sub / max(2 * n_pos, 1)),
            )
        )
        for (split, negset), part in sub.groupby(["split", "negset"], sort=False):
            same = (block["split"] == split) & (block["negset"] == negset)
            rows.append(
                {
                    "tf": str(tf),
                    "split": split,
                    "set": negset or POSITIVE_SET,
                    "n": len(part),
                    "n_dropped": int(block.loc[same, "_bad"].sum()),
                }
            )
    return out, pd.DataFrame(rows)


def chromosome_split(windows: list[TFWindows]) -> pd.DataFrame:
    """Which chromosomes the benchmark put on each side, checked for overlap.

    `D5` takes their split unchanged, so the one thing worth verifying is that it is a genuine
    partition — a chromosome appearing on both sides would mean a locus can straddle it.
    """
    seen: dict[str, set[str]] = defaultdict(set)
    for w in windows:
        for split in SPLITS:
            seen[split] |= set(np.unique(w.chrom[w.split == split]).tolist())  # noqa: PD011
    train, test = seen["Train"], seen["Test"]
    return pd.DataFrame(
        {
            "split": ["Train", "Test", "both"],
            "n_chromosomes": [len(train), len(test), len(train & test)],
            "chromosomes": [
                " ".join(sorted(train)),
                " ".join(sorted(test)),
                " ".join(sorted(train & test)),
            ],
        }
    )
