"""The numbers a model has to beat, per setting.

`GHT_PLAN.md` §9. Two settings and **a different field of competitors in each**, which is the
structure of the whole experiment:

| baseline | Setting 1 (held-out chromosomes) | Setting 2 (held-out TF) |
|---|---|---|
| PWM best-hit / sum-occupancy | yes — fitted on that TF's own motifs | **cannot run** |
| ArChIPelago | yes, by citation (`D11`) | **cannot run** |
| 1NN / 5NN motif transfer | (pointless — the TF's own motif is available) | yes |

A PWM and a random forest over PWM hits are *fitted per TF*. Asked about a protein they have
never seen they have nothing to compute. That is not a weakness of the baselines, it is the
capability gap the arm exists to demonstrate, and it is why Setting 1's goal is **parity** and
Setting 2's is the result.

## What the nearest-neighbour baseline copies, and why it is not the peaks

The PBM arm's `nn_lookup` copies the most identical training domain's *binary calls over a
shared 8-mer vocabulary*. That works because every protein is measured against the same 32,896
8-mers. Here the DNA axis is genomic loci and the chromosome split means a training TF has **no
measurement at any test-chromosome locus** — so a peak-copying lookup is not merely weak, it is
undefined.

So the transfer is of the neighbour's **motif**, not of its measurements: score a held-out TF's
windows with the nearest training TF's selected PWM. That is also the honest analogue, because
it is what the field actually does — CIS-BP infers a motif for an uncharacterised TF from the
most identical characterised DBD, and Weirauch et al. 2014 set the per-family identity
thresholds for exactly that. `k = 5` averages the top five neighbours' *standardised* scores
weighted by identity, for the reason the PBM arm gives: `k = 1` emits one neighbour's ranking and
`k = 5` emits a blend, and a model emitting a continuous score should be compared against both.

**Scores are standardised per (TF, PWM) before blending.** A log-odds best-hit has no common
scale across motifs of different width and information content, so averaging raw scores from five
neighbours would be a weighted vote dominated by whichever motif happens to have the largest
numbers. Ranking metrics are invariant to a per-TF monotone transform, so standardising costs the
`k = 1` baseline nothing and makes `k = 5` mean what it says.
"""

from __future__ import annotations

import re
import zipfile
from dataclasses import dataclass

import numpy as np
import pandas as pd
import torch

from snp2prot.evaluation import metrics
from snp2prot.ght import config

#: Where a PWM's provenance sits in its filename: `TF.CONSTRUCT@SOURCE@codename@tool@params.pwm`.
_SOURCE = re.compile(r"@([A-Za-z0-9.]+)@")


@dataclass(frozen=True)
class PWM:
    """One log-odds position weight matrix, `(width, 4)` in `ACGT` order."""

    name: str
    tf: str
    source: str
    matrix: np.ndarray

    @property
    def width(self) -> int:
        return self.matrix.shape[0]


def read_pwms(tfs, platform: str = "GHTS", path=None) -> dict[str, list[PWM]]:
    """Every PWM the release ships for each TF on one platform, deduplicated by name.

    `platform` is `GHTS` or `CHS`, and the choice is a leakage boundary the release made for us:
    "the PWM set for training CHS models does not include any GHTS-derived PWMs and vice versa"
    (their README). So the `GHTS` set carries HT-SELEX, PBM and SMiLE-seq motifs alongside the
    GHT-SELEX ones, and no ChIP-seq motif at all.

    **The GHT-SELEX-derived motifs in this set were discovered from the GHT-SELEX data itself,
    and the release does not say whether discovery was restricted to the training chromosomes.**
    If it was not, this baseline has seen test-chromosome sequence that the model has not, which
    makes Setting 1 parity a conservative claim rather than a flattering one. Recorded rather
    than corrected: re-deriving the motifs would abandon `D11`, which is what keeps our numbers
    comparable to theirs.
    """
    archive = zipfile.ZipFile(path or config.PWMS_ZIP)
    wanted = {str(t) for t in tfs}
    out: dict[str, list[PWM]] = {t: [] for t in wanted}
    seen: dict[str, set[str]] = {t: set() for t in wanted}
    for entry in archive.namelist():
        if entry.endswith("/") or not entry.startswith(f"{platform}/"):
            continue
        parts = entry.split("/")
        tf = parts[1]
        if tf not in wanted:
            continue
        name = parts[-1].removesuffix(".pwm")
        if name in seen[tf]:
            continue
        seen[tf].add(name)
        text = archive.read(entry).decode().splitlines()
        rows = [line.split() for line in text[1:] if line.strip()]
        matrix = np.array([[float(v) for v in row] for row in rows], dtype=np.float32)
        if matrix.ndim != 2 or matrix.shape[1] != 4:
            raise ValueError(f"{entry}: expected a (width, 4) matrix, got {matrix.shape}")
        found = _SOURCE.search(name)
        out[tf].append(PWM(name, tf, found.group(1) if found else "", matrix))
    return out


def scan(
    tokens: np.ndarray,
    pwms: list[PWM],
    device: torch.device | None = None,
    chunk: int = 4096,
) -> tuple[np.ndarray, np.ndarray]:
    """`(n, m)` best-hit and sum-occupancy scores for `n` windows against `m` PWMs.

    A PWM best-hit is a convolution followed by a maximum over positions and strands — literally
    `GenomicDNAEncoder` with one fixed filter and no learning — so it is computed the same way,
    on the GPU, grouped by width because a convolution needs one kernel size at a time.

    Sum-occupancy is `log sum exp` over the same positions rather than the max: the standard
    "total binding" aggregation, which counts several weak sites where the best hit sees only the
    strongest. Both are returned because which one suits a TF is an empirical question and
    choosing on the training chromosomes is free.
    """
    device = device or torch.device("cuda" if torch.cuda.is_available() else "cpu")
    best = np.full((len(tokens), len(pwms)), -np.inf, dtype=np.float32)
    total = np.full((len(tokens), len(pwms)), -np.inf, dtype=np.float32)
    by_width: dict[int, list[int]] = {}
    for i, pwm in enumerate(pwms):
        by_width.setdefault(pwm.width, []).append(i)

    source = torch.from_numpy(tokens)
    for indices in by_width.values():
        # (m, 4, w): the conv weight is the PWM transposed; the reverse complement is the same
        # matrix flipped along both axes, which is why only one kernel is built.
        stack = np.stack([pwms[i].matrix for i in indices])
        weight = torch.from_numpy(stack).permute(0, 2, 1).contiguous().to(device)
        reverse = torch.flip(weight, dims=[1, 2]).contiguous()
        for start in range(0, len(tokens), chunk):
            block = source[start : start + chunk].to(device).long()
            one_hot = torch.nn.functional.one_hot(block, 4).float().transpose(1, 2)
            forward = torch.nn.functional.conv1d(one_hot, weight)
            backward = torch.nn.functional.conv1d(one_hot, reverse)
            both = torch.cat([forward, backward], dim=-1)
            best[start : start + chunk, indices] = both.amax(-1).float().cpu().numpy()
            total[start : start + chunk, indices] = (
                torch.logsumexp(both, dim=-1).float().cpu().numpy()
            )
    return best, total


def select(scores: np.ndarray, labels: np.ndarray) -> tuple[int, float]:
    """The column with the highest auPRC against `labels`, and that auPRC.

    Selection happens on the **training** chromosomes and the chosen column is then read on the
    test chromosomes, which is what makes this a fitted per-TF baseline rather than an oracle.
    """
    if scores.size == 0:
        return -1, float("nan")
    values = [metrics.average_precision(labels, scores[:, j]) for j in range(scores.shape[1])]
    best = int(np.nanargmax(values))
    return best, float(values[best])


def standardise(scores: np.ndarray) -> np.ndarray:
    """Zero mean, unit variance per column, so scores from different PWMs can be averaged.

    A no-op for any ranking metric applied to a single column, which is why `k = 1` is unaffected.
    """
    mean = scores.mean(axis=0, keepdims=True)
    deviation = scores.std(axis=0, keepdims=True)
    return (scores - mean) / np.where(deviation > 0, deviation, 1.0)


def identity_matrix(panel: pd.DataFrame) -> np.ndarray:
    """`(n, n)` percent identity between the panel's domains, under the overlap guard.

    The guard is not optional: with free terminal gaps the aligner can park most of both
    sequences in gaps that cost nothing and report a handful of edits over the sliver that
    survives (`snp2prot.distances`, `T25`). Under a cross-family holdout — which every Setting 2
    fold is, mostly — that degenerate case is the typical one rather than the exception.
    """
    from snp2prot import distances as distances_module
    from snp2prot import thresholds

    distances = distances_module.build(list(panel.dbd_seq), processes=1)
    identity = distances.identity()
    guarded = distances.comparable(float(thresholds.load()["cluster"]["min_overlap"]))
    return np.where(guarded & np.isfinite(identity), identity, -np.inf)


def neighbours(
    identity: np.ndarray, test_rows, train_rows, k: int
) -> tuple[np.ndarray, np.ndarray]:
    """The `k` most identical training panel rows for each test row, best first.

    A stable sort rather than `argpartition`, for the reason `nn_lookup.choose` gives: ties in
    identity are common and `argpartition` promises no order among equal values, so a rule like
    "ties break towards the lower row" would simply be untrue.
    """
    block = identity[np.ix_(np.asarray(test_rows), np.asarray(train_rows))]
    k = min(k, block.shape[1])
    top = np.argsort(-block, axis=1, kind="stable")[:, :k]
    values = np.take_along_axis(block, top, axis=1)
    empty = ~np.isfinite(values)
    return np.where(empty, -1, top), np.where(empty, np.nan, values)


# --- reading a trained encoder as a set of motifs ---------------------------------------------
#
# `GHT_PLAN.md` §5 argues that the genomic DNA encoder is the learned generalisation of the PWM
# baseline. These three functions are what makes that checkable rather than rhetorical, and they
# live here rather than in the script that reports them because the figure script needs them too.

#: Columns two motifs must share before a similarity is reported. Below this a "match" is an
#: artefact of the search over offsets rather than a statement about the motif.
MIN_OVERLAP = 6


def centre(matrix: np.ndarray) -> np.ndarray:
    """Subtract each column's mean.

    **A `Conv1d` weight over a one-hot encoding is a position weight matrix** — `weight[c, b, j]`
    is what base `b` at offset `j` contributes to channel `c`, which is what a log-odds column is.
    The one freedom it leaves undetermined is an additive constant per column: a one-hot column
    always sums to 1, so adding a constant shifts every window's score equally and cannot change a
    ranking. Centring removes it, and makes a learned filter and a published PWM comparable.
    """
    return matrix - matrix.mean(axis=1, keepdims=True)


def best_alignment(
    a: np.ndarray, b: np.ndarray, min_overlap: int = MIN_OVERLAP
) -> tuple[float, bool, int]:
    """`(similarity, b_was_reverse_complemented, offset)` for the best alignment of `b` to `a`.

    **Which strand won is worth returning, not just how well.** The encoder scans both strands and
    takes the maximum, so a learned filter's orientation is arbitrary: a filter that has learned
    `CTTTGAT` and a published motif written `ATCAAAG` are the same motif. Reporting only the score
    would leave a reader of the figure to work that out from the letters.
    """
    a = centre(np.asarray(a, dtype=np.float64))
    raw = np.asarray(b, dtype=np.float64)
    best = (-1.0, False, 0)
    for flipped, candidate in ((False, raw), (True, np.flip(raw, axis=(0, 1)))):
        candidate = centre(candidate)
        for offset in range(-len(candidate) + min_overlap, len(a) - min_overlap + 1):
            lo = max(0, offset)
            hi = min(len(a), offset + len(candidate))
            if hi - lo < min_overlap:
                continue
            left = a[lo:hi].ravel()
            right = candidate[lo - offset : hi - offset].ravel()
            denominator = np.linalg.norm(left) * np.linalg.norm(right)
            if denominator > 0:
                score = float(left @ right / denominator)
                if score > best[0]:
                    best = (score, flipped, offset)
    return best


def motif_similarity(a: np.ndarray, b: np.ndarray, min_overlap: int = MIN_OVERLAP) -> float:
    """Best Pearson correlation between two `(width, 4)` matrices, over offsets and both strands.

    The ungapped core of what TOMTOM does, implemented here rather than adding a MEME-suite
    dependency for thirty lines. Both strands because the encoder scans both.
    """
    return best_alignment(a, b, min_overlap)[0]


def learned_filters(checkpoint, which: str = "best") -> dict[int, np.ndarray]:
    """`width -> (channels, width, 4)` learned filters, read out of a saved GHT checkpoint."""
    state = torch.load(checkpoint, map_location="cpu", weights_only=False)
    key = {"best": "state_dict", "final": "final_state_dict"}[which]
    out: dict[int, np.ndarray] = {}
    for name, tensor in state[key].items():
        if name.startswith("dna.convolutions.") and name.endswith(".weight"):
            # (channels, 4, width) -> (channels, width, 4), the orientation a PWM is written in.
            block = tensor.numpy().transpose(0, 2, 1)
            out[block.shape[1]] = block
    return out
