"""E-scores computed from raw universal-PBM scans, for sources that publish no E-scores.

Every parsed source so far handed us an 8-mer table someone else had already scored. Kock
et al. 2024 does not: GEO `GSE233827` deposits GenePix `.gpr` scanner output and nothing
else, and the paper's processed deliverable is a `upbm` affinity estimate on a per-allele
log-intensity scale that no fixed cutoff can turn into our labels
(`reports/kock2024_screening.md` §6). This module closes that gap by computing the statistic
the corpus is already built on — the Berger et al. 2006 E-score — from the scans directly.

The inputs are all public: the `.gpr` files, and the probe design for the Bulyk-lab
"all 10-mer" 8x60K array, kept at `data/external/pbm_design/` with a `PROVENANCE.md` row.

**What it reproduces, and what it does not.** Measured against the `BAR15A` E-scores we
already hold, recomputing from Barrera's own raw scans in that same GEO series reproduces the
*ranking* — Spearman 0.83-0.91 over the 8-mers a protein actually binds — and a rank-matched
positive set agrees at Jaccard 0.60-0.77, which is the level of a within-source replicate
(0.70-0.72, `reports/overlap.md`). It does **not** reproduce the *values*: our E-scores are
compressed at the top, so `E >= 0.45` selects 0-32 8-mers where UniPROBE's own values select
129-204. Numbers and the full validation are in `reports/escore_recomputation.md`. Admitting a
source through this path is therefore still a threshold decision (`TODO.md` `T19`, `T20`) —
but one parameter on the existing scale, not a second score type.

Three implementation decisions are load-bearing and each cost a measurement to find:

1. **Join the design to a scan by `(Column, Row)`, never by probe ID.** The two ID sets match
   exactly, once the `.gpr`'s `_v_` is normalised away, and they name different probes: 35 of
   62,976 grid positions agree. E-scores joined by ID correlate with UniPROBE's at
   rho = -0.003 and look perfectly clean.
2. **Keep saturated spots.** They are the brightest probes on the array, so dropping them
   decapitates exactly the ranking the E-score is made of: on one ARX array 260 of the top
   1,000 probes are saturated, and masking them cost the top 8-mer `CTAATTAG` 12 of its 18
   probes. The E-score is a rank statistic, so a ceiling-clipped spot still carries the one
   fact that matters — it is at the top.
3. **Normalise Cy3 as a residual, not a ratio.** Dividing the protein channel by the raw Cy3
   signal removes the sequence-driven part of the double-stranding signal, which is also the
   part that carries binding: it drops agreement from 0.489 to 0.123. Regressing Cy3 on
   sequence composition first and dividing only by the *unexplained* part removes the
   spot-specific defect and leaves the biology: 0.595 to 0.621.
"""

from __future__ import annotations

import gzip
import itertools
import re
from collections.abc import Iterable
from pathlib import Path

import numpy as np
import pandas as pd

from snp2prot.config import EXTERNAL_DIR

#: The constant primer at the 3' end of every de Bruijn probe. What precedes it is the
#: variable region the 8-mers are read from.
PRIMER = "GTCTGTGTTCCGTTGTCCGTGCTG"
#: Length of that variable region: 60-mer probe minus the 24 nt primer.
VARIABLE = 36
#: Probe-ID prefix of the de Bruijn probes. Everything else on the slide is a control
#: (`Cbf_*`, `DarkCorner`, `GE_BrightCorner`) and carries no sequence.
DE_BRUIJN = "dBr_"

DESIGN_FILE = EXTERNAL_DIR / "pbm_design" / "GPL34105_pbm_8x60k_v1.txt.gz"

_COMPLEMENT = str.maketrans("ACGT", "TGCA")
_NUCLEOTIDES = "ACGT"
_DINUCLEOTIDES = ["".join(p) for p in itertools.product(_NUCLEOTIDES, repeat=2)]


def revcomp(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def canonical_8mer(kmer: str) -> str:
    """The strand-independent name of an 8-mer: whichever of it and its reverse complement
    sorts first. This is the same convention UniPROBE's 8-mer tables use, which is why a
    recomputed table joins to `dna_seq` directly with no remapping."""
    return min(kmer, revcomp(kmer))


def load_design(path: str | Path = DESIGN_FILE) -> pd.DataFrame:
    """`Column`, `Row`, `seq` for every de Bruijn probe, `seq` being the variable region."""
    rows: list[tuple[int, int, str]] = []
    with gzip.open(path, "rt") as fh:
        next(fh)
        for line in fh:
            column, row, probe_id, seq = line.rstrip("\n").split("\t")
            if probe_id.startswith(DE_BRUIJN) and seq.endswith(PRIMER):
                rows.append((int(column), int(row), seq[:VARIABLE]))
    return pd.DataFrame(rows, columns=["Column", "Row", "seq"])


def load_gpr(path: str | Path) -> pd.DataFrame:
    """One GenePix results file as `Column, Row, fg, bg, sig, sat, flag`.

    The channel is whatever the scan used and differs between the two files of a pair: the
    protein scan is `F488` (Alexa-488 anti-GST) and the double-stranding control is `F1` or
    `F532` (Cy3), so the column names are sniffed rather than assumed.
    """
    opener = gzip.open if str(path).endswith(".gz") else open
    with opener(path, "rt") as fh:
        for i, line in enumerate(fh):
            if line.startswith('"Block"'):
                header = [c.strip('"') for c in line.rstrip("\n").split("\t")]
                skip = i
                break
        else:
            raise ValueError(f"{path}: no GenePix data block found")

    frame = pd.read_csv(
        path, sep="\t", skiprows=skip + 1, names=header, quotechar='"', low_memory=False
    )
    channel = next(
        (m.group(1) for m in (re.fullmatch(r"F(\d+) Median", c) for c in header) if m), None
    )
    if channel is None:
        raise ValueError(f"{path}: no 'F<channel> Median' column")

    foreground = pd.to_numeric(frame[f"F{channel} Median"], errors="coerce")
    background = pd.to_numeric(frame[f"B{channel} Median"], errors="coerce")
    return pd.DataFrame(
        {
            "Column": frame["Column"].astype(int),
            "Row": frame["Row"].astype(int),
            "fg": foreground,
            "bg": background,
            "sig": foreground - background,
            "sat": pd.to_numeric(frame[f"F{channel} % Sat."], errors="coerce"),
            "flag": pd.to_numeric(frame["Flags"], errors="coerce"),
        }
    )


def probe_signal(gpr: pd.DataFrame, design: pd.DataFrame) -> pd.DataFrame:
    """Background-subtracted signal per de Bruijn probe, joined by grid position.

    Only two things are dropped: spots GenePix flagged bad (`Flags < 0`, typically a handful)
    and spots whose background exceeds their foreground. Saturated spots are deliberately
    kept — see this module's docstring.
    """
    joined = design.merge(gpr, on=["Column", "Row"], how="inner")
    kept = joined[(joined.flag >= 0) & (joined.sig > 0)]
    return kept[["Column", "Row", "seq", "sig"]].reset_index(drop=True)


def _sequence_features(seqs: Iterable[str]) -> np.ndarray:
    s = pd.Series(list(seqs))
    counts = [s.str.count(n) for n in _NUCLEOTIDES]
    counts += [s.str.count(f"(?={d})") for d in _DINUCLEOTIDES]
    return np.column_stack([np.ones(len(s))] + [c.to_numpy(float) for c in counts])


def cy3_residual(cy3: pd.DataFrame, design: pd.DataFrame) -> pd.DataFrame:
    """Per-spot double-stranding defect: observed Cy3 over what its sequence predicts.

    The Cy3 control measures how much double-stranded DNA a spot carries, and that is partly
    a property of the sequence (mono- and dinucleotide composition explains 53-74% of the
    variance) and partly a property of the spot. Only the second part is an artifact. The
    returned `resid` is in log2 units, so dividing the protein signal by `2**resid` removes
    the defect and leaves the sequence effect where it belongs.
    """
    probes = probe_signal(cy3, design)
    features = _sequence_features(probes.seq)
    observed = np.log2(probes.sig.to_numpy(float))
    coefficients, *_ = np.linalg.lstsq(features, observed, rcond=None)
    return probes.assign(resid=observed - features @ coefficients)[["Column", "Row", "resid"]]


def normalised_array(
    protein_gpr: pd.DataFrame, cy3_gpr: pd.DataFrame | None, design: pd.DataFrame
) -> pd.DataFrame:
    """One array's probes with the adjusted intensity the E-score ranks on."""
    probes = probe_signal(protein_gpr, design)
    if cy3_gpr is None:
        return probes
    residual = cy3_residual(cy3_gpr, design)
    merged = probes.merge(residual, on=["Column", "Row"], how="inner")
    return merged.assign(sig=merged.sig / np.exp2(merged.resid))[["Column", "Row", "seq", "sig"]]


def tail_width(probes: pd.DataFrame) -> float:
    """How far the bright end runs above the median, in log2 units.

    Kock's own processing filters replicates on "upper-tail probe intensity width", and it
    does separate the arrays that reproduce from the ones that do not: on the HOXD13 series
    one slide gave 1.1-1.7 for all eight chambers and reproduced badly, the other gave
    3.5-4.1 and reproduced well. **It is not a QC gate on its own**, because a variant that
    has lost DNA binding has a narrow tail for the honest reason: three dead `ARX` variants
    measure 1.4-2.0 on arrays that are fine. Report it; do not filter on it blindly.
    """
    return float(np.log2(probes.sig.quantile(0.999)) - np.log2(probes.sig.median()))


def combine_arrays(arrays: list[pd.DataFrame]) -> pd.DataFrame:
    """Median-centre each array in log space, then average probe by probe.

    Replicates are scanned at different laser powers, so their raw scales differ by a
    constant that has to come out before averaging. Probes missing from any array are
    dropped, which is a few hundred of 41,944.
    """
    if len(arrays) == 1:
        return arrays[0][["seq", "sig"]]
    merged = None
    for i, part in enumerate(arrays):
        centred = np.log2(part.sig) - np.log2(part.sig).median()
        piece = part[["Column", "Row", "seq"]].assign(**{f"lg{i}": centred})
        merged = piece if merged is None else merged.merge(piece, on=["Column", "Row", "seq"])
    columns = [c for c in merged.columns if c.startswith("lg")]
    return merged.assign(sig=np.exp2(merged[columns].mean(axis=1)))[["seq", "sig"]]


def escores(probes: pd.DataFrame) -> pd.Series:
    """The Berger et al. 2006 E-score for every 8-mer on the array.

    Rank the probes by adjusted intensity, keep the top half, and score each 8-mer by the
    Mann-Whitney statistic of the probes containing it (or its reverse complement) against
    the probes that do not, shifted so that the result runs -0.5 to +0.5. Every 8-mer sits on
    16-47 probes of this design, median 37, so the statistic is well determined even for the
    palindromes, which appear least often.
    """
    ranked = probes.sort_values("sig", ascending=False).reset_index(drop=True)
    half = len(ranked) // 2
    ranks_by_8mer: dict[str, list[int]] = {}
    for rank, seq in enumerate(ranked.seq.iloc[:half], start=1):
        seen = set()
        for i in range(len(seq) - 7):
            kmer = canonical_8mer(seq[i : i + 8])
            if kmer not in seen:
                seen.add(kmer)
                ranks_by_8mer.setdefault(kmer, []).append(rank)

    scored: dict[str, float] = {}
    for kmer, ranks in ranks_by_8mer.items():
        n_fg = len(ranks)
        n_bg = half - n_fg
        if n_bg <= 0:
            continue
        u = n_fg * n_bg + n_fg * (n_fg + 1) / 2 - sum(ranks)
        scored[kmer] = u / (n_fg * n_bg) - 0.5
    return pd.Series(scored, name="escore").sort_index()


def score_allele(
    arrays: list[tuple[str | Path, str | Path | None]], design: pd.DataFrame | None = None
) -> tuple[pd.Series, list[float]]:
    """E-scores for one protein from its replicate arrays, plus each array's tail width.

    `arrays` pairs each protein scan with the double-stranding Cy3 scan of the same slide and
    chamber — matched on `(slide, chamber)`, never on date, because the control is often
    scanned days before the binding reaction.
    """
    design = load_design() if design is None else design
    normalised, widths = [], []
    for protein_path, cy3_path in arrays:
        one = normalised_array(
            load_gpr(protein_path), load_gpr(cy3_path) if cy3_path else None, design
        )
        widths.append(tail_width(one))
        normalised.append(one)
    return escores(combine_arrays(normalised)), widths
