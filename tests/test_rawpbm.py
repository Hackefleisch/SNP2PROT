"""Recomputing E-scores from raw scans.

The three cases that matter are the three that cost a measurement to find on real data:
the join key, saturated spots, and the form of the Cy3 correction. Each is asserted here on
a synthetic array small enough to reason about.
"""

from __future__ import annotations

import gzip
import random
import zlib

import numpy as np
import pandas as pd
import pytest

from snp2prot import rawpbm

MOTIF = "ACGTACGT"
GRID_COLUMNS, GRID_ROWS = 40, 40  # 1,600 probes, enough for a stable rank statistic


@pytest.fixture(scope="module")
def design() -> pd.DataFrame:
    """A fake slide: a quarter of the probes carry `MOTIF`, the rest are random."""
    rng = random.Random(20260818)
    rows = []
    for column in range(1, GRID_COLUMNS + 1):
        for row in range(1, GRID_ROWS + 1):
            seq = "".join(rng.choice("ACGT") for _ in range(rawpbm.VARIABLE))
            if (column + row) % 4 == 0:
                at = rng.randrange(0, rawpbm.VARIABLE - 8)
                seq = seq[:at] + MOTIF + seq[at + 8 :]
            rows.append((column, row, seq))
    return pd.DataFrame(rows, columns=["Column", "Row", "seq"])


def write_gpr(path, design: pd.DataFrame, signal, channel: str = "488", **overrides) -> None:
    """A GenePix results file with the preamble the reader has to skip."""
    columns = [
        "Block", "Column", "Row", "Name", "ID", "X", "Y", "Dia.",
        f"F{channel} Median", f"B{channel} Median", f"F{channel} % Sat.", "Flags",
    ]  # fmt: skip
    flags = overrides.get("flags", {})
    saturated = overrides.get("saturated", {})
    background = overrides.get("background", 100)
    with gzip.open(path, "wt") as fh:
        fh.write("ATF\t1.0\n7\t12\n")
        for line in ('"Type=GenePix Results 3"', '"Wavelengths=488"', '"Creator=test"'):
            fh.write(line + "\n")
        fh.write("\t".join(f'"{c}"' for c in columns) + "\n")
        for i, r in enumerate(design.itertuples()):
            key = (r.Column, r.Row)
            fh.write(
                "\t".join(
                    str(v)
                    for v in (
                        1,
                        r.Column,
                        r.Row,
                        f'"dBr_{i}"',
                        f'"dBr_{i}_v_Jan07"',
                        0,
                        0,
                        30,
                        int(signal(r.seq, key) + background),
                        background,
                        saturated.get(key, 0),
                        flags.get(key, 0),
                    )
                )  # fmt: skip
                + "\n"
            )


def flat(_seq, _key):
    return 1000.0


def motif_driven(seq, _key):
    """Bright where the motif is, plus a per-probe jitter so that no two probes tie: with
    tied intensities the top half is decided by row order, and nothing downstream should
    depend on that."""
    base = 8000.0 if MOTIF in seq or rawpbm.revcomp(MOTIF) in seq else 1000.0
    return base + zlib.crc32(seq.encode()) % 400


def test_canonical_8mer_is_the_uniprobe_convention():
    assert rawpbm.canonical_8mer("TTTTTTTT") == "AAAAAAAA"
    assert rawpbm.canonical_8mer("ACGTACGT") == "ACGTACGT"  # palindrome names itself


def test_the_reader_sniffs_whichever_channel_the_scan_used(tmp_path, design):
    protein = tmp_path / "a488.gpr.gz"
    control = tmp_path / "cy3.gpr.gz"
    write_gpr(protein, design, flat, channel="488")
    write_gpr(control, design, flat, channel="1")
    for path in (protein, control):
        g = rawpbm.load_gpr(path)
        assert len(g) == len(design)
        assert (g.sig == 1000).all()


def test_a_planted_motif_comes_out_on_top(tmp_path, design):
    path = tmp_path / "a488.gpr.gz"
    write_gpr(path, design, motif_driven)
    e = rawpbm.escores(rawpbm.probe_signal(rawpbm.load_gpr(path), design))
    assert e[rawpbm.canonical_8mer(MOTIF)] == pytest.approx(e.max())
    assert e.max() > 0.45


def test_saturated_spots_are_kept_because_they_are_the_brightest(tmp_path, design):
    """Masking them decapitates the ranking: on a real ARX array 260 of the top 1,000 probes
    are saturated, and the top 8-mer loses 12 of its 18 probes."""
    path = tmp_path / "a488.gpr.gz"
    bright = {
        (r.Column, r.Row): 100 for r in design.itertuples() if MOTIF in r.seq
    }  # every motif probe flagged as saturated
    write_gpr(path, design, motif_driven, saturated=bright)
    probes = rawpbm.probe_signal(rawpbm.load_gpr(path), design)
    assert len(probes) == len(design)
    e = rawpbm.escores(probes)
    assert e[rawpbm.canonical_8mer(MOTIF)] == pytest.approx(e.max())


def test_flagged_spots_are_dropped(tmp_path, design):
    path = tmp_path / "a488.gpr.gz"
    bad = {(r.Column, r.Row): -100 for r in design.head(7).itertuples()}
    write_gpr(path, design, flat, flags=bad)
    assert len(rawpbm.probe_signal(rawpbm.load_gpr(path), design)) == len(design) - 7


def test_the_design_joins_by_position_not_by_probe_id(tmp_path, design):
    """The .gpr's IDs are a red herring — here they are shuffled against the grid, exactly
    as the real files are, and the answer must not change."""
    path = tmp_path / "a488.gpr.gz"
    write_gpr(path, design, motif_driven)
    scrambled = design.sample(frac=1.0, random_state=0)  # same rows, different order
    e = rawpbm.escores(rawpbm.probe_signal(rawpbm.load_gpr(path), scrambled))
    straight = rawpbm.escores(rawpbm.probe_signal(rawpbm.load_gpr(path), design))
    both = pd.concat([e.rename("a"), straight.rename("b")], axis=1).dropna()
    assert np.corrcoef(both.a.rank(), both.b.rank())[0, 1] > 0.99
    assert set(e.nlargest(20).index) == set(straight.nlargest(20).index)
    assert e[rawpbm.canonical_8mer(MOTIF)] == pytest.approx(e.max())


def test_cy3_correction_lifts_a_spot_that_double_stranded_badly(tmp_path, design):
    protein = tmp_path / "a488.gpr.gz"
    control = tmp_path / "cy3.gpr.gz"
    victim = (1, 1)

    def half_dsdna(_seq, key):
        return 500.0 if key == victim else 1000.0

    write_gpr(protein, design, half_dsdna)
    write_gpr(control, design, half_dsdna, channel="1")
    adjusted = rawpbm.normalised_array(rawpbm.load_gpr(protein), rawpbm.load_gpr(control), design)
    corrected = adjusted.merge(
        pd.DataFrame([victim], columns=["Column", "Row"]), on=["Column", "Row"]
    )
    assert corrected.sig.iloc[0] == pytest.approx(adjusted.sig.median(), rel=0.05)


def test_replicates_are_centred_before_averaging(tmp_path, design):
    """Two scans of one array at different laser powers must not average as raw numbers."""
    dim = design.assign(sig=1000.0)[["Column", "Row", "seq", "sig"]]
    bright = design.assign(sig=8000.0)[["Column", "Row", "seq", "sig"]]
    dim.loc[dim.index[:10], "sig"] = 4000.0
    bright.loc[bright.index[:10], "sig"] = 32000.0
    combined = rawpbm.combine_arrays([dim, bright])
    top = combined.nlargest(10, "sig")
    assert set(top.seq) == set(design.seq.iloc[:10])
    assert np.isclose(top.sig.min() / combined.sig.median(), 4.0, rtol=0.01)


def test_tail_width_measures_the_bright_end(tmp_path, design):
    flat_probes = design.assign(sig=1000.0)
    assert rawpbm.tail_width(flat_probes) == pytest.approx(0.0, abs=1e-9)
    peaked = design.assign(sig=1000.0)
    peaked.loc[peaked.index[:20], "sig"] = 16000.0
    assert rawpbm.tail_width(peaked) > 3.0
