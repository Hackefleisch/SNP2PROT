"""The GHT arm's invariants: the encoders, the splits, and the things a fold must not do.

Not the numbers a randomly initialised model produces. What matters is that the DNA tower is
strand-symmetric and translation-invariant *by construction*, that the protein tower never reads
a pad position, that a fold's training and test rows cannot overlap, and that a fold's digest
covers both holdout axes rather than only the protein.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

torch = pytest.importorskip("torch")

from snp2prot.ght import splits as ght_splits  # noqa: E402
from snp2prot.ght.baselines import PWM, neighbours, scan, select, standardise  # noqa: E402
from snp2prot.ght.data import WindowCorpus  # noqa: E402
from snp2prot.ght.model import GenomicTwoTower, genomic_loss  # noqa: E402
from snp2prot.ght.training import fold_rows  # noqa: E402
from snp2prot.models.encoders import (  # noqa: E402
    GenomicDNAEncoder,
    ResidueProteinTower,
    reverse_complement,
)


def _tokens(sequences: list[str]) -> torch.Tensor:
    table = {"A": 0, "C": 1, "G": 2, "T": 3}
    return torch.tensor([[table[b] for b in s] for s in sequences], dtype=torch.long)


# --- the DNA tower ---------------------------------------------------------------------------


def test_the_genomic_encoder_is_strand_symmetric_by_construction():
    """A window and its reverse complement are the same piece of DNA. Taking the maximum over
    the concatenation of both strands makes that an identity rather than something learned."""
    encoder = GenomicDNAEncoder(width=32, widths=(4, 6), channels_per_width=8).eval()
    tokens = _tokens(["ACGTTGCAAACCGGTTACGT", "TTTTAAAACCCCGGGGATAT"])
    with torch.no_grad():
        forward = encoder(tokens)
        backward = encoder(reverse_complement(tokens))
    assert torch.allclose(forward, backward, atol=1e-6)


def test_the_genomic_encoder_cannot_see_where_a_motif_sits():
    """Global max pooling makes positional overfitting structurally impossible: the encoder never
    sees a position, so a motif moved within the window must give the same pooled vector
    (`GHT_PLAN.md` §5). The two windows below are the same length with the same flanking base, so
    they present exactly the same multiset of k-mers and only the offset differs."""
    encoder = GenomicDNAEncoder(width=32, widths=(4,), channels_per_width=8).eval()
    motif = "ACGTGGTCAA"
    left = _tokens(["A" * 10 + motif + "A" * 20])
    right = _tokens(["A" * 20 + motif + "A" * 10])
    with torch.no_grad():
        assert torch.allclose(encoder.pooled(left), encoder.pooled(right), atol=1e-6)


def test_the_genomic_encoder_emits_unit_vectors():
    encoder = GenomicDNAEncoder(width=16).eval()
    with torch.no_grad():
        out = encoder(torch.randint(0, 4, (7, 301)))
    assert torch.allclose(out.norm(dim=-1), torch.ones(7), atol=1e-5)


# --- the per-residue protein tower -----------------------------------------------------------


def test_the_residue_tower_never_reads_a_pad_position():
    """A maximum taken over padding silently invents a motif. Padding a batch with wildly
    different values must change nothing about the shorter sequence's output."""
    tower = ResidueProteinTower(8, 16, reduced=4, widths=(3,), channels_per_width=4).eval()
    residues = torch.randn(1, 12, 8)
    mask = torch.zeros(1, 12, dtype=torch.bool)
    mask[0, :5] = True
    padded = residues.clone()
    padded[0, 5:] = 1000.0
    with torch.no_grad():
        assert torch.allclose(tower(residues, mask), tower(padded, mask), atol=1e-5)


def test_the_residue_tower_handles_a_domain_shorter_than_its_widest_kernel():
    tower = ResidueProteinTower(8, 16, reduced=4, widths=(3, 15), channels_per_width=4).eval()
    mask = torch.zeros(1, 4, dtype=torch.bool)
    mask[0, :4] = True
    with torch.no_grad():
        out = tower(torch.randn(1, 4, 8), mask)
    assert torch.isfinite(out).all()


# --- the model and the loss ------------------------------------------------------------------


def test_mu_zero_is_exactly_plain_bce():
    """The knob mirrors `bce_weight` in the PBM arm: the headline setting must be the earlier
    objective *exactly*, not approximately, or the sweep measures two things at once."""
    model = GenomicTwoTower(protein_features=8, width=16)
    proteins = torch.randn(6, 8)
    tokens = torch.randint(0, 4, (6, 40))
    labels = torch.tensor([1, 0, 1, 0, 1, 0])
    group = torch.tensor([0, 0, 0, 1, 1, 1])
    cosine = model.cosine(proteins, tokens)
    logits = model.calibrate(cosine)
    plain = torch.nn.functional.binary_cross_entropy_with_logits(logits, labels.float())
    total, parts = genomic_loss(logits, labels, cosine, group, model.logit_scale, 0.0)
    assert torch.allclose(total, plain)
    assert parts["infonce"] == 0.0


def test_the_calibration_prior_starts_at_the_training_base_rate():
    model = GenomicTwoTower(protein_features=8, width=16)
    model.set_calibration_prior(0.25)
    assert model.calibration_bias.item() == pytest.approx(np.log(0.25 / 0.75), abs=1e-6)


def test_the_model_has_no_null_anchor():
    """`GHT_PLAN.md` §7: with explicit labelled negatives there is no row to normalise, so the
    anchor is absent rather than present and untrained — a parameter that exists and is never
    used is one more thing a later reader has to rule out."""
    assert "null" not in dict(GenomicTwoTower(protein_features=8, width=16).named_parameters())


# --- splits ----------------------------------------------------------------------------------


def _corpus() -> WindowCorpus:
    rng = np.random.default_rng(0)
    tfs = np.array(["TF1", "TF2", "TF3", "TF4"])
    rows = []
    for i in range(len(tfs)):
        for chrom, split in (("chr1", "Train"), ("chr3", "Train"), ("chr2", "Test")):
            for label, negset in ((1, ""), (0, "shades"), (0, "random")):
                for _ in range(5):
                    rows.append((i, label, split, negset, chrom))
    index, label, split, negset, chrom = (np.array(c) for c in zip(*rows, strict=True))
    return WindowCorpus(
        tfs=tfs,
        tokens=rng.integers(0, 4, (len(rows), 20)).astype(np.uint8),
        tf_index=index.astype(np.int32),
        label=label.astype(np.int8),
        split=split.astype(np.str_),
        negset=negset.astype(np.str_),
        chrom=chrom.astype(np.str_),
        ratio=np.full(len(tfs), 10.0, dtype=np.float32),
    )


def test_a_fold_never_puts_the_same_window_on_both_sides():
    corpus = _corpus()
    fold = ght_splits.chromosome_fold(
        list(corpus.tfs), ["chr1", "chr3"], ["chr2"], corpus.window_counts("Train"), 0.4, 1
    )
    rows = fold_rows(corpus, fold)
    assert len(rows["train"]) and len(rows["validation"]) and len(rows["test"])
    assert not np.intersect1d(rows["train"], rows["test"]).size
    assert not np.intersect1d(rows["train"], rows["validation"]).size


def test_a_tf_fold_holds_out_both_axes():
    """Setting 2 tests held-out TFs on held-out chromosomes. Holding out only the protein would
    leave a test window at a locus some training TF was trained on."""
    corpus = _corpus()
    folds = ght_splits._tf_folds(
        "G1", list(corpus.tfs), np.arange(4), ["chr1", "chr3"], ["chr2"], 2, 1, 0.25, ""
    )
    for fold in folds:
        assert set(fold.test_tfs).isdisjoint(fold.train_tfs)
        assert set(fold.test_chromosomes).isdisjoint(fold.fitting_chromosomes)
        rows = fold_rows(corpus, fold)
        assert not np.intersect1d(rows["train"], rows["test"]).size


def test_the_digest_changes_when_the_chromosome_set_does():
    """A TF list alone does not pin a split: the same held-out TFs scored on the *training*
    chromosomes would be a different and much easier experiment."""
    base = ght_splits.GHTFold("G1", "fold-0", ("TF1",), ("TF2",), ("chr2",), ("chr1",))
    moved = ght_splits.GHTFold("G1", "fold-0", ("TF1",), ("TF2",), ("chr1",), ("chr2",))
    assert base.digest() != moved.digest()


def test_the_validation_carve_never_splits_a_chromosome():
    corpus = _corpus()
    fold = ght_splits.chromosome_fold(
        list(corpus.tfs), ["chr1", "chr3"], ["chr2"], corpus.window_counts("Train"), 0.4, 1
    )
    assert set(fold.validation_chromosomes).isdisjoint(fold.fitting_chromosomes)
    assert set(fold.validation_chromosomes) | set(fold.fitting_chromosomes) == {"chr1", "chr3"}


# --- baselines -------------------------------------------------------------------------------


def test_a_pwm_best_hit_finds_its_own_motif():
    """The claim `GHT_PLAN.md` §5 rests on: a PWM best-hit is one fixed filter plus a global
    max, so scanning a window that contains the motif must beat one that does not."""
    matrix = np.full((4, 4), -5.0, dtype=np.float32)
    for position, base in enumerate([0, 1, 2, 3]):  # ACGT
        matrix[position, base] = 5.0
    tokens = np.array([[0, 1, 2, 3, 0, 0, 0, 0], [0, 0, 0, 0, 0, 0, 0, 0]], dtype=np.uint8)
    best, total = scan(tokens, [PWM("m", "TF", "src", matrix)], device=torch.device("cpu"))
    assert best[0, 0] > best[1, 0]
    assert total[0, 0] > total[1, 0]


def test_selection_picks_the_column_that_separates_the_labels():
    labels = np.array([1, 1, 0, 0])
    scores = np.column_stack([[0.0, 0.0, 1.0, 1.0], [1.0, 1.0, 0.0, 0.0]])
    index, aupr = select(scores, labels)
    assert index == 1
    assert aupr == pytest.approx(1.0)


def test_standardising_leaves_a_single_column_ranking_unchanged():
    """`k = 1` must be unaffected by the rescaling that makes `k = 5` meaningful."""
    scores = np.array([[3.0], [1.0], [2.0]])
    assert list(np.argsort(standardise(scores)[:, 0])) == list(np.argsort(scores[:, 0]))


def test_neighbours_ignores_a_candidate_that_fails_the_overlap_guard():
    identity = np.array([[1.0, 0.9, -np.inf], [0.9, 1.0, 0.4], [-np.inf, 0.4, 1.0]])
    picks, values = neighbours(identity, [0], [1, 2], k=2)
    assert picks[0][0] == 0  # position 0 of train_rows, i.e. domain 1
    assert values[0][0] == pytest.approx(0.9)
    assert np.isnan(values[0][1])


def test_the_window_corpus_selects_positives_plus_one_negative_set():
    corpus = _corpus()
    rows = corpus.rows(("TF1",), ("chr2",), "shades")
    assert set(np.unique(corpus.negset[rows])) == {"", "shades"}
    assert set(np.unique(corpus.chrom[rows])) == {"chr2"}
    assert set(np.unique(corpus.tf_index[rows])) == {0}


def test_the_panel_report_and_windows_agree_on_the_tf_set():
    """A guard against the two Phase A artifacts drifting: every panel TF must have windows."""
    from snp2prot.ght import config

    if not config.PANEL_TABLE.exists():
        pytest.skip("no panel built")
    table = pd.read_parquet(config.PANEL_TABLE)
    missing = [t for t in table.tf if not config.window_file(t).exists()]
    assert not missing, f"{len(missing)} panel TFs have no windows, e.g. {missing[:3]}"


# --- the protein ablations -------------------------------------------------------------------


def test_the_constant_ablation_gives_every_protein_the_same_vector():
    """The control the whole arm is read against: structurally protein-blind, by construction."""
    from snp2prot.ght.data import ablate

    vectors = np.arange(12, dtype=np.float32).reshape(4, 3)
    blind = ablate(vectors, "constant")
    assert np.allclose(blind, blind[0])
    assert np.allclose(blind[0], vectors.mean(axis=0))


def test_the_shuffled_ablation_is_a_derangement():
    """A permutation can leave a TF holding its own vector, and at n = 33 the expected number of
    such fixed points is exactly 1 — which would quietly make one protein un-ablated."""
    from snp2prot.ght.data import ablate

    vectors = np.arange(99, dtype=np.float32).reshape(33, 3)
    for seed in range(5):
        moved = ablate(vectors, "shuffled", seed=seed)
        assert not any(np.allclose(moved[i], vectors[i]) for i in range(len(vectors)))
        assert sorted(moved[:, 0].tolist()) == sorted(vectors[:, 0].tolist())


def test_an_unknown_protein_mode_is_an_error_not_a_pass_through():
    from snp2prot.ght.data import ablate

    with pytest.raises(ValueError):
        ablate(np.zeros((2, 2), dtype=np.float32), "definitely-not-a-mode")


# --- filters as motifs -----------------------------------------------------------------------


def test_a_filter_matches_itself_and_its_reverse_complement():
    """The comparison has to be strand-agnostic, because the encoder scans both strands."""
    from snp2prot.ght.baselines import motif_similarity as similarity

    rng = np.random.default_rng(0)
    motif = rng.normal(size=(10, 4))
    assert similarity(motif, motif) == pytest.approx(1.0, abs=1e-9)
    assert similarity(motif, np.flip(motif, axis=(0, 1))) == pytest.approx(1.0, abs=1e-9)


def test_adding_a_constant_to_a_column_changes_no_similarity():
    """A one-hot column sums to 1, so a per-column offset shifts every window's score equally and
    cannot change a ranking. The comparison must be blind to it or it compares an artefact."""
    from snp2prot.ght.baselines import motif_similarity as similarity

    rng = np.random.default_rng(1)
    a, b = rng.normal(size=(9, 4)), rng.normal(size=(9, 4))
    shifted = b + rng.normal(size=(9, 1))
    assert similarity(a, b) == pytest.approx(similarity(a, shifted), abs=1e-9)


# --- the co-binding overlap search ------------------------------------------------------------


def test_overlap_finds_a_window_on_the_same_grid_and_misses_one_a_window_away():
    """`aliens` and `positives` are fixed-width windows, so "overlapping" is |s - p| < width."""
    from scripts_check_ght_cobinding import overlaps  # type: ignore

    peaks = {"chr1": np.array([1000, 5000], dtype=np.int64)}
    chroms = np.array(["chr1", "chr1", "chr1", "chr2"])
    # 1300 is 300 bases away and two 301 bp windows still touch; 1301 is 301 and they do not.
    starts = np.array([1000, 1300, 1301, 1000], dtype=np.int64)
    assert list(overlaps(chroms, starts, peaks)) == [True, True, False, False]
    # chr2 has no peaks at all, so no window on it can overlap one.
    assert not overlaps(np.array(["chr2"]), np.array([1000], dtype=np.int64), peaks).any()


def test_overlap_on_a_chromosome_with_no_peaks_is_all_false():
    from scripts_check_ght_cobinding import overlaps  # type: ignore

    hit = overlaps(np.array(["chrX", "chrX"]), np.array([1, 2], dtype=np.int64), {})
    assert not hit.any()


def test_each_protein_configuration_gets_its_own_checkpoint_name():
    """The grid trains the same (arm, fold, seed) five ways. Without a distinguishing component
    each would overwrite the last, and the directory would say nothing about which produced what."""
    from snp2prot.ght.config import checkpoint_file, run_tag

    names = {
        checkpoint_file("A1", "C1", "all", 1, run_tag(mode, tower)).name
        for mode in ("real", "constant", "shuffled", "pbm_frozen")
        for tower in ("pooled", "residue")
    }
    assert len(names) == 8
    # The default configuration keeps the plain name the rest of the arm refers to.
    assert checkpoint_file("A1", "C1", "all", 1, run_tag()).name == "ght_A1_C1_all_seed1.pt"


def test_grouping_the_protein_axis_changes_no_score():
    """Each distinct protein is embedded once per pass rather than once per window — an 8,192-row
    batch would otherwise gather an (8192, 219, 1280) tensor for the per-residue tower, which is
    larger than the card. The optimisation must be exactly that and nothing else."""
    group = torch.tensor([0, 0, 1, 1, 2, 2, 2])
    tokens = torch.randint(0, 4, (7, 40))

    pooled = GenomicTwoTower(protein_features=8, width=16).eval()
    vectors = torch.randn(3, 8)
    with torch.no_grad():
        assert torch.allclose(
            pooled.cosine_grouped(vectors, group, tokens),
            pooled.cosine(vectors[group], tokens),
            atol=1e-6,
        )

    residue = GenomicTwoTower(protein_features=8, width=16, protein_tower="residue").eval()
    per_residue = torch.randn(3, 12, 8)
    mask = torch.zeros(3, 12, dtype=torch.bool)
    for i, length in enumerate((12, 5, 9)):
        mask[i, :length] = True
    with torch.no_grad():
        assert torch.allclose(
            residue.cosine_grouped(per_residue, group, tokens, mask),
            residue.cosine(per_residue[group], tokens, mask[group]),
            atol=1e-5,
        )
