"""`metadata.organism` — recovering the organism UniPROBE's template failed to substitute."""

from __future__ import annotations

import pytest

from snp2prot.metadata import organism

HEADER = ">sp|P46581|CND1_CAEEL Protein cnd-1 OS=Caenorhabditis elegans OX=6239 GN=cnd-1 PE=2 SV=2"


@pytest.fixture
def cache(tmp_path, monkeypatch):
    monkeypatch.setattr(organism, "CACHE", tmp_path)
    organism.from_cached_fasta.cache_clear()
    yield tmp_path
    organism.from_cached_fasta.cache_clear()


def test_a_real_organism_is_passed_through_untouched(cache):
    assert organism.resolve("Mus musculus", "P00001") == "Mus musculus"


def test_placeholder_resolves_from_the_cached_header(cache):
    (cache / "P46581.fasta").write_text(HEADER + "\nMSTNMDVSSF\n")
    assert organism.resolve("$species", "P46581") == "Caenorhabditis elegans"


def test_missing_accession_gives_nothing_rather_than_a_guess(cache):
    assert organism.resolve("$species", "Q99999") == ""


@pytest.mark.parametrize(
    "raw", ["Chimera", "N/A", "PBM CONSTRUCTS", "None Available", "unknown", "  "]
)
def test_non_organisms_are_stored_empty_not_as_a_word(cache, raw):
    assert organism.resolve(raw, "") == ""


def test_a_placeholder_accession_is_not_looked_up(cache):
    """ROG18A writes `None Available` in the accession field too; it is not a file name."""
    (cache / "None Available.fasta").write_text(HEADER + "\nM\n")
    assert organism.resolve("Chimera", "None Available") == ""


@pytest.mark.parametrize(
    ("raw", "expected"),
    [
        ("C. elegans", "Caenorhabditis elegans"),
        ("c. ELEGANS", "Caenorhabditis elegans"),
        ("Acyrtosiphon pisum", "Acyrthosiphon pisum"),
        ("Malus x domestica", "Malus domestica"),
    ],
)
def test_known_variant_spellings_map_to_the_binomial(cache, raw, expected):
    assert organism.resolve(raw, "") == expected


@pytest.mark.parametrize(
    "raw",
    [
        "Mus musculus",
        # A genuinely unnamed species, not a formatting variant: normalising it would
        # assert a species that has not been assigned.
        "Sarsia sp. Long Island Sound",
        # Correctly named in three words.
        "Acanthamoeba polyphaga mimivirus",
        # A taxonomic revision, not a spelling disagreement -- left as deposited.
        "Physcomitrella patens",
    ],
)
def test_conventional_and_out_of_scope_names_are_untouched(cache, raw):
    assert organism.resolve(raw, "") == raw


def test_internal_whitespace_is_collapsed(cache):
    assert organism.resolve("  Mus   musculus ", "") == "Mus musculus"


def test_the_deposit_wins_over_uniprot_when_it_names_an_organism(cache):
    """UniProt is a fallback for a missing name, not an override for a present one."""
    (cache / "P46581.fasta").write_text(HEADER + "\nM\n")
    assert organism.resolve("Caenorhabditis briggsae", "P46581") == "Caenorhabditis briggsae"


def test_zero_byte_file_is_a_negative_cache_not_an_organism(cache):
    """`build_protein_table.py` writes an empty file for an accession UniProt does not serve."""
    (cache / "Q17588.fasta").write_text("")
    assert organism.resolve("$species", "Q17588") == ""


def test_header_without_an_os_field(cache):
    (cache / "X00001.fasta").write_text(">tr|X00001|X_NOTHING some protein\nMSTN\n")
    assert organism.resolve("$species", "X00001") == ""


@pytest.mark.parametrize("nxt", ["OX=6239", "GN=cnd-1", "PE=2", "SV=2"])
def test_os_field_stops_at_whichever_key_follows_it(cache, nxt):
    (cache / "A00001.fasta").write_text(f">sp|A00001|A_X p OS=Danio rerio {nxt} rest\nM\n")
    assert organism.resolve("$species", "A00001") == "Danio rerio"
