"""Smoke tests — the package imports and the directory contract holds."""

from snp2prot import __version__
from snp2prot.config import (
    CONFIG_DIR,
    PROJECT_ROOT,
    RAW_DIR,
    TESTSET_DIR,
    interim_table,
    raw_dir,
)
from snp2prot.utils.seed import set_seed


def test_version():
    assert __version__


def test_project_root_is_repo_root():
    assert (PROJECT_ROOT / "pyproject.toml").exists()


def test_expected_directories_exist():
    for p in (RAW_DIR, CONFIG_DIR, TESTSET_DIR):
        assert p.is_dir(), p


def test_source_path_helpers():
    assert raw_dir("BAR15A") == RAW_DIR / "BAR15A"
    assert interim_table("BAR15A").name == "BAR15A.parquet"


def test_set_seed_runs():
    set_seed(0)
