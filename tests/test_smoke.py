"""Smoke tests — the package imports and paths resolve."""

from snp2prot import __version__
from snp2prot.config import PROJECT_ROOT
from snp2prot.utils.seed import set_seed


def test_version():
    assert __version__


def test_project_root_is_repo_root():
    assert (PROJECT_ROOT / "pyproject.toml").exists()


def test_set_seed_runs():
    set_seed(0)
