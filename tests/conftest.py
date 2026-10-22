"""Shared fixtures. Tests that need the real data skip cleanly without it."""
from __future__ import annotations

import pathlib
import sys

import pytest

SRC = pathlib.Path(__file__).resolve().parents[1] / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from phmimic.core import load_config, read_fasta  # noqa: E402
from phmimic.scan import Proteome  # noqa: E402


@pytest.fixture(scope="session")
def cfg():
    return load_config()


@pytest.fixture(scope="session")
def proteome_records(cfg):
    path = cfg.path_of("data.proteome")
    if not path.exists():
        pytest.skip("human proteome not present")
    return read_fasta(path)


@pytest.fixture(scope="session")
def proteome(proteome_records):
    return Proteome(proteome_records)


@pytest.fixture(scope="session")
def reports(cfg):
    return cfg.reports
