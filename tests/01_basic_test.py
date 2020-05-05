"""Test of basic PDB2PQR functionality"""
import logging
import importlib
from pathlib import Path
import pytest
import common


_LOGGER = logging.getLogger(__name__)
TEST_PQRS = list(common.DATA_DIR.glob("*.pqr"))


@pytest.mark.parametrize("module", ["pdb2pqr", "propka"], ids=str)
def test_import(module):
    """Test module import"""
    importlib.import_module(module)


@pytest.mark.parametrize("pqr", TEST_PQRS, ids=str)
def test_pqr_compare(pqr):
    """Test comparison of PQRs."""
    common.compare_pqr(pqr, pqr)
