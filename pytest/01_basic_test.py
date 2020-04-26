"""Test of basic PDB2PQR functionality"""
import pytest
import pathlib
import logging
import common


_LOGGER = logging.getLogger(__name__)


def test_import():
    """Test module import"""
    import pdb2pqr


def test_compare():
    """Test comparison test"""
    for pqr in common.DATA_DIR.glob("*.pqr"):
        _LOGGER.info("Testing comparison for %s" % pqr)
        common.compare_pqr(pqr, pqr)

