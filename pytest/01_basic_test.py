"""Test of basic PDB2PQR functionality"""
import pytest
import pathlib
import logging
import common


_LOGGER = logging.getLogger(__name__)
TEST_PQRs = common.DATA_DIR.glob("*.pqr")


def test_import():
    """Test module import"""
    import pdb2pqr

def pqr_ids(pqr):
    return str(pqr)


@pytest.mark.parametrize("pqr", TEST_PQRs, ids=pqr_ids)
def test_compare(pqr):
    _LOGGER.info("Testing comparison for %s" % pqr)
    common.compare_pqr(pqr, pqr)

