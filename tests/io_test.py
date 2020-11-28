"""Tests of I/O functions."""
import logging
from difflib import Differ
from pathlib import Path
import pytest
from pdb2pqr.io import read_pqr, read_dx, write_cube, read_qcd


_LOGGER = logging.getLogger(__name__)
DATA_DIR = Path("tests/data")
PQR_LIST = list(DATA_DIR.glob("**/*.pqr"))


@pytest.mark.parametrize("input_pqr", PQR_LIST, ids=str)
def test_read_pqr(input_pqr):
    """Test that :func:`read_pqr` doesn't raise an error.

    Doesn't test functionality since that is implicit in several other tests
    that parse generated PQR output.

    :param input_pqr:  path to PQR file to test
    :type input_pqr:  str
    """
    with open(input_pqr, "rt") as pqr_file:
        read_pqr(pqr_file)


def test_read_qcd():
    """Test that :func:`read_pqr` doesn't raise an error.

    Doesn't test functionality.
    """
    qcd_path = DATA_DIR / "dummy.qcd"
    with open(qcd_path, "rt") as qcd_file:
        read_qcd(qcd_file)


def test_dx2cube(tmp_path):
    """Test conversion of OpenDX files to Cube files."""
    pqr_path = DATA_DIR / "dx2cube.pqr"
    dx_path = DATA_DIR / "dx2cube.dx"
    cube_gen = tmp_path / "test.cube"
    cube_test = DATA_DIR / "dx2cube.cube"
    _LOGGER.info(f"Reading PQR from {pqr_path}...")
    with open(pqr_path, "rt") as pqr_file:
        atom_list = read_pqr(pqr_file)
    _LOGGER.info(f"Reading DX from {dx_path}...")
    with open(dx_path, "rt") as dx_file:
        dx_dict = read_dx(dx_file)
    _LOGGER.info(f"Writing Cube to {cube_gen}...")
    with open(cube_gen, "wt") as cube_file:
        write_cube(cube_file, dx_dict, atom_list)
    _LOGGER.info(f"Reading this cube from {cube_gen}...")
    this_lines = [line.strip() for line in open(cube_gen, "rt")]
    _LOGGER.info(f"Reading test cube from {cube_test}...")
    test_lines = [line.strip() for line in open(cube_test, "rt")]
    differ = Differ()
    differences = [
        line
        for line in differ.compare(this_lines, test_lines)
        if line[0] != " "
    ]

    if differences:
        for diff in differences:
            _LOGGER.error(f"Found difference:  {diff}")
        raise ValueError()
    _LOGGER.info("No differences found in output")
