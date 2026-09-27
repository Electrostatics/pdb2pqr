"""Tests of I/O functions."""

import logging
from argparse import Namespace
from difflib import Differ
from pathlib import Path

import pytest

from pdb2pqr import main
from pdb2pqr.io import read_dx, read_pqr, read_qcd, write_cube
from pdb2pqr.structures import Atom

_LOGGER = logging.getLogger(__name__)
DATA_DIR = Path("tests/data")
PQR_LIST = list(DATA_DIR.glob("**/*.pqr"))


def _make_atom(serial=1, chain_id="A", res_seq=1, x=0.0):
    atom = Atom(type_="ATOM")
    atom.serial = serial
    atom.name = "N"
    atom.res_name = "GLY"
    atom.chain_id = chain_id
    atom.res_seq = res_seq
    atom.ins_code = ""
    atom.x = x
    atom.y = atom.z = 0.0
    atom.occupancy = 1.0
    atom.temp_factor = 0.0
    atom.seg_id = ""
    atom.element = "N"
    atom.charge = ""
    return atom


@pytest.mark.parametrize(
    "atoms, reason",
    [
        ([_make_atom(chain_id="AA")], "chain"),
        ([_make_atom(res_seq=10000)], "residue"),
        ([_make_atom(res_seq=-1000)], "residue"),
        ([_make_atom()] * 100000, "count"),
        ([_make_atom(x=10000.0)], "coordinate"),
        ([_make_atom(x=-1000.0)], "coordinate"),
        ([_make_atom(x=float("inf"))], "coordinate"),
    ],
)
def test_print_pdb_rejects_unrepresentable_atoms(tmp_path, atoms, reason):
    """PDB output fails before creating a file when fields would overflow."""
    output = tmp_path / "output.pdb"
    args = Namespace(pdb_output=output, keep_chain=True)

    with pytest.raises(RuntimeError, match=reason):
        main.print_pdb(
            args=args,
            atomlist=atoms,
            header_lines="",
            missing_lines=[],
            is_cif=False,
        )

    assert not output.exists()


def test_print_pdb_writes_representable_atoms(tmp_path):
    """Representable atoms continue to use fixed-column PDB output."""
    output = tmp_path / "output.pdb"
    args = Namespace(pdb_output=output, keep_chain=True)

    main.print_pdb(
        args=args,
        atomlist=[_make_atom(x=9999.999), _make_atom(x=-999.999)],
        header_lines="",
        missing_lines=[],
        is_cif=False,
    )

    assert output.read_text().startswith("ATOM      1  N   GLY A   1")


def test_print_pdb_renumbers_input_serials(tmp_path):
    """Large input serials are safe because output atoms are reserialized."""
    output = tmp_path / "output.pdb"
    args = Namespace(pdb_output=output, keep_chain=True)

    main.print_pdb(
        args=args,
        atomlist=[_make_atom(serial=100000)],
        header_lines="",
        missing_lines=[],
        is_cif=False,
    )

    assert output.read_text().startswith("ATOM      1")


@pytest.mark.parametrize("input_pqr", PQR_LIST, ids=str)
def test_read_pqr(input_pqr):
    """Test that :func:`read_pqr` doesn't raise an error.

    Doesn't test functionality since that is implicit in several other tests
    that parse generated PQR output.

    :param input_pqr:  path to PQR file to test
    :type input_pqr:  str
    """
    with open(input_pqr) as pqr_file:
        read_pqr(pqr_file)


def test_read_qcd():
    """Test that :func:`read_pqr` doesn't raise an error.

    Doesn't test functionality.
    """
    qcd_path = DATA_DIR / "dummy.qcd"
    with open(qcd_path) as qcd_file:
        read_qcd(qcd_file)


def test_dx2cube(tmp_path):
    """Test conversion of OpenDX files to Cube files."""
    pqr_path = DATA_DIR / "dx2cube.pqr"
    dx_path = DATA_DIR / "dx2cube.dx"
    cube_gen = tmp_path / "test.cube"
    cube_test = DATA_DIR / "dx2cube.cube"
    _LOGGER.info(f"Reading PQR from {pqr_path}...")
    with open(pqr_path) as pqr_file:
        atom_list = read_pqr(pqr_file)
    _LOGGER.info(f"Reading DX from {dx_path}...")
    with open(dx_path) as dx_file:
        dx_dict = read_dx(dx_file)
    _LOGGER.info(f"Writing Cube to {cube_gen}...")
    with open(cube_gen, "w") as cube_file:
        write_cube(cube_file, dx_dict, atom_list)
    _LOGGER.info(f"Reading this cube from {cube_gen}...")
    this_lines = [line.strip() for line in open(cube_gen)]
    _LOGGER.info(f"Reading test cube from {cube_test}...")
    test_lines = [line.strip() for line in open(cube_test)]
    differ = Differ()
    differences = [
        line
        for line in differ.compare(this_lines, test_lines)
        if line[0] != " "
    ]

    if differences:
        for diff in differences:
            _LOGGER.error(f"Found difference:  {diff}")
        raise ValueError
    _LOGGER.info("No differences found in output")
