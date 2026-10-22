"""Tests of mmCIF reading in :mod:`pdb2pqr.cif`."""

from io import StringIO

import pdbx
import pytest

from pdb2pqr import cif, pdb

ATOM_SITE_HEADER = """data_test
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.label_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.auth_seq_id
_atom_site.auth_comp_id
_atom_site.auth_asym_id
_atom_site.auth_atom_id
_atom_site.pdbx_PDB_model_num
"""


def read_atoms(rows, header=ATOM_SITE_HEADER):
    """Read ATOM/HETATM records from ``_atom_site`` rows."""
    (block,) = pdbx.load(StringIO(header + "\n".join(rows) + "\n"))
    records, _ = cif.atom_site(block)
    return [
        record
        for record in records
        if isinstance(record, (pdb.ATOM, pdb.HETATM))
    ]


def test_insertion_codes_are_preserved():
    """Residues 100 and 100A remain distinct."""
    atoms = read_atoms(
        [
            "ATOM 1 C CA . ALA A 1 ? 1.0 2.0 3.0 1.00 10.00 100 ALA A CA 1",
            "ATOM 2 C CA . GLY A 2 A 4.0 5.0 6.0 1.00 10.00 100 GLY A CA 1",
        ]
    )
    assert [(atom.res_seq, atom.ins_code) for atom in atoms] == [
        (100, ""),
        (100, "A"),
    ]


def test_values_beyond_pdb_columns_are_preserved():
    """Values that do not fit PDB columns are read, not dropped."""
    atoms = read_atoms(
        [
            "HETATM 100000 C C1 . LIG1 BA 1 . -1234.5678 10000.0 0.0 1.00 "
            "10.00 10000 LIG1 AA C101 1",
        ]
    )
    assert len(atoms) == 1
    atom = atoms[0]
    assert isinstance(atom, pdb.HETATM)
    assert atom.serial == 100000
    assert atom.chain_id == "AA"
    assert atom.res_seq == 10000
    assert atom.res_name == "LIG1"
    assert atom.name == "C101"
    assert (atom.x, atom.y, atom.z) == (-1234.568, 10000.0, 0.0)


def test_author_identity_is_preferred():
    """Author items are used when present, as in a PDB-format file."""
    atoms = read_atoms(
        ["HETATM 1 O O . HOH B . ? 1.0 2.0 3.0 1.00 10.00 201 HOH A O 1"]
    )
    atom = atoms[0]
    assert (atom.chain_id, atom.res_seq, atom.res_name, atom.name) == (
        "A",
        201,
        "HOH",
        "O",
    )


def test_label_identity_is_used_without_author_items():
    """Label items are used when author items are absent or null."""
    header = ATOM_SITE_HEADER.replace("_atom_site.auth_asym_id\n", "")
    atoms = read_atoms(
        ["ATOM 1 C CA . ALA B 7 ? 1.0 2.0 3.0 1.00 10.00 ? ? ? 1"],
        header=header,
    )
    atom = atoms[0]
    assert (atom.chain_id, atom.res_seq, atom.res_name, atom.name) == (
        "B",
        7,
        "ALA",
        "CA",
    )


def test_null_values_are_blank():
    """CIF nulls are stored as empty strings, not "?" or ".", or None."""
    atoms = read_atoms(
        ["ATOM 1 C CA ? ALA A 1 . 1.0 2.0 3.0 1.00 10.00 1 ALA ? CA 1"]
    )
    atom = atoms[0]
    assert (atom.alt_loc, atom.ins_code, atom.chain_id) == ("", "", "A")


def test_unreadable_row_is_an_error():
    """A row without coordinates raises an error instead of being dropped."""
    with pytest.raises(ValueError, match="atom id 1"):
        read_atoms(
            ["ATOM 1 C CA . ALA A 1 ? ? 2.0 3.0 1.00 10.00 1 ALA A CA 1"]
        )


def test_multiple_models_preserve_identity():
    """Atoms in each model keep their identity."""
    (block,) = pdbx.load(
        StringIO(
            ATOM_SITE_HEADER
            + "ATOM 1 C CA . ALA A 1 A 1.0 2.0 3.0 1.00 10.00 100 ALA AA CA 1\n"
            + "ATOM 2 C CA . ALA A 1 A 1.5 2.0 3.0 1.00 10.00 100 ALA AA CA 2\n"
        )
    )
    records, _ = cif.atom_site(block)
    assert [type(record).__name__ for record in records] == [
        "MODEL",
        "ATOM",
        "ENDMDL",
        "MODEL",
        "ATOM",
        "ENDMDL",
    ]
    assert [
        (record.chain_id, record.res_seq, record.ins_code, record.x)
        for record in records
        if isinstance(record, pdb.ATOM)
    ] == [("AA", 100, "A", 1.0), ("AA", 100, "A", 1.5)]
