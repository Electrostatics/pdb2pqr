"""Tests for ligand functionality."""

import logging
from math import isclose
from pathlib import Path

import common
import pandas as pd
import pytest
from ligand_results import (
    CHARGES_1HPX,
    FORMAL_CHARGE_RESULTS,
    PARAMETER_RESULTS,
    RING_RESULTS,
    TORSION_RESULTS,
)
from numpy.testing import assert_almost_equal

from pdb2pqr.ligand.mol2 import Mol2Molecule

_LOGGER = logging.getLogger(__name__)


ALL_LIGANDS = set(TORSION_RESULTS) | set(RING_RESULTS)
ALL_LIGANDS |= {
    "1HPX-ligand.mol2",
    "1QBS-ligand.mol2",
    "1US0-ligand.mol2",
    "adp.mol2",
    "acetate.mol2",
}
ALL_LIGANDS = sorted(ALL_LIGANDS)


def test_peoe_charges():
    """Specifically test PEOE charges."""
    ligand = Mol2Molecule()
    mol2_path = Path("tests/data/1HPX-ligand.mol2")
    with open(mol2_path) as mol2_file:
        ligand.read(mol2_file)
    # WARNING: Nathan? The old_total_charge is never used
    # old_total_charge = sum(
    #    atom.formal_charge for atom in ligand.atoms.values()
    # )
    ligand.assign_parameters()
    new_total_charge = 0
    test_charges = []
    for atom in ligand.atoms.values():
        new_total_charge += atom.charge
        test_charges.append(atom.charge)
    assert_almost_equal(test_charges, CHARGES_1HPX, verbose=True)


@pytest.mark.parametrize("input_mol2", ALL_LIGANDS)
def test_assign_parameters(input_mol2):
    """Tests to break basic ligand functionality."""
    ligand = Mol2Molecule()
    mol2_path = Path("tests/data") / input_mol2
    with open(mol2_path) as mol2_file:
        ligand.read(mol2_file)
    old_total_charge = sum(
        atom.formal_charge for atom in ligand.atoms.values()
    )
    ligand.assign_parameters()
    new_total_charge = 0
    test_results = []
    for atom in ligand.atoms.values():
        test_row = {
            "name": atom.name,
            "type": atom.type,
            "charge": atom.charge,
            "radius": atom.radius,
        }
        test_results.append(test_row)
        new_total_charge += atom.charge
    _LOGGER.debug(f"Test results: {test_results}")
    test_results = pd.DataFrame(test_results)
    test_results = test_results.set_index("name")
    # _LOGGER.debug(f"Test results:\n{test_results.to_string()}")
    _LOGGER.info(
        f"Total charge: {old_total_charge:5.2f} -> {new_total_charge:5.2f}"
    )
    expected_results = pd.DataFrame(PARAMETER_RESULTS[input_mol2])
    expected_results = expected_results.set_index("name")
    assert_almost_equal(
        test_results["charge"].to_numpy(),
        expected_results["charge"].to_numpy(),
    )
    assert_almost_equal(
        test_results["radius"].to_numpy(),
        expected_results["radius"].to_numpy(),
        verbose=True,
    )


@pytest.mark.parametrize("input_mol2", ALL_LIGANDS)
def test_formal_charge(input_mol2):
    """Testing formal charge calculation."""
    ligand = Mol2Molecule()
    mol2_path = Path("tests/data") / input_mol2
    with open(mol2_path) as mol2_file:
        ligand.read(mol2_file)
    expected_results = FORMAL_CHARGE_RESULTS[input_mol2]
    errors = []
    for iatom, atom in enumerate(ligand.atoms.values()):
        try:
            expected_charge = expected_results[iatom]
        except IndexError as e:
            err = f"Missing result for {atom.name}, {atom.type}, "
            err += f"{atom.formal_charge}"
            raise IndexError(err) from e
        if not isclose(atom.formal_charge, expected_charge):
            err = (
                f"Atom {atom.name} {atom.type} with bond order "
                f"{atom.bond_order}: expected {expected_charge}, "
                f"got {atom.formal_charge}"
            )
            errors.append(err)
    if errors:
        err = "Errors in test values:\n"
        err += "\n".join(errors)
        raise ValueError(err)


@pytest.mark.parametrize("input_mol2", TORSION_RESULTS)
def test_torsions(input_mol2):
    """Test assignment of torsion angles."""
    ligand = Mol2Molecule()
    mol2_path = Path("tests/data") / input_mol2
    with open(mol2_path) as mol2_file:
        ligand.read(mol2_file)
    test = set()
    # Only test heavy-atom torsions
    for torsion in ligand.torsions:
        has_hydrogen = any(atom.startswith("H") for atom in torsion)
        if not has_hydrogen:
            test.add(torsion)
    try:
        expected = TORSION_RESULTS[input_mol2]
        diff = test ^ expected
        if len(diff) > 0:
            err = (
                f"Torsion test failed for {input_mol2}:\n"
                f"Got: {sorted(test)}\n"
                f"Expected: {sorted(expected)}\n"
                f"Difference: {sorted(diff)}"
            )
            raise ValueError(err)
    except KeyError as e:
        err = f"No results for {input_mol2}: {sorted(ligand.torsions)}"
        raise KeyError(err) from e


@pytest.mark.parametrize("input_mol2", ALL_LIGANDS)
def test_rings(input_mol2):
    """Test assignment of torsion angles."""
    ligand = Mol2Molecule()
    mol2_path = Path("tests/data") / input_mol2
    with open(mol2_path) as mol2_file:
        ligand.read(mol2_file)
    test = ligand.rings
    try:
        benchmark = RING_RESULTS[input_mol2]
    except KeyError as e:
        err = f"Missing expected results for {input_mol2}: {test}"
        raise KeyError(err) from e
    diff = test ^ benchmark
    if len(diff) > 0:
        err = (
            f"Ring test failed for {input_mol2}:\n"
            f"Got: {sorted(test)}\n"
            f"Expected: {sorted(benchmark)}\n"
            f"Difference: {sorted(diff)}"
        )
        raise ValueError(err)
    for atom_name in ligand.atoms:
        atom = ligand.atoms[atom_name]
        if atom.num_rings > 0:
            str_ = f"{atom.num_rings} rings: {atom}"
            _LOGGER.debug(str_)


@pytest.mark.parametrize("input_pdb", ["1HPX"], ids=str)
def test_ligand_biomolecule(input_pdb, tmp_path):
    """PROPKA non-regression tests on biomolecules without ligands."""
    input_pdb = Path(input_pdb)
    ligand = Path("tests/data") / f"{input_pdb.stem}-ligand.mol2"
    args = f"--log-level=INFO --ff=AMBER --drop-water --ligand={ligand}"
    output_pqr = Path(input_pdb).stem + ".pqr"
    _LOGGER.debug(f"Running test in {tmp_path}")
    common.run_pdb2pqr_for_tests(
        args=args,
        input_pdb=input_pdb,
        output_pqr=output_pqr,
        tmp_path=tmp_path,
    )
