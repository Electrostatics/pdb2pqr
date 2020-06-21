"""Tests for ligand functionality."""
import logging
import random
from pathlib import Path
import pytest
from pdb2pqr.ligand import parameterize
import common


_LOGGER = logging.getLogger(__name__)
_LOGGER.warning("Need functional and regression test coverage for --ligand")


@pytest.mark.parametrize("input_mol2", [
    "1HPX-ligand.mol2", "1QBS-ligand.mol2", "1US0-ligand.mol2", "adp.mol2"])
def test_parameterization(input_mol2):
    """Testing basic aspects of code breaking."""
    _LOGGER.warning("Ideally, this would be a regression test.")
    ligand = parameterize.ParameterizedMolecule()
    mol2_path = Path("tests/data") / input_mol2
    with open(mol2_path, "rt") as mol2_file:
        ligand.read(mol2_file)
        for atom in ligand.atoms.values():
            atom.charge = random.uniform(-1, 1)
            atom.old_charge = atom.charge
        ligand.update(ligand)
        for atom in ligand.atoms.values():
            fmt = "{a!s} -- {a.old_charge:5.2f} -> {a.charge:5.2f}"
            _LOGGER.info(fmt.format(a=atom))


TORSION_RESULTS = {
    "ethanol.mol2": {
        ('CAA', 'CAB', 'OAC', 'HAF'), ('HAF', 'OAC', 'CAB', 'CAA')},
    "glycerol.mol2": {
        ('CAA', 'CAB', 'CAC', 'OAF'), ('CAA', 'CAB', 'OAE', 'HAG'),
        ('CAB', 'CAA', 'OAD', 'HAF'), ('CAB', 'CAC', 'OAF', 'HAH'),
        ('CAC', 'CAB', 'CAA', 'OAD'), ('CAC', 'CAB', 'OAE', 'HAG'),
        ('HAF', 'OAD', 'CAA', 'CAB'), ('HAG', 'OAE', 'CAB', 'CAA'),
        ('HAG', 'OAE', 'CAB', 'CAC'), ('HAH', 'OAF', 'CAC', 'CAB'),
        ('OAD', 'CAA', 'CAB', 'CAC'), ('OAD', 'CAA', 'CAB', 'OAE'),
        ('OAE', 'CAB', 'CAA', 'OAD'), ('OAE', 'CAB', 'CAC', 'OAF'),
        ('OAF', 'CAC', 'CAB', 'CAA'), ('OAF', 'CAC', 'CAB', 'OAE')},
    "cyclohexane.mol2": {
        ('CAA', 'CAB', 'CAC', 'CAF'), ('CAA', 'CAD', 'CAE', 'CAF'),
        ('CAB', 'CAA', 'CAD', 'CAE'), ('CAB', 'CAC', 'CAF', 'CAE'),
        ('CAC', 'CAB', 'CAA', 'CAD'), ('CAC', 'CAF', 'CAE', 'CAD'),
        ('CAD', 'CAA', 'CAB', 'CAC'), ('CAD', 'CAE', 'CAF', 'CAC'),
        ('CAE', 'CAD', 'CAA', 'CAB'), ('CAE', 'CAF', 'CAC', 'CAB'),
        ('CAF', 'CAC', 'CAB', 'CAA'), ('CAF', 'CAE', 'CAD', 'CAA')}
}

@pytest.mark.parametrize("input_mol2", [
    "cyclohexane.mol2", "ethanol.mol2", "glycerol.mol2"])
def test_torsions(input_mol2):
    """Test assignment of torsion angles."""
    ligand = parameterize.ParameterizedMolecule()
    mol2_path = Path("tests/data") / input_mol2
    with open(mol2_path, "rt") as mol2_file:
        ligand.read(mol2_file)
    try:
        benchmark = TORSION_RESULTS[input_mol2]
        diff = ligand.torsions ^ benchmark
        if len(diff) > 0:
            err = "Torsion test failed for %s: %s" % (
                input_mol2, sorted(list(diff)))
            raise ValueError(err)
    except KeyError:
        _LOGGER.warning(
            "Skipping torsion test for %s: %s", input_mol2,
            sorted(list(ligand.torsions)))


RING_RESULTS = {
    "ethanol.mol2": set(),
    "glycerol.mol2": set(),
    "cyclohexane.mol2": {('CAA', 'CAD', 'CAE', 'CAF', 'CAC', 'CAB')},
    "naphthalene.mol2": {
        ('CAA', 'CAB', 'CAC', 'CAH', 'CAG', 'CAF'),
        ('CAC', 'CAH', 'CAI', 'CAJ', 'CAE', 'CAD')},
    "anthracene.mol2": {
        ('CAC', 'CAJ', 'CAK', 'CAL', 'CAE', 'CAD'),
        ('CAE', 'CAL', 'CAM', 'CAN', 'CAG', 'CAF'),
        ('CAA', 'CAB', 'CAC', 'CAJ', 'CAI', 'CAH')},
    "phenalene.mol2": {
        ('CAA', 'CAE', 'CAF', 'CAG', 'CAC', 'CAB'),
        ('CAD', 'CAI', 'CAJ', 'CAK', 'CAF', 'CAE'),
        ('CAF', 'CAG', 'CAH', 'CAM', 'CAL', 'CAK')}
}


@pytest.mark.parametrize("input_mol2", ["phenalene.mol2"])
def test_bad_rings(input_mol2):
    """Test assignment of torsion angles."""
    ligand = parameterize.ParameterizedMolecule()
    mol2_path = Path("tests/data") / input_mol2
    with open(mol2_path, "rt") as mol2_file:
        ligand.read(mol2_file)
    benchmark = RING_RESULTS[input_mol2]
    with pytest.raises(ValueError) as err:
        diff = ligand.rings ^ benchmark
        if len(diff) > 0:
            err = "Ring test failed for %s: %s" % (
                input_mol2, sorted(list(diff)))
            raise ValueError(err)
    err = "Known bond detection failure for %s: %s" % (input_mol2, err)
    _LOGGER.error(err)


@pytest.mark.parametrize("input_mol2", [
    "cyclohexane.mol2", "ethanol.mol2", "glycerol.mol2", "anthracene.mol2",
    "naphthalene.mol2"])
def test_rings(input_mol2):
    """Test assignment of torsion angles."""
    ligand = parameterize.ParameterizedMolecule()
    mol2_path = Path("tests/data") / input_mol2
    with open(mol2_path, "rt") as mol2_file:
        ligand.read(mol2_file)
    benchmark = RING_RESULTS[input_mol2]
    diff = ligand.rings ^ benchmark
    if len(diff) > 0:
        err = "Ring test failed for %s: %s" % (
            input_mol2, sorted(list(diff)))
        raise ValueError(err)
    for atom_name in ligand.atoms:
        atom = ligand.atoms[atom_name]
        if atom.num_rings > 0:
            str_ = "%d rings: %s" % (atom.num_rings, atom)
            _LOGGER.debug(str_)


BOND_RESULTS = {
    "cyclohexane.mol2": 6 * ["single"],
    "ethanol.mol2": [
        "single", "single", None, "single", "single"],
    "glycerol.mol2": [
        None, "single", "single", "single", "single", None, "single", None],
    "acetylcholine.mol2": [
        "single", "double", "single", "single", "single", "single", "single",
        "single", "single"],
    "cyanide.mol2": [
        "triple", "single"],
    "pyrrole.mol2": [
        "aromatic", "aromatic", "aromatic", "aromatic", "aromatic", None],
    "fatty-acid.mol2": [
        "double", "double", "single", "single", "single", "single", "double",
        "single", "single", "single", "single"],
    "trimethylamine.mol2": ["single", None, "single", "single"],
    "naphthalene.mol2": 11 * ["aromatic"]
}


@pytest.mark.parametrize("input_mol2", ["fatty-acid.mol2", "pyrrole.mol2"])
def test_bad_bonds(input_mol2):
    """Test known failure of detected bond types."""
    ligand = parameterize.ParameterizedMolecule()
    mol2_path = Path("tests/data") / input_mol2
    with open(mol2_path, "rt") as mol2_file:
        ligand.read(mol2_file)
    results = BOND_RESULTS[input_mol2]
    for ibond, bond in enumerate(ligand.bonds):
        try:
            if bond.bond_order != results[ibond]:
                err = "Incorrect order for %s. Got %s, expected %s" % (
                    str(bond), bond.bond_order, results[ibond])
                err = "Known bond detection failure for %s: %s" % (
                    input_mol2, err)
                _LOGGER.error(err)
        except IndexError:
            err = "Add test for %s -- %s (%s)" % (
                input_mol2, str(bond), bond.bond_order)
            raise IndexError(err)


@pytest.mark.parametrize("input_mol2", [
    "cyclohexane.mol2", "ethanol.mol2", "glycerol.mol2", "acetylcholine.mol2",
    "cyanide.mol2", "trimethylamine.mol2", "naphthalene.mol2"])
def test_bonds(input_mol2):
    """Test detection of bond types."""
    ligand = parameterize.ParameterizedMolecule()
    mol2_path = Path("tests/data") / input_mol2
    with open(mol2_path, "rt") as mol2_file:
        ligand.read(mol2_file)
    results = BOND_RESULTS[input_mol2]
    for ibond, bond in enumerate(ligand.bonds):
        try:
            if bond.bond_order != results[ibond]:
                err = "Incorrect order for %s. Got %s, expected %s" % (
                    str(bond), bond.bond_order, results[ibond])
                raise ValueError(err)
        except IndexError:
            err = "Add test for %s -- %s (%s)" % (
                input_mol2, str(bond), bond.bond_order)
            raise IndexError(err)


@pytest.mark.parametrize("input_pdb", ["1HPX", "1QBS", "1US0"], ids=str)
def test_ligand(input_pdb, tmp_path):
    """PROPKA non-regression tests on proteins without ligands."""
    ligand = Path("tests/data") / ("%s-ligand.mol2" % input_pdb)
    args = "--log-level=INFO --ff=AMBER --drop-water --ligand=%s" % ligand
    output_pqr = Path(input_pdb).stem + ".pqr"
    common.run_pdb2pqr(
        args=args, input_pdb=input_pdb, output_pqr=output_pqr,
        tmp_path=tmp_path)


# @pytest.mark.parametrize(
#     "args, input_pdb, input_mol2, output_pqr",
#     [
#         pytest.param(
#             "--log-level=INFO --ff=AMBER",
#             "1HPX",
#             common.DATA_DIR / "1HPX-ligand.mol2",
#             "output.pqr",
#             id="1HPX-ligand AMBER"
#         ),
#         pytest.param(
#             "--log-level=INFO --ff=AMBER",
#             common.DATA_DIR / "1QBS.pdb",
#             common.DATA_DIR / "1QBS-ligand.mol2",
#             "output.pqr",
#             id="1QBS-ligand AMBER"
#         ),
#         pytest.param(
#             "--log-level=INFO --ff=AMBER",
#             common.DATA_DIR / "1US0.pdb",
#             common.DATA_DIR / "1US0-ligand.mol2",
#             "output.pqr",
#             id="1US0-ligand AMBER"
#         ),
#     ]
# )
# def test_ligand(args, input_pdb, input_mol2, output_pqr, tmp_path):
#     """Test ligand handling."""
#     args_ = "{args} --ligand={ligand}".format(args=args, ligand=input_mol2)
#     run_pdb2pqr(args_, input_pdb, output_pqr, tmp_path)
#     _LOGGER.warning("This test needs better checking to avoid silent failure.")