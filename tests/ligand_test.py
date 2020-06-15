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
        torsions = set()
        for name, atom in ligand.atoms.items():
            torsions |= set(atom.torsions)
        try:
            benchmark = TORSION_RESULTS[input_mol2]
            diff = torsions ^ benchmark
            if len(diff) > 0:
                err = "Torsion test failed for %s: %s" % (
                    input_mol2, sorted(list(diff)))
                raise ValueError(err)
        except KeyError:
            _LOGGER.warning(
                "Skipping torsions for %s: %s", input_mol2,
                sorted(list(torsions)))


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