"""Basic tests to see if the code raises exceptions."""
import logging
from pathlib import Path
import pytest
from pdb2pqr.ligand import mol2
import common


_LOGGER = logging.getLogger(__name__)


_LOGGER.error("Need functional and regression test coverage for --userff")
_LOGGER.error("Need functional and regression test coverage for --usernames")
_LOGGER.error("Need functional and regression test coverage for --ligand")
_LOGGER.error("Need functional and regression test coverage for --apbs-input")


@pytest.mark.parametrize("input_pdb", ["1K1I", "1AFS", "1FAS", "5DV8", "5D8V"], ids=str)
def test_basic_apo(input_pdb, tmp_path):
    """Basic non-regression tests on proteins without ligands."""
    args = "--log-level=INFO --ff=AMBER --drop-water"
    output_pqr = Path(input_pdb).stem + ".pqr"
    common.run_pdb2pqr(args=args, input_pdb=input_pdb, output_pqr=output_pqr,
                        tmp_path=tmp_path)


@pytest.mark.parametrize("input_pdb", ["1K1I", "1AFS", "1FAS", "5DV8", "5D8V"], ids=str)
def test_propka_apo(input_pdb, tmp_path):
    """PROPKA non-regression tests on proteins without ligands."""
    args = "--log-level=INFO --ff=AMBER --drop-water --titration-state-method=propka"
    output_pqr = Path(input_pdb).stem + ".pqr"
    common.run_pdb2pqr(args=args, input_pdb=input_pdb, output_pqr=output_pqr,
                        tmp_path=tmp_path)


@pytest.mark.parametrize("input_mol2", [
    "1HPX-ligand.mol2", "1QBS-ligand.mol2", "1US0-ligand.mol2", "adp.mol2"])
def test_ligand_read(input_mol2):
    """Testing basic aspects of code breaking."""
    ligand = mol2.Mol2Molecule()
    mol2_path = Path("tests/data") / input_mol2
    with open(mol2_path, "rt") as mol2_file:
        ligand.read(mol2_file)


# @pytest.mark.parametrize("input_pdb", ["1K1I", "1FAS"], ids=str)
# def test_propka_apo(input_pdb, tmp_path):
#     """PROPKA titration of proteins without ligands."""
#     args = "--log-level=INFO --ff=AMBER --drop-water --titration-state-method=propka"
#     output_pqr = Path(input_pdb).stem + ".pqr"
#     run_pdb2pqr(args, input_pdb, output_pqr, tmp_path)


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
