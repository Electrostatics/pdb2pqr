"""Basic tests to see if the code raises exceptions."""
import logging
from pathlib import Path
import pytest
import common


_LOGGER = logging.getLogger(__name__)


_LOGGER.warning("Need functional and regression test coverage for --userff")
_LOGGER.warning("Need functional and regression test coverage for --usernames")
_LOGGER.warning(
    "Need functional and regression test coverage for --apbs-input"
)


@pytest.mark.parametrize(
    "input_pdb",
    ["4UN3", "1K1I", "1AFS", "1FAS", "5DV8", "5D8V", "1E7G"],
    ids=str,
)
def test_basic_pdb(input_pdb, tmp_path):
    """Non-regression tests on PDB-format biomolecules without ligands."""
    args = "--log-level=INFO --ff=AMBER --drop-water --apbs-input=apbs.in"
    output_pqr = Path(input_pdb).stem + ".pqr"
    common.run_pdb2pqr(
        args=args,
        input_pdb=input_pdb,
        output_pqr=output_pqr,
        tmp_path=tmp_path,
    )


@pytest.mark.parametrize("input_pdb", ["1FAS"], ids=str)
def test_basic_cif(input_pdb, tmp_path):
    """Non-regression tests on CIF-format biomolecules without ligands."""
    args = "--log-level=INFO --ff=AMBER --drop-water --apbs-input=apbs.in"
    output_pqr = Path(input_pdb).stem + ".cif"
    common.run_pdb2pqr(
        args=args,
        input_pdb=input_pdb,
        output_pqr=output_pqr,
        tmp_path=tmp_path,
    )


@pytest.mark.parametrize("input_pdb", ["4E8M"], ids=str)
def test_nucleic_only(input_pdb, tmp_path):
    """Non-regression tests on structures that contain only nucleic acids."""
    args = "--log-level=INFO --ff=AMBER --drop-water --apbs-input=apbs.in"
    output_pqr = Path(input_pdb).stem + ".pqr"
    common.run_pdb2pqr(
        args=args,
        input_pdb=input_pdb,
        output_pqr=output_pqr,
        tmp_path=tmp_path,
    )


@pytest.mark.parametrize("input_pdb", ["1EJG", "3U7T"], ids=str)
@pytest.mark.xfail
def test_broken_backbone(input_pdb, tmp_path):
    """Test graceful failure of optimization with missing backbone atoms."""
    args = "--log-level=INFO --ff=AMBER --drop-water"
    output_pqr = Path(input_pdb).stem + ".pqr"
    common.run_pdb2pqr(
        args=args,
        input_pdb=input_pdb,
        output_pqr=output_pqr,
        tmp_path=tmp_path,
    )
