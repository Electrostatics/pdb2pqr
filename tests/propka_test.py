"""Tests for PROPKA functionality."""
import logging
from pathlib import Path
import pytest
import common


_LOGGER = logging.getLogger(__name__)


@pytest.mark.parametrize("input_pdb", ["1K1I", "1AFS", "1FAS", "5DV8", "5D8V"], ids=str)
def test_propka_apo(input_pdb, tmp_path):
    """PROPKA non-regression tests on proteins without ligands."""
    args = "--log-level=INFO --ff=AMBER --drop-water --titration-state-method=propka"
    output_pqr = Path(input_pdb).stem + ".pqr"
    common.run_pdb2pqr(args=args, input_pdb=input_pdb, output_pqr=output_pqr,
                        tmp_path=tmp_path)


# @pytest.mark.parametrize("input_pdb", ["1K1I", "1FAS"], ids=str)
# def test_propka_apo(input_pdb, tmp_path):
#     """PROPKA titration of proteins without ligands."""
#     args = "--log-level=INFO --ff=AMBER --drop-water --titration-state-method=propka"
#     output_pqr = Path(input_pdb).stem + ".pqr"
#     run_pdb2pqr(args, input_pdb, output_pqr, tmp_path)

