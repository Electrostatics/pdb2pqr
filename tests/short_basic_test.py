"""Basic tests to see if the code raises exceptions."""
import logging
from pathlib import Path
import pytest
import common


_LOGGER = logging.getLogger(__name__)


_LOGGER.warning("Need functional and regression test coverage for --userff")
_LOGGER.warning("Need functional and regression test coverage for --usernames")
_LOGGER.warning("Need functional and regression test coverage for --apbs-input")


@pytest.mark.parametrize("input_pdb", ["1K1I", "1AFS", "1FAS", "5DV8", "5D8V"], ids=str)
def test_basic_apo(input_pdb, tmp_path):
    """Basic non-regression tests on proteins without ligands."""
    args = "--log-level=INFO --ff=AMBER --drop-water"
    output_pqr = Path(input_pdb).stem + ".pqr"
    common.run_pdb2pqr(args=args, input_pdb=input_pdb, output_pqr=output_pqr,
                        tmp_path=tmp_path)
