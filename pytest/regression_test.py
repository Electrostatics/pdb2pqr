"""Regression tests for PDB2PQR behavior."""
import pytest
import logging
import common
from pdb2pqr import cli
from pdb2pqr import main
from pathlib import Path


_LOGGER = logging.getLogger(__name__)
PARSER = cli.build_parser()


@pytest.mark.parametrize(
    "args, input_pdb, output_pqr, expected_pqr",
    [
        pytest.param(
            "--log-level=INFO --ff=AMBER",
            common.DATA_DIR / "1AFS.pdb",
            "output.pqr",
            common.DATA_DIR / "1AFS_ff=AMBER.pqr",
            id="1AFS local"
        ),
        pytest.param(
            "--log-level=INFO --ff=AMBER",
            "1AFS",
            "output.pqr",
            common.DATA_DIR / "1AFS_ff=AMBER.pqr",
            id="1AFS remote"
        )
    ]
)
def test_fast(args, input_pdb, output_pqr, expected_pqr, tmp_path):
    """Basic code to run 1AFS."""
    arg_str = args + " {inp} {out}"
    output_pqr = tmp_path / output_pqr
    _LOGGER.debug("Writing output to %s", output_pqr)
    arg_str = arg_str.format(inp=input_pdb, out=output_pqr)
    args = PARSER.parse_args(arg_str.split())
    main(args)
    common.compare_pqr(output_pqr, expected_pqr)

