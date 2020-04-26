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
            id="1AFS basic local"
        ),
        pytest.param(
            "--log-level=INFO --ff=AMBER",
            "1AFS",
            "output.pqr",
            common.DATA_DIR / "1AFS_ff=AMBER.pqr",
            id="1AFS basic remote"
        )
    ]
)
def test_basic(args, input_pdb, output_pqr, expected_pqr, tmp_path):
    """Basic code to run 1AFS."""
    arg_str = args + " {inp} {out}"
    output_pqr = tmp_path / output_pqr
    _LOGGER.debug("Writing output to %s", output_pqr)
    arg_str = arg_str.format(inp=input_pdb, out=output_pqr)
    args = PARSER.parse_args(arg_str.split())
    main(args)
    common.compare_pqr(output_pqr, expected_pqr)


@pytest.mark.parametrize(
    "args, input_pdb, output_pqr, expected_pqr",
    [
        pytest.param(
            "--log-level=INFO --whitespace --ff=AMBER",
            common.DATA_DIR / "1AFS.pdb",
            "output.pqr",
            common.DATA_DIR / "1AFS_whitespace_ff=AMBER.pqr",
            id="1AFS whitespace AMBER"
        ),
        pytest.param(
            "--log-level=INFO --whitespace --ff=CHARMM",
            common.DATA_DIR / "1AFS.pdb",
            "output.pqr",
            common.DATA_DIR / "1AFS_whitespace_ff=CHARMM.pqr",
            id="1AFS whitespace CHARMM"
        ),
        pytest.param(
            "--log-level=INFO --whitespace --ff=PARSE",
            common.DATA_DIR / "1AFS.pdb",
            "output.pqr",
            common.DATA_DIR / "1AFS_whitespace_ff=PARSE.pqr",
            id="1AFS whitespace PARSE"
        ),
        pytest.param(
            "--log-level=INFO --whitespace --ff=PEOEPB",
            common.DATA_DIR / "1AFS.pdb",
            "output.pqr",
            common.DATA_DIR / "1AFS_whitespace_ff=PEOEPB.pqr",
            id="1AFS whitespace PEOEPB"
        ),
        pytest.param(
            "--log-level=INFO --whitespace --ff=SWANSON",
            common.DATA_DIR / "1AFS.pdb",
            "output.pqr",
            common.DATA_DIR / "1AFS_whitespace_ff=SWANSON.pqr",
            id="1AFS whitespace SWANSON"
        ),
        pytest.param(
            "--log-level=INFO --whitespace --ff=TYL06",
            common.DATA_DIR / "1AFS.pdb",
            "output.pqr",
            common.DATA_DIR / "1AFS_whitespace_ff=TYL06.pqr",
            id="1AFS whitespace TYL06"
        )
    ]
)
def test_whitespace(args, input_pdb, output_pqr, expected_pqr, tmp_path):
    """Basic code to run 1AFS with --whitespace."""
    arg_str = args + " {inp} {out}"
    output_pqr = tmp_path / output_pqr
    _LOGGER.debug("Writing output to %s", output_pqr)
    arg_str = arg_str.format(inp=input_pdb, out=output_pqr)
    args = PARSER.parse_args(arg_str.split())
    main(args)
    common.compare_pqr(output_pqr, expected_pqr)


@pytest.mark.slow
def test_slow():
    _LOGGER.error("Need to add slow tests")