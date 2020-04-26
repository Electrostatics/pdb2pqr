"""Regression tests for PDB2PQR behavior."""
import pytest
import logging
import common
from pdb2pqr import cli
from pdb2pqr import main
from pathlib import Path


_LOGGER = logging.getLogger(__name__)
PARSER = cli.build_parser()


def run(input_pdb, output_pqr, ref_pqr, arg_str):
    """Basic code to run 1AFS."""
    arg_str = arg_str + " {inp} {out}"
    arg_str = arg_str.format(inp=input_pdb, out=output_pqr)
    args = PARSER.parse_args(arg_str.split())
    main(args)
    common.compare_pqr(output_pqr, ref_pqr)


def test_1AFS_local(tmpdir):
    """Test basic options for 1AFS using local copy of PDB."""
    ref_pqr = common.DATA_DIR / "1AFS_ff=AMBER.pqr"
    output_pqr = Path(tmpdir) / "output.pqr"
    input_pdb = common.DATA_DIR / "1AFS.pdb"
    arg_str = "--log-level=INFO --ff=AMBER"
    run(input_pdb, output_pqr, ref_pqr, arg_str)


def test_1AFS_remote(tmpdir):
    """Test basic options for 1AFS using remote copy of PDB."""
    ref_pqr = common.DATA_DIR / "1AFS_ff=AMBER.pqr"
    output_pqr = Path(tmpdir) / "output.pqr"
    input_pdb = "1AFS"
    arg_str = "--log-level=INFO --ff=AMBER"
    run(input_pdb, output_pqr, ref_pqr, arg_str)
