"""Regression tests for PDB2PQR behavior."""
import pytest
import logging
import common
from pdb2pqr import cli
from pdb2pqr import main
from pathlib import Path


_LOGGER = logging.getLogger(__name__)
PARSER = cli.build_parser()


def run_pdb2pqr(args, input_pdb, output_pqr, expected_pqr, tmp_path_):
    """Basic code for invoking PDB2PQR."""
    arg_str = args + " {inp} {out}"
    output_pqr = tmp_path_ / output_pqr
    _LOGGER.debug("Writing output to %s", output_pqr)
    arg_str = arg_str.format(inp=input_pdb, out=output_pqr)
    args = PARSER.parse_args(arg_str.split())
    main(args)
    common.compare_pqr(output_pqr, expected_pqr)


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
    run_pdb2pqr(args, input_pdb, output_pqr, expected_pqr, tmp_path)


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
def test_forcefields(args, input_pdb, output_pqr, expected_pqr, tmp_path):
    """Basic code to run 1AFS with --whitespace for different forcefields."""
    run_pdb2pqr(args, input_pdb, output_pqr, expected_pqr, tmp_path)


@pytest.mark.parametrize(
    "args, input_pdb, output_pqr, expected_pqr",
    [
        pytest.param(
            "--log-level=INFO --whitespace --clean",
            common.DATA_DIR / "1AFS.pdb",
            "output.pqr",
            common.DATA_DIR / "1AFS_clean_whitespace.pqr",
            id="1AFS whitespace clean"
        ),
        pytest.param(
            "--log-level=INFO --whitespace --assign-only --ff=AMBER",
            common.DATA_DIR / "1A1P.pdb",
            "output.pqr",
            common.DATA_DIR / '1A1P_assign-only_whitespace_ff=AMBER.pqr',
            id="1A1P assign-only whitespace AMBER"
        ),
        pytest.param(
            "--log-level=INFO --whitespace --nodebump --noopt --ff=AMBER",
            common.DATA_DIR / "1A1P.pdb",
            "output.pqr",
            common.DATA_DIR / '1AFS_nodebump_noopt_whitespace_ff=AMBER.pqr',
            id="1A1P nodebump noopt whitespace AMBER"
        ),
        pytest.param(
            "--ff=AMBER --chain --whitespace",
            common.DATA_DIR / "1AFS.pdb",
            "output.pqr",
            common.DATA_DIR / '1AFS_chain_whitespace_ff=AMBER.pqr',
            id="1AFS chain whitespace AMBER"
        ),
        pytest.param(
            "--ff=PARSE --neutraln --neutralc --whitespace",
            common.DATA_DIR / "1AFS.pdb",
            "output.pqr",
            common.DATA_DIR / '1AFS_neutralc_neutraln_whitespace_ff=PARSE.pqr',
            id="1AFS neutralc neutraln PARSE"
        ),
        pytest.param(
            "--ff=AMBER --drop-water",
            common.DATA_DIR / "1AFS.pdb",
            "output.pqr",
            common.DATA_DIR / '1AFS_drop-water_ff=AMBER.pqr',
            id="1AFS drop-water AMBER"
        ),
        pytest.param(
            "--ff=AMBER --apbs-input --whitespace --include-header",
            common.DATA_DIR / "1AFS.pdb",
            "output.pqr",
            common.DATA_DIR / '1AFS_apbs-input_include-header_whitespace_ff=AMBER.pqr',
            id="1AFS apbs-input include-header whitespace AMBER"
        ),
        pytest.param(
            "--userff={ff} --usernames={names} --whitespace".format(
                ff=common.DATA_DIR / "CUSTOM-FF.DAT",
                names=common.DATA_DIR / "CUSTOM.names"),
            common.DATA_DIR / "1AFS.pdb",
            "output.pqr",
            common.DATA_DIR / '1AFS_userff_usernames_whitespace.pqr',
            id="1AFS userff usernames whitespace"
        )
    ]
)
def test_other_options(args, input_pdb, output_pqr, expected_pqr, tmp_path):
    """Basic code to run 1AFS with --whitespace."""
    run_pdb2pqr(args, input_pdb, output_pqr, expected_pqr, tmp_path)


@pytest.mark.slow
def test_slow():
    _LOGGER.warning("Need to add/mark slow tests")