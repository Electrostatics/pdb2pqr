"""Basic PDB2PQR tests.

Basic tests do not include ligand support or titration state assignment.
"""
import logging
import random
import pytest
import common


_LOGGER = logging.getLogger(__name__)


# This is the value passed to the --log-level option (not being tested separately)
LOG_LEVEL = "INFO"
# These are PDB files that should pose few problems with reconstruction
GOOD_PDBS = ["1K1I", "1AFS", "1FAS", "5DV8", "5D8V"]
# We are missing some tests
_LOGGER.error("Need --assign-only test coverage.")
_LOGGER.error("Need --userff test coverage.")
_LOGGER.error("Need --usernames test coverage.")
_LOGGER.error("Need --ligand test coverage")
_LOGGER.error("Need --titration-state-method test coverage")
_LOGGER.error("Need --with-ph test coverage")
_LOGGER.error("Need --pdb2pka-out test coverage")
_LOGGER.error("Need --pdb2pka-resume test coverage")
_LOGGER.error("Need --pdie test coverage")
_LOGGER.error("Need --sdie coverage")
_LOGGER.error("Need --pairene coverage")
_LOGGER.error("Need --propka-reference test coverage")
_LOGGER.error("Need --apbs-input test coverage")
_LOGGER.error("Need to implement regression tests")


# Combinations of force field options
FF_COMBOS = common.generate_ff_combinations()
_LOGGER.info("Found %d force field option combinations", len(FF_COMBOS))
# These are options that speed up execution when added to list
FAST_OPTIONS = ["--clean", "--nodebump", "--noopt"]
# These are other options to test
SIMPLE_OPTIONS = ["--whitespace", "--neutraln", "--neutralc", "--drop-water",
                  "--include-header"] + FAST_OPTIONS
SIMPLE_COMBOS = common.generate_combinations(SIMPLE_OPTIONS)
_LOGGER.info("Found %d simple option combinations", len(SIMPLE_COMBOS))


# Generate an "exhaustive" list of options
OPTION_LIST = []
for ff_opts in FF_COMBOS:
    for simp_opts in SIMPLE_COMBOS:
        if "--ff=PARSE" not in ff_opts:
            opts = ff_opts + list(simp_opts)
        else:
            opts = ff_opts + list(set(simp_opts) - {"--neutralc", "--neutraln"})
        OPTION_LIST += [opts]
_LOGGER.info("Found %d option combinations", len(OPTION_LIST))


# Add PDBs to test list
_LOGGER.info("Testing with small group of PDBs: %s", GOOD_PDBS)
PARAM_LIST = []
for pdb in GOOD_PDBS:
    for opt in OPTION_LIST:
        marks = []
        if len(set(opt) & set(FAST_OPTIONS)) > 0:
            marks.append(pytest.mark.fast)
        args = " ".join(opt)
        id_ = "%s %s" % (pdb, args)
        PARAM_LIST += [{"args": args, "pdb": common.DATA_DIR/pdb, "id": id_,
                        "marks": marks}]
_LOGGER.info("Generated %d total tests.", len(PARAM_LIST))


# Generate test subsets
_LOGGER.info("Generating test subsets")
for param in random.sample(PARAM_LIST, 10):
    param["marks"].append(pytest.mark.random_10)
for param in random.sample(PARAM_LIST, 100):
    param["marks"].append(pytest.mark.random_100)
for param in random.sample(PARAM_LIST, 1000):
    param["marks"].append(pytest.mark.random_1000)
for iparam, param in enumerate(PARAM_LIST):
    if iparam % 11 == 0:
        param["marks"].append(pytest.mark.every_11)
    if iparam % 101 == 0:
        param["marks"].append(pytest.mark.every_101)
    if iparam % 1009 == 0:
        param["marks"].append(pytest.mark.every_1009)
    if iparam % 10007 == 0:
        param["marks"].append(pytest.mark.every_10007)


# Generate parameter list
_LOGGER.info("Generating parameter list")
EXHAUST_PARAMS = [pytest.param(p["args"], p["pdb"],
                               id=p["id"], marks=p["marks"]) for p in PARAM_LIST]


@pytest.mark.parametrize("args, input_pdb", EXHAUST_PARAMS)
def test_exhaustive(args, input_pdb, tmp_path):
    """Run an exhaustive test of PDB2PQR functionality.

    See error messages above for functionality we're not testing.
    """
    common.run_pdb2pqr(args=args, input_pdb=input_pdb, tmp_path=tmp_path)
