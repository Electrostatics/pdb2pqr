"""Basic tests of simple core functionality."""

from pathlib import Path

import common
import pytest

# fmt: off
#: Protein-nucleic acid complexes
PROTEIN_NUCLEIC_SET = {"4UN3"}
#: Proteins without nucleic acids
PROTEIN_ONLY_SET = {"1K1I", "1FAS", "5D8V"}
#: DNA without protein
DNA_ONLY_SET = {"1NAJ", "7BNA"}
#: RNA without protein
RNA_ONLY_SET = {"1AQO", "4E8M"}
#: Nucleic acids without proteins
NUCLEIC_ONLY_SET = DNA_ONLY_SET | RNA_ONLY_SET
#: Basic test set
SHORT_SET = PROTEIN_NUCLEIC_SET | PROTEIN_ONLY_SET | NUCLEIC_ONLY_SET
#: Long test set of known good structures
LONG_PROTEIN_SET = {
    "5DV8", "1AFS", "1E7G", "3H0X", "1DUP", "4RWU", "3JTE", "3CP1", "1TG0",
    "1USM", "2OO2", "2J6B", "2SFA", "3US6", "3CSR", "3V1A", "3LAA", "3LS0",
    "2OP6", "2NWD", "3ONJ", "1Z0P", "4YPC", "3CQT", "1YQB", "1OX3", "2P5K",
    "1R7J", "1H75", "1YIB", "4Q2Q", "4GQM", "1UG4", "4TTW", "2JIC", "4GS3",
    "3I1G", "4PGR", "1HOE", "2CKX", "1HUF", "1BM8", "2G7O", "1OA4", "3RJS",
    "1TEN", "2HC8", "3KZD", "3FH2", "1X91", "3T8R", "3B7H", "4N6T", "1J2A",
    "3QC7", "4S11", "3NGP", "2IGD", "2NYC", "1NXO", "4OSN", "3SO2", "2J70",
    "3ZSL", "2IVY", "2NSN", "4IL7", "1NH9", "2EBB", "2REY", "2JCP", "2CMP",
    "3ZHI", "2YXF", "1TUD", "4ZQA", "1X3O", "4O7Q", "3I35", "3CE7", "2CJJ",
    "2QR3", "1TGN", "3KT2", "1YZM", "1FAS", "2X4J", "1TIG", "1L2P", "4EO0",
    "2IBL", "1JQ0", "2QVO", "4G3O", "2ZQE", "2FD4", "1DF4", "2PMR", "2PII",
    "3HRO", "2D8E", "2PCY", "4GMQ", "1J27", "4CVD", "1OPC", "1P1L", "1I2T",
    "4NPN", "2J9V", "2CWR", "1ZX6", "4L9E", "1RJ1", "3ISU", "3RFI", "4LJ1",
    "3TSV", "3HRN", "2HP7", "2SGA", "3C97", "1R5Q", "1NKD", "2IWN", "1TQG",
    "3QE1", "1X6J", "1PZC", "2VWR", "1PDO", "2BYG", "1ZZK", "2YGS", "2ERW",
    "4HTI", "3LLB", "2FI9", "1P9I", "4I2T", "2GUS", "3T1S", "1V2Z", "2QBV",
    "3ZZP", "1PHT", "4HS2", "1USE", "5EE2", "4YKD", "4CFI", "1WOU", "3RFF",
    "1NOA", "2CGQ", "2TGI", "2R6Q", "1IS5", "4CRH", "1B0X", "2FWG", "2HDZ",
    "3WFV", "4AXT", "2ERA", "4ZMK", "1ZVA", "3S0A", "2PTV", "2FRG", "4POY",
    "1PGX", "1LN4", "1SJV", "2ESK", "2PND", "1XAK", "3NX6", "1TVQ", "3RVC",
    "1RIS", "1AIL", "3IDW", "1KFN", "3LMO", "1OPS", "2H2C", "3DVI", "2GZV",
    "2XXC", "1ZLM", "2OCH",
}
LONG_NUCLEIC_SET = {"5V0O"}
LONG_SET = LONG_PROTEIN_SET | LONG_NUCLEIC_SET
#: Tests that should fail (broken backbones)
BROKEN_SET = {"1EJG", "3U7T", "4MGP", "2V75"}
#: Tests for naming residues with different protonation states
NAMING_TESTS = []
for ff in ["CHARMM", "AMBER", "PARSE", "SWANSON", "TYL06"]:
    for pdb in ["1AJJ", "1BX8"]:
        for pH in [2, 7, 14]:
            options = {"pdb": pdb, "ff": ff, "pH": pH}
            if ff == "CHARMM":
                NAMING_TESTS.append(pytest.param(options, marks=pytest.mark.xfail(reason="CHARMM force field is broken!")))
            else:
                NAMING_TESTS.append(pytest.param(options))
# fmt: on


@pytest.mark.parametrize("input_pdb", list(SHORT_SET), ids=str)
def test_short_pdb(input_pdb, tmp_path):
    """Non-regression tests on short list of PDB-format biomolecules."""
    args = "--log-level=INFO --ff=AMBER --drop-water --apbs-input=apbs.in"
    output_pqr = Path(input_pdb).stem + ".pqr"
    common.run_pdb2pqr_for_tests(
        args=args,
        input_pdb=input_pdb,
        output_pqr=output_pqr,
        tmp_path=tmp_path,
    )


@pytest.mark.parametrize("input_pdb", list(SHORT_SET), ids=str)
def test_basic_cif(input_pdb, tmp_path):
    """Non-regression tests on short list of CIF-format biomolecules."""
    args = "--log-level=INFO --ff=AMBER --drop-water --apbs-input=apbs.in"
    output_pqr = Path(input_pdb).stem + ".cif"
    common.run_pdb2pqr_for_tests(
        args=args,
        input_pdb=input_pdb,
        output_pqr=output_pqr,
        tmp_path=tmp_path,
    )


@pytest.mark.parametrize("input_pdb", list(LONG_SET), ids=str)
@pytest.mark.long_test
def test_long_pdb(input_pdb, tmp_path):
    """Non-regression tests on short list of PDB-format biomolecules."""
    args = "--log-level=INFO --ff=AMBER --drop-water --apbs-input=apbs.in"
    output_pqr = Path(input_pdb).stem + ".pqr"
    common.run_pdb2pqr_for_tests(
        args=args,
        input_pdb=input_pdb,
        output_pqr=output_pqr,
        tmp_path=tmp_path,
    )


@pytest.mark.parametrize("input_pdb", list(BROKEN_SET), ids=str)
@pytest.mark.xfail
def test_broken_backbone(input_pdb, tmp_path):
    """Test graceful failure of optimization with missing backbone atoms."""
    args = "--log-level=INFO --ff=AMBER --drop-water"
    output_pqr = Path(input_pdb).stem + ".pqr"
    common.run_pdb2pqr_for_tests(
        args=args,
        input_pdb=input_pdb,
        output_pqr=output_pqr,
        tmp_path=tmp_path,
    )


@pytest.mark.parametrize(
    "input_pdb, expected_pqr",
    [pytest.param("cterm_hid.pdb", "cterm_hid_out.pqr", id="C-terminal HID")],
)
def test_protonated_terminals(input_pdb, expected_pqr, tmp_path):
    """Tests for terminal residue protonation."""
    args = "--log-level=INFO --ff=AMBER --ffout AMBER"
    output_pqr = Path(input_pdb).stem + ".pqr"
    common.run_pdb2pqr_for_tests(
        args=args,
        input_pdb=common.DATA_DIR / input_pdb,
        output_pqr=output_pqr,
        expected_pqr=common.DATA_DIR / expected_pqr,
        tmp_path=tmp_path,
        compare_resnames=True,
    )


@pytest.mark.parametrize(
    "input_pdb, expected_pqr",
    [
        pytest.param(
            "5vav_cyclic_peptide.pdb",
            "5vav_cyclic_peptide_out.pqr",
            id="Cyclic peptide",
        )
    ],
)
def test_cyclic_peptide(input_pdb, expected_pqr, tmp_path):
    """Tests for cyclic peptide protonation."""
    args = "--log-level=INFO --ff=AMBER --ffout AMBER"
    output_pqr = Path(input_pdb).stem + ".pqr"
    common.run_pdb2pqr_for_tests(
        args=args,
        input_pdb=common.DATA_DIR / input_pdb,
        output_pqr=output_pqr,
        expected_pqr=common.DATA_DIR / expected_pqr,
        tmp_path=tmp_path,
        compare_resnames=True,
    )


@pytest.mark.parametrize("naming_test", NAMING_TESTS, ids=str)
def test_ph_naming(naming_test, tmp_path):
    """Non-regression tests on naming schemes at different pH values."""
    args = "--log-level=INFO --ff=AMBER --drop-water --apbs-input=apbs.in"
    input_pdb = naming_test["pdb"] + ".pdb"
    output_pqr = Path(input_pdb).stem + ".pqr"
    expected_pqr = (
        f"{naming_test['pdb']}_pH{naming_test['pH']}_{naming_test['ff']}.pqr"
    )
    args = (
        f"--log-level=INFO --ff={naming_test['ff']} --ffout={naming_test['ff']} "
        f"--drop-water --whitespace --with-ph={naming_test['pH']} "
        f"--titration-state-method=propka"
    )
    common.run_pdb2pqr_for_tests(
        args=args,
        input_pdb=common.DATA_DIR / input_pdb,
        output_pqr=output_pqr,
        expected_pqr=common.DATA_DIR / expected_pqr,
        tmp_path=tmp_path,
        compare_resnames=True,
    )
