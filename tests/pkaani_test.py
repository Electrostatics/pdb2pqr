"""Tests for PKAANI functionality."""

import sys
from pathlib import Path

import common
import pytest

PKAANI_TEST_DIR = Path("tests/pkaani_tests")


@pytest.mark.skipif(
    sys.version_info >= (3, 13),
    reason="pkaani dependencies (torchani/torch) do not support Python 3.13 yet",
)
@pytest.mark.parametrize(
    "input_pdb",
    [
        str(PKAANI_TEST_DIR / "1brs_clean.pdb"),
        str(PKAANI_TEST_DIR / "1zy8_clean.pdb"),
        str(PKAANI_TEST_DIR / "6oge_de.pdb"),
        "1FCC.pdb",
    ],
    ids=str,
)
def test_pkaani_apo(input_pdb, tmp_path):
    """PKAANI non-regression tests on biomolecules without ligands."""
    args = (
        "--log-level=INFO --ff=AMBER --drop-water "
        "--titration-state-method=pkaani"
    )
    output_pqr = Path(input_pdb).stem + ".pqr"
    common.run_pdb2pqr_for_tests(
        args=args,
        input_pdb=input_pdb,
        output_pqr=output_pqr,
        tmp_path=tmp_path,
    )


@pytest.mark.skipif(
    sys.version_info >= (3, 13),
    reason="pkaani dependencies (torchani/torch) do not support Python 3.13 yet",
)
@pytest.mark.parametrize(
    "input_pdb",
    [
        str(PKAANI_TEST_DIR / "1brs_clean.pdb"),
        str(PKAANI_TEST_DIR / "1zy8_clean.pdb"),
        str(PKAANI_TEST_DIR / "6oge_de.pdb"),
        "1FCC.pdb",
    ],
    ids=str,
)
def test_pkaani_pka(input_pdb):
    """PKAANI non-regression tests for pKa values on biomolecules without ligands."""
    pdb_filename = Path(input_pdb).name
    output_csv = PKAANI_TEST_DIR / f"{pdb_filename}_pka.csv"
    common.run_pkaani_for_tests(
        input_pdb=input_pdb, compare_file=output_csv, pH=7.4
    )
