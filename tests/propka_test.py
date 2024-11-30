"""Tests for PROPKA functionality."""

from pathlib import Path

import common
import pytest


@pytest.mark.parametrize(
    "input_pdb", ["1K1I", "1AFS", "1FAS", "5DV8", "5D8V"], ids=str
)
def test_propka_apo(input_pdb, tmp_path):
    """PROPKA non-regression tests on biomolecules without ligands."""
    args = (
        "--log-level=INFO --ff=AMBER --drop-water "
        "--titration-state-method=propka"
    )
    output_pqr = Path(input_pdb).stem + ".pqr"
    common.run_pdb2pqr_for_tests(
        args=args,
        input_pdb=input_pdb,
        output_pqr=output_pqr,
        tmp_path=tmp_path,
    )


@pytest.mark.parametrize(
    "input_pdb", ["1K1I", "1AFS", "1FAS", "5DV8", "5D8V"], ids=str
)
def test_propka_pka(input_pdb):
    """PROPKA non-regression tests for pKa values on biomolecules without ligands."""
    output_csv = Path("tests/data") / f"{input_pdb}_pka.csv"
    common.run_propka_for_tests(
        input_pdb=input_pdb, compare_file=output_csv, pH=7.4
    )
