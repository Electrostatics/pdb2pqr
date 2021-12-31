"""Basic tests of simple core functionality."""
import logging
from pathlib import Path
import pytest
import common


_LOGGER = logging.getLogger(__name__)


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
#: Long test set, primarily protein-based
LONG_PROTEIN_SET = {"5DV8", "1AFS", "1E7G"}
LONG_NUCLEIC_SET = {"5V0O"}
LONG_SET = LONG_PROTEIN_SET | LONG_NUCLEIC_SET
#: Tests that should fail (broken backbones)
BROKEN_SET = {"1EJG", "3U7T"}


_LOGGER.warning("Need functional and regression test coverage for --userff")
_LOGGER.warning("Need functional and regression test coverage for --usernames")
_LOGGER.warning(
    "Need functional and regression test coverage for --apbs-input"
)


@pytest.mark.parametrize("input_pdb", list(SHORT_SET), ids=str)
def test_short_pdb(input_pdb, tmp_path):
    """Non-regression tests on short list of PDB-format biomolecules."""
    args = "--log-level=INFO --ff=AMBER --drop-water --apbs-input=apbs.in"
    output_pqr = Path(input_pdb).stem + ".pqr"
    common.run_pdb2pqr(
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
    common.run_pdb2pqr(
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
    common.run_pdb2pqr(
        args=args,
        input_pdb=input_pdb,
        output_pqr=output_pqr,
        tmp_path=tmp_path,
    )
