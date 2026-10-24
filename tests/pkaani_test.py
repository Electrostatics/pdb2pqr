"""Tests for PKAANI functionality."""

import subprocess
import sys
from pathlib import Path

import common
import pytest
from torch.utils._filelock import FileLock

PKAANI_TEST_DIR = Path("tests/pkaani_tests")


def test_torch_filelock_excludes_other_processes(tmp_path):
    """Torch's lock wrapper preserves exclusion and releases after errors."""
    lock_path = tmp_path / "torch.lock"
    command = [
        sys.executable,
        "-c",
        (
            "import sys\n"
            "from filelock import FileLock, Timeout\n"
            "try:\n"
            "    with FileLock(sys.argv[1], timeout=0.1):\n"
            "        pass\n"
            "except Timeout:\n"
            "    sys.exit(23)\n"
        ),
        str(lock_path),
    ]
    lock = FileLock(str(lock_path), timeout=1)
    with pytest.raises(ValueError, match="release on error"):
        with lock:
            assert lock.is_locked
            contender = subprocess.run(
                command, capture_output=True, text=True, timeout=30
            )
            assert contender.returncode == 23, contender.stderr
            raise ValueError("release on error")

    assert not lock.is_locked
    released = subprocess.run(
        command, capture_output=True, text=True, timeout=30
    )
    assert released.returncode == 0, released.stderr


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
