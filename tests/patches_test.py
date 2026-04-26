from pathlib import Path

from pdb2pqr.definitions import Definition


def test_dt5_patch_applies_to_dt_reference() -> None:
    dat_dir = Path(__file__).resolve().parents[1] / "pdb2pqr" / "dat"

    with (
        (dat_dir / "AA.xml").open() as aa_file,
        (dat_dir / "NA.xml").open() as na_file,
        (dat_dir / "PATCHES.xml").open() as patch_file,
    ):
        definition = Definition(aa_file, na_file, patch_file)

    assert "DT" in definition.map
    assert "DT5" in definition.map
    assert "P" not in definition.map["DT5"].map
    assert "O1P" not in definition.map["DT5"].map
    assert "H5T" in definition.map["DT5"].map
