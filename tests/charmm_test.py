"""Tests for CHARMM parameter assignment."""

import math

import common
import pytest


@pytest.mark.parametrize("ffout", [False, True], ids=["ff", "ffout"])
def test_rna_2prime_hydrogens(ffout, tmp_path):
    """The O2' and C2' hydrogens of RNA get their own CHARMM charges (#496).

    CHARMM names the O2' hydrogen H2' and the C2' hydrogen H2'', the
    opposite of the NA.xml names, so the charges were swapped.
    """
    args = "--ff=CHARMM --whitespace" + (" --ffout=CHARMM" if ffout else "")
    common.run_pdb2pqr_for_tests(
        args=args,
        input_pdb=common.DATA_DIR / "1AQO_3nt_noH.pdb",
        output_pqr="rna.pqr",
        tmp_path=tmp_path,
    )
    atoms = {}
    with open(tmp_path / "rna.pqr", encoding="utf-8") as pqr_file:
        for line in pqr_file:
            if line.startswith("ATOM"):
                fields = line.split()
                coords = tuple(map(float, fields[-5:-2]))
                atoms[(fields[-6], fields[2])] = (coords, float(fields[-2]))
    hydrogens = [
        (resid, name) for resid, name in atoms if name in ("H2'", "H2''")
    ]
    assert len(hydrogens) == 6
    for resid, name in hydrogens:
        coords, charge = atoms[(resid, name)]
        to_o2 = math.dist(coords, atoms[(resid, "O2'")][0])
        to_c2 = math.dist(coords, atoms[(resid, "C2'")][0])
        expected = 0.43 if to_o2 < to_c2 else 0.09
        assert charge == pytest.approx(expected), (resid, name)
