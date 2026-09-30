"""Tests for APBS input generation method selection.

Covers the ``mg-para``/``mg-auto`` paths in :class:`inputgen.Elec` that
were left untested when the ``Psize.getSmallest()`` removal broke the
parallel-focusing caller (see issue #492 and PR #474).
"""

from pdb2pqr.inputgen import Elec
from pdb2pqr.psize import Psize


def make_size(ngrid, nsmall=(33, 33, 33)):
    """Build a Psize with hand-set grid attributes (no structure needed)."""
    size = Psize()
    size.ngrid = list(ngrid)
    size.nsmall = list(nsmall)
    size.proc_grid = [2.0, 2.0, 2.0]
    size.coarse_length = [10.0, 10.0, 10.0]
    size.fine_length = [8.0, 8.0, 8.0]
    return size


def gmem(ngrid):
    """Mirror the memory estimate in Elec.__init__."""
    return 200.0 * ngrid[0] * ngrid[1] * ngrid[2] / 1024.0 / 1024.0


def test_auto_selects_para_above_ceiling():
    size = make_size((161, 161, 161))
    assert gmem(size.ngrid) > size.gmemceil
    elec = Elec("mol.pqr", size, "", asyncflag=False)
    assert elec.method == "mg-para"


def test_para_uses_nsmall_for_dime():
    size = make_size((161, 161, 161), nsmall=(97, 97, 97))
    elec = Elec("mol.pqr", size, "mg-para", asyncflag=False)
    assert elec.dime == [97, 97, 97]


def test_para_uses_proc_grid_for_pdime():
    size = make_size((161, 161, 161))
    elec = Elec("mol.pqr", size, "mg-para", asyncflag=False)
    assert elec.pdime == [2.0, 2.0, 2.0]


def test_explicit_para_kept_for_small_grid():
    size = make_size((65, 65, 65))
    assert gmem(size.ngrid) < size.gmemceil
    elec = Elec("mol.pqr", size, "mg-para", asyncflag=False)
    assert elec.method == "mg-para"
    assert elec.dime == [33, 33, 33]


def test_auto_keeps_mg_auto_below_ceiling():
    size = make_size((65, 65, 65))
    assert gmem(size.ngrid) < size.gmemceil
    elec = Elec("mol.pqr", size, "", asyncflag=False)
    assert elec.method == "mg-auto"
    assert elec.dime == [65, 65, 65]


def test_serialized_para_has_parallel_directives():
    size = make_size((161, 161, 161))
    text = str(Elec("mol.pqr", size, "", asyncflag=False))
    assert "mg-para" in text
    assert "pdime 2 2 2" in text
    assert "ofrac 0.1" in text
    assert "dime 33 33 33" in text


def test_serialized_auto_has_focusing_directives():
    size = make_size((65, 65, 65))
    text = str(Elec("mol.pqr", size, "", asyncflag=False))
    assert "mg-auto" in text
    assert "cglen" in text
    assert "fglen" in text
    assert "pdime" not in text
