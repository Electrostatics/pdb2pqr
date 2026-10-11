"""Tests for biomolecule behavior."""

from types import SimpleNamespace

import pytest

from pdb2pqr import aa, biomolecule, main


def _make_asp(ins_code):
    residue = object.__new__(aa.ASP)
    residue.name = "ASP"
    residue.res_seq = 100
    residue.ins_code = ins_code
    residue.chain_id = "A"
    residue.is_n_term = False
    residue.is_c_term = False
    return residue


def test_apply_pka_values_distinguishes_insertion_codes():
    """Residues differing only by insertion code receive distinct pKas."""
    plain = _make_asp("")
    inserted = _make_asp("A")
    molecule = object.__new__(biomolecule.Biomolecule)
    molecule.residues = [plain, inserted]
    applied = []
    molecule.apply_patch = lambda patchname, residue: applied.append(
        (patchname, residue.ins_code)
    )
    pka_values = {
        biomolecule.pka_key("ASP", 100, "A", ""): 6.0,
        biomolecule.pka_key("ASP", 100, "A", "A"): 8.0,
    }

    molecule.apply_pka_values("parse", 7.0, pka_values)

    assert applied == [("ASH", "A")]
    assert pka_values == {}


def test_pka_key_preserves_legacy_format_without_insertion_code():
    """Blank and CIF-null insertion codes keep the historical key format."""
    expected = "ASP 100 A"
    assert biomolecule.pka_key("ASP", 100, "A") == expected
    assert biomolecule.pka_key("ASP", 100, "A", " ") == expected
    assert biomolecule.pka_key("ASP", 100, "A", ".") == expected
    assert biomolecule.pka_key("ASP", 100, "A", "?") == expected


def test_pkaani_rejects_residue_insertion_codes():
    """pKa-ANI fails before collapsing insertion-coded residue identities."""
    molecule = SimpleNamespace(residues=[_make_asp("A")])

    with pytest.raises(RuntimeError, match="cannot distinguish"):
        main.run_pkaani(SimpleNamespace(), molecule)
