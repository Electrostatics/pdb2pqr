"""Tests for biomolecule behavior."""

import pytest

from pdb2pqr import aa, biomolecule


def _make_asp(ins_code):
    residue = object.__new__(aa.ASP)
    residue.name = "ASP"
    residue.res_seq = 100
    residue.ins_code = ins_code
    residue.chain_id = "A"
    residue.is_n_term = False
    residue.is_c_term = False
    return residue


def test_set_reference_distance_matches_shortest_path():
    """Intra-residue bonds are symmetric and CA BFS matches shortest_path."""
    defn = biomolecule.io.get_definitions()
    for res in defn.map.values():
        if "CA" in res.map:
            for name, atom in res.map.items():
                for bond in atom.bonds:
                    if bond in res.map:
                        assert name in res.map[bond].bonds

    pdb_path = biomolecule.io.Path(__file__).resolve().parent / "data/1AFS.pdb"
    pdblist, _ = biomolecule.io.get_molecule(pdb_path)
    mol = biomolecule.Biomolecule(pdblist, defn)
    mol.set_termini()
    mol.update_bonds()
    mol.add_hydrogens()
    mol.set_reference_distance()
    for residue in mol.residues:
        if not isinstance(residue, aa.Amino):
            continue
        caatom = residue.get_atom("CA")
        map_ = {atom: atom.bonds for atom in residue.atoms}
        for atom in residue.atoms:
            for nbr in atom.bonds:
                assert atom in nbr.bonds
            if atom.is_backbone or (residue.is_c_term and atom.name == "HO"):
                continue
            if residue.is_n_term and atom.name in ("H3", "H2"):
                continue
            path = biomolecule.util.shortest_path(map_, atom, caatom)
            assert path is not None
            assert atom.refdistance == len(path) - 1

    target_atom = next(
        atom
        for residue in mol.residues
        if isinstance(residue, aa.Amino)
        for atom in residue.atoms
        if not atom.is_backbone
        and not (residue.is_c_term and atom.name == "HO")
        and not (residue.is_n_term and atom.name in ("H3", "H2"))
    )
    for nbr in list(target_atom.bonds):
        nbr.bonds.remove(target_atom)
    target_atom.bonds.clear()
    with pytest.raises(ValueError, match="Found gap in biomolecule structure"):
        mol.set_reference_distance()


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
