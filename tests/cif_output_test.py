"""Tests for large-assembly CIF support and output-format selection.

Covers the extensions added for structures that exceed the fixed PDB columns
(> 99999 atoms, multi-character chain ids, residue numbers outside
[-999, 9999]):

* the limit-free CIF reader (:func:`pdb2pqr.cif.read_cif`),
* the ``exceeds_pdb_limits`` gate and the mmCIF / free-format PQR writers
  (:mod:`pdb2pqr.io`, :class:`pdb2pqr.structures.Atom`),
* output-format resolution and the PROPKA / pKa-ANI overflow guards
  (:mod:`pdb2pqr.main`),
* ``psize`` parsing of free-format (whitespace-delimited) PQR.
"""

import argparse

import pytest

from pdb2pqr import cif, io, pdb, psize
from pdb2pqr.structures import Atom

# ``main`` pulls in propka at import time; keep the reader/writer tests usable
# even if propka is unavailable by importing it lazily.
try:
    from pdb2pqr import main

    HAVE_MAIN = True
except Exception:  # pragma: no cover - depends on optional propka install
    main = None
    HAVE_MAIN = False

needs_main = pytest.mark.skipif(not HAVE_MAIN, reason="requires propka/main")


def make_atom(
    serial=1,
    chain="A",
    res_seq=1,
    name="N",
    res_name="THR",
    x=1.234,
    y=5.678,
    z=9.012,
    charge=-0.32,
    radius=1.5,
):
    """Build a fully-populated in-memory :class:`Atom` for writer tests."""
    atom = Atom(type_="ATOM")
    atom.serial = serial
    atom.chain_id = chain
    atom.res_seq = res_seq
    atom.name = name
    atom.res_name = res_name
    atom.x, atom.y, atom.z = x, y, z
    atom.ffcharge = charge
    atom.radius = radius
    atom.occupancy = 1.0
    atom.temp_factor = 0.0
    atom.element = name[0]
    # Real atoms carry empty-string (not None) insertion code / alt loc; the
    # fixed-column writer renders a None ins_code as the literal "None".
    atom.ins_code = ""
    atom.alt_loc = ""
    return atom


def make_args(output_format="auto", keep_chain=False):
    return argparse.Namespace(
        output_format=output_format, keep_chain=keep_chain
    )


# ---------------------------------------------------------------------------
# exceeds_pdb_limits
# ---------------------------------------------------------------------------


def test_exceeds_pdb_limits_fits():
    exceeds, reasons = io.exceeds_pdb_limits([make_atom()])
    assert exceeds is False
    assert reasons == []


@pytest.mark.parametrize(
    "atoms, needle",
    [
        ([make_atom(chain="AA")], "chain"),
        ([make_atom(serial=100001)], "99999"),
        ([make_atom(res_seq=10001)], "9999"),
        ([make_atom(res_seq=-1000)], "-999"),
        # count > 99999 without building 100k distinct objects
        ([make_atom()] * 100001, "99999"),
    ],
)
def test_exceeds_pdb_limits_overflow(atoms, needle):
    exceeds, reasons = io.exceeds_pdb_limits(atoms)
    assert exceeds is True
    assert any(needle in reason for reason in reasons)


# ---------------------------------------------------------------------------
# Atom writers
# ---------------------------------------------------------------------------


def test_free_pqr_string_untruncated_fields():
    atom = make_atom(serial=100001, chain="A", res_seq=10001)
    tokens = atom.get_free_pqr_string(include_chain=False).split()
    # ATOM serial name resname resseq x y z charge radius  (no chain)
    assert tokens[0] == "ATOM"
    assert int(tokens[1]) == 100001
    assert int(tokens[4]) == 10001  # residue number, untruncated
    assert float(tokens[-2]) == pytest.approx(-0.32, abs=1e-4)
    assert float(tokens[-1]) == pytest.approx(1.5, abs=1e-4)


def test_free_pqr_string_chain_toggle():
    atom = make_atom(chain="A", res_seq=7)
    with_chain = atom.get_free_pqr_string(include_chain=True).split()
    without = atom.get_free_pqr_string(include_chain=False).split()
    assert "A" in with_chain
    assert len(with_chain) == len(without) + 1


def test_cif_atom_dict_auth_label_equal_and_untruncated():
    row = make_atom(
        serial=123456, chain="AA", res_seq=10001
    ).get_cif_atom_dict()
    # auth_* and label_* must agree so the PROPKA round-trip matches
    assert row["auth_asym_id"] == row["label_asym_id"] == "AA"
    assert row["auth_seq_id"] == row["label_seq_id"] == "10001"
    assert row["auth_comp_id"] == row["label_comp_id"]
    assert row["auth_atom_id"] == row["label_atom_id"]
    assert row["id"] == "123456"  # untruncated serial


# ---------------------------------------------------------------------------
# _cif_token quoting
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("name", ["O5'", "C2'", "H5''", "N", "CA", "AA"])
def test_cif_token_leaves_bare_values_bare(name):
    # A value with no whitespace and no reserved leading character is valid
    # bare -- including nucleic-acid atom names whose apostrophe is embedded
    # (NOT leading).  It must not be quoted.
    assert io._cif_token(name) == name


def test_cif_token_quotes_with_nonconflicting_delimiter():
    # Whitespace forces quoting; single quotes are fine when the value has none.
    assert io._cif_token("a b") == "'a b'"
    # A value *starting* with a quote must pick the other delimiter rather than
    # emit a broken ''x' token.
    assert io._cif_token("'x") == '"' + "'x" + '"'
    # Leading ';' opens a text field at line start, so it must be quoted.
    assert io._cif_token(";x") == "';x'"


def test_cif_token_empty_is_null():
    assert io._cif_token("") == "."


def test_cif_prime_atom_name_roundtrips_via_gemmi():
    # End-to-end: a bare O5' in the loop must read back unchanged.
    gemmi = pytest.importorskip("gemmi")
    atoms = [make_atom(name="O5'", res_name="G")]
    text = "".join(io.print_biomolecule_atoms_cif(atoms))
    block = gemmi.cif.read_string(text).sole_block()
    assert list(block.find_values("_atom_site.label_atom_id")) == ["O5'"]


# ---------------------------------------------------------------------------
# io writers (loop emitters)
# ---------------------------------------------------------------------------


def test_print_free_omits_multichar_chain(caplog):
    atoms = [
        make_atom(chain="AA", res_seq=1),
        make_atom(chain="AA", res_seq=2),
    ]
    import logging

    with caplog.at_level(logging.WARNING):
        lines = io.print_biomolecule_atoms_free(atoms, keep_chain=True)
    atom_lines = [ln for ln in lines if ln.startswith("ATOM")]
    # multi-char chain -> chain omitted -> 10 tokens (no chain field)
    assert all(len(ln.split()) == 10 for ln in atom_lines)
    assert any("chain" in rec.message.lower() for rec in caplog.records)


def test_print_free_keeps_single_char_chain():
    atoms = [make_atom(chain="A")]
    lines = io.print_biomolecule_atoms_free(atoms, keep_chain=True)
    atom_lines = [ln for ln in lines if ln.startswith("ATOM")]
    assert atom_lines
    assert atom_lines[0].split()[4] == "A"


def test_print_cif_reparses_with_gemmi():
    gemmi = pytest.importorskip("gemmi")
    atoms = [
        make_atom(serial=100001, chain="AA", res_seq=10001, name="N"),
        make_atom(serial=100002, chain="AA", res_seq=10001, name="CA"),
    ]
    text = "".join(io.print_biomolecule_atoms_cif(atoms))
    block = gemmi.cif.read_string(text).sole_block()
    ids = block.find_values("_atom_site.id")
    auth = list(block.find_values("_atom_site.auth_asym_id"))
    label = list(block.find_values("_atom_site.label_asym_id"))
    assert len(ids) == 2
    assert auth == label == ["AA", "AA"]
    # PQR charge/radius are carried as dedicated columns and round-trip
    assert "pdb2pqr_charge" in io.CIF_ATOM_SITE_COLUMNS
    charges = list(block.find_values("_atom_site.pdb2pqr_charge"))
    assert float(charges[0]) == pytest.approx(-0.32, abs=1e-4)


# ---------------------------------------------------------------------------
# resolve_output_format
# ---------------------------------------------------------------------------


@needs_main
@pytest.mark.parametrize(
    "fmt, atoms, expected",
    [
        ("auto", [make_atom()], "pdb"),
        ("auto", [make_atom(chain="AA")], "free"),
        ("pdb", [make_atom()], "pdb"),
        ("free", [make_atom()], "free"),
        ("cif", [make_atom()], "cif"),
        ("cif", [make_atom(chain="AA")], "cif"),
    ],
)
def test_resolve_output_format(fmt, atoms, expected):
    assert main.resolve_output_format(make_args(fmt), atoms) == expected


@needs_main
def test_resolve_output_format_pdb_overflow_errors():
    with pytest.raises(RuntimeError, match="PDB-format"):
        main.resolve_output_format(make_args("pdb"), [make_atom(chain="AA")])


# ---------------------------------------------------------------------------
# CIF reader: limit-free + robustness
# ---------------------------------------------------------------------------

_BIG_ATOM_SITE = """data_big
loop_
_atom_site.group_PDB
_atom_site.id
_atom_site.type_symbol
_atom_site.label_atom_id
_atom_site.label_alt_id
_atom_site.label_comp_id
_atom_site.label_asym_id
_atom_site.auth_seq_id
_atom_site.pdbx_PDB_ins_code
_atom_site.Cartn_x
_atom_site.Cartn_y
_atom_site.Cartn_z
_atom_site.occupancy
_atom_site.B_iso_or_equiv
_atom_site.pdbx_PDB_model_num
ATOM 100001 N N . THR AA 10001 ? 1.234 5.678 9.012 1.00 0.00 1
ATOM 100002 C CA . THR AA 10001 ? 2.234 6.678 10.012 1.00 0.00 1
"""

# struct_conn with a multi-character chain used to crash pdb.SSBOND
# (int(line[17:21]) on shifted columns). It must now be ignored entirely.
_STRUCT_CONN_TRAP = """
loop_
_struct_conn.id
_struct_conn.conn_type_id
_struct_conn.ptnr1_auth_comp_id
_struct_conn.ptnr1_auth_asym_id
_struct_conn.ptnr1_auth_seq_id
_struct_conn.pdbx_ptnr1_PDB_ins_code
_struct_conn.ptnr2_auth_comp_id
_struct_conn.ptnr2_auth_asym_id
_struct_conn.ptnr2_auth_seq_id
_struct_conn.pdbx_ptnr2_PDB_ins_code
_struct_conn.ptnr1_symmetry
_struct_conn.ptnr2_symmetry
_struct_conn.pdbx_dist_value
disulf1 disulf CYS AAA 10001 ? CYS AAA 10055 ? 1_555 1_555 2.03
"""


def _read_cif_text(tmp_path, text, name="test.cif"):
    path = tmp_path / name
    path.write_text(text)
    with open(path) as handle:
        return cif.read_cif(handle)


def test_reader_limit_free(tmp_path):
    pdblist, _ = _read_cif_text(tmp_path, _BIG_ATOM_SITE)
    atoms = [a for a in pdblist if isinstance(a, (pdb.ATOM, pdb.HETATM))]
    assert len(atoms) == 2
    first = atoms[0]
    assert first.serial == 100001  # not truncated to 5 columns
    assert first.chain_id == "AA"  # not truncated to 1 column
    assert first.res_seq == 10001  # not truncated to 4 columns


def test_reader_ignores_struct_conn_trap(tmp_path):
    # Regression: a multi-character chain in struct_conn previously crashed
    # pdb.SSBOND parsing. The connectivity records are now not parsed at all.
    pdblist, _ = _read_cif_text(tmp_path, _BIG_ATOM_SITE + _STRUCT_CONN_TRAP)
    atoms = [a for a in pdblist if isinstance(a, (pdb.ATOM, pdb.HETATM))]
    assert len(atoms) == 2
    assert not any(isinstance(rec, pdb.SSBOND) for rec in pdblist)


def test_reader_missing_header_categories_nonfatal(tmp_path):
    # atom_site-only file: header categories are absent; read must not raise.
    pdblist, _ = _read_cif_text(tmp_path, _BIG_ATOM_SITE)
    assert any(isinstance(a, (pdb.ATOM, pdb.HETATM)) for a in pdblist)


# ---------------------------------------------------------------------------
# psize: free-format PQR parsing
# ---------------------------------------------------------------------------


def test_psize_parses_free_format_untruncated():
    atoms = [
        make_atom(
            serial=100001, res_seq=10001, x=10.0, y=0.0, z=0.0, radius=2.0
        ),
        make_atom(
            serial=100002, res_seq=10001, x=-10.0, y=0.0, z=0.0, radius=2.0
        ),
    ]
    text = "".join(io.print_biomolecule_atoms_free(atoms, keep_chain=False))
    size = psize.Psize()
    size.parse_string(text)
    # x extent spans -10-2 .. 10+2
    assert size.maxlen[0] == pytest.approx(12.0, abs=1e-3)
    assert size.minlen[0] == pytest.approx(-12.0, abs=1e-3)


def test_psize_fixed_and_free_agree():
    atoms = [
        make_atom(serial=1, res_seq=1, x=3.0, y=-4.0, z=5.0, radius=1.5),
        make_atom(serial=2, res_seq=2, x=-6.0, y=7.0, z=-8.0, radius=1.2),
    ]
    fixed = "".join(io.print_biomolecule_atoms(atoms, chainflag=False))
    free = "".join(io.print_biomolecule_atoms_free(atoms, keep_chain=False))
    size_fixed = psize.Psize()
    size_fixed.parse_string(fixed)
    size_free = psize.Psize()
    size_free.parse_string(free)
    assert size_fixed.maxlen == pytest.approx(size_free.maxlen)
    assert size_fixed.minlen == pytest.approx(size_free.minlen)


@pytest.mark.parametrize(
    "x, y, z",
    [
        (-22.340, 15.180, -188.097),  # value <= -100 fuses: "15.180-188.097"
        (15.180, 1188.097, 2.104),  # value >= 1000 fuses: "15.1801188.097"
    ],
)
def test_psize_parses_colliding_fixed_columns(x, y, z):
    # Adjacent 8-char fixed-column coordinate fields collide with no separator;
    # psize must read them by column, not by whitespace (regression for 1NH9).
    atom = make_atom(x=x, y=y, z=z, radius=1.5)
    text = "".join(io.print_biomolecule_atoms([atom], chainflag=False))
    # Confirm the collision is real: a clean ATOM PQR line has 10 whitespace
    # tokens (ATOM serial name res resseq x y z charge radius); a fused
    # coordinate pair drops that to 9.
    atom_line = next(ln for ln in text.splitlines() if ln.startswith("ATOM"))
    assert len(atom_line.split()) == 9
    size = psize.Psize()
    size.parse_string(text)
    assert size.gotatom == 1
    assert size.maxlen[0] == pytest.approx(x + 1.5, abs=1e-3)
    assert size.minlen[1] == pytest.approx(y - 1.5, abs=1e-3)
    assert size.maxlen[2] == pytest.approx(z + 1.5, abs=1e-3)


# ---------------------------------------------------------------------------
# Overflow guards on the titration bridges
# ---------------------------------------------------------------------------


class _BioStub:
    def __init__(self, atoms):
        self.atoms = atoms


@needs_main
def test_pkaani_rejects_oversized_structure():
    bio = _BioStub([make_atom(chain="AA")])
    with pytest.raises(RuntimeError, match="pKa-ANI"):
        main.run_pkaani(make_args(), bio)


@needs_main
def test_propka_requires_cif_capable_when_oversized(monkeypatch):
    # Simulate an older/stock propka without mmCIF support.
    monkeypatch.delattr(main.pk_in, "read_mmcif", raising=False)
    bio = _BioStub([make_atom(chain="AA")])
    args = argparse.Namespace(parameters=None, keep_chain=False)
    with pytest.raises(RuntimeError, match="CIF-capable"):
        main.run_propka(args, bio)
