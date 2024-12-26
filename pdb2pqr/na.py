"""Nucleic Acid Structures for PDB2PQR

This module contains the base nucleic acid structures for pdb2pqr.

.. codeauthor::  Todd Dolinsky
.. codeauthor::  Nathan Baker
"""

from . import residue
from . import structures as struct
from .definitions import DefinitionResidue
from .pdb import ATOM, HETATM


class Nucleic(residue.Residue):
    """This class provides standard features of the nucleic acids listed below."""

    def __init__(self, atoms: list[ATOM | HETATM], ref: DefinitionResidue):
        sample_atom = atoms[-1]

        self.atoms: list[struct.Atom] = []
        self.name = sample_atom.res_name
        self.chain_id = sample_atom.chain_id
        self.res_seq = sample_atom.res_seq
        self.ins_code = sample_atom.ins_code

        self.ffname = self.name
        self.map: dict[str, struct.Atom] = {}
        self.dihedrals = []
        self.patches = []
        self.is3term = 0
        self.is5term = 0
        self.is_c_term = 0
        self.is_n_term = 0
        self.missing = []
        self.reference = ref

        # Create each atom
        for atom in atoms:
            if atom.name in ref.altnames:  # Rename atoms
                atom.name = ref.altnames[atom.name]

            if atom.name not in self.map:
                atom_ = struct.Atom(atom, "ATOM", self)
                self.add_atom(atom_)

    def create_atom(
        self,
        atomname: str,
        newcoords: list[float] | tuple[float, float, float],
    ):
        """Create an atom.

        Overrides the generic residue's create_atom().

        .. todo:: This code is duplicated in several places.

        :param atomname:  the name of the atom to add
        :type atomname:  str
        :param newcoords:  the coordinates of the atom
        :type newcoords:  [(float, float, float)]
        """
        oldatom = self.atoms[0]
        newatom = struct.Atom(oldatom, "ATOM", self)
        newatom.x = newcoords[0]
        newatom.y = newcoords[1]
        newatom.z = newcoords[2]
        newatom.name = atomname
        newatom.element = atomname[0]
        newatom.occupancy = 1.0
        newatom.temp_factor = 0.0
        newatom.added = 1
        self.add_atom(newatom)

    def add_atom(self, atom: struct.Atom):
        """Add existing atom to system.

        Override the existing add_atom - include the link to the reference
        object.

        :param atom:  atom to add to system.
        :type atom:  Atom
        """
        self.atoms.append(atom)
        atomname = atom.name
        self.map[atomname] = atom
        try:
            atom.reference = self.reference.map[atomname]
            for bond in atom.reference.bonds:
                if self.has_atom(bond):
                    bondatom = self.map[bond]
                    if bondatom not in atom.bonds:
                        atom.bonds.append(bondatom)
                    if atom not in bondatom.bonds:
                        bondatom.bonds.append(atom)
        except KeyError:
            atom.reference = None

    def add_dihedral_angle(self, value):
        """Add the value to the list of chi angles.

        :param value:  dihedral angle to add to list (in degrees)
        :type value:  float
        """
        self.dihedrals.append(value)

    def set_state(self):
        """Adds the termini for all inherited objects."""
        if self.is5term:
            self.ffname = self.ffname + "5"
        if self.is3term:
            self.ffname = self.ffname + "3"


class ADE(Nucleic):
    """Adenosine class."""

    def __init__(self, atoms, ref: DefinitionResidue):
        """Initialize residue.

        :param atoms:  add atoms to residue
        :type atoms:  [Atom]
        """
        Nucleic.__init__(self, atoms, ref)
        self.reference = ref

    def letter_code(self):
        """Return letter code for nucleic acid.

        :return:  letter code for nucleic acid
        :rtype:  str
        """
        return "A"

    def set_state(self):
        """Set ribo- vs. deoxyribo- state of this residue."""
        self.ffname = "RA" if self.has_atom("O2'") else "DA"
        Nucleic.set_state(self)


DA = ADE
RA = ADE
A = ADE


class CYT(Nucleic):
    """Cytidine class"""

    def __init__(self, atoms, ref: DefinitionResidue):
        """Initialize residue.

        :param atoms:  add atoms to residue
        :type atoms:  [Atom]
        """
        Nucleic.__init__(self, atoms, ref)
        self.reference = ref

    def letter_code(self):
        """Return letter code for nucleic acid.

        :return:  letter code for nucleic acid
        :rtype:  str
        """
        return "C"

    def set_state(self):
        """Set ribo- vs. deoxyribo- state of this residue."""
        self.ffname = "RC" if self.has_atom("O2'") else "DC"
        Nucleic.set_state(self)


DC = CYT
RC = CYT
C = CYT


class GUA(Nucleic):
    """Guanosine class"""

    def __init__(self, atoms, ref: DefinitionResidue):
        """Initialize residue.

        :param atoms:  add atoms to residue
        :type atoms:  [Atom]
        """
        Nucleic.__init__(self, atoms, ref)
        self.reference = ref

    def letter_code(self):
        """Return letter code for nucleic acid.

        :return:  letter code for nucleic acid
        :rtype:  str
        """
        return "G"

    def set_state(self):
        """Set ribo- vs. deoxyribo- state of this residue."""
        self.ffname = "RG" if self.has_atom("O2'") else "DG"
        Nucleic.set_state(self)


DG = GUA
RG = GUA
G = GUA


class THY(Nucleic):
    """Thymine class"""

    def __init__(self, atoms, ref: DefinitionResidue):
        """Initialize residue.

        :param atoms:  add atoms to residue
        :type atoms:  [Atom]
        """
        Nucleic.__init__(self, atoms, ref)
        self.reference = ref

    def letter_code(self):
        """Return letter code for nucleic acid.

        :return:  letter code for nucleic acid
        :rtype:  str
        """
        return "T"

    def set_state(self):
        """Set ribo- vs. deoxyribo- state of this residue."""
        self.ffname = "DT"
        Nucleic.set_state(self)


DT = THY


class URA(Nucleic):
    """Uridine class"""

    def __init__(self, atoms, ref: DefinitionResidue):
        """Initialize residue.

        :param atoms:  add atoms to residue
        :type atoms:  [Atom]
        """
        Nucleic.__init__(self, atoms, ref)
        self.reference = ref

    def letter_code(self):
        """Return letter code for nucleic acid.

        :return:  letter code for nucleic acid
        :rtype:  str
        """
        return "U"

    def set_state(self):
        """Set ribo- vs. deoxyribo- state of this residue."""
        self.ffname = "RU"
        Nucleic.set_state(self)


RU = URA
U = URA
