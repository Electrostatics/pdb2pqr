"""Nucleic Acid Structures for PDB2PQR

This module contains the base nucleic acid structures for pdb2pqr.

Author:  Todd Dolinsky
"""
from .structures import Residue, Atom


class Nucleic(Residue):
    """This class provides standard features of the nucleic acids listed
    below.
    """
    def __init__(self, atoms, ref):
        sample_atom = atoms[-1]

        self.atoms = []
        self.name = sample_atom.res_name
        self.chain_id = sample_atom.chain_id
        self.res_seq = sample_atom.res_seq
        self.ins_code = sample_atom.ins_code

        self.ffname = self.name
        self.map = {}
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
            if atom.name in ref.altnames: # Rename atoms
                atom.name = ref.altnames[atom.name]

            if atom.name not in self.map:
                atom_ = Atom(atom, "ATOM", self)
                self.add_atom(atom_)

    # TODO - here's some more code that's duplicated all over the place.
    def create_atom(self, atomname, newcoords):
        """Create an atom.  Overrides the generic residue's create_atom().

        Args:
            atomname:  The name of the atom to add (string)
            newcoords: The coordinates of the atom (list)
        """
        oldatom = self.atoms[0]
        newatom = Atom(oldatom, "ATOM", self)
        newatom.x = newcoords[0]
        newatom.y = newcoords[1]
        newatom.z = newcoords[2]
        newatom.name = atomname
        newatom.occupancy = 1.0
        newatom.temp_factor = 0.0
        newatom.added = 1
        self.add_atom(newatom)

    def add_atom(self, atom):
        """Override the existing add_atom - include the link to the reference
        object.
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
        """Add the value to the list of chiangles."""
        self.dihedrals.append(value)

    def set_state(self):
        """Adds the termini for all inherited objects."""
        if self.is5term:
            self.ffname = self.ffname + "5"
        if self.is3term:
            self.ffname = self.ffname + "3"


class ADE(Nucleic):
    """Adenosine class"""

    def __init__(self, atoms, ref):
        Nucleic.__init__(self, atoms, ref)
        self.reference = ref

    def letter_code(self):
        return 'A'

    def set_state(self):
        if self.has_atom("O2'"):
            self.ffname = "RA"
        else:
            self.ffname = "DA"
        Nucleic.set_state(self)


class CYT(Nucleic):
    """Cytidine class"""

    def __init__(self, atoms, ref):
        Nucleic.__init__(self, atoms, ref)
        self.reference = ref

    def letter_code(self):
        return 'C'

    def set_state(self):
        if self.has_atom("O2'"):
            self.ffname = "RC"
        else:
            self.ffname = "DC"
        Nucleic.set_state(self)


class GUA(Nucleic):
    """Guanosine class"""

    def __init__(self, atoms, ref):
        Nucleic.__init__(self, atoms, ref)
        self.reference = ref

    def letter_code(self):
        return 'G'

    def set_state(self):
        if self.has_atom("O2'"):
            self.ffname = "RG"
        else:
            self.ffname = "DG"
        Nucleic.set_state(self)


class THY(Nucleic):
    """Thymine class"""

    def __init__(self, atoms, ref):
        Nucleic.__init__(self, atoms, ref)
        self.reference = ref

    def letter_code(self):
        return 'T'

    def set_state(self):
        self.ffname = "DT"
        Nucleic.set_state(self)


class URA(Nucleic):
    """Uridine class"""

    def __init__(self, atoms, ref):
        Nucleic.__init__(self, atoms, ref)
        self.reference = ref

    def letter_code(self):
        return 'U'

    def set_state(self):
        self.ffname = "RU"
        Nucleic.set_state(self)


class RA(ADE):
    """Ribo-ADE"""
    pass


class RC(CYT):
    """Ribo-CYT"""
    pass


class RG(GUA):
    """Ribo-GUA"""
    pass


class DT(THY):
    """Deoxyribo-THY"""
    pass


class RU(URA):
    """Ribo-URA"""
    pass
