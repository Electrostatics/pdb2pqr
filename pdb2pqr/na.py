"""Nucleic Acid Structures for PDB2PQR

This module contains the base nucleic acid structures for pdb2pqr.

Author:  Todd Dolinsky
"""
import string
from .structures import Residue


class Nucleic(Residue):
    """This class provides standard features of the nucleic acids listed
    below.
    """
    def __init__(self, atoms, ref):
        sample_atom = atoms[-1]
        
        self.atoms = []
        self.name = sample_atom.resName
        self.chainID = sample_atom.chainID
        self.resSeq = sample_atom.resSeq
        self.iCode = sample_atom.iCode

        self.ffname = self.name
        self.map = {}
        self.dihedrals = []
        self.patches = []
        self.is3term = 0
        self.is5term = 0
        self.isCterm = 0
        self.isNterm = 0
        self.missing = []
        self.reference = ref
     
        # Create each atom

        for a in atoms:
            if a.name in ref.altnames: # Rename atoms
                a.name = ref.altnames[a.name]

            if a.name not in self.map:
                atom = Atom(a, "ATOM", self)
                self.add_atom(atom)

    def create_atom(self, atomname, newcoords):
        """Create an atom.  Overrides the generic residue's create_atom().

        Args:
            atomname:  The name of the atom to add (string)
            newcoords: The coordinates of the atom (list)
        """
        oldatom = self.atoms[0]
        newatom = Atom(oldatom, "ATOM", self)
        newatom.set("x",newcoords[0])
        newatom.set("y",newcoords[1])
        newatom.set("z",newcoords[2])
        newatom.set("name", atomname)
        newatom.set("occupancy",1.00)
        newatom.set("tempFactor",0.00)
        newatom.added = 1
        self.add_atom(newatom) 

    def add_atom(self, atom):
        """Override the existing add_atom - include the link to the reference
        object.
        """
        self.atoms.append(atom)
        atomname = atom.get("name")
        self.map[atomname] = atom
        try:
            atom.reference = self.reference.map[atomname]
            for bond in atom.reference.bonds:
                if self.has_atom(bond):
                    bondatom = self.map[bond]
                    if bondatom not in atom.bonds: atom.bonds.append(bondatom)
                    if atom not in bondatom.bonds: bondatom.bonds.append(atom)
        except KeyError:
            atom.reference = None

    def add_dihedral_angle(self, value):
        """Add the value to the list of chiangles."""
        self.dihedrals.append(value)

    def set_state(self):
        """Adds the termini for all inherited objects."""
        if self.is5term: self.ffname = self.ffname + "5"
        if self.is3term: self.ffname = self.ffname + "3"
 

class ADE(Nucleic):
    """Adenosine class"""

    def __init__(self, atoms, ref):
        Nucleic.__init__(self, atoms, ref)
        self.reference = ref
        
    def letterCode(self):
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
    
    def letterCode(self):
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
        
    def letterCode(self):
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
        
    def letterCode(self):
        return 'T'
        
    def set_state(self):
        self.ffname = "DT"
        Nucleic.set_state(self)


class URA(Nucleic):
    """Uridine class"""

    def __init__(self, atoms, ref):
        Nucleic.__init__(self, atoms, ref)
        self.reference = ref
        
    def letterCode(self):
        return 'U'
        
    def set_state(self):
        self.ffname = "RU"
        Nucleic.set_state(self)


class RA(ADE):
    pass


class RC(CYT):
    pass


class RG(GUA):
    pass


class DT(THY):
    pass


class RU(URA):
    pass
