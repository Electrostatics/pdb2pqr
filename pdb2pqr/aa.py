"""Amino Acid Structures for PDB2PQR

This module contains the base amino acid structures for pdb2pqr.

Author:  Todd Dolinsky
"""
import logging
import string
# TODO - remove import * statement
# from .structures import Residue
from .structures import *


_LOGGER = logging.getLogger(__name__)


class Amino(Residue):
    """Amino class
    
    This class provides standard features of the amino acids
    """
    def __init__(self, atoms, ref):
        """Constructor

        Args:
            atoms:  A list of Atom objects to be stored in this class (list)
            ref:  The reference object for the amino acid.  Used to convert
            from the alternate naming scheme to the main naming scheme.
        """
        sampleAtom = atoms[-1]

        self.atoms = []
        self.name = sampleAtom.resName
        self.chainID = sampleAtom.chainID
        self.resSeq = sampleAtom.resSeq
        self.iCode = sampleAtom.iCode

        self.ffname = self.name
        self.map = {}
        self.dihedrals = []
        self.patches = []
        self.peptideC = None
        self.peptideN = None
        self.isNterm = 0
        self.isCterm = 0
        self.is5term = 0
        self.is3term = 0
        self.missing = []
        self.reference = ref
        self.fixed = 0
        self.stateboolean = {}

        # Create each atom
        for a in atoms:
            if a.name in ref.altnames: # Rename atoms
                a.name = ref.altnames[a.name]
            if a.name not in self.map:
                atom = Atom(a, "ATOM", self)
                self.addAtom(atom)
            else:
                _LOGGER.debug("Ignoring atom %s", a.name)

    def createAtom(self, atomname, newcoords):
        """Create an atom.  Override the generic residue's version of
        createAtom().
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
        self.addAtom(newatom)

    def addAtom(self, atom):
        """Override the existing addAtom - include the link to the reference
        object
        """
        self.atoms.append(atom)
        atomname = atom.get("name")
        self.map[atomname] = atom
        try:
            atom.reference = self.reference.map[atomname]
            for bond in atom.reference.bonds:
                if self.hasAtom(bond):
                    bondatom = self.map[bond]
                    if bondatom not in atom.bonds: atom.bonds.append(bondatom)
                    if atom not in bondatom.bonds: bondatom.bonds.append(atom)
        except KeyError:
            _LOGGER.debug("Skipping atom reference for %s", atomname)
            atom.reference = None

    def addDihedralAngle(self, value):
        """Add the value to the list of chiangles"""
        self.dihedrals.append(value)

    def setState(self):
        """Set the name to use for the forcefield based on the current state.
        Uses N* and C* for termini.
        """
        if self.isNterm:
            if "NEUTRAL-NTERM" in self.patches:
                self.ffname = "NEUTRAL-N%s" % self.ffname
            else:
                self.ffname = "N%s" % self.ffname
        elif self.isCterm:
            if "NEUTRAL-CTERM" in self.patches:
                self.ffname = "NEUTRAL-C%s" % self.ffname
            else:
                self.ffname = "C%s" % self.ffname
        return


class ALA(Amino):
    """Alanine class"""

    def __init__(self, atoms, ref):
        Amino.__init__(self, atoms, ref)
        self.reference = ref

    def letterCode(self):
        return 'A'


class ARG(Amino):
    """Arginine class"""

    def __init__(self, atoms, ref):
        """
            Initialize the class

            Parameters
                atoms:      A list of Atom objects to be stored in this class
                            (list)
        """
        Amino.__init__(self, atoms, ref)
        self.reference = ref

    def letterCode(self):
        return 'R'

    def setState(self):
        """
           Set the name to use for the forcefield based on the current
           state.
        """
        if "AR0" in self.patches or self.name == "AR0":
            self.ffname = "AR0"
        Amino.setState(self)


class ASN(Amino):
    """Asparagine class"""

    def __init__(self, atoms, ref):
        Amino.__init__(self, atoms, ref)
        self.reference = ref

    def letterCode(self):
        return 'N'


class ASP(Amino):
    """Aspartic Acid class"""

    def __init__(self, atoms, ref):
        Amino.__init__(self, atoms, ref)
        self.reference = ref

    def letterCode(self):
        return 'D'

    def setState(self):
        """Set the name to use for the forcefield based on the current state."""
        if "ASH" in self.patches or self.name == "ASH":
            self.ffname = "ASH"
        Amino.setState(self)


class CYS(Amino):
    """Cysteine class"""

    def __init__(self, atoms, ref):
        Amino.__init__(self, atoms, ref)
        self.reference = ref
        self.SSbonded = 0
        self.SSbondedpartner = None

    def letterCode(self):
        return 'C'

    def setState(self):
        """Set the state of the CYS object.
        If SS-bonded, use CYX.  If negatively charged, use CYM.  If HG is not
        present, use CYX.
        """
        if "CYX" in self.patches or self.name == "CYX":
            self.ffname = "CYX"
        elif self.SSbonded:
            self.ffname = "CYX"
        elif "CYM" in self.patches or self.name == "CYM":
            self.ffname = "CYM"
        elif not self.hasAtom("HG"):
            self.ffname = "CYX"
        Amino.setState(self)


class GLN(Amino):
    """Glutamine class"""

    def __init__(self, atoms, ref):
        Amino.__init__(self, atoms, ref)
        self.reference = ref

    def letterCode(self):
        return 'Q'


class GLU(Amino):
    """Glutamic Acid class"""

    def __init__(self, atoms, ref):
        Amino.__init__(self, atoms, ref)
        self.reference = ref

    def letterCode(self):
        return 'E'

    def setState(self):
        """Set the name to use for the forcefield based on the current state."""
        if "GLH" in self.patches or self.name == "GLH":
            self.ffname = "GLH"
        Amino.setState(self)


class GLY(Amino):
    """Glycine class"""

    def __init__(self, atoms, ref):
        Amino.__init__(self, atoms, ref)
        self.reference = ref

    def letterCode(self):
        return 'G'


class HIS(Amino):
    """Histidine class"""

    def __init__(self, atoms, ref):
        Amino.__init__(self, atoms, ref)
        self.reference = ref

    def letterCode(self):
        return 'H'

    def setState(self):
        """Histidines are a special case due to the presence of several
        different forms.  This function sets all non- positive incarnations
        of HIS to neutral HIS by checking to see if optimization removed
        hacceptor or hdonor flags.  Otherwise HID is used as the default.
        """
        if "HIP" not in self.patches and self.name not in ["HIP", "HSP"]:
            if self.getAtom("ND1").hdonor and not self.getAtom("ND1").hacceptor:
                if self.hasAtom("HE2"):
                    self.removeAtom("HE2")
            elif self.getAtom("NE2").hdonor and not self.getAtom("NE2").hacceptor:
                if self.hasAtom("HD1"):
                    self.removeAtom("HD1")
            elif self.getAtom("ND1").hacceptor and not self.getAtom("ND1").hdonor:
                if self.hasAtom("HD1"):
                    self.removeAtom("HD1")
            else: # Default to HID
                if self.hasAtom("HE2"):
                    self.removeAtom("HE2")

        if self.hasAtom("HD1") and self.hasAtom("HE2"):
            self.ffname = "HIP"
        elif self.hasAtom("HD1"):
            self.ffname = "HID"
        elif self.hasAtom("HE2"):
            self.ffname = "HIE"
        else:
            errstr = "Invalid type for %s! Missing both HD1 and HE2 atoms. If you receive this error while using the --assign-only option you can only resolve it by adding HD1, HE2 or both to this residue." % str(self)
            raise TypeError(errstr)
        Amino.setState(self)


class ILE(Amino):
    """Isoleucine class"""

    def __init__(self, atoms, ref):
        Amino.__init__(self, atoms, ref)
        self.reference = ref

    def letterCode(self):
        return 'I'


class LEU(Amino):
    """Leucine class"""

    def __init__(self, atoms, ref):
        Amino.__init__(self, atoms, ref)
        self.reference = ref

    def letterCode(self):
        return 'L'


class LYS(Amino):
    """Lysine class"""

    def __init__(self, atoms, ref):
        Amino.__init__(self, atoms, ref)
        self.reference = ref

    def letterCode(self):
        return 'K'

    def setState(self):
        """Determine if this is LYN or not"""
        if "LYN" in self.patches or self.name == "LYN":
            self.ffname = "LYN"
        Amino.setState(self)


class MET(Amino):
    """Methionine class"""

    def __init__(self, atoms, ref):
        Amino.__init__(self, atoms, ref)
        self.reference = ref

    def letterCode(self):
        return 'M'


class PHE(Amino):
    """Phenylalanine class"""

    def __init__(self, atoms, ref):
        Amino.__init__(self, atoms, ref)
        self.reference = ref

    def letterCode(self):
        return 'F'


class PRO(Amino):
    """Proline class"""

    def __init__(self, atoms, ref):
        Amino.__init__(self, atoms, ref)
        self.reference = ref

    def letterCode(self):
        return 'P'

    def setState(self):
        """Set the name to use for the forcefield based on the current state.
        Uses N* and C* for termini.
        """
        if self.isNterm:
            self.ffname = "N%s" % self.ffname
        elif self.isCterm:
            if "NEUTRAL-CTERM" in self.patches:
                self.ffname = "NEUTRAL-C%s" % self.ffname
            else:
                self.ffname = "C%s" % self.ffname


class SER(Amino):
    """Serine class"""

    def __init__(self, atoms, ref):
        Amino.__init__(self, atoms, ref)
        self.reference = ref

    def letterCode(self):
        return 'S'


class THR(Amino):
    """Threonine class"""

    def __init__(self, atoms, ref):
        Amino.__init__(self, atoms, ref)
        self.reference = ref

    def letterCode(self):
        return 'T'


class TRP(Amino):
    """Tryptophan class """

    def __init__(self, atoms, ref):
        Amino.__init__(self, atoms, ref)
        self.reference = ref

    def letterCode(self):
        return 'W'


class TYR(Amino):
    """Tyrosine class"""

    def __init__(self, atoms, ref):
        Amino.__init__(self, atoms, ref)
        self.reference = ref

    def letterCode(self):
        return 'Y'

    def setState(self):
        """See if the TYR is negative or not"""
        if "TYM" in self.patches or self.name == "TYM":
            self.ffname = "TYM"
        Amino.setState(self)


class VAL(Amino):
    """Valine class"""

    def __init__(self, atoms, ref):
        Amino.__init__(self, atoms, ref)
        self.reference = ref

    def letterCode(self):
        return 'V'


class WAT(Residue):
    """Water class"""
    water_residue_names = ['HOH', 'WAT']

    def __init__(self, atoms, ref):
        sampleAtom = atoms[-1]

        self.atoms = []
        self.name = sampleAtom.resName
        self.chainID = sampleAtom.chainID
        self.resSeq = sampleAtom.resSeq
        self.iCode = sampleAtom.iCode

        self.fixed = 0
        self.ffname = "WAT"
        self.map = {}
        self.reference = ref

        # Create each atom
        for a in atoms:
            if a.name in ref.altnames: # Rename atoms
                a.name = ref.altnames[a.name]

            atom = Atom(a, "HETATM", self)
            atomname = atom.get("name")
            if atomname not in self.map:
                self.addAtom(atom)
            else: # Don't add duplicate atom with altLoc field
                oldatom = self.getAtom(atomname)
                oldatom.set("altLoc","")

    def createAtom(self, atomname, newcoords):
        """Create a water atom.  Note the HETATM field."""
        oldatom = self.atoms[0]
        newatom = Atom(oldatom, "HETATM", self)
        newatom.set("x",newcoords[0])
        newatom.set("y",newcoords[1])
        newatom.set("z",newcoords[2])
        newatom.set("name", atomname)
        newatom.set("occupancy",1.00)
        newatom.set("tempFactor",0.00)
        newatom.added = 1
        self.addAtom(newatom)

    def addAtom(self, atom):
        """Override the existing addAtom - include the link to the reference
        object.
        """
        self.atoms.append(atom)
        atomname = atom.get("name")
        self.map[atomname] = atom
        try:
            atom.reference = self.reference.map[atomname]
            for bond in atom.reference.bonds:
                if self.hasAtom(bond):
                    bondatom = self.map[bond]
                    if bondatom not in atom.bonds: atom.bonds.append(bondatom)
                    if atom not in bondatom.bonds: bondatom.bonds.append(atom)
        except KeyError:
            _LOGGER.debug("Ignoring reference for WAT atom %s", atomname)
            atom.reference = None


class LIG(Residue):
    """Generic ligand class"""

    def __init__(self, atoms, ref):
        sampleAtom = atoms[-1]

        self.atoms = []
        self.name = sampleAtom.resName
        self.chainID = sampleAtom.chainID
        self.resSeq = sampleAtom.resSeq
        self.iCode = sampleAtom.iCode

        self.fixed = 0
        # TODO - why is the ffname "WAT"?
        self.ffname = "WAT"
        self.map = {}
        self.reference = ref

        self.isNterm = 0
        self.isCterm = 0

        # Create each atom
        for a in atoms:
            if a.name in ref.altnames: # Rename atoms
                a.name = ref.altnames[a.name]

            atom = Atom(a, "HETATM", self)
            atomname = atom.get("name")
            if atomname not in self.map:
                self.addAtom(atom)
            else: # Don't add duplicate atom with altLoc field
                oldatom = self.getAtom(atomname)
                oldatom.set("altLoc","")

    def createAtom(self, atomname, newcoords):
        oldatom = self.atoms[0]
        newatom = Atom(oldatom, "HETATM", self)
        newatom.set("x",newcoords[0])
        newatom.set("y",newcoords[1])
        newatom.set("z",newcoords[2])
        newatom.set("name", atomname)
        newatom.set("occupancy",1.00)
        newatom.set("tempFactor",0.00)
        newatom.added = 1
        self.addAtom(newatom)

    def addAtom(self, atom):
        self.atoms.append(atom)
        atomname = atom.get("name")
        self.map[atomname] = atom
        try:
            atom.reference = self.reference.map[atomname]
            for bond in atom.reference.bonds:
                if self.hasAtom(bond):
                    bondatom = self.map[bond]
                    if bondatom not in atom.bonds:
                        atom.bonds.append(bondatom)
                    if atom not in bondatom.bonds:
                        bondatom.bonds.append(atom)
        except KeyError:
            _LOGGER.debug("Ignoring atom reference for ligand %s", atomname)
            atom.reference = None
