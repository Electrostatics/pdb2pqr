"""Amino Acid Structures for PDB2PQR

This module contains the base amino acid structures for pdb2pqr.

Author:  Todd Dolinsky
"""
import logging
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
        # TODO - why isn't super().__init__() called?
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
        self.peptide_c = None
        self.peptide_n = None
        self.is_n_term = 0
        self.is_c_term = 0
        self.is5term = 0
        self.is3term = 0
        self.missing = []
        self.reference = ref
        self.fixed = 0
        self.stateboolean = {}

        # Create each atom
        for atom_ in atoms:
            if atom_.name in ref.altnames: # Rename atoms
                atom_.name = ref.altnames[atom_.name]
            if atom_.name not in self.map:
                atom = Atom(atom_, "ATOM", self)
                self.add_atom(atom)
            else:
                _LOGGER.debug("Ignoring atom %s", atom_.name)

    def create_atom(self, atomname, newcoords):
        """Create an atom.  Override the generic residue's version of
        create_atom().
        """
        # TODO - OK to add a default type=ATOM argument like superclass?
        oldatom = self.atoms[0]
        newatom = Atom(oldatom, "ATOM", self)
        newatom.set("x", newcoords[0])
        newatom.set("y", newcoords[1])
        newatom.set("z", newcoords[2])
        newatom.set("name", atomname)
        newatom.set("occupancy", 1.00)
        newatom.set("temp_factor", 0.00)
        newatom.added = 1
        self.add_atom(newatom)

    def add_atom(self, atom):
        """Override the existing add_atom - include the link to the reference
        object
        """
        self.atoms.append(atom)
        atomname = atom.get("name")
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
            _LOGGER.debug("Skipping atom reference for %s", atomname)
            atom.reference = None

    def add_dihedral_angle(self, value):
        """Add the value to the list of chiangles"""
        self.dihedrals.append(value)

    def set_state(self):
        """Set the name to use for the forcefield based on the current state.
        Uses N* and C* for termini.
        """
        if self.is_n_term:
            if "NEUTRAL-NTERM" in self.patches:
                self.ffname = "NEUTRAL-N%s" % self.ffname
            else:
                self.ffname = "N%s" % self.ffname
        elif self.is_c_term:
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

    def set_state(self):
        """
           Set the name to use for the forcefield based on the current
           state.
        """
        if "AR0" in self.patches or self.name == "AR0":
            self.ffname = "AR0"
        Amino.set_state(self)


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

    def set_state(self):
        """Set the name to use for the forcefield based on the current state."""
        if "ASH" in self.patches or self.name == "ASH":
            self.ffname = "ASH"
        Amino.set_state(self)


class CYS(Amino):
    """Cysteine class"""

    def __init__(self, atoms, ref):
        Amino.__init__(self, atoms, ref)
        self.reference = ref
        self.ss_bonded = 0
        self.ss_bonded_partner = None

    def letterCode(self):
        return 'C'

    def set_state(self):
        """Set the state of the CYS object.
        If SS-bonded, use CYX.  If negatively charged, use CYM.  If HG is not
        present, use CYX.
        """
        if "CYX" in self.patches or self.name == "CYX":
            self.ffname = "CYX"
        elif self.ss_bonded:
            self.ffname = "CYX"
        elif "CYM" in self.patches or self.name == "CYM":
            self.ffname = "CYM"
        elif not self.has_atom("HG"):
            self.ffname = "CYX"
        Amino.set_state(self)


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

    def set_state(self):
        """Set the name to use for the forcefield based on the current state."""
        if "GLH" in self.patches or self.name == "GLH":
            self.ffname = "GLH"
        Amino.set_state(self)


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

    def set_state(self):
        """Histidines are a special case due to the presence of several
        different forms.  This function sets all non- positive incarnations
        of HIS to neutral HIS by checking to see if optimization removed
        hacceptor or hdonor flags.  Otherwise HID is used as the default.
        """
        if "HIP" not in self.patches and self.name not in ["HIP", "HSP"]:
            if self.get_atom("ND1").hdonor and not self.get_atom("ND1").hacceptor:
                if self.has_atom("HE2"):
                    self.removeAtom("HE2")
            elif self.get_atom("NE2").hdonor and not self.get_atom("NE2").hacceptor:
                if self.has_atom("HD1"):
                    self.removeAtom("HD1")
            elif self.get_atom("ND1").hacceptor and not self.get_atom("ND1").hdonor:
                if self.has_atom("HD1"):
                    self.removeAtom("HD1")
            else: # Default to HID
                if self.has_atom("HE2"):
                    self.removeAtom("HE2")

        if self.has_atom("HD1") and self.has_atom("HE2"):
            self.ffname = "HIP"
        elif self.has_atom("HD1"):
            self.ffname = "HID"
        elif self.has_atom("HE2"):
            self.ffname = "HIE"
        else:
            errstr = ("Invalid type for %s! Missing both HD1 and HE2 atoms. If "
                      "you receive this error while using the --assign-only "
                      "option you can only resolve it by adding HD1, HE2 or "
                      "both to this residue.") % str(self)
            raise TypeError(errstr)
        Amino.set_state(self)


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

    def set_state(self):
        """Determine if this is LYN or not"""
        if "LYN" in self.patches or self.name == "LYN":
            self.ffname = "LYN"
        Amino.set_state(self)


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

    def set_state(self):
        """Set the name to use for the forcefield based on the current state.
        Uses N* and C* for termini.
        """
        if self.is_n_term:
            self.ffname = "N%s" % self.ffname
        elif self.is_c_term:
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

    def set_state(self):
        """See if the TYR is negative or not"""
        if "TYM" in self.patches or self.name == "TYM":
            self.ffname = "TYM"
        Amino.set_state(self)


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
        sample_atom = atoms[-1]

        self.atoms = []
        self.name = sample_atom.res_name
        self.chain_id = sample_atom.chain_id
        self.res_seq = sample_atom.res_seq
        self.ins_code = sample_atom.ins_code

        self.fixed = 0
        self.ffname = "WAT"
        self.map = {}
        self.reference = ref

        # Create each atom
        for atom_ in atoms:
            if atom_.name in ref.altnames: # Rename atoms
                atom_.name = ref.altnames[atom_.name]

            atom = Atom(atom_, "HETATM", self)
            atomname = atom.get("name")
            if atomname not in self.map:
                self.add_atom(atom)
            else: # Don't add duplicate atom with alt_loc field
                oldatom = self.get_atom(atomname)
                oldatom.set("alt_loc", "")

    def create_atom(self, atomname, newcoords):
        """Create a water atom.  Note the HETATM field."""
        oldatom = self.atoms[0]
        newatom = Atom(oldatom, "HETATM", self)
        newatom.set("x", newcoords[0])
        newatom.set("y", newcoords[1])
        newatom.set("z", newcoords[2])
        newatom.set("name", atomname)
        newatom.set("occupancy", 1.00)
        newatom.set("temp_factor", 0.00)
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
                    if bondatom not in atom.bonds:
                        atom.bonds.append(bondatom)
                    if atom not in bondatom.bonds:
                        bondatom.bonds.append(atom)
        except KeyError:
            _LOGGER.debug("Ignoring reference for WAT atom %s", atomname)
            atom.reference = None


class LIG(Residue):
    """Generic ligand class"""

    def __init__(self, atoms, ref):
        sample_atom = atoms[-1]

        self.atoms = []
        self.name = sample_atom.res_name
        self.chain_id = sample_atom.chain_id
        self.res_seq = sample_atom.res_seq
        self.ins_code = sample_atom.ins_code

        self.fixed = 0
        # TODO - why is the ffname "WAT"?
        self.ffname = "WAT"
        self.map = {}
        self.reference = ref

        self.is_n_term = 0
        self.is_c_term = 0

        # Create each atom
        for atom_ in atoms:
            if atom_.name in ref.altnames: # Rename atoms
                atom_.name = ref.altnames[atom_.name]

            atom = Atom(atom_, "HETATM", self)
            atomname = atom.get("name")
            if atomname not in self.map:
                self.add_atom(atom)
            else: # Don't add duplicate atom with alt_loc field
                oldatom = self.get_atom(atomname)
                oldatom.set("alt_loc", "")

    def create_atom(self, atomname, newcoords):
        oldatom = self.atoms[0]
        newatom = Atom(oldatom, "HETATM", self)
        newatom.set("x", newcoords[0])
        newatom.set("y", newcoords[1])
        newatom.set("z", newcoords[2])
        newatom.set("name", atomname)
        newatom.set("occupancy", 1.00)
        newatom.set("temp_factor", 0.00)
        newatom.added = 1
        self.add_atom(newatom)

    def add_atom(self, atom):
        self.atoms.append(atom)
        atomname = atom.get("name")
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
            _LOGGER.debug("Ignoring atom reference for ligand %s", atomname)
            atom.reference = None
