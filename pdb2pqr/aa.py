"""Amino Acid Structures for PDB2PQR

This module contains the base amino acid structures for pdb2pqr.

.. codeauthor:: Todd Dolinsky
.. codeauthor:: Nathan Baker
"""

import logging

from . import quatfit as quat
from . import residue
from . import structures as struct
from . import utilities as util
from .definitions import DefinitionResidue
from .pdb import ATOM, HETATM

_LOGGER = logging.getLogger(__name__)


class Amino(residue.Residue):
    """Amino acid class

    This class provides standard features of the amino acids.
    """

    def __init__(self, atoms: list[ATOM | HETATM], ref: DefinitionResidue):
        """Initialize object.

        .. todo:: need to see whether :func:`super().__init__()` should be
            called

        :param atoms:  A list of :class:`Atom` objects to be stored in this
            object
        :type atoms:  [Atom]
        :param ref:  The reference object for the amino acid.
            Used to convert from the alternate naming scheme to the main naming
            scheme.
        :type ref:  DefinitionResidue
        """
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
            if atom_.name in ref.altnames:  # Rename atoms
                atom_.name = ref.altnames[atom_.name]
            if atom_.name not in self.map:
                atom = struct.Atom(atom_, "ATOM", self)
                self.add_atom(atom)
            else:
                _LOGGER.debug(f"Ignoring atom {atom_.name}")

    def create_atom(
        self,
        atomname: str,
        newcoords: list[float] | tuple[float, float, float],
    ):
        """Create an atom.

        .. todo:: Determine why this is different than superclass method.

        Override the generic residue's version of create_atom().

        :param atomname:  name of atom
        :type atomname:  str
        :param newcoords:  new coordinates for atom
        :type newcoords:  [float, float, float]
        """
        oldatom = self.atoms[0]
        newatom = struct.Atom(oldatom, "ATOM", self)
        newatom.x = newcoords[0]
        newatom.y = newcoords[1]
        newatom.z = newcoords[2]
        newatom.element = atomname[0]
        newatom.name = atomname
        newatom.occupancy = 1.00
        newatom.temp_factor = 0.00
        newatom.added = 1
        self.add_atom(newatom)

    def add_atom(self, atom: struct.Atom):
        """Add atom to residue.

        Override the existing add_atom; include the link to the reference
        object.

        :param atom:  atom to add
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
            _LOGGER.debug(f"Skipping atom reference for {atomname}")
            atom.reference = None

    def add_dihedral_angle(self, value):
        """Add a dihedral angle to the residue list.

        :param value:  dihedral angle (in degrees) to add
        :type value:  float
        """
        self.dihedrals.append(value)

    def set_state(self):
        """Set the name to use for the forcefield based on the current state.

        Uses ``N*`` and ``C*`` for termini.
        """
        if self.is_n_term:
            if "NEUTRAL-NTERM" in self.patches:
                self.ffname = f"NEUTRAL-N{self.ffname}"
            else:
                self.ffname = f"N{self.ffname}"
        elif self.is_c_term:
            if "NEUTRAL-CTERM" in self.patches:
                self.ffname = f"NEUTRAL-C{self.ffname}"
            else:
                self.ffname = f"C{self.ffname}"

    def rebuild_tetrahedral(self, atomname):
        """Rebuild a tetrahedral hydrogen group.

        This is necessary due to the shortcomings of the quatfit routine -
        given a tetrahedral geometry and two existing hydrogens, the quatfit
        routines have two potential solutions.
        This function uses basic tetrahedral geometry to fix this issue.

        :param atomname:  the atom name to add
        :type atomname:  str
        :return:  indication of whether this was successful
        :rtype:  bool
        """
        hcount = 0
        nextatomname = None
        atomref = self.reference.map.get(atomname)
        if atomref is None:
            return False
        bondname = atomref.bonds[0]
        # Return if the bonded atom does not exist
        if not self.has_atom(bondname):
            return False
        # This group is tetrahedral if bondatom has 4 bonds,
        #  3 of which are hydrogens
        for bond in self.reference.map[bondname].bonds:
            if bond.startswith("H"):
                hcount += 1
            elif bond not in ("C-1", "N+1"):
                nextatomname = bond
        # Check if this is a tetrahedral group
        if hcount != 3 or nextatomname is None:
            return False
        # Now rebuild according to the tetrahedral geometry
        bondatom = self.get_atom(bondname)
        nextatom = self.get_atom(nextatomname)
        numbonds = len(bondatom.bonds)
        if numbonds == 1:
            # Place according to two atoms
            coords = [bondatom.coords, nextatom.coords]
            refcoords = [
                self.reference.map[bondname].coords,
                self.reference.map[nextatomname].coords,
            ]
            refatomcoords = atomref.coords
            newcoords = quat.find_coordinates(
                2, coords, refcoords, refatomcoords
            )
            self.create_atom(atomname, newcoords)
            # For LEU and ILE residues only: make sure the Hydrogens are in
            # staggered conformation instead of eclipsed.
            if isinstance(self, LEU):
                hcoords = newcoords
                cbatom = self.get_atom("CB")
                ang = util.dihedral(
                    cbatom.coords, nextatom.coords, bondatom.coords, hcoords
                )
                diffangle = 60 - ang
                self.rotate_tetrahedral(nextatom, bondatom, diffangle)
            elif isinstance(self, ILE):
                hcoords = newcoords
                cg1atom = self.get_atom("CG1")
                cbatom = self.get_atom("CB")
                if bondatom.name == "CD1" or bondatom.name != "CG2":
                    ang = util.dihedral(
                        cbatom.coords,
                        nextatom.coords,
                        bondatom.coords,
                        hcoords,
                    )
                else:
                    ang = util.dihedral(
                        cg1atom.coords,
                        nextatom.coords,
                        bondatom.coords,
                        hcoords,
                    )
                diffangle = 60 - ang
                self.rotate_tetrahedral(nextatom, bondatom, diffangle)
            return True
        elif numbonds == 2:
            # Get the single hydrogen coordinates
            hatom = None
            for bond in bondatom.reference.bonds:
                if self.has_atom(bond) and bond.startswith("H"):
                    hatom = self.get_atom(bond)
                    break
            # Use the existing hydrogen and rotate about the bond
            self.rotate_tetrahedral(nextatom, bondatom, 120)
            newcoords = hatom.coords
            self.rotate_tetrahedral(nextatom, bondatom, -120)
            self.create_atom(atomname, newcoords)
            return True
        elif numbonds == 3:
            # Find the one spot the atom can be
            hatoms = [
                self.get_atom(bond)
                for bond in bondatom.reference.bonds
                if self.has_atom(bond) and bond.startswith("H")
            ]

            # If this is more than two something is wrong
            if len(hatoms) != 2:
                return False
            # Use the existing hydrogen and rotate about the bond
            self.rotate_tetrahedral(nextatom, bondatom, 120)
            newcoords1 = hatoms[0].coords
            self.rotate_tetrahedral(nextatom, bondatom, 120)
            newcoords2 = hatoms[0].coords
            self.rotate_tetrahedral(nextatom, bondatom, 120)
            # Determine which one hatoms[1] is not in
            if util.distance(hatoms[1].coords, newcoords1) > 0.1:
                self.create_atom(atomname, newcoords1)
            else:
                self.create_atom(atomname, newcoords2)
            return True
        return False


class ALA(Amino):
    """Alanine class."""

    def __init__(self, atoms, ref: DefinitionResidue):
        """Initialize object.

        :param atoms:  A list of :class:`Atom` objects to be stored in this
            object
        :type atoms:  [Atom]
        :param ref:  The reference object for the amino acid.
            Used to convert from the alternate naming scheme to the main naming
            scheme.
        :type ref:  DefinitionResidue
        """
        Amino.__init__(self, atoms, ref)
        self.reference = ref

    def letter_code(self) -> str:
        """Return letter code for amino acid.

        :return:  amino acid 1-letter code
        :rtype:  str
        """
        return "A"


class ARG(Amino):
    """Arginine class."""

    def __init__(self, atoms, ref: DefinitionResidue):
        """Initialize object.

        :param atoms:  A list of :class:`Atom` objects to be stored in this
            object
        :type atoms:  [Atom]
        :param ref:  The reference object for the amino acid.
            Used to convert from the alternate naming scheme to the main naming
            scheme.
        :type ref:  DefinitionResidue
        """
        Amino.__init__(self, atoms, ref)
        self.reference = ref

    def letter_code(self):
        """Return letter code for amino acid.

        :return:  amino acid 1-letter code
        :rtype:  str
        """
        return "R"

    def set_state(self):
        """Set forcefield name based on current titration state."""
        if "AR0" in self.patches or self.name == "AR0":
            self.ffname = "AR0"
        Amino.set_state(self)


class ASN(Amino):
    """Asparagine class."""

    def __init__(self, atoms, ref: DefinitionResidue):
        """Initialize object.

        :param atoms:  A list of :class:`Atom` objects to be stored in this
            object
        :type atoms:  [Atom]
        :param ref:  The reference object for the amino acid.
            Used to convert from the alternate naming scheme to the main naming
            scheme.
        :type ref:  DefinitionResidue
        """
        Amino.__init__(self, atoms, ref)
        self.reference = ref

    def letter_code(self):
        """Return letter code for amino acid.

        :return:  amino acid 1-letter code
        :rtype:  str
        """
        return "N"


class ASP(Amino):
    """Aspartic acid class."""

    def __init__(self, atoms, ref: DefinitionResidue):
        """Initialize object.

        :param atoms:  A list of :class:`Atom` objects to be stored in this
            object
        :type atoms:  [Atom]
        :param ref:  The reference object for the amino acid.
            Used to convert from the alternate naming scheme to the main naming
            scheme.
        :type ref:  DefinitionResidue
        """
        Amino.__init__(self, atoms, ref)
        self.reference = ref

    def letter_code(self):
        """Return letter code for amino acid.

        :return:  amino acid 1-letter code
        :rtype:  str
        """
        return "D"

    def set_state(self):
        """Set forcefield name based on current titration state."""
        if "ASH" in self.patches or self.name == "ASH":
            self.ffname = "ASH"
        Amino.set_state(self)


class CYS(Amino):
    """Cysteine class."""

    def __init__(self, atoms, ref: DefinitionResidue):
        """Initialize object.

        :param atoms:  A list of :class:`Atom` objects to be stored in this
            object
        :type atoms:  [Atom]
        :param ref:  The reference object for the amino acid.
            Used to convert from the alternate naming scheme to the main naming
            scheme.
        :type ref:  DefinitionResidue
        """
        Amino.__init__(self, atoms, ref)
        self.reference = ref
        self.ss_bonded = 0
        self.ss_bonded_partner = None

    def letter_code(self):
        """Return letter code for amino acid.

        :return:  amino acid 1-letter code
        :rtype:  str
        """
        return "C"

    def set_state(self):
        """Set forcefield name based on current state.

        If SS-bonded, use CYX.  If negatively charged, use CYM.  If HG is not
        present, use CYX.
        """
        if "CYX" in self.patches or self.name == "CYX" or self.ss_bonded:
            self.ffname = "CYX"
        elif "CYM" in self.patches or self.name == "CYM":
            self.ffname = "CYM"
        elif not self.has_atom("HG"):
            self.ffname = "CYX"
        Amino.set_state(self)


class GLN(Amino):
    """Glutamine class."""

    def __init__(self, atoms, ref: DefinitionResidue):
        """Initialize object.

        :param atoms:  A list of :class:`Atom` objects to be stored in this
            object
        :type atoms:  [Atom]
        :param ref:  The reference object for the amino acid.
            Used to convert from the alternate naming scheme to the main naming
            scheme.
        :type ref:  DefinitionResidue
        """
        Amino.__init__(self, atoms, ref)
        self.reference = ref

    def letter_code(self):
        """Return letter code for amino acid.

        :return:  amino acid 1-letter code
        :rtype:  str
        """
        return "Q"


class GLU(Amino):
    """Glutamic acid class."""

    def __init__(self, atoms, ref: DefinitionResidue):
        """Initialize object.

        :param atoms:  A list of :class:`Atom` objects to be stored in this
            object
        :type atoms:  [Atom]
        :param ref:  The reference object for the amino acid.
            Used to convert from the alternate naming scheme to the main naming
            scheme.
        :type ref:  DefinitionResidue
        """
        Amino.__init__(self, atoms, ref)
        self.reference = ref

    def letter_code(self):
        """Return letter code for amino acid.

        :return:  amino acid 1-letter code
        :rtype:  str
        """
        return "E"

    def set_state(self):
        """Set forcefield name based on current titration state."""
        if "GLH" in self.patches or self.name == "GLH":
            self.ffname = "GLH"
        Amino.set_state(self)


class GLY(Amino):
    """Glycine class."""

    def __init__(self, atoms, ref: DefinitionResidue):
        """Initialize object.

        :param atoms:  A list of :class:`Atom` objects to be stored in this
            object
        :type atoms:  [Atom]
        :param ref:  The reference object for the amino acid.
            Used to convert from the alternate naming scheme to the main naming
            scheme.
        :type ref:  DefinitionResidue
        """
        Amino.__init__(self, atoms, ref)
        self.reference = ref

    def letter_code(self):
        """Return letter code for amino acid.

        :return:  amino acid 1-letter code
        :rtype:  str
        """
        return "G"


class HIS(Amino):
    """Histidine class."""

    def __init__(self, atoms, ref: DefinitionResidue):
        """Initialize object.

        :param atoms:  A list of :class:`Atom` objects to be stored in this
            object
        :type atoms:  [Atom]
        :param ref:  The reference object for the amino acid.
            Used to convert from the alternate naming scheme to the main naming
            scheme.
        :type ref:  DefinitionResidue
        """
        Amino.__init__(self, atoms, ref)
        self.reference = ref

    def letter_code(self):
        """Return letter code for amino acid.

        :return:  amino acid 1-letter code
        :rtype:  str
        """
        return "H"

    def set_state(self):
        """Set forcefield name based on current titration state.

        Histidines are a special case due to the presence of several different
        forms.
        This function sets all neutral forms of HIS to neutral HIS by checking
        to see if optimization removed :makevar:`hacceptor` or
        :makevar:`hdonor` flags.
        Otherwise HID is used as the default.
        """
        errstr = ""
        if "HIP" not in self.patches and self.name not in ["HIP", "HSP"]:
            if (
                self.get_atom("ND1").hdonor
                and not self.get_atom("ND1").hacceptor
            ):
                if self.has_atom("HE2"):
                    self.remove_atom("HE2")
            elif (
                self.get_atom("NE2").hdonor
                and not self.get_atom("NE2").hacceptor
            ) or (
                self.get_atom("ND1").hacceptor
                and not self.get_atom("ND1").hdonor
            ):
                if self.has_atom("HD1"):
                    self.remove_atom("HD1")
            else:  # Default to HID
                if self.has_atom("HE2"):
                    self.remove_atom("HE2")
        if self.has_atom("HD1") and self.has_atom("HE2"):
            self.ffname = "HIP"
        elif self.has_atom("HD1"):
            self.ffname = "HID"
        elif self.has_atom("HE2"):
            self.ffname = "HIE"
        else:
            errstr = (
                f"Invalid type for {self!s}! Missing both HD1 and HE2 "
                "atoms. If you receive this error while using the "
                "--assign-only option you can only resolve it by adding HD1, "
                "HE2 or both to this residue."
            )
            raise TypeError(errstr)
        Amino.set_state(self)


class ILE(Amino):
    """Isoleucine class."""

    def __init__(self, atoms, ref: DefinitionResidue):
        """Initialize object.

        :param atoms:  A list of :class:`Atom` objects to be stored in this
            object
        :type atoms:  [Atom]
        :param ref:  The reference object for the amino acid.
            Used to convert from the alternate naming scheme to the main naming
            scheme.
        :type ref:  DefinitionResidue
        """
        Amino.__init__(self, atoms, ref)
        self.reference = ref

    def letter_code(self):
        """Return letter code for amino acid.

        :return:  amino acid 1-letter code
        :rtype:  str
        """
        return "I"


class LEU(Amino):
    """Leucine class."""

    def __init__(self, atoms, ref: DefinitionResidue):
        """Initialize object.

        :param atoms:  A list of :class:`Atom` objects to be stored in this
            object
        :type atoms:  [Atom]
        :param ref:  The reference object for the amino acid.
            Used to convert from the alternate naming scheme to the main naming
            scheme.
        :type ref:  DefinitionResidue
        """
        Amino.__init__(self, atoms, ref)
        self.reference = ref

    def letter_code(self):
        """Return letter code for amino acid.

        :return:  amino acid 1-letter code
        :rtype:  str
        """
        return "L"


class LYS(Amino):
    """Lysine class."""

    def __init__(self, atoms, ref: DefinitionResidue):
        """Initialize object.

        :param atoms:  A list of :class:`Atom` objects to be stored in this
            object
        :type atoms:  [Atom]
        :param ref:  The reference object for the amino acid.
            Used to convert from the alternate naming scheme to the main naming
            scheme.
        :type ref:  DefinitionResidue
        """
        Amino.__init__(self, atoms, ref)
        self.reference = ref

    def letter_code(self):
        """Return letter code for amino acid.

        :return:  amino acid 1-letter code
        :rtype:  str
        """
        return "K"

    def set_state(self):
        """Set forcefield name based on current titration state."""
        if "LYN" in self.patches or self.name == "LYN":
            self.ffname = "LYN"
        Amino.set_state(self)


class MET(Amino):
    """Methionine class."""

    def __init__(self, atoms, ref: DefinitionResidue):
        """Initialize object.

        :param atoms:  A list of :class:`Atom` objects to be stored in this
            object
        :type atoms:  [Atom]
        :param ref:  The reference object for the amino acid.
            Used to convert from the alternate naming scheme to the main naming
            scheme.
        :type ref:  DefinitionResidue
        """
        Amino.__init__(self, atoms, ref)
        self.reference = ref

    def letter_code(self):
        """Return letter code for amino acid.

        :return:  amino acid 1-letter code
        :rtype:  str
        """
        return "M"


class PHE(Amino):
    """Phenylalanine class."""

    def __init__(self, atoms, ref: DefinitionResidue):
        """Initialize object.

        :param atoms:  A list of :class:`Atom` objects to be stored in this
            object
        :type atoms:  [Atom]
        :param ref:  The reference object for the amino acid.
            Used to convert from the alternate naming scheme to the main naming
            scheme.
        :type ref:  DefinitionResidue
        """
        Amino.__init__(self, atoms, ref)
        self.reference = ref

    def letter_code(self):
        """Return letter code for amino acid.

        :return:  amino acid 1-letter code
        :rtype:  str
        """
        return "F"


class PRO(Amino):
    """Proline class."""

    def __init__(self, atoms, ref: DefinitionResidue):
        """Initialize object.

        :param atoms:  A list of :class:`Atom` objects to be stored in this
            object
        :type atoms:  [Atom]
        :param ref:  The reference object for the amino acid.
            Used to convert from the alternate naming scheme to the main naming
            scheme.
        :type ref:  DefinitionResidue
        """
        Amino.__init__(self, atoms, ref)
        self.reference = ref

    def letter_code(self):
        """Return letter code for amino acid.

        :return:  amino acid 1-letter code
        :rtype:  str
        """
        return "P"

    def set_state(self):
        """Set forcefield name based on the current state.

        Uses ``N*`` and ``C*`` for termini.
        """
        if self.is_n_term:
            self.ffname = f"N{self.ffname}"
        elif self.is_c_term:
            if "NEUTRAL-CTERM" in self.patches:
                self.ffname = f"NEUTRAL-C{self.ffname}"
            else:
                self.ffname = f"C{self.ffname}"


class SER(Amino):
    """Serine class."""

    def __init__(self, atoms, ref: DefinitionResidue):
        """Initialize object.

        :param atoms:  A list of :class:`Atom` objects to be stored in this
            object
        :type atoms:  [Atom]
        :param ref:  The reference object for the amino acid.
            Used to convert from the alternate naming scheme to the main naming
            scheme.
        :type ref:  DefinitionResidue
        """
        Amino.__init__(self, atoms, ref)
        self.reference = ref

    def letter_code(self):
        """Return letter code for amino acid.

        :return:  amino acid 1-letter code
        :rtype:  str
        """
        return "S"


class THR(Amino):
    """Threonine class."""

    def __init__(self, atoms, ref: DefinitionResidue):
        """Initialize object.

        :param atoms:  A list of :class:`Atom` objects to be stored in this
            object
        :type atoms:  [Atom]
        :param ref:  The reference object for the amino acid.
            Used to convert from the alternate naming scheme to the main naming
            scheme.
        :type ref:  DefinitionResidue
        """
        Amino.__init__(self, atoms, ref)
        self.reference = ref

    def letter_code(self):
        """Return letter code for amino acid.

        :return:  amino acid 1-letter code
        :rtype:  str
        """
        return "T"


class TRP(Amino):
    """Tryptophan class."""

    def __init__(self, atoms, ref: DefinitionResidue):
        """Initialize object.

        :param atoms:  A list of :class:`Atom` objects to be stored in this
            object
        :type atoms:  [Atom]
        :param ref:  The reference object for the amino acid.
            Used to convert from the alternate naming scheme to the main naming
            scheme.
        :type ref:  DefinitionResidue
        """
        Amino.__init__(self, atoms, ref)
        self.reference = ref

    def letter_code(self):
        """Return letter code for amino acid.

        :return:  amino acid 1-letter code
        :rtype:  str
        """
        return "W"


class TYR(Amino):
    """Tyrosine class."""

    def __init__(self, atoms, ref: DefinitionResidue):
        """Initialize object.

        :param atoms:  A list of :class:`Atom` objects to be stored in this
            object
        :type atoms:  [Atom]
        :param ref:  The reference object for the amino acid.
            Used to convert from the alternate naming scheme to the main naming
            scheme.
        :type ref:  DefinitionResidue
        """
        Amino.__init__(self, atoms, ref)
        self.reference = ref

    def letter_code(self):
        """Return letter code for amino acid.

        :return:  amino acid 1-letter code
        :rtype:  str
        """
        return "Y"

    def set_state(self):
        """Set forcefield name based on current titration state."""
        if "TYM" in self.patches or self.name == "TYM":
            self.ffname = "TYM"
        Amino.set_state(self)


class VAL(Amino):
    """Valine class."""

    def __init__(self, atoms, ref: DefinitionResidue):
        """Initialize object.

        :param atoms:  A list of :class:`Atom` objects to be stored in this
            object
        :type atoms:  [Atom]
        :param ref:  The reference object for the amino acid.
            Used to convert from the alternate naming scheme to the main naming
            scheme.
        :type ref:  DefinitionResidue
        """
        Amino.__init__(self, atoms, ref)
        self.reference = ref

    def letter_code(self):
        """Return letter code for amino acid.

        :return:  amino acid 1-letter code
        :rtype:  str
        """
        return "V"


class WAT(residue.Residue):
    """Water class.

    .. todo:: Why is water in the amino acid module?

    """

    water_residue_names = ["HOH", "WAT"]

    def __init__(self, atoms, ref: DefinitionResidue):
        """Initialize object.

        :param atoms:  A list of :class:`Atom` objects to be stored in this
            object
        :type atoms:  [Atom]
        :param ref:  The reference object for the residue.
            Used to convert from the alternate naming scheme to the main naming
            scheme.
        :type ref:  DefinitionResidue
        """
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
            if atom_.name in ref.altnames:  # Rename atoms
                atom_.name = ref.altnames[atom_.name]
            atom = struct.Atom(atom_, "HETATM", self)
            atomname = atom.name
            if atomname not in self.map:
                self.add_atom(atom)
            else:  # Don't add duplicate atom with alt_loc field
                oldatom = self.get_atom(atomname)
                oldatom.alt_loc = ""

    def create_atom(self, atomname, newcoords):
        """Create a water atom.

        Note the HETATM field.

        .. todo::  There is a huge amount of duplicated code in this module.

        :param atomname:  name of atom to be added
        :type atomname:  str
        :param newcoords:  coordinates for new atom
        :type newcoords:  [float, float, float]
        """
        oldatom = self.atoms[0]
        newatom = struct.Atom(oldatom, "HETATM", self)
        newatom.x = newcoords[0]
        newatom.y = newcoords[1]
        newatom.z = newcoords[2]
        newatom.name = atomname
        newatom.element = atomname[0]
        newatom.occupancy = 1.00
        newatom.temp_factor = 0.00
        newatom.added = 1
        self.add_atom(newatom)

    def add_atom(self, atom):
        """Add an atom to the residue.

        Override the existing add_atom - include the link to the reference
        object.

        :param atom:  add atom to residue
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
            _LOGGER.debug(f"Ignoring reference for WAT atom {atomname}")
            atom.reference = None


class LIG(residue.Residue):
    """Generic ligand class."""

    def __init__(self, atoms, ref):
        """Initialize this object.

        .. todo:: why is the force field name "WAT" for this?

        :param atoms:  A list of :class:`Atom` objects to be stored in this
            object
        :type atoms:  [Atom]
        :param ref:  The reference object for the residue.
            Used to convert from the alternate naming scheme to the main naming
            scheme.
        :type ref:  DefinitionResidue
        """
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
        self.is_n_term = 0
        self.is_c_term = 0
        # Create each atom
        for atom_ in atoms:
            if atom_.name in ref.altnames:  # Rename atoms
                atom_.name = ref.altnames[atom_.name]
            atom = struct.Atom(atom_, "HETATM", self)
            atomname = atom.name
            if atomname not in self.map:
                self.add_atom(atom)
            else:  # Don't add duplicate atom with alt_loc field
                oldatom = self.get_atom(atomname)
                oldatom.alt_loc = ""

    def create_atom(self, atomname, newcoords):
        """Create a ligand atom.

        :param atomname:  name of atom to be added
        :type atomname:  str
        :param newcoords:  coordinates for new atom
        :type newcoords:  [float, float, float]
        """
        oldatom = self.atoms[0]
        newatom = struct.Atom(oldatom, "HETATM", self)
        newatom.x = newcoords[0]
        newatom.y = newcoords[1]
        newatom.z = newcoords[2]
        newatom.name = atomname
        newatom.element = atomname[0]
        newatom.occupancy = 1.00
        newatom.temp_factor = 0.00
        newatom.added = 1
        self.add_atom(newatom)

    def add_atom(self, atom):
        """Add an atom to the residue.

        Override the existing add_atom - include the link to the reference
        object.

        :param atom:  add atom to residue
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
            _LOGGER.debug(f"Ignoring atom reference for ligand {atomname}")
            atom.reference = None
