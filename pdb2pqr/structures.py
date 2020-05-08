"""Simple biomolecular structures

This module contains the simpler structure objects used in PDB2PQR and their
associated methods.

Author: Todd Dolinsky
"""
from . import pdb
from .config import BACKBONE


class Chain(object):
    """Chain class

    The chain class contains information about each chain within a given
    Protein object.
    """

    def __init__(self, chain_id):
        """Initialize the class

        Args:
            chain_id: The chain_id for this chain as denoted in the PDB file (string)
        """
        self.chain_id = chain_id
        self.residues = []
        self.name = None

    def add_residue(self, residue):
        """Add a residue to the chain

        Args:
            residue: The residue to be added (Residue)
        """
        self.residues.append(residue)

    def renumber_residues(self):
        """Renumber Atoms based on actual Residue number and not PDB res_seq"""
        count = 1
        for residue in self.residues:
            residue.set_res_seq(count)
            count += 1

    @property
    def atoms(self):
        """Return a list of Atom objects contained in this chain

        Returns
            atomlist: List of Atom objects (list)
        """
        atomlist = []
        for residue in self.residues:
            my_list = residue.atoms
            for atom in my_list:
                atomlist.append(atom)
        return atomlist

    def __str__(self):
        output = []
        for residue in self.residues:
            output.append(residue.letter_code())
        return ''.join(output)


class Atom(pdb.ATOM):
    """Class Atom

    The Atom class inherits from the ATOM object in pdb.py.  It is used
    for adding fields not found in the pdb that may be useful for analysis.
    Also simplifies code by combining ATOM and HETATM objects into a
    single class.
    """

    def __init__(self, atom, type_, residue):
        """Initialize the new Atom object by using the old object.

        Args:
            atom:  The original ATOM object (ATOM)
            type:  Either ATOM or HETATM (string)
            residue:  A pointer back to the parent residue object (Residue)
        """
        if type_ == "ATOM" or type_ == "HETATM":
            self.type = type_
        else:
            raise ValueError("Invalid atom type %s (Atom Class IN structures.py)!" % type_)
        self.serial = atom.serial
        self.name = atom.name
        self.alt_loc = atom.alt_loc
        self.res_name = atom.res_name
        self.chain_id = atom.chain_id
        self.res_seq = atom.res_seq
        self.ins_code = atom.ins_code
        self.x = atom.x
        self.y = atom.y
        self.z = atom.z
        self.occupancy = atom.occupancy
        self.temp_factor = atom.temp_factor
        self.seg_id = atom.seg_id
        self.element = atom.element
        self.charge = atom.charge
        self.bonds = []
        self.reference = None
        self.residue = residue
        self.radius = None
        self.ffcharge = None
        self.hdonor = 0
        self.hacceptor = 0
        self.cell = None
        self.added = 0
        self.optimizeable = 0
        self.refdistance = 0
        self.id = None
        try:
            self.mol2charge = atom.mol2charge
        except AttributeError:
            self.mol2charge = None

    def get_common_string_rep(self, chainflag=False):
        """Returns a string of the common column of the new atom type.

        Uses the ATOM string output but changes the first field to either by ATOM
        or HETATM as necessary. This is used to create the output for pqr and pdb files.

        Returns
            outstr: String with ATOM/HETATM field set appropriately
        """
        outstr = ""
        tstr = self.type
        outstr += str.ljust(tstr, 6)[:6]
        tstr = "%d" % self.serial
        outstr += str.rjust(tstr, 5)[:5]
        outstr += " "
        tstr = self.name
        if len(tstr) == 4 or len(tstr.strip("FLIP")) == 4:
            outstr += str.ljust(tstr, 4)[:4]
        else:
            outstr += " " + str.ljust(tstr, 3)[:3]

        tstr = self.res_name
        if len(tstr) == 4:
            outstr += str.ljust(tstr, 4)[:4]
        else:
            outstr += " " + str.ljust(tstr, 3)[:3]

        outstr += " "
        if chainflag:
            tstr = self.chain_id
        else:
            tstr = ''
        outstr += str.ljust(tstr, 1)[:1]
        tstr = "%d" % self.res_seq
        outstr += str.rjust(tstr, 4)[:4]
        if self.ins_code != "":
            outstr += "%s   " % self.ins_code
        else:
            outstr += "    "
        tstr = "%8.3f" % self.x
        outstr += str.ljust(tstr, 8)[:8]
        tstr = "%8.3f" % self.y
        outstr += str.ljust(tstr, 8)[:8]
        tstr = "%8.3f" % self.z
        outstr += str.ljust(tstr, 8)[:8]
        return outstr

    def __str__(self):
        """Returns a string of the new atom type.

        Uses the ATOM string output but changes the first field to either by
        ATOM or HETATM as necessary.
        This is used to create the output for pqr files!

        Returns
            str: String with ATOM/HETATM field set appropriately
        """
        return self.get_pqr_string()

    def get_pqr_string(self, chainflag=False):
        """Returns a string of the new atom type.

        Uses the ATOM string output but changes the first field to either by
        ATOM or HETATM as necessary. This is used to create the output for pqr
        files!

        Returns
            str: String with ATOM/HETATM field set appropriately
        """
        outstr = self.get_common_string_rep(chainflag=chainflag)
        if self.ffcharge != None:
            ffcharge = "%.4f" % self.ffcharge
        else:
            ffcharge = "0.0000"
        outstr += str.rjust(ffcharge, 8)[:8]
        if self.radius != None:
            ffradius = "%.4f" % self.radius
        else:
            ffradius = "0.0000"
        outstr += str.rjust(ffradius, 7)[:7]
        return outstr

    def get_pdb_string(self):
        """Returns a string of the new atom type.

        Uses the ATOM string output but changes the first field to either by
        ATOM or HETATM as necessary. This is for the pdb representation of the
        atom. The propka30 module depends on this being correct.

        Returns
            str: String with ATOM/HETATM field set appropriately
        """
        outstr = self.get_common_string_rep(chainflag=True)

        tstr = "%6.2f" % self.occupancy
        outstr += str.ljust(tstr, 6)[:6]
        tstr = "%6.2f" % self.temp_factor
        outstr += str.rjust(tstr, 6)[:6]
        #padding between temp factor and seg_id
        outstr += ' ' * 7
        tstr = self.seg_id
        outstr += str.ljust(tstr, 4)[:4]
        tstr = self.element
        outstr += str.ljust(tstr, 2)[:2]
        tstr = str(self.charge)
        outstr += str.ljust(tstr, 2)[:2]
        return outstr

    @property
    def coords(self):
        """Return the x,y,z coordinates of the atom in list form

        TODO - this should be converted to numpy

        Returns
            List of the coordinates (list)
        """
        return [self.x, self.y, self.z]

    def add_bond(self, bondedatom):
        """Add a bond to the list of bonds

        Args:
            bondedatom: The atom to bond to (Atom)
        """
        self.bonds.append(bondedatom)

    @property
    def is_hydrogen(self):
        """Is this atom a Hydrogen atom?"""
        return self.name[0] == "H"

    @property
    def is_backbone(self):
        """Return true if atom name is in backbone, otherwise false"""
        return self.name in BACKBONE

    @property
    def has_reference(self):
        """Determine if the atom object has a reference object or not.
            All known atoms should have reference objects.
        """
        return self.reference is not None
