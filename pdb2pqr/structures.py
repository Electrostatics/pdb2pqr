"""Simple biomolecular structures.

This module contains the simpler structure objects used in PDB2PQR and their
associated methods.

.. codeauthor:: Todd Dolinsky
.. codeauthor:: Nathan Baker
"""

from typing import Self

from .config import BACKBONE
from .pdb import ATOM, HETATM


class Chain:
    """Chain class

    The chain class contains information about each chain within a given
    :class:`Biomolecule` object.
    """

    def __init__(self, chain_id: str):
        """Initialize the class.

        :param chain_id:  ID for this chain as denoted in the PDB
        :type chain_id:  str
        """
        self.chain_id = chain_id
        self.residues = []
        self.name = None

    def add_residue(self, residue):
        """Add a residue to the chain

        :param residue:  residue to be added
        :type residue:  Residue
        """
        self.residues.append(residue)

    def renumber_residues(self):
        """Renumber atoms.

        Renumber based on actual residue number and not PDB :makevar:`res_seq`
        """
        for count, residue in enumerate(self.residues, start=1):
            residue.set_res_seq(count)

    @property
    def atoms(self):
        """Return a list of Atom objects contained in this chain.

        :return: list of Atom objects
        :rtype: [Atom]
        """
        atomlist = []
        for residue in self.residues:
            my_list = residue.atoms
            for atom in my_list:
                atomlist.append(atom)
        return atomlist

    def __str__(self):
        output = [residue.letter_code() for residue in self.residues]
        return "".join(output)


class Atom:
    """Represent an atom.

    The Atom class inherits from the :class:`ATOM` object in :mod:`pdb`.
    This class used for adding fields not found in the PDB that may be useful
    for analysis.
    This class also simplifies code by combining :class:`ATOM` and
    :class:`HETATM` objects into a single class.
    """

    def __init__(
        self,
        atom: ATOM | HETATM | Self | None = None,
        type_="ATOM",
        residue=None,
    ):
        """Initialize the new Atom object by using the old object.

        :param atom:  the original ATOM object (could be None)
        :type atom:  ATOM
        :param type_:  either ATOM or HETATM
        :type type_:  str
        :param residue:  a pointer back to the parent residue object (could be
            None)
        :type residue:  Residue
        """
        self.type = None
        self.serial = None
        self.name = None
        self.alt_loc = None
        self.res_name = None
        self.chain_id = None
        self.res_seq = None
        self.ins_code = None
        self.x = None
        self.y = None
        self.z = None
        self.occupancy = None
        self.temp_factor = None
        self.seg_id = None
        self.element = None
        self.charge = None
        self.bonds = []
        self.reference = None
        self.residue = None
        self.radius = None
        self.ffcharge = None
        self.hdonor = 0
        self.hacceptor = 0
        self.cell = None
        self.added = 0
        self.optimizeable = 0
        self.refdistance = 0
        self.id = None
        self.mol2charge = None
        if type_ in ["ATOM", "HETATM"]:
            self.type = type_
        else:
            err = f"Invalid atom type {type_} (Atom Class IN structures.py)!"
            raise ValueError(err)
        if atom is not None:
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
            self.residue = residue
            if isinstance(atom, ATOM):
                # ATOM class doesn't have mol2charge
                self.mol2charge = None
            else:
                self.mol2charge = atom.mol2charge

    @classmethod
    def from_pqr_line(cls, line):
        """Create an atom from a PQR line.

        :param cls:  class for classmethod
        :type cls:  Atom
        :param line:  PQR line
        :type line:  str
        :returns:  new atom or None (for REMARK and similar lines)
        :rtype:  Atom
        :raises ValueError:  for problems parsing
        """
        atom = cls()
        words = [w.strip() for w in line.split()]
        token = words.pop(0)
        if token in [
            "REMARK",
            "TER",
            "END",
            "HEADER",
            "TITLE",
            "COMPND",
            "SOURCE",
            "KEYWDS",
            "EXPDTA",
            "AUTHOR",
            "REVDAT",
            "JRNL",
        ]:
            return None
        if token in ["ATOM", "HETATM"]:
            atom.type = token
        elif token[:4] == "ATOM":
            atom.type = "ATOM"
            words = [token[4:], *words]
        elif token[:6] == "HETATM":
            atom.type = "HETATM"
            words = [token[6:], *words]
        else:
            err = f"Unable to parse line: {line}"
            raise ValueError(err)
        atom.serial = int(words.pop(0))
        atom.name = words.pop(0)
        atom.res_name = words.pop(0)
        token = words.pop(0)
        try:
            atom.res_seq = int(token)
        except ValueError:
            atom.chain_id = token
            atom.res_seq = int(words.pop(0))
        token = words.pop(0)
        try:
            atom.x = float(token)
        except ValueError:
            atom.ins_code = token
            atom.x = float(words.pop(0))
        atom.y = float(words.pop(0))
        atom.z = float(words.pop(0))
        atom.charge = float(words.pop(0))
        atom.radius = float(words.pop(0))
        return atom

    @classmethod
    def from_qcd_line(cls, line, atom_serial):
        """Create an atom from a QCD (UHBD QCARD format) line.

        :param Atom cls:  class for classmethod
        :param str line:  PQR line
        :param int atom_serial:  atom serial number
        :returns:  new atom or None (for REMARK and similar lines)
        :rtype:  Atom
        :raises ValueError:  for problems parsing
        """
        atom = cls()
        words = [w.strip() for w in line.split()]
        token = words.pop(0)
        if token in [
            "REMARK",
            "TER",
            "END",
            "HEADER",
            "TITLE",
            "COMPND",
            "SOURCE",
            "KEYWDS",
            "EXPDTA",
            "AUTHOR",
            "REVDAT",
            "JRNL",
        ]:
            return None
        if token in ["ATOM", "HETATM"]:
            atom.type = token
        elif token[:4] == "ATOM":
            atom.type = "ATOM"
            words = [token[4:], *words]
        elif token[:6] == "HETATM":
            atom.type = "HETATM"
            words = [token[6:], *words]
        else:
            err = f"Unable to parse line: {line}"
            raise ValueError(err)
        atom.serial = int(atom_serial)
        atom.res_seq = int(words.pop(0))
        atom.res_name = words.pop(0)
        atom.name = words.pop(0)
        atom.x = float(words.pop(0))
        atom.y = float(words.pop(0))
        atom.z = float(words.pop(0))
        atom.charge = float(words.pop(0))
        atom.radius = float(words.pop(0))
        return atom

    def get_common_string_rep(self, chainflag=False):
        """Returns a string of the common column of the new atom type.

        Uses the :class:`ATOM` string output but changes the first field to
        either be ``ATOM`` or ``HETATM`` as necessary.
        This is used to create the output for PQR and PDB files.

        :return:  string with ATOM/HETATM field set appropriately
        :rtype:  str
        """
        outstr = ""
        tstr = self.type
        outstr += str.ljust(tstr, 6)[:6]
        tstr = f"{self.serial:d}"
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
        tstr = self.chain_id if chainflag else ""
        outstr += str.ljust(tstr, 1)[:1]
        tstr = f"{self.res_seq:d}"
        outstr += str.rjust(tstr, 4)[:4]
        outstr += f"{self.ins_code}   " if self.ins_code != "" else "    "
        tstr = f"{self.x:8.3f}"
        outstr += str.ljust(tstr, 8)[:8]
        tstr = f"{self.y:8.3f}"
        outstr += str.ljust(tstr, 8)[:8]
        tstr = f"{self.z:8.3f}"
        outstr += str.ljust(tstr, 8)[:8]
        return outstr

    def __str__(self):
        return self.get_pqr_string()

    def get_pqr_string(self, chainflag=False):
        """Returns a string of the atom type.

        Uses the :class:`ATOM` string output but changes the first field to
        either be ``ATOM`` or ``HETATM`` as necessary.
        This is used to create the output for PQR files.

        :return:  string with ATOM/HETATM field set appropriately
        :rtype:  str
        """
        outstr = self.get_common_string_rep(chainflag=chainflag)
        ffcharge = (
            f"{self.ffcharge:.4f}" if self.ffcharge is not None else "0.0000"
        )
        outstr += str.rjust(ffcharge, 8)[:8]
        ffradius = (
            f"{self.radius:.4f}" if self.radius is not None else "0.0000"
        )
        outstr += str.rjust(ffradius, 7)[:7]
        return outstr

    def get_pdb_string(self):
        """Returns a string of the atom type.

        Uses the :class:`ATOM` string output but changes the first field to
        either be ``ATOM`` or ``HETATM`` as necessary.
        This is for the PDB representation of the atom.
        The :mod:`propka` module depends on this being correct.

        :return:  string with ATOM/HETATM field set appropriately
        :rtype:  str
        """
        outstr = self.get_common_string_rep(chainflag=True)
        outstr += f"{self.occupancy:>6.2f}{self.temp_factor:>6.2f}      "
        outstr += f"{self.seg_id:4.4s}{self.element:>2.2s}{self.charge:2.2s}"
        return outstr

    @property
    def coords(self):
        """Return the x,y,z coordinates of the atom.

        .. todo:
           All atom coordinates should be converted to :mod:`numpy` arrays

        :return:  list of the coordinates
        :rtype:  [float, float, float]
        """
        return [self.x, self.y, self.z]

    def add_bond(self, bondedatom):
        """Add a bond to the list of bonds.

        :param bondedatom:  the atom to bond to
        :type bondedatom:  ATOM
        """
        self.bonds.append(bondedatom)

    @property
    def is_hydrogen(self):
        """Is this atom a Hydrogen atom?

        :return:  whether this atom is a hydrogen
        :rtype:  bool
        """
        return self.name[0] == "H"

    @property
    def is_backbone(self):
        """Return True if atom name is in backbone, otherwise False.

        :return:  whether atom is in backbone
        :rtype:  bool
        """
        return self.name in BACKBONE

    @property
    def has_reference(self):
        """Determine if the object has a reference object or not.

        All known atoms should have reference objects.

        :return:  whether atom has reference object
        :rtype:  bool
        """
        return self.reference is not None
