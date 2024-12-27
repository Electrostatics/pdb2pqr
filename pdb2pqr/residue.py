"""Biomolecular residue class.

.. codeauthor:: Todd Dolinsky
.. codeauthor:: Nathan Baker
"""

import logging

from . import pdb, structures
from . import quatfit as quat
from . import utilities as util

_LOGGER = logging.getLogger(__name__)


class Residue:
    """Residue class

    .. todo::
       Should this class have a member variable for dihedrals? Pylint
       complains!

    The residue class contains a list of Atom objects associated with that
    residue and other helper functions.
    """

    def __init__(self, atoms: list[pdb.ATOM | pdb.HETATM]):
        """Initialize the class

        :param atoms:  list of atom-like (:class:`HETATM` or :class:`ATOM`)
            objects to be stored
        :type atoms:  list
        """
        sample_atom = atoms[-1]
        self.atoms: list[structures.Atom] = []
        self.name = sample_atom.res_name
        self.chain_id = sample_atom.chain_id
        self.res_seq = sample_atom.res_seq
        self.ins_code = sample_atom.ins_code
        self.map: dict[str, structures.Atom] = {}
        self.naname = None
        self.reference = None
        self.is_n_term = None
        self.is_c_term = None
        self.dihedrals = []
        atomclass = ""
        for atom in atoms:
            if isinstance(atom, pdb.ATOM):
                atomclass = "ATOM"
            elif isinstance(atom, pdb.HETATM):
                atomclass = "HETATM"
            atom = structures.Atom(atom, atomclass, self)
            atomname = atom.name
            if atomname not in self.map:
                self.add_atom(atom)
            else:  # Don't add duplicate atom
                oldatom = self.get_atom(atomname)
                oldatom.alt_loc = ""
        if self.name == "HOH":
            self.name = "WAT"
            for atom in self.atoms:
                atom.res_name = "WAT"

    def __str__(self):
        return f"{self.name} {self.chain_id} {self.res_seq}{self.ins_code}"

    def get_moveable_names(self, pivot):
        """Return all atom names that are further away than the pivot atom.

        :param residue:  the residue to use
        :type residue:  Residue
        :param pivot:  the pivot atomname
        :type pivot:  str
        :return:  names of atoms further away than pivot atom
        :rtype:  [str]
        """
        refdist = self.get_atom(pivot).refdistance
        return [atom.name for atom in self.atoms if atom.refdistance > refdist]

    def update_terminus_status(self):
        """Update the :makevar:`is_n_terms` and :makevar:`is_c_term` flags."""
        # If Nterm then update counter of hydrogens
        if self.is_n_term:
            count = 0
            atoms = ["H", "H2", "H3"]
            for atom in atoms:
                for atom2 in self.atoms:
                    atomname = atom2.name
                    if atom == atomname:
                        count += 1
            self.is_n_term = count
        # If Cterm then update counter
        if self.is_c_term:
            self.is_c_term = None
            for atom in self.atoms:
                atomname = atom.name
                if atomname == "HO":
                    self.is_c_term = 2
                    break
            if not self.is_c_term:
                self.is_c_term = 1

    def set_res_seq(self, value):
        """Change the residue sequence number.

        Set the atom field :makevar:`res_seq` and change the residue's
        information.
        The :makevar:`icode` field is no longer useful.

        :param value:  the new value of :makevar:`res_seq`
        :type value:  int
        """
        self.ins_code = ""
        self.res_seq = value
        for atom in self.atoms:
            atom.res_seq = value

    def set_chain_id(self, value):
        """Set the chain ID.

        :param value:  new :makevar:`chain_id` value
        :type value:  str
        """
        self.chain_id = value
        for atom in self.atoms:
            atom.chain_id = value

    def add_atom(self, atom):
        """Add the atom object to the residue.

        :param atom: atom-like object, e.g., :class:`HETATM` or :class:`ATOM`
        """
        self.atoms.append(atom)
        self.map[atom.name] = atom

    def remove_atom(self, atomname):
        """Remove an atom from the residue object.

        :param atomname:  the name of the atom to be removed
        :type atomname:  str
        """
        # Delete the atom from the map
        atom = self.map[atomname]
        bonds = atom.bonds
        del self.map[atomname]
        # Delete the atom from the list
        self.atoms.remove(atom)
        # Delete all instances of the atom as a bond
        for bondatom in bonds:
            if atom in bondatom.bonds:
                bondatom.bonds.remove(atom)
        del atom

    def rename_atom(self, oldname, newname):
        """Rename an atom to a new name.

        :param oldname:  old atom name
        :type oldname:  str
        :param newname:  new atom name
        :type newname:  str
        """
        atom = self.map[oldname]
        atom.name = newname
        self.map[newname] = atom
        del self.map[oldname]

    def get_atom(self, name: str) -> structures.Atom | None:
        """Retrieve a residue atom based on its name.

        :param resname:  name of the residue to retrieve
        :type resname:  str
        :return:  residue
        :rtype:  structures.Atom | None
        """
        return self.map.get(name)

    def has_atom(self, name):
        """Return True if atom in residue.

        :param name:  atom name in question
        :type name:  str
        :return:  True if atom in residue
        :rtype:  bool
        """
        return name in self.map

    @property
    def charge(self):
        """Get the total charge of the residue.

        In order to get rid of floating point rounding error, do a string
        transformation.

        Returns:
            charge: The charge of the residue (float)
        """
        charge = (atom.ffcharge for atom in self.atoms if atom.ffcharge)
        charge = sum(charge)
        charge = float(f"{charge:.4f}")
        return charge

    def rename_residue(self, name):
        """Rename the residue.

        :param name:  the new name of the residue
        :type name:  str
        """
        self.name = name
        for atom in self.atoms:
            atom.res_name = name

    @classmethod
    def rotate_tetrahedral(cls, atom1, atom2, angle):
        """Rotate about the atom1-atom2 bond by a given angle.

        All atoms connected to atom2 will rotate.

        :param atom1:  first atom of the bond to rotate about
        :type atom1:  structures.Atom
        :param atom2:  second atom of the bond to rotate about
        :type atom2:  structures.Atom
        :param angle:  degrees to rotate
        :type angle:  float
        """
        moveatoms = []
        movecoords = []
        initcoords = util.subtract(atom2.coords, atom1.coords)
        # Determine which atoms to rotate
        for atom in atom2.bonds:
            if atom == atom1:
                continue
            moveatoms.append(atom)
            movecoords.append(util.subtract(atom.coords, atom1.coords))
        newcoords = quat.qchichange(initcoords, movecoords, angle)
        for iatom, atom in enumerate(moveatoms):
            atom.x = newcoords[iatom][0] + atom1.x
            atom.y = newcoords[iatom][1] + atom1.y
            atom.z = newcoords[iatom][2] + atom1.z

    def pick_dihedral_angle(self, conflict_names, oldnum=None):
        """Choose an angle number to use in debumping.

        Instead of simply picking a random chiangle, this function
        uses a more intelligent method to improve efficiency.
        The algorithm uses the names of the conflicting atoms
        within the residue to determine which angle number
        has the best chance of fixing the problem(s). The method
        also insures that the same chiangle will not be run twice
        in a row.

        :param residue:  residue that is being debumped
        :type residue:  Residue
        :param conflict_names: list of atom names that are currently
            conflicts
        :type conflict_names: [str]
        :param oldnum:  old dihedral angle number
        :type oldnum:  int
        :return:  new dihedral angle number
        :rtype:  int
        """
        bestnum = -1
        best = 0
        ilist = list(range(len(self.dihedrals)))
        if oldnum is not None and oldnum >= 0 and ilist:
            del ilist[oldnum]
            test_dihedral_indices = ilist[oldnum:] + ilist[:oldnum]
        else:
            test_dihedral_indices = ilist
        for i in test_dihedral_indices:
            if i == oldnum:
                continue
            if self.dihedrals[i] is None:
                continue
            score = 0
            atomnames = self.reference.dihedrals[i].split()
            pivot = atomnames[2]
            moveablenames = self.get_moveable_names(pivot)
            if conflict_names == moveablenames:
                return i
            for name in conflict_names:
                if name in moveablenames:
                    score += 1
                    if score > best:
                        best = score
                        bestnum = i
        return bestnum

    def set_donors_acceptors(self):
        """Set the donors and acceptors within the residue."""
        if self.reference is None:
            return
        for atom in self.atoms:
            atomname = atom.name
            atom.hdonor = False
            atom.hacceptor = False
            if atomname.startswith("N"):
                bonded = 0
                for bondedatom in atom.bonds:
                    if bondedatom.is_hydrogen:
                        atom.hdonor = True
                        bonded = 1
                        break
                if not bonded and self.reference.name == "HIS":
                    atom.hacceptor = True

            elif atomname.startswith("O") or (
                atomname.startswith("S") and self.reference.name == "CYS"
            ):
                atom.hacceptor = True
                for bondedatom in atom.bonds:
                    if bondedatom.is_hydrogen:
                        atom.hdonor = True
                        break

    def reorder(self):
        """Reorder the atoms to start with N, CA, C, O if they exist."""
        templist = []
        if self.has_atom("N"):
            templist.append(self.get_atom("N"))
        if self.has_atom("CA"):
            templist.append(self.get_atom("CA"))
        if self.has_atom("C"):
            templist.append(self.get_atom("C"))
        if self.has_atom("O"):
            templist.append(self.get_atom("O"))
        # Add remaining atoms
        for atom in self.atoms:
            if atom.name not in ["N", "CA", "C", "O"]:
                templist.append(atom)
        # Change the list pointer
        self.atoms = templist[:]

    def letter_code(self) -> str:
        """Default letter code for residue.

        :return:  letter code for residue
        :rtype:  str
        """
        return "X"
