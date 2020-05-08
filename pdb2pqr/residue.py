"""Biomolecular residue class.

Author:  Todd Dolinsky
"""
from . import pdb
from . import structures
from . import utilities as util
from . import quatfit as quat


class Residue(object):
    """Residue class

    TODO - move this class to a separate file

    The residue class contains a list of Atom objects associated with that
    residue and other helper functions.
    """

    def __init__(self, atoms):
        """Initialize the class

        Args:
            atoms:  list of Atom objects to be stored in this class (list)
        """
        sample_atom = atoms[-1]
        self.atoms = []
        self.name = sample_atom.res_name
        self.chain_id = sample_atom.chain_id
        self.res_seq = sample_atom.res_seq
        self.ins_code = sample_atom.ins_code
        self.map = {}
        self.naname = None
        self.reference = None
        self.is_n_term = None
        self.is_c_term = None

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
            else: # Don't add duplicate atom
                oldatom = self.get_atom(atomname)
                oldatom.alt_loc = ""

        if self.name == "HOH":
            self.name = "WAT"
            for atom in self.atoms:
                atom.res_name = "WAT"

    def __str__(self):
        text = "%s %s %i%s" % (self.name, self.chain_id, self.res_seq, self.ins_code)
        return text

    def get_moveable_names(self, pivot):
        """Return all atomnames that are further away than the pivot atom.

        Parameters
            residue:  The residue to use
            pivot:    The pivot atomname
        """
        movenames = []
        refdist = self.get_atom(pivot).refdistance
        for atom in self.atoms:
            if atom.refdistance > refdist:
                movenames.append(atom.name)
        return movenames

    def update_terminus_status(self):
        """Update the is_n_terms and is_c_term flags"""
        # If Nterm then update counter of hydrogens
        if self.is_n_term:
            count = 0
            atoms = ['H', 'H2', 'H3']
            for atom in atoms:
                for atom2 in self.atoms:
                    atomname = atom2.name
                    if atom == atomname:
                        count = count+1
            self.is_n_term = count

        # If Cterm then update counter
        if self.is_c_term:
            self.is_c_term = None
            for atom in self.atoms:
                atomname = atom.name
                if atomname == 'HO':
                    self.is_c_term = 2
                    break
            if not self.is_c_term:
                self.is_c_term = 1

    def set_res_seq(self, value):
        """Set the atom field res_seq to a certain value and change the
        residue's information.  The icode field is no longer useful.

        Args:
            value:  The new value of res_seq (int)
        """
        self.ins_code = ""
        self.res_seq = value
        for atom in self.atoms:
            atom.res_seq = value

    def set_chain_id(self, value):
        """Set the chain_id field to a certain value"""
        self.chain_id = value
        for atom in self.atoms:
            atom.chain_id = value

    def add_atom(self, atom):
        """Add the atom object to the residue.

        Args:
            atom: The object to be added (ATOM)
        """
        self.atoms.append(atom)
        self.map[atom.name] = atom

    def remove_atom(self, atomname):
        """Remove an atom from the residue object.

        Args:
            atomname: The name of the atom to be removed (string)
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
        """Rename an atom to a new name

        Args:
            oldname: The old atom name (string)
            newname: The new atom name (string)
        """
        atom = self.map[oldname]
        atom.name = newname
        self.map[newname] = atom
        del self.map[oldname]

    def create_atom(self, name, newcoords, type_):
        """Add a new atom object to the residue.

        Uses an atom currently in the residue to seed the new atom object,
        then replaces the coordinates and name accordingly.

        Args:
            name:  The name of the new atom (string)
            newcoords:  The x,y,z coordinates of the new atom (list)
            type:  The type of atom, ATOM or HETATM
        """
        oldatom = self.atoms[0]
        newatom = structures.Atom(oldatom, type_, self)
        newatom.x = newcoords[0]
        newatom.y = newcoords[1]
        newatom.z = newcoords[2]
        newatom.name = name
        newatom.occupancy = 1.00
        newatom.temp_factor = 0.00
        self.add_atom(newatom)

    def get_atom(self, name):
        """Retrieve an atom from the mapping

        Args:
            resname: The name of the residue to retrieve (string)
        """
        return self.map.get(name)

    def has_atom(self, name):
        """Return True if atom in residue."""
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
        charge = float("%.4f" % charge)
        return charge

    def rename_residue(self, name):
        """Rename a given residue

        Args:
            name:  The new name of the residue
        """
        self.name = name
        for atom in self.atoms:
            atom.res_name = name

    @classmethod
    def rotate_tetrahedral(cls, atom1, atom2, angle):
        """Rotate about the atom1-atom2 bond by a given angle All atoms connected
        to atom2 will rotate.

        Args:
            atom1:  The first atom of the bond to rotate about (atom)
            atom2:  The second atom of the bond to rotate about (atom)
            angle:  The number of degrees to rotate (float)
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
            x = newcoords[iatom][0] + atom1.x
            y = newcoords[iatom][1] + atom1.y
            z = newcoords[iatom][2] + atom1.z
            atom.x = x
            atom.y = y
            atom.z = z

    def set_donors_acceptors(self):
        """Set the donors and acceptors within the residue"""
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

            elif atomname.startswith("O") or \
                 (atomname.startswith("S") and self.reference.name == "CYS"):
                atom.hacceptor = True
                for bondedatom in atom.bonds:
                    if bondedatom.is_hydrogen:
                        atom.hdonor = True
                        break

    def reorder(self):
        """Reorder the atoms to start with N, CA, C, O if they exist"""
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

    @classmethod
    def letter_code(cls):
        """Letter code for residue."""
        return 'X'
