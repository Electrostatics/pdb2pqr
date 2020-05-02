"""Structures for PDB2PQR

This module contains the structure objects used in PDB2PQR and their associated
methods.

Author: Todd Dolinsky
"""
import string
from .pdb import ATOM, HETATM
from .utilities import subtract
from .quatfit import qchichange


# TODO - why is the backbone defined here?
BACKBONE = ["N", "CA", "C", "O", "O2", "HA", "HN", "H", "tN"]


class Chain:
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

    def addResidue(self, residue):
        """Add a residue to the chain

        Args:
            residue: The residue to be added (Residue)
        """
        self.residues.append(residue)

    def renumberResidues(self):
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
            myList = residue.atoms
            for atom in myList:
                atomlist.append(atom)
        return atomlist

    def __str__(self):
        output = []
        for residue in self.residues:
            output.append(residue.letterCode())
        return ''.join(output)


class Residue:
    """Residue class

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

        atomclass = ""
        for a in atoms:
            if isinstance(a, ATOM):
                atomclass = "ATOM"
            elif isinstance(a, HETATM):
                atomclass = "HETATM"
            atom = Atom(a, atomclass, self)
            atomname = atom.name
            if atomname not in self.map:
                self.add_atom(atom)
            else: # Don't add duplicate atom
                oldatom = self.get_atom(atomname)
                oldatom.set("alt_loc", "")

        if self.name == "HOH":
            self.name = "WAT"
            for atom in self.atoms:
                atom.set("res_name", "WAT")

    def __str__(self):
        text = "%s %s %i%s" % (self.name, self.chain_id, self.res_seq, self.ins_code)
        return text

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
        atom.set("name", newname)
        self.map[newname] = atom
        del self.map[oldname]

    def create_atom(self, name, newcoords, type):
        """Add a new atom object to the residue.
        
        Uses an atom currently in the residue to seed the new atom object,
        then replaces the coordinates and name accordingly.

        Args:
            name:  The name of the new atom (string)
            newcoords:  The x,y,z coordinates of the new atom (list)
            type:  The type of atom, ATOM or HETATM
        """
        oldatom = self.atoms[0]
        newatom = Atom(oldatom, type, self)
        newatom.set("x", newcoords[0])
        newatom.set("y", newcoords[1])
        newatom.set("z", newcoords[2])
        newatom.set("name", name)
        newatom.set("occupancy", 1.00)
        newatom.set("temp_factor", 0.00)
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

    def rotate_tetrahedral(self, atom1, atom2, angle):
        """
            Rotate about the atom1-atom2 bond by a given angle
            All atoms connected to atom2 will rotate.

            Parameters:
                atom1:  The first atom of the bond to rotate about (atom)
                atom2:  The second atom of the bond to rotate about (atom)
                angle:  The number of degrees to rotate (float)
        """
        moveatoms = []
        movecoords = []
        initcoords = subtract(atom2.getCoords(), atom1.getCoords())

        # Determine which atoms to rotate
        for atom in atom2.bonds:
            if atom == atom1: continue
            moveatoms.append(atom)
            movecoords.append(subtract(atom.getCoords(), atom1.getCoords()))

        newcoords = qchichange(initcoords, movecoords, angle)
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
            resname = self.name
            atom.hdonor = False
            atom.hacceptor = False

            if atomname.startswith("N"):
                bonded = 0
                for bondedatom in atom.bonds:
                    if bondedatom.isHydrogen():
                        atom.hdonor = True
                        bonded = 1
                        break
                if not bonded and self.reference.name == "HIS":
                    atom.hacceptor = True

            elif atomname.startswith("O") or \
                 (atomname.startswith("S") and self.reference.name == "CYS"):
                atom.hacceptor = True
                for bondedatom in atom.bonds:
                    if bondedatom.isHydrogen():
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

    def letterCode(self):
        return 'X'


class Atom(ATOM):
    """Class Atom

        The Atom class inherits off the ATOM object in pdb.py.  It is used
        for adding fields not found in the pdb that may be useful for analysis.
        Also simplifies code by combining ATOM and HETATM objects into a
        single class.
    """

    def __init__(self, atom, type, residue):
        """Initialize the new Atom object by using the old object.

        Args:
            atom:  The original ATOM object (ATOM)
            type:  Either ATOM or HETATM (string)
            residue:  A pointer back to the parent residue object (Residue)
        """
        if type == "ATOM" or type == "HETATM":
            self.type = type
        else:
            raise PDBInternalError("Invalid atom type %s (Atom Class IN structures.py)!" % type)
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


    def getCommonStringRep(self, chainflag=False):
        """
            Returns a string of the common column of the new atom type.
            Uses the ATOM string output but changes the first field
            to either by ATOM or HETATM as necessary.

            This is used to create the output for pqr and pdb files! No touchy!

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
        """
            Returns a string of the new atom type.  Uses the ATOM string
            output but changes the first field to either by ATOM or
            HETATM as necessary.

            This is used to create the output for pqr files! No touchy!

            Returns
                str: String with ATOM/HETATM field set appropriately
        """
        return self.getPQRString()


    def getPQRString(self, chainflag=False):
        """
            Returns a string of the new atom type.  Uses the ATOM string
            output but changes the first field to either by ATOM or
            HETATM as necessary.

            This is used to create the output for pqr files! No touchy!

            Returns
                str: String with ATOM/HETATM field set appropriately
        """
        outstr = self.getCommonStringRep(chainflag=chainflag)

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


    def getPDBString(self):
        """
            Returns a string of the new atom type.  Uses the ATOM string
            output but changes the first field to either by ATOM or
            HETATM as necessary.

            This is for the pdb representation of the atom. The propka30 module
            depends on this being correct. No touchy!

            Returns
                str: String with ATOM/HETATM field set appropriately
        """
        outstr = self.getCommonStringRep(chainflag=True)

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

    #TODO: What? Why? Isn't this Python?
    #Are we really doing attribute access based
    # on dynamic names? A search of the code says no.
    def get(self, name):
        """
            Get a member of the Atom class

            Parameters
                name:       The name of the member (string)
            Possible Values
                type:       The type of Atom (either ATOM or HETATM)
                serial:     Atom serial number
                name:       Atom name
                alt_loc:     Alternate location
                res_name:    Residue name
                chain_id:    Chain identifier
                res_seq:     Residue sequence number
                ins_code:      Code for insertion of residues
                x:          Orthogonal coordinates for X in Angstroms.
                y:          Orthogonal coordinates for Y in Angstroms.
                z:          Orthogonal coordinates for Z in Angstroms.
                occupancy:  Occupancy
                temp_factor: Temperature Factor
                seg_id:      Segment identifier
                element:    Element symbol
                charge:     Charge on the atom
                bonds:      The bonds associated with the atom
                interbonds: The intrabonds associated with the atom
                extrabonds: The extrabonds assocaited with the atom
                residue:    The parent residue of the atom
                radius:     The radius of the atom
                ffcharge:   The forcefield charge on the atom
                hdonor:     Whether the atom is a hydrogen donor
                hacceptor:  Whether the atom is a hydrogen acceptor
            Returns
                item:       The value of the member
        """
        try:
            item = getattr(self, name)
            return item
        except AttributeError:
            message = "Unable to access object \"%s\" in class Atom" % name
            raise PDBInternalError(message)

    def set(self, name, value):
        """
            Set a member of the Atom class

            Parameters
                name:       The name of the member (string)
                value:      The value to set the member to
            Possible Values
                type:       The type of Atom (either ATOM or HETATM)
                serial:     Atom serial number
                name:       Atom name
                alt_loc:     Alternate location
                res_name:    Residue name
                chain_id:    Chain identifier
                res_seq:     Residue sequence number
                ins_code:      Code for insertion of residues
                x:          Orthogonal coordinates for X in Angstroms.
                y:          Orthogonal coordinates for Y in Angstroms.
                z:          Orthogonal coordinates for Z in Angstroms.
                occupancy:  Occupancy
                temp_factor: Temperature Factor
                seg_id:      Segment identifier
                element:    Element symbol
                charge:     Charge on the atom
                residue:    The parent residue of the atom
                radius:     The radius of the atom
                ffcharge:   The forcefield charge on the atom
                hdonor:     Whether the atom is a hydrogen donor
                hacceptor:  Whether the atom is a hydrogen acceptor
        """
        try:
            setattr(self, name, value)
        except AttributeError:
            message = "Unable to set object \"%s\" in class Atom" % name
            raise PDBInternalError(message)

    def getCoords(self):
        """
            Return the x,y,z coordinates of the atom in list form

            Returns
                List of the coordinates (list)
        """
        return [self.x, self.y, self.z]

    def addBond(self, bondedatom):
        """
            Add a bond to the list of bonds

            Parameters:
                bondedatom: The atom to bond to (Atom)
        """
        self.bonds.append(bondedatom)

    def isHydrogen(self):
        """
            Is this atom a Hydrogen atom?

            Returns
                value: 1 if Atom is a Hydrogen, 0 otherwise
        """
        return self.name[0] == "H"

    def isBackbone(self):
        """
            Return true if atom name is in backbone, otherwise false

            Returns
                state: 1 if true, 0 if false
        """
        return self.name in BACKBONE

    def hasReference(self):
        """
            Determine if the atom object has a reference object or not.
            All known atoms should have reference objects.

            Returns
                True if atom has a reference object, False otherwise.
        """

        return self.reference != None
