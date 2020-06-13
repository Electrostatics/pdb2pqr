"""Support molecules in Tripos MOL2 format.

    For further information look at (web page exists: 25 August 2005):
    http://www.tripos.com/index.php?family=modules,SimplePage,,,&page=sup_mol2&s=0
"""
import logging
import copy
from ..pdb import HETATM


_LOGGER = logging.getLogger(__name__)


class Mol2Bond(object):
    """Bonding of MOL2 files"""
    def __init__(self, frm, to, type_, bond_id=0):
        self.bond_to_self = to # bond to this atom
        self.bond_from_self = frm # bond from atom
        self.type = type_ # 1=single, 2=double, ar=aromatic
        self.bond_id = bond_id # bond_id


class Mol2Molecule(object):
    """Tripos MOL2 molecule"""
    def __init__(self):
        self.l_atoms = [] # all atoms of class <ATOM>
        self.l_bonds = [] # all bonds of class <BOND>
        self.l_pdb_atoms = [] # PDB-like list of all atoms
        self.serial = None
        self.name = None
        self.res_name = None
        self.res_seq = None
        self.x = None
        self.y = None
        self.z = None

    def read(self, file_):
        """Routines for reading MOL2 file"""
        data = file_.read()
        data = data.replace("\r\n", "\n")
        data = data.replace("\r", "\n")

        # ATOM section
        start = data.find("@<TRIPOS>ATOM")
        stop = data.find("@<TRIPOS>BOND")

        # Do some error checking
        if start == -1:
            raise ValueError("Unable to find '@<TRIPOS>ATOM' in MOL2 file!")
        elif stop == -1:
            raise ValueError("Unable to find '@<TRIPOS>BOND' in MOL2 file!")

        atoms = data[start+14:stop-2].split("\n")
        # BOND section
        start = data.find("@<TRIPOS>BOND")
        stop = data.find("@<TRIPOS>SUBSTRUCTURE")

        # More error checking
        if stop == -1:
            raise ValueError("Unable to find '@<TRIPOS>SUBSTRUCTURE' in MOL2 file!")

        bonds = data[start+14:stop-1].split("\n")
        self.parse_atoms(atoms)
        self.parse_bonds(bonds)
        self.create_bonded_atoms()

    def parse_atoms(self, atom_list):
        """For parsing @<TRIPOS>ATOM"""
        for atom_line in atom_list:
            separated_atom_line = atom_line.split()

            # Special handling for blank lines
            if len(separated_atom_line) == 0:
                continue

            # Error checking
            if len(separated_atom_line) < 8:
                raise ValueError("Bad atom entry in MOL2 file: %s" % atom_line)

            fake_record = "HETATM"
            fake_chain = " L"

            try:
                mol2pdb = '%s%5i%5s%4s%2s%4i    %8.3f%8.3f%8.3f' % \
                    (fake_record, int(separated_atom_line[0]),
                     separated_atom_line[1], separated_atom_line[7][:4],
                     fake_chain, int(separated_atom_line[6]),
                     float(separated_atom_line[2]), float(separated_atom_line[3]),
                     float(separated_atom_line[4]))

            except ValueError:
                raise ValueError("Bad atom entry in MOL2 file: %s" % atom_line)

            this_atom = HETATM(mol2pdb, separated_atom_line[5], [], [])
            if len(separated_atom_line) > 8:
                charge = separated_atom_line[8]
                try:
                    this_atom.mol2charge = float(charge)
                except TypeError:
                    _LOGGER.warning('Warning. Non-float charge (%s) in mol2 file.', charge)
                    this_atom.mol2charge = None
            self.l_pdb_atoms.append(mol2pdb)
            self.l_atoms.append(this_atom)

    def parse_bonds(self, bond_list):
        """For parsing @<TRIPOS>BOND"""
        for bond_line in bond_list:
            separated_bond_line = bond_line.split()
            # Special handling for blank lines
            if len(separated_bond_line) == 0:
                continue
            if len(separated_bond_line) < 4:
                raise ValueError("Bad bond entry in MOL2 file: %s" % bond_line)
            try:
                this_bond = Mol2Bond(
                    int(separated_bond_line[1]), # bond frm
                    int(separated_bond_line[2]), # bond to
                    separated_bond_line[3],      # bond type
                    int(separated_bond_line[0])  # bond id
                    )
            except ValueError:
                raise ValueError("Bad bond entry in MOL2 file: %s" % bond_line)
            self.l_bonds.append(this_bond)

    def create_bonded_atoms(self):
        """Creates for each atom a list of the bonded Atoms

        This becomes one attribute of MOL2ATOM!
        """
        for bond in self.l_bonds:
            self.l_atoms[bond.bond_from_self-1].l_bonded_atoms\
                .append(self.l_atoms[bond.bond_to_self-1])

            self.l_atoms[bond.bond_to_self-1].l_bonded_atoms\
                .append(self.l_atoms[bond.bond_from_self-1])

            atbond = copy.deepcopy(bond)
            atbond.other_atom = self.l_atoms[bond.bond_to_self-1]
            self.l_atoms[bond.bond_from_self-1].l_bonds.append(atbond)

            atbond = copy.deepcopy(bond)
            atbond.other_atom = self.l_atoms[bond.bond_from_self-1]
            self.l_atoms[bond.bond_to_self-1].l_bonds.append(atbond)

    def create_pdb_line_from_mol2(self):
        """Generate PDB line from MOL2."""
        raise NotImplementedError("TODO - FIX THIS CODE")
        # fake_type = "HETATM"
        # rstr = "%s%5i%5s%4s%2s%5s   %8.3f%8.3f%8.3f\n" % (fake_type, self.serial,
        #                                                   self.name, self.res_name, ' L',
        #                                                   self.res_seq, self.x, self.y, self.z)
        # return rstr
