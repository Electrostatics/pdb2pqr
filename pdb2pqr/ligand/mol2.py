"""Support molecules in Tripos MOL2 format.

    For further information look at (web page exists: 25 August 2005):
    http://www.tripos.com/index.php?family=modules,SimplePage,,,&page=sup_mol2&s=0
"""
import logging
from collections import OrderedDict
import numpy


_LOGGER = logging.getLogger(__name__)


class Mol2Bond:
    """MOL2 molecule bonds."""
    def __init__(self, atom1, atom2, bond_type, bond_id=0):
        """Initialize bond.

        Args:
            atom1:  name of first atom in bond
            atom2:  name of second atom in bond
            bond_type:  type of bond:  1 (single), 2 (double), or ar (aromatic)
            bond_id:  integer ID of bond
        """
        self.atoms = (atom1, atom2)
        self.type = bond_type
        self.bond_id = int(bond_id)

    def __str__(self):
        fmt = "{b.atoms[0]:s} {b.type:s}-bonded to {b.atoms[1]:s}"
        return fmt.format(b=self)


class Mol2Atom:
    """MOL2 molecule atoms."""
    def __init__(self):
        self.serial = None
        self.name = None
        self.alt_loc = None
        self.res_name = None
        self.chain_id = None
        self.res_seq = None
        self.x = None
        self.y = None
        self.z = None
        self.atom_type = None
        self.radius = None
        self.is_c_term = False
        self.is_n_term = False
        self.mol2charge = None
        self.occupancy = 0.00
        self.temp_factor = 0.00
        self.seg_id = None
        self.charge = None
        self.formal_charge = None
        self.radius = None
        self.bonded_atoms = []
        self.bonds = []
        self.torsions = []
        self.rings = []
        # Terms for calculating atom electronegativity
        self.poly_terms = None
        # Atom electronegativity
        self.chi = None
        # Atom charge change during equilibration
        self.delta_charge = None

    def __str__(self):
        """Generate PDB line from MOL2."""
        pdb_fmt = (
            "HETATM{a.serial:5d}{a.name:>5s}{a.res_name:>4s} L"
            "{a.res_seq!s:>5s}   {a.x:8.3f}{a.y:8.3f}{a.z:8.3f}"
        )
        return pdb_fmt.format(a=self)

    @property
    def coords(self):
        """Return coordinates as numpy vector."""
        return numpy.array([self.x, self.y, self.z])

    @property
    def bonded_atom_names(self):
        """Return a list of bonded atom names."""
        return [a.name for a in self.bonded_atoms]

    @property
    def num_bonded_heavy(self):
        """Return the number of heavy atoms bonded to this atom."""
        return len([a for a in self.bonded_atoms if a.atom_type != "H"])

    @property
    def num_bonded_hydrogen(self):
        """Return the number of hydrogen atoms bonded to this atom."""
        return len([a for a in self.bonded_atoms if a.atom_type == "H"])

    @property
    def element(self):
        """Return a string with the element for this atom (uppercase)."""
        return self.atom_type.split(".")[0].upper()


class Mol2Molecule:
    """Tripos MOL2 molecule"""
    def __init__(self):
        self.atoms = OrderedDict()
        self.bonds = []
        self.torsions = set()
        self.rings = set()
        self.serial = None
        self.name = None
        self.res_name = None
        self.res_seq = None

    def find_atom_torsions(self, start_atom):
        """Set the torsion angles that start with this atom (name).

        Args:
            start_atom:  starting atom name
        Returns:
            list of 4-tuples containing atom names comprising torsions
        """
        torsions = []
        for bonded1 in self.atoms[start_atom].bonded_atom_names:
            for bonded2 in self.atoms[bonded1].bonded_atom_names:
                if bonded2 == start_atom:
                    continue
                for end_atom in self.atoms[bonded2].bonded_atom_names:
                    if end_atom == bonded1:
                        continue
                    torsions.append((start_atom, bonded1, bonded2, end_atom))
        return torsions

    def set_torsions(self):
        """Set all torsions in molecule."""
        for atom_name, atom in self.atoms.items():
            atom.torsions = self.find_atom_torsions(atom_name)
            for torsion in atom.torsions:
                self.torsions.add(torsion)

    def set_rings(self):
        """Set all rings in molecule.

        Like many things, this was borrowed from StackOverflow:
        https://stackoverflow.com/questions/12367801/finding-all-cycles-in-undirected-graphs
        """
        for bond in self.bonds:
            _LOGGER.error(str(bond))
            raise NotImplementedError()

    def read(self, mol2_file):
        """Routines for reading MOL2 file.

        Args:
            mol2_file:  file-like object with MOL2 data.
        """
        mol2_file = self.parse_atoms(mol2_file)
        mol2_file = self.parse_bonds(mol2_file)

    def parse_atoms(self, mol2_file):
        """Parse @<TRIPOS>ATOM section of file.

        Args:
            mol2_file:  file-like object with MOL2 data.
        Returns:
            file object advanced to bonds section
        Raises:
            ValueError for bad MOL2 ATOM lines
            TypeError for bad charge entries
        """
        # Skip material before atoms section
        for line in mol2_file:
            if "@<TRIPOS>ATOM" in line:
                break
            _LOGGER.debug("Skipping: %s", line.strip())
        duplicates = set()
        for line in mol2_file:
            line = line.strip()
            if not line:
                continue
            if "@<TRIPOS>BOND" in line:
                break
            words = line.split()
            if len(words) < 8:
                err = "Bad entry in MOL2 file: %s" % line
                raise ValueError(err)
            atom = Mol2Atom()
            atom.name = words[1]
            atom.atom_type = words[5]
            atom.chain_id = "L"
            try:
                atom.serial = int(words[0])
                atom.res_name = words[7][:4]
                atom.res_seq = int(words[6])
                atom.x = float(words[2])
                atom.y = float(words[3])
                atom.z = float(words[4])
            except ValueError as exc:
                err = "Error (%s) parsing atom line: %s" % (exc, line)
                raise ValueError(err)
            if len(line) > 8:
                try:
                    atom.mol2charge = float(words[8])
                except TypeError:
                    err = "Unable to parse %s as charge in atom line: %s" % (
                        words[8], line)
                    _LOGGER.warning(err)
            if atom.name in self.atoms:
                duplicates.add(atom.name)
            else:
                self.atoms[atom.name] = atom
        if len(duplicates) > 0:
            raise KeyError(
                "Found duplicate atoms names in MOL2 file: %s" % duplicates)
        return mol2_file

    def parse_bonds(self, mol2_file):
        """Parse @<TRIPOS>BOND section of file.

        Atoms must already have been parsed.
        Also sets up torsions and rings.

        Args:
            mol2_file:  file-like object with MOL2 data.
        Returns:
            file object advanced to SUBSTRUCTURE section
        """
        atom_names = list(self.atoms.keys())
        for line in mol2_file:
            line = line.strip()
            if not line:
                continue
            if "@<TRIPOS>SUBSTRUCTURE" in line:
                break
            words = line.split()
            if len(words) < 4:
                err = "Bond line too short: %s" % line
                raise ValueError(err)
            bond_type = words[3]
            bond_id = int(words[0])
            atom_id1 = int(words[1])
            atom_id2 = int(words[2])
            atom_name1 = atom_names[atom_id1-1]
            atom1 = self.atoms[atom_name1]
            atom_name2 = atom_names[atom_id2-1]
            atom2 = self.atoms[atom_name2]
            bond = Mol2Bond(
                atom1=atom_name1, atom2=atom_name2, bond_type=bond_type,
                bond_id=bond_id)
            atom1.bonds.append(bond)
            atom1.bonded_atom_names.append(atom_name2)
            atom1.bonded_atoms.append(atom2)
            atom2.bonds.append(bond)
            atom2.bonded_atom_names.append(atom_name1)
            atom2.bonded_atoms.append(atom1)
            self.bonds.append(bond)
        self.set_torsions()
        # self.set_rings()
        return mol2_file
