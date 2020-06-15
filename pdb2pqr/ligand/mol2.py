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
    def __init__(self, bond_from, bond_to, bond_type, bond_id=0):
        """Initialize bond.

        Args:
            bond_from:  bond from this atom
            bond_to:  bond to this atom
            bond_type:  type of bond:  1 (single), 2 (double), or ar (aromatic)
            bond_id:  integer ID of bond
        """
        self.bond_to_self = bond_to
        self.bond_from_self = bond_from
        self.type = bond_type
        self.bond_id = bond_id


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
        self.bonded_atoms = None
        self.torsions = None
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
        self.serial = None
        self.name = None
        self.res_name = None
        self.res_seq = None
        self.x = None
        self.y = None
        self.z = None

    def get_torsions(self, start_atom):
        """Get the torsion angles that start with this atom (name).

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

    def read(self, mol2_file):
        """Routines for reading MOL2 file.

        Args:
            mol2_file:  file-like object with MOL2 data.
        """
        mol2_file = self.parse_atoms(mol2_file)
        mol2_file = self.parse_bonds(mol2_file)
        self.create_bonded_atoms()

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
            raise KeyError("Found duplicate atoms names in MOL2 file: %s" %
                duplicates)
        return mol2_file

    def parse_bonds(self, mol2_file):
        """Parse @<TRIPOS>BOND section of file.

        Args:
            mol2_file:  file-like object with MOL2 data.
        Returns:
            file object advanced to bonds section
        Raises:
            ValueError for problems parsing bond information
        """
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
            try:
                bond_from = int(words[1])
                bond_to = int(words[2])
                bond_id = int(words[0])
                bond = Mol2Bond(
                    bond_from=bond_from, bond_to=bond_to, bond_type=bond_type,
                    bond_id=bond_id)
            except ValueError as exc:
                err = "Got error (%s) when parsing bond line: %s" % (exc, line)
                raise ValueError(err)
            self.bonds.append(bond)
        return mol2_file

    def create_bonded_atoms(self):
        """Create a list of bonded atoms and torsions for each atom."""
        atom_names = list(self.atoms.keys())
        for bond in self.bonds:
            from_atom_name = atom_names[bond.bond_from_self-1]
            from_atom = self.atoms[from_atom_name]
            to_atom_name = atom_names[bond.bond_to_self-1]
            to_atom = self.atoms[to_atom_name]
            if from_atom.bonded_atoms is None:
                from_atom.bonded_atoms = [to_atom]
            else:
                from_atom.bonded_atoms.append(to_atom)
            if to_atom.bonded_atoms is None:
                to_atom.bonded_atoms = [from_atom]
            else:
                to_atom.bonded_atoms.append(from_atom)
        for atom_name, atom in self.atoms.items():
            atom.torsions = self.get_torsions(atom_name)
