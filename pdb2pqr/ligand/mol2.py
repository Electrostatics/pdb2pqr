"""Support molecules in Tripos MOL2 format.

For further information look at (web page exists: 25 August 2005):
http://www.tripos.com/index.php?family=modules,SimplePage,,,&page=sup_mol2&s=0
"""

import logging
from collections import OrderedDict
from itertools import combinations

from numpy import array
from numpy.linalg import norm

from . import NONBONDED_BY_TYPE, RADII, VALENCE_BY_ELEMENT, peoe

_LOGGER = logging.getLogger(__name__)


# These are the allowed bond types
BOND_TYPES = {"single", "double", "triple", "aromatic"}
# This is the maximum deviation from an ideal bond distance
BOND_DIST = 2.0


class Mol2Bond:
    """MOL2 molecule bonds."""

    def __init__(self, atom1, atom2, bond_type, bond_id=0):
        """Initialize bond.

        :param atom1:  name of first atom in bond
        :type atom1:  str
        :param atom2:  name of second atom in bond
        :type atom2:  str
        :param bond_type:  type of bond:  1 (single), 2 (double), or ar
            (aromatic)
        :type bond_type:  int
        :param bond_id:  integer ID of bond
        :type bond_id:  int
        """
        self.atoms = (atom1, atom2)
        self.bond_id = int(bond_id)
        if bond_type in BOND_TYPES:
            self.type = bond_type
        else:
            err = f"Unknown bond type: {bond_type}"
            raise ValueError(err)

    @property
    def atom_names(self):
        """Get atom names in bond.

        :return:  tuple with names of atoms in bond.
        :rtype:  (str, str)
        """
        return (self.atoms[0].name, self.atoms[1].name)

    @property
    def length(self):
        """Get bond length.

        :return:  bond length
        :rtype:  float
        """
        return self.atoms[0].distance(self.atoms[1])

    def __str__(self):
        mol2 = f"{self.atoms[0].name:s} {self.type:s}-bonded to "
        mol2 += f"{self.atoms[1].name:s}"
        return mol2


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
        self.type = None
        self.radius = None
        self.is_c_term = False
        self.is_n_term = False
        self.mol2charge = None
        self.occupancy = 0.00
        self.temp_factor = 0.00
        self.seg_id = None
        self.charge = None
        self.num_rings = 0
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

    def distance(self, other):
        """Get distance between two atoms.

        :param other:  other atom for distance measurement
        :type other:  Mol2Atom
        :return:  distance
        :rtype:  float
        """
        return norm(other.coords - self.coords)

    def __str__(self):
        """Generate PDB line from MOL2."""
        mol2 = (
            f"HETATM{self.serial:5d}{self.name:>5s}{self.res_name:>4s} L"
            f"{self.res_seq!s:>5s}   {self.x:8.3f}{self.y:8.3f}{self.z:8.3f}"
        )
        return mol2

    def assign_radius(self, primary_dict, secondary_dict):
        """Assign radius to atom.

        .. todo::
           It seems inconsistent that this function pulls radii from a
           dictionary and the biomolecule routines use force field files.

        :param primary_dict:  primary dictionary of radii indexed by atom
            type or element
        :type primary_dict:  dict
        :param secondary_dict:  backup dictionary for radii not found in
            primary dictionary
        :type secondary_dict:  dict
        """
        radius = None
        for rdict in [primary_dict, secondary_dict]:
            if radius is not None:
                break
            for key in [self.type, self.element]:
                if key in rdict:
                    radius = rdict[key]
                    break
        if radius is not None:
            self.radius = radius
        else:
            err = (
                f"Unable to find radius parameter for self of type "
                f"{self.type} in radius dictionary: {primary_dict}"
            )
            raise KeyError(err)

    @property
    def coords(self):
        """Coordinates.

        :return:  coordinates
        :rtype:  numpy.ndarray
        """
        return array([self.x, self.y, self.z])

    @property
    def bonded_atom_names(self):
        """Bonded atom names.

        :return:  bonded atom names
        :rtype:  list
        """
        return [a.name for a in self.bonded_atoms]

    @property
    def num_bonded_heavy(self):
        """Number of heavy atoms bonded to this atom.

        :return:  number of heavy atoms
        :rtype: int
        """
        return len([a for a in self.bonded_atoms if a.type != "H"])

    @property
    def num_bonded_hydrogen(self):
        """Number of hydrogen atoms bonded to this atom.

        :return:  number of hydrogen atoms
        :rtype:  int
        """
        return len([a for a in self.bonded_atoms if a.type == "H"])

    @property
    def element(self):
        """Element for this atom (uppercase).

        :return:  element for this atom
        :rtype:  str
        """
        return self.type.split(".")[0].upper()

    @property
    def bond_order(self):
        """Total number of electrons in bonds with other atoms.

        :return:  total number of electrons in bonds with other atoms
        :rtype:  int
        """
        order = 0
        num_aromatic = 0
        for bond in self.bonds:
            if bond.type == "single":
                order += 1
            elif bond.type == "double":
                order += 2
            elif bond.type == "triple":
                order += 3
            elif bond.type == "aromatic":
                num_aromatic += 1
            else:
                err = f"Unknown bond type: {bond.type}"
                raise ValueError(err)
        if num_aromatic > 0:
            order = order + num_aromatic + 1
        return order

    @property
    def formal_charge(self):
        """Formal charge for this atom

        :return:  formal charge for this atom
        :rtype:  int
        """
        element = self.type.split(".")[0]
        valence = VALENCE_BY_ELEMENT[element]
        nonbonded = NONBONDED_BY_TYPE[self.type]
        bond_order = self.bond_order
        formal_charge = valence - nonbonded - bond_order
        if (
            (self.type in ["N.pl3", "N.am"])
            and (bond_order == 3)
            and (formal_charge != 0)
        ):
            # Planar nitrogen bond orders are not always correct in MOL2
            _LOGGER.warning("Correcting planar/amide bond order.")
            formal_charge = 0
        elif (
            (self.type in ["N.ar"])
            and (bond_order == 4)
            and (formal_charge != 0)
        ):
            # Aromatic nitrogen bond orders are not always correct in MOL2
            _LOGGER.warning("Correcting aromatic nitrogen bond order.")
            formal_charge = 0
        elif (
            (self.type in ["C.ar"])
            and (bond_order == 5)
            and (formal_charge != 0)
        ):
            # Aromatic carbon bond orders are not always correct in MOL2
            _LOGGER.warning("Correcting aromatic carbon bond order.")
            formal_charge = 0
        elif (
            (self.type in ["O.co2"])
            and (bond_order == 1)
            and (formal_charge != -0.5)
        ):
            # CO2 bond orders are hardly ever set correctly in MOL2
            formal_charge = -0.5
        elif (
            (self.type in ["C.2"])
            and (bond_order == 5)
            and (formal_charge == -1)
        ):
            # CO2 bond orders are hardly ever set correctly in MOL2
            formal_charge = 0
        elif (
            (self.type in ["N.3"])
            and (bond_order == 4)
            and (formal_charge == -1)
        ):
            # Tetravalent nitrogen atom types are sometimes wrong in MOL2
            _LOGGER.warning("Correcting ammonium atom type.")
            formal_charge = 1
        elif (
            (self.type in ["O.3"])
            and (bond_order == 1)
            and (formal_charge == 1)
        ):
            # Phosphate groups are sometimes confused in MOL2
            # Assign negative charge to first O.3 with bond order 1
            # attached to phosphorous
            elements = [a.type[0] for a in self.bonds[0].atoms]
            p_atom = self.bonds[0].atoms[elements.index("P")]
            _LOGGER.warning("Correcting phosphate bond order.")
            o_atoms = []
            for bond in p_atom.bonds:
                for atom in bond.atoms:
                    if atom.type[0] == "O" and atom.bond_order == 1:
                        o_atoms.append(atom.name)
            formal_charge = -1 if o_atoms.index(self.name) == 0 else 0
        return formal_charge


class Mol2Molecule:
    """Tripos MOL2 molecule."""

    def __init__(self):
        self.atoms = OrderedDict()
        self.bonds = []
        self.torsions = set()
        self.rings = set()
        self.serial = None
        self.name = None
        self.res_name = None
        self.res_seq = None

    def assign_parameters(
        self, primary_dict=RADII["zap9"], secondary_dict=RADII["bondi"]
    ):
        """Assign charges and radii to atoms in molecule.

        Args:
            primary_dict:  primary dictionary of radii indexed by atom type or
                           element
            secondary_dict:  backup dictionary for radii not found in primary
                             dictionary
        """
        self.assign_radii(primary_dict, secondary_dict)
        self.assign_charges()

    def assign_radii(self, primary_dict, secondary_dict):
        """Assign radii to atoms in molecule.

        :param primary_dict:  primary dictionary of radii indexed by atom
            type or element
        :type primary_dict:  dict
        :param secondary_dict:  backup dictionary for radii not found in
            primary dictionary
        :type secondary_dict:  dict
        """
        for atom in self.atoms.values():
            atom.assign_radius(primary_dict, secondary_dict)

    def assign_charges(self):
        """Assign charges to atoms in molecule."""
        for atom in self.atoms.values():
            atom.charge = atom.formal_charge
        peoe.equilibrate(self.atoms.values())

    def find_atom_torsions(self, start_atom):
        """Set the torsion angles that start with this atom (name).

        :param start_atom:  starting atom name
        :type start_atom:  str
        :return: list of 4-tuples containing atom names comprising torsions
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

    @staticmethod
    def rotate_to_smallest(path):
        """Rotate cycle path so that it begins with the smallest node.

        This was borrowed from StackOverflow: https://j.mp/2AHaukj

        :param path:  list of atom names
        :type path:  list of str
        :return:  rotated path
        :rtype:  list of str
        """
        n = path.index(min(path))
        return path[n:] + path[:n]

    def find_new_rings(self, path, rings, level=0):
        """Find new rings in molecule.

        This was borrowed from StackOverflow: https://j.mp/2AHaukj

        :param path:  list of atom names
        :type path:  list of str
        :param rings:  current list of rings
        :type rings:  list of str
        :param level:  recursion level
        :type level:  int
        :return:  new list of rings
        :rtype:  int
        """
        start_node = path[0]
        next_node = None
        sub_path = []
        for bond in self.bonds:
            atom1 = bond.atoms[0].name
            atom2 = bond.atoms[1].name
            if start_node in (atom1, atom2):
                next_node = atom2 if atom1 == start_node else atom1
                if next_node not in path:
                    sub_path = [next_node, *path]
                    rings = self.find_new_rings(sub_path, rings, level + 1)
                elif len(path) > 2 and next_node == path[-1]:
                    path_ = self.rotate_to_smallest(path)
                    inv_path = tuple(self.rotate_to_smallest(path_[::-1]))
                    path_ = tuple(path_)
                    if (path_ not in rings) and (inv_path not in rings):
                        rings.add(tuple(path_))
        return rings

    def set_rings(self):
        """Set all rings in molecule.

        This was borrowed from StackOverflow: https://j.mp/2AHaukj
        """
        self.rings = set()
        rings = set()
        # Generate all rings
        for bond in self.bonds:
            for atom_name in bond.atom_names:
                rings = self.find_new_rings([atom_name], rings)
        # Prune rings that are products of other rings
        # TODO - testing on molecules like phenalene shows that this is broken
        ring_sets = []
        for i in range(2, len(rings) + 1):
            for combo in combinations(rings, i):
                ring_set = set().union(*combo)
                ring_sets.append(ring_set)
        for ring in rings:
            ring_set = set(ring)
            if ring_set in ring_sets:
                _LOGGER.debug(f"Fused ring: {ring}")
            else:
                _LOGGER.debug(f"Unfused ring: {ring}")
                self.rings.add(ring)
        for ring in self.rings:
            for atom in ring:
                self.atoms[atom].num_rings += 1

    def read(self, mol2_file):
        """Routines for reading MOL2 file.

        :param mol2_file:  file-like object with MOL2 data
        """
        mol2_file = self.parse_atoms(mol2_file)
        mol2_file = self.parse_bonds(mol2_file)

    def parse_atoms(self, mol2_file):
        """Parse @<TRIPOS>ATOM section of file.

        :param mol2_file:  file-like object with MOL2 data
        :return:  file-like object advanced to bonds section
        :raises ValueError:  for bad MOL2 ATOM lines
        :raises TypeError:  for bad charge entries
        """
        # Skip material before atoms section
        for line in mol2_file:
            if "@<TRIPOS>ATOM" in line:
                break
            _LOGGER.debug(f"Skipping: {line.strip()}")
        duplicates = set()
        for line in mol2_file:
            line = line.strip()
            if not line:
                continue
            if "@<TRIPOS>BOND" in line:
                break
            words = line.split()
            if len(words) < 8:
                err = f"Bad entry in MOL2 file: {line}"
                raise ValueError(err)
            atom = Mol2Atom()
            atom.name = words[1]
            atom_type = words[5]
            type_parts = atom_type.split(".")
            type_parts[0] = type_parts[0].capitalize()
            if len(type_parts) == 2:
                type_parts[1] = type_parts[1].lower()
            elif len(type_parts) > 2:
                err = f"Invalid atom type: {atom_type}"
                raise ValueError(err)
            atom.type = ".".join(type_parts)
            atom.chain_id = "L"
            try:
                atom.serial = int(words[0])
                atom.res_name = words[7][:4]
                atom.res_seq = int(words[6])
                atom.x = float(words[2])
                atom.y = float(words[3])
                atom.z = float(words[4])
            except ValueError as exc:
                err = f"Error ({exc}) parsing atom line: {line}"
                raise ValueError(err)
            if len(line) > 8:
                try:
                    atom.mol2charge = float(words[8])
                except TypeError:
                    err = f"Unable to parse {words[8]} as charge in atom "
                    err += f"line: {line}"
                    _LOGGER.warning(err)
            if atom.name in self.atoms:
                duplicates.add(atom.name)
            else:
                self.atoms[atom.name] = atom
        if duplicates:
            raise KeyError(
                f"Found duplicate atoms names in MOL2 file: {duplicates}"
            )
        return mol2_file

    def parse_bonds(self, mol2_file):
        """Parse @<TRIPOS>BOND section of file.

        Atoms must already have been parsed.
        Also sets up torsions and rings.

        :param mol2_file:  file-like object with MOL2 data
        :return:  file-like object advanced to SUBSTRUCTURE section
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
                err = f"Bond line too short: {line}"
                raise ValueError(err)
            bond_type = words[3]
            if bond_type == "1":
                bond_type = "single"
            elif bond_type == "2":
                bond_type = "double"
            elif bond_type == "3":
                bond_type = "triple"
            elif bond_type == "am":
                raise NotImplementedError(
                    "PDB2PQR does not currently support the amide (am) bond "
                    "type."
                )
            elif bond_type == "ar":
                bond_type = "aromatic"
            elif bond_type == "du":
                raise NotImplementedError(
                    "PDB2PQR does not currently support the dummy (du) bond "
                    "type."
                )
            elif bond_type == "un":
                raise NotImplementedError(
                    "PDB2PQR does not currently support the unknown (un) bond "
                    "type."
                )
            elif bond_type == "nc":
                raise NotImplementedError(
                    "PDB2PQR does not currently support the not-connected "
                    "(nc) bond type."
                )
            else:
                err = f"Unknown bond type: {bond_type}"
                raise ValueError(err)
            bond_id = int(words[0])
            atom_id1 = int(words[1])
            atom_id2 = int(words[2])
            atom_name1 = atom_names[atom_id1 - 1]
            atom1 = self.atoms[atom_name1]
            atom_name2 = atom_names[atom_id2 - 1]
            atom2 = self.atoms[atom_name2]
            bond = Mol2Bond(
                atom1=atom1, atom2=atom2, bond_type=bond_type, bond_id=bond_id
            )
            atom1.bonds.append(bond)
            atom1.bonded_atom_names.append(atom_name2)
            atom1.bonded_atoms.append(atom2)
            atom2.bonds.append(bond)
            atom2.bonded_atom_names.append(atom_name1)
            atom2.bonded_atoms.append(atom1)
            self.bonds.append(bond)
        self.set_torsions()
        self.set_rings()
        return mol2_file
