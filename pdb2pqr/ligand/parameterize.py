"""Calculating and assigning ligand charges and radii."""
import logging
from .mol2 import Mol2Molecule
from .peoe import equilibrate
from . import PARSE_RADII
from ..forcefield import ForcefieldAtom
from ..forcefield import ForcefieldResidue


_LOGGER = logging.getLogger(__name__)


class ParameterizedMolecule(Mol2Molecule):
    """Ligand with charge and radius assignments."""

    def __init__(self):
        super().__init__()
        self.ligand_properties = {}

    def update(self, ligand):
        """Update self with latest version of ligand (if needed).

        Args:
            ligand:  latest version of ligand
        """
        prev_atom_names = set(self.ligand_properties)
        curr_atom_names = {a.name for a in ligand.atoms}
        if len(prev_atom_names ^ curr_atom_names) > 0:
            for atom in ligand.atoms:
                atom.formal_charge = 0.0
            self.reparameterize(ligand)

    def read(self, mol2_file):
        """Routines for reading MOL2 file.

        Args:
            mol2_file:  file-like object with MOL2 data.
        """
        super().read(mol2_file)
        atom_names = set()
        duplicates = set()
        for atom in self.atoms:
            if atom.name in atom_names:
                duplicates.add(atom.name)
            else:
                atom_names.add(atom.name)
        if len(duplicates) > 0:
            err = "Found duplicate atom names: %s" % duplicates
            raise KeyError(err)

    def reparameterize(self, ligand):
        """Reassign parameters given new ligand.

        Args:
            ligand:  latest version of ligand
        """
        self.ligand_properties = {}
        for atom in ligand.atoms:
            atom.charge = atom.formal_charge
        ligand.atoms = equilibrate(ligand.atoms)
        for atom in ligand.atoms:
            elem = atom.atom_type.split(".")[0].upper()
            charge = atom.charge
            try:
                radius = PARSE_RADII[elem]
                atom.radius = radius
            except KeyError:
                raise KeyError(
                    "Unable to assign radius for element %s in atom %s" % (
                        elem, atom))
            self.ligand_properties[atom.name] = {
                "charge": charge, "radius": radius}
