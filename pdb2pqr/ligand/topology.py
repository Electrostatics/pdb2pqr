"""Ligand topology classes."""

import logging

_LOGGER = logging.getLogger(__name__)


class Topology:
    """Ligand topology class."""

    def __init__(self, molecule):
        """Initialize with molecule.

        :param molecule:  Mol2Molecule object
        :type molecule:  Mol2Molecule
        """
        self.atom_dict = {}
        for atom in molecule.atoms:
            self.atom_dict[atom.name] = atom
            raise NotImplementedError
