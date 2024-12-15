"""Routines for biomolecule optimization.

.. codeauthor::  Jens Erik Nielsen
.. codeauthor::  Todd Dolinsky
.. codeauthor::  Yong Huang
.. codeauthor::  Nathan Baker
"""

import logging

from . import aa, cells, io
from . import quatfit as quat
from . import utilities as util
from .config import (
    BUMP_HEAVY_SIZE,
    BUMP_HYDROGEN_SIZE,
    CELL_SIZE,
    DEBUMP_ANGLE_STEP_SIZE,
    DEBUMP_ANGLE_STEPS,
    DEBUMP_ANGLE_TEST_COUNT,
    SMALL_NUMBER,
)

_LOGGER = logging.getLogger(__name__)
_LOGGER.addFilter(io.DuplicateFilter())


class Debump:
    """Grab bag of random stuff that apparently didn't fit elsewhere.

    .. todo::
        This class needs to be susbtantially refactored in to multiple classes
        with clear responsibilities.
    """

    def __init__(self, biomolecule, definition=None):
        """Initialize the Debump object.

        :param biomolecule:  the biomolecule to debump
        :type biomolecule:  Biomolecule
        :param definition:  topology definition file
        :type definition:  Definition
        """
        self.biomolecule = biomolecule
        self.definition = definition
        self.aadef = None
        self.cells = {}
        if definition is not None:
            self.aadef = definition.getAA()
            self.nadef = definition.getNA()

    def get_bump_score(self, residue):
        """Get a bump score for a residue.

        :param residue:  residue with bumping to evaluate
        :type residue:  Residue
        :return:  bump score
        :rtype:  float
        """
        # Do some setup
        self.cells = cells.Cells(CELL_SIZE)
        self.cells.assign_cells(self.biomolecule)
        self.biomolecule.calculate_dihedral_angles()
        self.biomolecule.set_donors_acceptors()
        self.biomolecule.update_internal_bonds()
        self.biomolecule.set_reference_distance()
        bumpscore = 0.0
        if not isinstance(residue, aa.Amino):
            return 0.0
        # Initialize variables
        for atom in residue.atoms:
            atomname = atom.name
            if atomname[0] != "H":
                continue
            bumpscore += self.get_bump_score_atom(atom)
        return bumpscore

    def get_bump_score_atom(self, atom):
        """Find nearby atoms for conflict-checking.

        Uses neighboring cells to compare atoms rather than an all versus all
        O(n^2) algorithm, which saves a great deal of time.  There are several
        instances where we ignore potential conflicts; these include
        donor/acceptor pairs, atoms in the same residue, and bonded CYS
        bridges.

        :param atom:  find nearby atoms to this atom
        :type atom:  Atom
        :return:  a bump score sum((dist-cutoff)**20 for all nearby atoms
        :rtype:  float
        """
        # Initialize some variables
        residue = atom.residue
        atom_size = BUMP_HYDROGEN_SIZE if atom.is_hydrogen else BUMP_HEAVY_SIZE
        # Get atoms from nearby cells
        closeatoms = self.cells.get_near_cells(atom)
        # Loop through and see if any are within the cutoff
        bumpscore = 0.0
        for closeatom in closeatoms:
            closeresidue = closeatom.residue
            if closeresidue == residue and (
                closeatom in atom.bonds or atom in closeatom.bonds
            ):
                continue
            if not isinstance(closeresidue, aa.Amino):
                continue
            if (
                isinstance(residue, aa.CYS)
                and residue.ss_bonded_partner == closeatom
            ):
                continue
            if (
                atom.is_hydrogen
                and len(atom.bonds) != 0
                and atom.bonds[0].hdonor
                and closeatom.hacceptor
            ):
                continue
            if (
                closeatom.is_hydrogen
                and len(closeatom.bonds) != 0
                and closeatom.bonds[0].hdonor
                and atom.hacceptor
            ):
                continue
            dist = util.distance(atom.coords, closeatom.coords)
            other_size = (
                BUMP_HYDROGEN_SIZE
                if closeatom.is_hydrogen
                else BUMP_HEAVY_SIZE
            )
            cutoff = atom_size + other_size
            if dist < cutoff:
                bumpscore += 1000.0
        _LOGGER.debug(f"BUMPSCORE {bumpscore!s}")
        return bumpscore

    def debump_biomolecule(self):
        """Minimize bump score for molecule.

        Make sure that none of the added atoms were rebuilt on top of existing
        atoms. See each called function for more information.

        :raises ValueError:  if missing (backbone) atoms are encountered
        """
        # Do some setup
        self.cells = cells.Cells(CELL_SIZE)
        self.cells.assign_cells(self.biomolecule)
        self.biomolecule.calculate_dihedral_angles()
        self.biomolecule.set_donors_acceptors()
        self.biomolecule.update_internal_bonds()
        try:
            self.biomolecule.set_reference_distance()
        except ValueError as err:
            err = f"Biomolecular structure is incomplete:  {err}"
            raise ValueError(err)
        # Determine which residues to debump
        for residue in self.biomolecule.residues:
            if not isinstance(residue, aa.Amino):
                continue
            # Initialize variables
            conflict_names = self.find_residue_conflicts(
                residue, write_conflict_info=True
            )
            if not conflict_names:
                continue
            # Otherwise debump the residue
            _LOGGER.debug(f"Starting to debump {residue}...")
            _LOGGER.debug(
                f"Debumping cutoffs: {BUMP_HEAVY_SIZE * 2:2.1f} for "
                f"heavy-heavy, {BUMP_HYDROGEN_SIZE + BUMP_HEAVY_SIZE:2.1f} "
                f"for hydrogen-heavy, and {BUMP_HYDROGEN_SIZE * 2:2.1f} "
                f"for hydrogen-hydrogen."
            )
            if self.debump_residue(residue, conflict_names):
                _LOGGER.debug("Debumping Successful!")
            else:
                text = f"WARNING: Unable to debump {residue}"
                _LOGGER.warning(text)
        _LOGGER.debug("Done checking if we must debump any residues.")

    def find_residue_conflicts(self, residue, *, write_conflict_info=False):
        """Find conflicts between residues.

        :param residue:  residue to check
        :type residue:  Residue
        :param write_conflict_info:  write verbose output about conflict
        :type write_conflict_info:  bool
        :return: list of conflicts
        :rtype: [str]
        """
        conflict_names = []
        for atom in residue.atoms:
            atomname = atom.name
            if not atom.added:
                continue
            if atomname == "H":
                continue
            if atom.optimizeable:
                continue
            nearatoms = self.find_nearby_atoms(atom)
            # If something is too close, we must debump the residue
            if nearatoms != {}:
                conflict_names.append(atomname)
                if write_conflict_info:
                    for repatom in nearatoms:
                        _LOGGER.debug(
                            f"{residue} {atomname} is too close to "
                            f"{repatom.residue} {repatom.name}"
                        )
        return conflict_names

    def score_dihedral_angle(self, residue, anglenum):
        """Assign score to dihedral angle.

        :param residue:  residue with dihedral angles to score
        :type residue:  Residue
        :param anglenum:  specific dihedral angle index
        :type anglenum:  int
        :return:  score for dihedral angle
        :rtype:  float
        """
        score = 0
        atomnames = residue.reference.dihedrals[anglenum].split()
        pivot = atomnames[2]
        moveablenames = residue.get_moveable_names(pivot)
        for name in moveablenames:
            nearatoms = self.find_nearby_atoms(residue.get_atom(name))
            for value in nearatoms.values():
                score += value
        return score

    def debump_residue(self, residue, conflict_names):
        """Debump a specific residue.

        Only should be called if the residue has been detected to have a
        conflict.
        If called, try to rotate about dihedral angles to resolve the conflict.

        :param residue:  the residue in question
        :type residue:  Residue
        :param conflict_names:  a list of atomnames that were rebuilt too close
            to other atoms
        :type conflict_names:  [str]
        :return: True if successful, False otherwise
        :rtype: bool
        """
        # Initialize some variables
        anglenum = -1
        curr_conflict_names = conflict_names
        # Try to find a workable solution
        for _ in range(DEBUMP_ANGLE_TEST_COUNT):
            anglenum = residue.pick_dihedral_angle(
                curr_conflict_names, anglenum
            )
            if anglenum == -1:
                return False
            _LOGGER.debug(
                f"Using dihedral angle number {anglenum} to debump "
                "the residue."
            )
            bestscore = self.score_dihedral_angle(residue, anglenum)
            found_improved = False
            bestangle = orig_angle = residue.dihedrals[anglenum]
            # Skip the first angle as it's already known.
            for i in range(1, DEBUMP_ANGLE_STEPS):
                newangle = orig_angle + (DEBUMP_ANGLE_STEP_SIZE * i)
                self.set_dihedral_angle(residue, anglenum, newangle)
                # Check for conflicts
                score = self.score_dihedral_angle(residue, anglenum)
                if score == 0:
                    if not self.find_residue_conflicts(residue):
                        _LOGGER.debug(
                            f"No conflicts found at angle {newangle!r}"
                        )
                        return True
                    else:
                        bestangle = newangle
                        found_improved = True
                        break
                # Set the best angle
                elif score < bestscore:
                    diff = abs(bestscore - score)
                    # Don't update if it's effectively a tie
                    if diff > SMALL_NUMBER:
                        bestscore = score
                        bestangle = newangle
                        found_improved = True
            self.set_dihedral_angle(residue, anglenum, bestangle)
            curr_conflict_names = self.find_residue_conflicts(residue)
            if found_improved:
                err = f"Best score of {bestscore} at angle {bestangle}."
                _LOGGER.debug(err)
                _LOGGER.debug(f"New conflict set: {curr_conflict_names}")
            else:
                _LOGGER.debug("No improvement found for this dihedral angle.")
        # If we're here, debumping was unsuccessful
        return False

    def get_closest_atom(self, atom):
        """Get the closest atom that does not form a donor/acceptor pair.

        Used to detect potential conflicts.

        .. note::  Cells must be set before using this function.

        :param atom:  the atom to test
        :type atom:  Atom
        :return:  the closest atom to the input atom that does not satisfy a
            donor/acceptor pair.
        :rtype:  Atom
        """
        # Initialize some variables
        bestdist = 999.99
        bestwatdist = 999.99
        bestatom = None
        bestwatatom = None
        residue = atom.residue
        # Get atoms from nearby cells
        closeatoms = self.cells.get_near_cells(atom)
        # Loop through and see which is the closest
        for closeatom in closeatoms:
            closeresidue = closeatom.residue
            if closeresidue == residue:
                continue
            if not isinstance(closeresidue, (aa.Amino, aa.WAT)):
                continue
            if (
                isinstance(residue, aa.CYS)
                and residue.ss_bonded_partner == closeatom
            ):
                continue
            # Also ignore if this is a donor/acceptor pair
            if (
                atom.is_hydrogen
                and atom.bonds[0].hdonor
                and closeatom.hacceptor
            ):
                continue
            if (
                closeatom.is_hydrogen
                and closeatom.bonds[0].hdonor
                and atom.hacceptor
            ):
                continue
            dist = util.distance(atom.coords, closeatom.coords)
            if isinstance(closeresidue, aa.WAT):
                if dist < bestwatdist:
                    bestwatdist = dist
                    bestwatatom = closeatom
            else:
                if dist < bestdist:
                    bestdist = dist
                    bestatom = closeatom
        if bestdist > bestwatdist:
            txt = (
                f"Skipped atom during water optimization: {bestwatatom.name} "
                f"in {bestwatatom.residue} skipped "
                f"when optimizing {atom.name} in {residue}"
            )
            _LOGGER.warning(txt)
        return bestatom

    def find_nearby_atoms(self, atom):
        """Find nearby atoms for conflict-checking.

        Uses neighboring cells to compare atoms rather than an all versus all
        O(n^2) algorithm, which saves a great deal of time.
        There are several instances where we ignore potential conflicts; these
        include donor/acceptor pairs, atoms in the same residue, and bonded CYS
        bridges.

        :param atom:  find nearby atoms to this atom
        :type atom:  Atom
        :return:  a dictionary of ``Atom too close`` to ``amount of overlap
            for that atom``
        :rtype:  dict
        """
        # Initialize some variables
        nearatoms = {}
        residue = atom.residue
        atom_size = BUMP_HYDROGEN_SIZE if atom.is_hydrogen else BUMP_HEAVY_SIZE
        # Get atoms from nearby cells
        closeatoms = self.cells.get_near_cells(atom)
        # Loop through and see if any are within the cutoff
        for closeatom in closeatoms:
            closeresidue = closeatom.residue
            if closeresidue == residue and (
                closeatom in atom.bonds or atom in closeatom.bonds
            ):
                continue
            if not isinstance(closeresidue, (aa.Amino, aa.WAT)):
                continue
            if (
                isinstance(residue, aa.CYS)
                and residue.ss_bonded_partner == closeatom
            ):
                continue
            # Also ignore if this is a donor/acceptor pair
            if (
                atom.is_hydrogen
                and len(atom.bonds) != 0
                and atom.bonds[0].hdonor
                and closeatom.hacceptor
            ):
                continue
            if (
                closeatom.is_hydrogen
                and len(closeatom.bonds) != 0
                and closeatom.bonds[0].hdonor
                and atom.hacceptor
            ):
                continue
            dist = util.distance(atom.coords, closeatom.coords)
            other_size = (
                BUMP_HYDROGEN_SIZE
                if closeatom.is_hydrogen
                else BUMP_HEAVY_SIZE
            )
            cutoff = atom_size + other_size
            if dist < cutoff:
                nearatoms[closeatom] = cutoff - dist
        return nearatoms

    def set_dihedral_angle(self, residue, anglenum, angle):
        """Rotate a residue about a given angle.

        Uses :mod:`quatfit` methods to perform the matrix mathematics.

        :param residue:  residue to rotate
        :type residue:  Residue
        :param anglenum:  specific dihedral angle index
        :param angle:  the angle to set the dihedral to
        :type angle:  float
        """
        coordlist = []
        initcoords = []
        movecoords = []
        pivot = ""
        oldangle = residue.dihedrals[anglenum]
        diff = angle - oldangle
        atomnames = residue.reference.dihedrals[anglenum].split()
        pivot = atomnames[2]
        for atomname in atomnames:
            if residue.has_atom(atomname):
                coordlist.append(residue.get_atom(atomname).coords)
            else:
                raise ValueError("Error occurred while trying to debump!")
        initcoords = util.subtract(coordlist[2], coordlist[1])
        moveablenames = residue.get_moveable_names(pivot)
        for name in moveablenames:
            atom = residue.get_atom(name)
            movecoords.append(util.subtract(atom.coords, coordlist[1]))
        newcoords = quat.qchichange(initcoords, movecoords, diff)
        for iatom, atom_name in enumerate(moveablenames):
            atom = residue.get_atom(atom_name)
            self.cells.remove_cell(atom)
            atom.x = newcoords[iatom][0] + coordlist[1][0]
            atom.y = newcoords[iatom][1] + coordlist[1][1]
            atom.z = newcoords[iatom][2] + coordlist[1][2]
            self.cells.add_cell(atom)
        # Set the new angle
        coordlist = []
        for atomname in atomnames:
            if residue.has_atom(atomname):
                coordlist.append(residue.get_atom(atomname).coords)
            else:
                raise ValueError("Error occurred while trying to debump!")
        dihed = util.dihedral(
            coordlist[0], coordlist[1], coordlist[2], coordlist[3]
        )
        residue.dihedrals[anglenum] = dihed
