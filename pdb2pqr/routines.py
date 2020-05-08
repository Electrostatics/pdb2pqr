"""Routines for PDB2PQR

This module contains debumping routines to optimize the biomolecule.

Authors:  Jens Erik Nielsen, Todd Dolinsky, Yong Huang
"""
import logging
from . import aa
from . import utilities as util
from . import io
from . import quatfit as quat
from . import cells
from .config import DEBUMP_ANGLE_STEP_SIZE, DEBUMP_ANGLE_STEPS, DEBUMP_ANGLE_TEST_COUNT
from .config import SMALL_NUMBER, CELL_SIZE, BUMP_HYDROGEN_SIZE, BUMP_HEAVY_SIZE


_LOGGER = logging.getLogger(__name__)
_LOGGER.addFilter(io.DuplicateFilter())


class Debump(object):
    """Grab bag of random stuff that apparently didn't fit elsewhere.

    TODO - needs to be susbtantially refactored in to multiple classes with clear
    responsibilities.
    """
    def __init__(self, protein, definition=None):
        """Initialize the Debump class.

        The class contains most of the main routines that run PDB2PQR

        Parameters
            protein:  The protein to run PDB2PQR on (Protein)
        """
        self.protein = protein
        self.definition = definition
        self.aadef = None
        self.cells = {}
        if definition is not None:
            self.aadef = definition.getAA()
            self.nadef = definition.getNA()

    def get_bump_score(self, residue):
        """Get an bump score for the current structure"""

        # Do some setup
        self.cells = cells.Cells(CELL_SIZE)
        self.cells.assign_cells(self.protein)

        self.protein.calculate_dihedral_angles()
        self.set_donors_acceptors()
        self.protein.update_internal_bonds()
        self.protein.set_reference_distance()
        bumpscore = 0.0
        if not isinstance(residue, aa.Amino):
            return 0.0

        # Initialize variables
        for atom in residue.atoms:
            atomname = atom.name
            if atomname[0] != "H":
                continue
            bumpscore = bumpscore + self.get_bump_score_atom(atom)
        return bumpscore

    def get_bump_score_atom(self, atom):
        """Find nearby atoms for conflict-checking.

        Uses neighboring cells to compare atoms rather than an all versus all
        O(n^2) algorithm, which saves a great deal of time.  There are several
        instances where we ignore potential conflicts; these include donor/acceptor
        pairs, atoms in the same residue, and bonded CYS bridges.

        Parameters
            atom:  Find nearby atoms to this atom (Atom)
        Returns
            bumpscore: a bump score sum((dist-cutoff)**20 for all near atoms
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
            if closeresidue == residue and (closeatom in atom.bonds or atom in closeatom.bonds):
                continue

            if not isinstance(closeresidue, aa.Amino):
                continue
            if isinstance(residue, aa.CYS):
                if residue.ss_bonded_partner == closeatom:
                    continue

            # Also ignore if this is a donor/acceptor pair
            pair_ignored = False
            if atom.is_hydrogen and len(atom.bonds) != 0 and atom.bonds[0].hdonor \
               and closeatom.hacceptor:
                continue

            if closeatom.is_hydrogen and len(closeatom.bonds) != 0 and closeatom.bonds[0].hdonor \
                   and atom.hacceptor:
                continue

            dist = util.distance(atom.coords, closeatom.coords)
            other_size = BUMP_HYDROGEN_SIZE if closeatom.is_hydrogen else BUMP_HEAVY_SIZE
            cutoff = atom_size + other_size
            if dist < cutoff:
                bumpscore = bumpscore + 1000.0
                if pair_ignored:
                    _LOGGER.debug('This bump is a donor/acceptor pair.')
        _LOGGER.debug('BUMPSCORE %s', str(bumpscore))
        return bumpscore

    def debump_protein(self):
        """Make sure that none of the added atoms were rebuilt on top of existing
        atoms.  See each called function for more information.
        """
        _LOGGER.info("Checking if we must debump any residues... ")

        # Do some setup
        self.cells = cells.Cells(CELL_SIZE)
        self.cells.assign_cells(self.protein)

        self.protein.calculate_dihedral_angles()
        self.protein.set_donors_acceptors()
        self.protein.update_internal_bonds()
        self.protein.set_reference_distance()

        # Determine which residues to debump
        for residue in self.protein.residues:
            if not isinstance(residue, aa.Amino):
                continue

            # Initialize variables
            conflict_names = self.find_residue_conflicts(residue, True)

            if not conflict_names:
                continue

            # Otherwise debump the residue
            _LOGGER.debug("Starting to debump %s...", residue)
            _LOGGER.debug(("Debumping cutoffs: %2.1f for heavy-heavy, %2.1f for "
                           "hydrogen-heavy, and %2.1f for hydrogen-hydrogen."),
                          BUMP_HEAVY_SIZE*2,
                          BUMP_HYDROGEN_SIZE+BUMP_HEAVY_SIZE,
                          BUMP_HYDROGEN_SIZE*2)
            if self.debump_residue(residue, conflict_names):
                _LOGGER.debug("Debumping Successful!")
            else:
                text = "WARNING: Unable to debump %s", residue
                _LOGGER.warning(text)

        _LOGGER.debug("Done checking if we must debump any residues.")

    def find_residue_conflicts(self, residue, write_conflict_info=False):
        """Find conflicts between residues."""
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
                        _LOGGER.debug("%s %s is too close to %s %s", residue,
                                      atomname, repatom.residue, repatom.name)
        return conflict_names

    def score_dihedral_angle(self, residue, anglenum):
        """Assign score to dihedral angle"""
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
        conflict. If called, try to rotate about dihedral angles to resolve the
        conflict.

        Parameters
            residue:  The residue in question
            conflict_names:  A list of atomnames that were rebuilt too close to other atoms
        Returns
            True if successful, False otherwise
        """
        # Initialize some variables
        anglenum = -1
        curr_conflict_names = conflict_names

        # Try to find a workable solution
        for _ in range(DEBUMP_ANGLE_TEST_COUNT):
            anglenum = residue.pick_dihedral_angle(curr_conflict_names, anglenum)
            if anglenum == -1:
                return False

            _LOGGER.debug("Using dihedral angle number %i to debump the residue.", anglenum)
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
                        _LOGGER.debug("No conflicts found at angle %s", repr(newangle))
                        return True
                    else:
                        bestangle = newangle
                        found_improved = True
                        break

                # Set the best angle
                elif score < bestscore:
                    diff = abs(bestscore - score)
                    #Don't update if it's effectively a tie
                    if diff > SMALL_NUMBER:
                        bestscore = score
                        bestangle = newangle
                        found_improved = True

            self.set_dihedral_angle(residue, anglenum, bestangle)
            curr_conflict_names = self.find_residue_conflicts(residue)

            if found_improved:
                err = "Best score of {best} at angle {angle}."
                err = err.format(best=repr(bestscore), angle=repr(bestangle))
                _LOGGER.debug(err)
                _LOGGER.debug("New conflict set: %s", str(curr_conflict_names))
            else:
                _LOGGER.debug("No improvement found for this dihedral angle.")

        # If we're here, debumping was unsuccessful
        return False

    def get_closest_atom(self, atom):
        """Get the closest atom that does not form a donor/acceptor pair.

        Used to detect potential conflicts.
        NOTE:  Cells must be set before using this function.

        Parameters
            atom:  The atom in question (Atom)
        Returns
            bestatom:  The closest atom to the input atom that does not satisfy
                       a donor/acceptor pair.
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
            if isinstance(residue, aa.CYS):
                if residue.ss_bonded_partner == closeatom:
                    continue

            # Also ignore if this is a donor/acceptor pair
            if atom.is_hydrogen and atom.bonds[0].hdonor and closeatom.hacceptor:
                continue
            if closeatom.is_hydrogen and closeatom.bonds[0].hdonor and atom.hacceptor:
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
            txt = ("Skipped atom during water optimization: %s in %s skipped "
                   "when optimizing %s in %s") % (bestwatatom.name,
                                                  bestwatatom.residue,
                                                  atom.name, residue)
            _LOGGER.warning(txt)
        return bestatom

    def find_nearby_atoms(self, atom):
        """Find nearby atoms for conflict-checking.

        Uses neighboring cells to compare atoms rather than an all versus all
        O(n^2) algorithm, which saves a great deal of time.  There are several
        instances where we ignore potential conflicts; these include donor/acceptor
        pairs, atoms in the same residue, and bonded CYS bridges.

        Parameters
            atom:  Find nearby atoms to this atom (Atom)
        Returns
            nearatoms:  A dictionary of <Atom too close> to <amount of overlap for that atom>.
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
            if closeresidue == residue and (closeatom in atom.bonds or atom in closeatom.bonds):
                continue

            if not isinstance(closeresidue, (aa.Amino, aa.WAT)):
                continue
            if isinstance(residue, aa.CYS) and residue.ss_bonded_partner == closeatom:
                continue

            # Also ignore if this is a donor/acceptor pair
            if (atom.is_hydrogen and len(atom.bonds) != 0 and \
                atom.bonds[0].hdonor and closeatom.hacceptor):
                continue

            if (closeatom.is_hydrogen and len(closeatom.bonds) != 0 and \
                closeatom.bonds[0].hdonor and atom.hacceptor):
                continue

            dist = util.distance(atom.coords, closeatom.coords)
            other_size = BUMP_HYDROGEN_SIZE if closeatom.is_hydrogen else BUMP_HEAVY_SIZE
            cutoff = atom_size + other_size
            if dist < cutoff:
                nearatoms[closeatom] = cutoff - dist

        return nearatoms

    def set_dihedral_angle(self, residue, anglenum, angle):
        """Rotate a residue about a given angle.

        Uses the quatfit methods to perform the matrix mathematics.

        Parameters
            residue:   The residue to rotate
            anglenum:  The number of the angle to rotate as listed in residue.dihedrals
            angle:     The desired angle.
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
            x = (newcoords[iatom][0] + coordlist[1][0])
            y = (newcoords[iatom][1] + coordlist[1][1])
            z = (newcoords[iatom][2] + coordlist[1][2])
            atom.x = x
            atom.y = y
            atom.z = z
            self.cells.add_cell(atom)

        # Set the new angle
        coordlist = []
        for atomname in atomnames:
            if residue.has_atom(atomname):
                coordlist.append(residue.get_atom(atomname).coords)
            else:
                raise ValueError("Error occurred while trying to debump!")

        dihed = util.dihedral(coordlist[0], coordlist[1], coordlist[2], coordlist[3])
        residue.dihedrals[anglenum] = dihed
