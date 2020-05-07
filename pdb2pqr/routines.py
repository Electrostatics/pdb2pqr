"""Routines for PDB2PQR

This module contains the protein object used in PDB2PQR and methods used to
correct, analyze, and optimize that protein.

Authors:  Jens Erik Nielsen, Todd Dolinsky, Yong Huang
"""
import logging
# TODO - replace os with pathlib
import os
import tempfile
from pprint import pformat
from .aa import Amino, PRO, WAT, CYS, LEU, ILE
from .na import Nucleic
from .utilities import distance, dihedral, shortest_path, subtract
from .io import DuplicateFilter
from .quatfit import find_coordinates, qchichange
from .cells import Cells
from .config import DEBUMP_ANGLE_STEP_SIZE, DEBUMP_ANGLE_STEPS, DEBUMP_ANGLE_TEST_COUNT
from .config import SMALL_NUMBER, CELL_SIZE, BUMP_HYDROGEN_SIZE, BUMP_HEAVY_SIZE
from .config import BONDED_SS_LIMIT, PEPTIDE_DIST, REPAIR_LIMIT


_LOGGER = logging.getLogger(__name__)
_LOGGER.addFilter(DuplicateFilter())


def drop_water(pdblist):
    """Drop waters from a list of PDB records.

    TODO - this module is already too long but this function fits better here.
    Other possible place would be utilities.

    Args:
        pdb_list:  list of PDB records as returned by io.get_molecule
    Returns:
        new list of PDB records with waters removed.
    """
    pdblist_new = []
    for record in pdblist:
        record_type = record.record_type()
        if record_type in ["HETATM", "ATOM", "SIGATM", "SEQADV"]:
            if record.res_name in WAT.water_residue_names:
                continue
        pdblist_new.append(record)
    return pdblist_new


class Routines(object):
    """Grab bag of random stuff that apparently didn't fit elsewhere.

    TODO - needs to be susbtantially refactored in to multiple classes with clear
    responsibilities.
    """
    def __init__(self, protein, definition=None):
        """Initialize the Routines class.

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
        self.cells = Cells(CELL_SIZE)
        self.cells.assign_cells(self.protein)

        self.protein.calculate_dihedral_angles()
        self.set_donors_acceptors()
        self.protein.update_internal_bonds()
        self.protein.set_reference_distance()
        bumpscore = 0.0
        if not isinstance(residue, Amino):
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

            if not isinstance(closeresidue, Amino):
                continue
            if isinstance(residue, CYS):
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

            dist = distance(atom.coords, closeatom.coords)
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
        self.cells = Cells(CELL_SIZE)
        self.cells.assign_cells(self.protein)

        self.protein.calculate_dihedral_angles()
        self.protein.set_donors_acceptors()
        self.protein.update_internal_bonds()
        self.protein.set_reference_distance()

        # Determine which residues to debump
        for residue in self.protein.residues:
            if not isinstance(residue, Amino):
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
            anglenum = self.pick_dihedral_angle(residue, curr_conflict_names, anglenum)
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
            if not isinstance(closeresidue, (Amino, WAT)):
                continue
            if isinstance(residue, CYS):
                if residue.ss_bonded_partner == closeatom:
                    continue

            # Also ignore if this is a donor/acceptor pair
            if atom.is_hydrogen and atom.bonds[0].hdonor and closeatom.hacceptor:
                continue
            if closeatom.is_hydrogen and closeatom.bonds[0].hdonor and atom.hacceptor:
                continue

            dist = distance(atom.coords, closeatom.coords)

            if isinstance(closeresidue, WAT):
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

            if not isinstance(closeresidue, (Amino, WAT)):
                continue
            if isinstance(residue, CYS) and residue.ss_bonded_partner == closeatom:
                continue

            # Also ignore if this is a donor/acceptor pair
            if (atom.is_hydrogen and len(atom.bonds) != 0 and \
                atom.bonds[0].hdonor and closeatom.hacceptor):
                continue

            if (closeatom.is_hydrogen and len(closeatom.bonds) != 0 and \
                closeatom.bonds[0].hdonor and atom.hacceptor):
                continue

            dist = distance(atom.coords, closeatom.coords)
            other_size = BUMP_HYDROGEN_SIZE if closeatom.is_hydrogen else BUMP_HEAVY_SIZE
            cutoff = atom_size + other_size
            if dist < cutoff:
                nearatoms[closeatom] = cutoff - dist

        return nearatoms

    def pick_dihedral_angle(self, residue, conflict_names, oldnum=None):
        """Choose an angle number to use in debumping

        Algorithm
            Instead of simply picking a random chiangle, this function
            uses a more intelligent method to improve efficiency.
            The algorithm uses the names of the conflicting atoms
            within the residue to determine which angle number
            has the best chance of fixing the problem(s). The method
            also insures that the same chiangle will not be run twice
            in a row.

        Parameters
            residue:  The residue that is being debumped (Residue)
            conflict_names: A list of atom names that are currently conflicts (list)
            oldnum:  The old dihedral angle number (int)
        Returns
            bestnum:  The new dihedral angle number (int)
        """
        bestnum = -1
        best = 0

        ilist = list(range(len(residue.dihedrals)))
        #Make sure our testing is done round robin.
        if oldnum is not None and oldnum >= 0 and len(ilist) > 0:
            del ilist[oldnum]
            test_dihedral_indices = ilist[oldnum:] + ilist[:oldnum]
        else:
            test_dihedral_indices = ilist

        for i in test_dihedral_indices:
            if i == oldnum:
                continue
            if residue.dihedrals[i] is None:
                continue

            score = 0
            atomnames = residue.reference.dihedrals[i].split()
            pivot = atomnames[2]

            moveablenames = residue.get_moveable_names(pivot)

            # If this pivot only moves the conflict atoms, pick it
            if conflict_names == moveablenames:
                return i

            # Otherwise find the pivot with the most matches
            for name in conflict_names:
                if name in moveablenames:
                    score += 1
                    if score > best:
                        best = score
                        bestnum = i

        # Return the best angle.  If none were found, return -1.
        return bestnum

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

        initcoords = subtract(coordlist[2], coordlist[1])

        moveablenames = residue.get_moveable_names(pivot)

        for name in moveablenames:
            atom = residue.get_atom(name)
            movecoords.append(subtract(atom.coords, coordlist[1]))

        newcoords = qchichange(initcoords, movecoords, diff)

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

        dihed = dihedral(coordlist[0], coordlist[1], coordlist[2], coordlist[3])
        residue.dihedrals[anglenum] = dihed

    def run_pdb2pka(self, ph, force_field, pdb_list, ligand, pdb2pka_params):
        """Run PDB2PKA"""
        # TODO - we are not ready to deal with PDB2PKA yet
        raise NotImplementedError("TODO - fix and re-enable PDB2PKA")
        # if force_field.lower() != 'parse':
        #     PDB2PKAError('PDB2PKA can only be run with the PARSE force field.')

        # _LOGGER.info("Running PDB2PKA and applying at pH %.2f... ", ph)
        # init_params = pdb2pka_params.copy()
        # init_params.pop('pairene')
        # init_params.pop('clean_output')
        # results = pka.pre_init(original_pdb_list=pdb_list, ff=force_field, ligand=ligand,
        #                        **init_params)
        # TODO - this is a messed-up variable unpacking:
        # output_dir, protein, routines, forcefield, apbs_setup, \
        #     ligand_titratable_groups, maps, sd = results
        # mypkaRoutines = pka_routines.pKaRoutines(protein, routines, forcefield,
        #                                          apbs_setup, output_dir, maps, sd,
        #                                          restart=pdb2pka_params.get('clean_output'),
        #                                          pairene=pdb2pka_params.get('pairene'))

        # _LOGGER.info('Doing full pKa calculation')
        # mypkaRoutines.runpKa()
        # pdb2pka_warnings = mypkaRoutines.warnings[:]
        # _LOGGER.warning(pdb2pka_warnings)

        # residue_ph = {}
        # for pka_residue_tuple, calc_ph in mypkaRoutines.ph_at_0_5.items():
        #     tit_type, chain_id, number_str = pka_residue_tuple
        #     if tit_type == 'NTR':
        #         tit_type = 'N+'
        #     elif tit_type == 'CTR':
        #         tit_type = 'C-'

        #     key = ' '.join([tit_type, number_str, chain_id])
        #     residue_ph[key] = calc_ph
        # pformat(residue_ph)
        # self.apply_pka_values(ff, ph, residue_ph)
        # _LOGGER.debug('Finished running PDB2PKA.')

    @classmethod
    def delete_propka_input(cls, file_name):
        """Delete input files created by PROPKA."""
        _, file_ = os.path.split(file_name)
        file_ = file_.replace('.pdb', '.propka_input')
        os.remove(file_)

    def run_propka_31(self, pka_options):
        """Run PROPKA 3.1 on the current protein, setting protonation states to
        the correct values. pH is set in pka_options

        Parameters
            pka_options: Options for propKa 3.1, including pH

        Returns
            pka_molecule: pKa's internal molecule object (including pKa's, etc)
            not_found:    dict of residues found in pka_molecule but not in PDB2PQR (with pKa)
        """
        # See https://github.com/jensengroup/propka-3.1/blob/master/scripts/propka31.py

        ph = pka_options.ph
        _LOGGER.info("Running propka 3.1 at pH %.2f... ", ph)

        # Initialize some variables
        pkadic = {}

        # Reorder the atoms in each residue to start with N - TONI is this necessary?
        for residue in self.protein.residues:
            residue.reorder()

        # TONI Make a string with all non-hydrogen atoms. Previously it was removing the "element"
        # column and hydrogens. This does not seem to be necessary in propKa 3.1 .
        with tempfile.NamedTemporaryFile(mode="w+", suffix=".pdb") as h_free_file:
            for atom in self.protein.atoms:
                if not atom.is_hydrogen:
                    atomtxt = atom.get_pdb_string()
                    h_free_file.write(atomtxt + '\n')

            # Run PropKa 3.1 -------------
            # Creating protein object. Annoyingly, at this stage propka generates a
            # *.propka_input file in PWD and does not delete it (irrespective of the original
            # .pdb location)
            pka_molecule = propka.molecular_container.Molecular_container(h_free_file.name,
                                                                          pka_options)

        # calculating pKa values for ionizable residues -
        pka_molecule.calculate_pka()

        ##  pka_molecule.write_pka()
        for grp in pka_molecule.conformations['AVR'].groups:
            key = str.strip('%s %s %s' % (grp.residue_type, grp.atom.resNumb, grp.atom.chain_id))
            pkadic[key] = grp.pka_value

        self.protein.pka_protein = pka_molecule
        return pkadic

    def run_propka(self, ph, force_field, options, version=30):
        """Run PROPKA on the current protein, setting protonation states to the correct values

        Parameters
            ph:  The desired pH of the system
            force_field:  The forcefield name to be used
            outname: The name of the PQR outfile
            options: Options to propka
            version: may be 30 or 31 (uses external propka 3.1)
        """
        _LOGGER.info("Running PROPKA v%d and applying at pH %.2f... ", version, ph)
        pkadic = self.run_propka_31(options)

        if len(pkadic) == 0:
            raise ValueError("PROPKA returned empty results!")

        # Now apply each pka to the appropriate residue
        self.apply_pka_values(force_field, ph, pkadic)
        _LOGGER.debug("Done running PROPKA")

    def apply_pka_values(self, force_field, ph, pkadic):
        """Apply calculated pKa values to assign titration states."""
        _LOGGER.info('Applying pKa values at a pH of %.2f:', ph)
        formatted_pkadict = pformat(pkadic)
        _LOGGER.debug("%s", formatted_pkadict)

        for residue in self.protein.residues:
            if not isinstance(residue, Amino):
                continue
            resname = residue.name
            resnum = residue.res_seq
            chain_id = residue.chain_id

            if residue.is_n_term:
                key = "N+ %i %s" % (resnum, chain_id)
                key = key.strip()
                if key in pkadic:
                    value = pkadic[key]
                    del pkadic[key]
                    if ph >= value:
                        if force_field in ["amber", "charmm", "tyl06", "peoepb", "swanson"]:
                            warn = ("N-terminal %s" % key, "neutral")
                            _LOGGER.warning(warn)
                        else:
                            self.protein.apply_patch("NEUTRAL-NTERM", residue)

            if residue.is_c_term:
                key = "C- %i %s" % (resnum, chain_id)
                key = key.strip()
                if key in pkadic:
                    value = pkadic[key]
                    del pkadic[key]
                    if ph < value:
                        if force_field in ["amber", "charmm", "tyl06", "peoepb", "swanson"]:
                            warn = ("C-terminal %s" % key, "neutral")
                            _LOGGER.warning(warn)
                        else:
                            self.protein.apply_patch("NEUTRAL-CTERM", residue)

            key = "%s %i %s" % (resname, resnum, chain_id)
            key = key.strip()
            if key in pkadic:
                value = pkadic[key]
                del pkadic[key]
                if resname == "ARG" and ph >= value:
                    if force_field == "parse":
                        self.protein.apply_patch("AR0", residue)
                        _LOGGER.warning(("Neutral arginines are very rare. Please "
                                         "double-check your setup."))
                    else:
                        warn = (key, "neutral")
                        _LOGGER.warning(warn)
                elif resname == "ASP" and ph < value:
                    if residue.is_c_term and force_field in ["amber", "tyl06", "swanson"]:
                        warn = (key, "Protonated at C-Terminal")
                        _LOGGER.warning(warn)
                    elif residue.is_n_term and force_field in ["amber", "tyl06", "swanson"]:
                        warn = (key, "Protonated at N-Terminal")
                        _LOGGER.warning(warn)
                    else:
                        self.protein.apply_patch("ASH", residue)
                elif resname == "CYS" and ph >= value:
                    if force_field == "charmm":
                        warn = (key, "negative")
                        _LOGGER.warning(warn)
                    else:
                        self.protein.apply_patch("CYM", residue)
                elif resname == "GLU" and ph < value:
                    if residue.is_c_term and force_field in ["amber", "tyl06", "swanson"]:
                        warn = (key, "Protonated at C-Terminal")
                        _LOGGER.warning(warn)
                    elif residue.is_n_term and force_field in ["amber", "tyl06", "swanson"]:
                        warn = (key, "Protonated at N-Terminal")
                        _LOGGER.warning(warn)
                    else:
                        self.protein.apply_patch("GLH", residue)
                elif resname == "HIS" and ph < value:
                    self.protein.apply_patch("HIP", residue)
                elif resname == "LYS" and ph >= value:
                    if force_field == "charmm":
                        warn = (key, "neutral")
                        _LOGGER.warning(warn)
                    elif force_field in ["amber", "tyl06", "swanson"] and residue.is_c_term:
                        warn = (key, "neutral at C-Terminal")
                        _LOGGER.warning(warn)
                    elif force_field == "tyl06" and residue.is_n_term:
                        warn = (key, "neutral at N-Terminal")
                        _LOGGER.warning(warn)
                    else:
                        self.protein.apply_patch("LYN", residue)
                elif resname == "TYR" and ph >= value:
                    if force_field in ["charmm", "amber", "tyl06", "peoepb", "swanson"]:
                        warn = (key, "negative")
                        _LOGGER.warning(warn)
                    else:
                        self.protein.apply_patch("TYM", residue)

        if len(pkadic) > 0:
            warn = ("PDB2PQR could not identify the following residues and residue "
                    "numbers as returned by PROPKA or PDB2PKA")
            _LOGGER.warning(warn)
            for item in pkadic:
                text = "             %s" % item
                _LOGGER.warning(text)

    def hold_residues(self, hlist):
        """Set the stateboolean dictionary to residues in hlist."""
        if not hlist:
            return

        for residue in self.protein.residues:
            reskey = (residue.res_seq, residue.chain_id, residue.ins_code)
            if reskey in hlist:
                hlist.remove(reskey)
                if isinstance(residue, Amino):
                    residue.stateboolean = {'FIXEDSTATE': False}
                    _LOGGER.debug("Setting residue {:s} as fixed.".format(str(residue)))
                else:
                    err = "Matched residue {:s} but not subclass of Amino."
                    _LOGGER.warning(err.format(str(residue)))

        if len(hlist) > 0:
            err = "The following fixed residues were not matched (possible internal error): {:s}."
            _LOGGER.warning(err.format(str(hlist)))


