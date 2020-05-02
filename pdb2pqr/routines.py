"""Routines for PDB2PQR

This module contains the protein object used in PDB2PQR and methods used to
correct, analyze, and optimize that protein.

Authors:  Jens Erik Nielsen, Todd Dolinsky, Yong Huang
"""
import copy
import string
import logging
# TODO - replace os with pathlib
import os
import tempfile
from pprint import pformat
from .aa import Amino, PRO, WAT, CYS, LEU, ILE
from .na import Nucleic
from .utilities import distance, getDihedral, shortestPath, subtract
from .utilities import DuplicateFilter
from .quatfit import find_coordinates, qchichange
from .structures import Chain
from . import propka
# TODO - PDB2PKA is still broken
# from .pdb2pka import pka


_LOGGER = logging.getLogger(__name__)
_LOGGER.addFilter(DuplicateFilter())


ANGLE_STEPS = 72
ANGLE_STEP_SIZE = float(360 // ANGLE_STEPS)
ANGLE_TEST_COUNT = 10
EPSILON = 0.0000001
CELL_SIZE = 2
BUMP_DIST = 2.0
BUMP_HDIST = 1.5
BUMP_HYDROGEN_SIZE = 0.5
BUMP_HEAVY_SIZE = 1.0
BONDED_SS_LIMIT = 2.5
PEPTIDE_DIST = 1.7
REPAIR_LIMIT = 10


# TODO -- seems like AAS should go in aa module.
AAS = ["ALA", "ARG", "ASH", "ASN", "ASP", "CYS", "CYM", "GLN", "GLU", "GLH", "GLY", \
       "HIS", "HID", "HIE", "HIP", "HSD", "HSE", "HSP", "ILE", "LEU", "LYS", "LYN", \
       "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "TYM", "VAL"]


# TODO -- seems like NAS should go in aa module.
NAS = ["A", "A5", "A3", "C", "C5", "C3", "G", "G5", "G3", "T", "T5", "T3", "U", \
       "U5", "U3", "RA", "RG", "RC", "RU", "DA", "DG", "DC", "DT"]


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
        if definition != None:
            self.aadef = definition.getAA()
            self.nadef = definition.getNA()

    def apply_name_scheme(self, forcefield):
        """Apply the naming scheme of the given forcefield to the atoms within the protein

        Parameters
            forcefield: The forcefield object (forcefield)
        """
        _LOGGER.info("Applying the naming scheme to the protein...")
        for residue in self.protein.get_residues():
            if isinstance(residue, (Amino, WAT, Nucleic)):
                resname = residue.ffname
            else:
                resname = residue.name

            for atom in residue.get_atoms():
                rname, aname = forcefield.get_names(resname, atom.name)
                if resname not in ['LIG', 'WAT', 'ACE', 'NME'] and rname != None:
                    try:
                        if (residue.is_n_term or residue.is_c_term) and rname != residue.name:
                            rname = residue.name
                    except AttributeError:
                        pass
                if aname != None and rname != None:
                    atom.res_name = rname
                    atom.name = aname

        _LOGGER.debug("Done applying naming scheme.")

    def apply_force_field(self, forcefield):
        """Apply the forcefield to the atoms within the protein

        Parameters
            forcefield: The forcefield object (forcefield)
        Returns
            hitlist:    A list of atoms that were found in the forcefield (list)
            misslist:   A list of atoms that were not found in the forcefield (list)
        """
        _LOGGER.info("Applying the forcefield to the protein...")
        misslist = []
        hitlist = []
        for residue in self.protein.get_residues():
            if isinstance(residue, (Amino, WAT, Nucleic)):
                resname = residue.ffname
            else:
                resname = residue.name

            # Apply the parameters
            for atom in residue.get_atoms():
                atomname = atom.get("name")
                charge, radius = forcefield.get_params(resname, atomname)
                if charge != None and radius != None:
                    atom.set("ffcharge", charge)
                    atom.set("radius", radius)
                    hitlist.append(atom)
                else:
                    misslist.append(atom)

        _LOGGER.debug("Done applying forcefield.")
        return hitlist, misslist

    def update_residue_types(self):
        """Find the type of residue as notated in the Amino Acid definition"""
        _LOGGER.info("Updating residue types... ")
        for chain in self.protein.get_chains():
            for residue in chain.get("residues"):
                # TODO - why are we setting residue types to numeric values?
                name = residue.get("name")
                if name in AAS:
                    residue.set("type", 1)
                elif name == "WAT":
                    residue.set("type", 3)
                elif name in NAS:
                    residue.set("type", 4)
                else: # Residue is a ligand or unknown
                    residue.set("type", 2)
        _LOGGER.debug("Done updating residue types.")

    def update_ss_bridges(self):
        """Check for SS-bridge partners, and if present, set appropriate partners."""
        _LOGGER.info("Updating SS bridges...")
        sg_partners = {}
        for residue in self.protein.get_residues():
            if isinstance(residue, CYS):
                atom = residue.get_atom("SG")
                if atom != None:
                    sg_partners[atom] = []

        for atom in sg_partners:
            for partner in sg_partners:
                if atom == partner or sg_partners[atom] != []:
                    continue
                dist = distance(atom.getCoords(), partner.getCoords())
                if dist < BONDED_SS_LIMIT:
                    sg_partners[atom].append(partner)
                    sg_partners[partner].append(atom)

        for atom in sg_partners:
            res1 = atom.get("residue")
            numpartners = len(sg_partners[atom])
            if numpartners == 1:
                partner = sg_partners[atom][0]
                res2 = partner.get("residue")
                res1.set("ss_bonded", 1)
                res1.set("ss_bonded_partner", partner)
                self.apply_patch("CYX", res1)
                _LOGGER.debug("%s - %s", res1, res2)
            elif numpartners > 1:
                error = "WARNING: %s has multiple potential " % res1
                error += "SS-bridge partners"
                _LOGGER.warn(error)
            elif numpartners == 0:
                _LOGGER.debug("%s is a free cysteine", res1)
        _LOGGER.debug("Done updating SS bridges.")

    def update_internal_bonds(self):
        """Update the internal bonding network using the reference objects in each atom."""
        for residue in self.protein.get_residues():
            if isinstance(residue, (Amino, WAT, Nucleic)):
                for atom in residue.get_atoms():
                    if not atom.hasReference():
                        continue
                    for bond in atom.reference.bonds:
                        if not residue.has_atom(bond):
                            continue
                        bondatom = residue.get_atom(bond)
                        if bondatom not in atom.bonds:
                            atom.addBond(bondatom)

    def update_bonds(self):
        """Update the bonding network of the protein.

        This happens in 3 steps:
            1.  Applying the PEPTIDE patch to all Amino residues so as to add
                reference for the N(i+1) and C(i-1) atoms
            2.  UpdateInternal_bonds for inter-residue linking
            3.  Set the links to the N(i+1) and C(i-1) atoms
        """
        # Apply the peptide patch
        for residue in self.protein.get_residues():
            if isinstance(residue, Amino):
                if residue.is_n_term or residue.is_c_term:
                    continue
                else:
                    self.apply_patch("PEPTIDE", residue)

        # Update all internal bonds
        self.update_internal_bonds()

        # Set the peptide bond pointers
        for chain in self.protein.get_chains():
            for i in range(chain.num_residues() - 1):
                res1 = chain.residues[i]
                res2 = chain.residues[i + 1]
                if not isinstance(res1, Amino) or not isinstance(res2, Amino):
                    continue
                atom1 = res1.get_atom("C")
                atom2 = res2.get_atom("N")

                if atom1 is not None:
                    res2.peptide_c = atom1
                if atom2 is not None:
                    res1.peptide_n = atom2
                if atom1 is None or atom2 is None:
                    continue

                if distance(atom1.getCoords(), atom2.getCoords()) > PEPTIDE_DIST:
                    text = "Gap in backbone detected between %s and %s!" % \
                           (res1, res2)
                    _LOGGER.warn(text)
                    res2.peptide_c = None
                    res1.peptide_n = None

    def apply_patch(self, patchname, residue):
        """Apply a patch to the given residue.

        This is one of the key functions in PDB2PQR.  A similar function appears
        in definitions.py - that version is needed for residue level subtitutions
        so certain protonation states (i.e. CYM, HSE) are detectatble on input.

        This version looks up the particular patch name in the patchmap stored
        in the protein, and then applies the various commands to the reference
        and actual residue structures.

        See the inline comments for a more detailed explanation.

        Parameters
            patchname:  The name of the patch (string)
            residue:    The residue to apply the patch to (residue)
        """
        if patchname not in self.protein.patchmap:
            raise KeyError("Unable to find patch %s!" % patchname)

        _LOGGER.debug('PATCH INFO: %s patched with %s', residue, patchname)

        if patchname == "PEPTIDE":
            newreference = residue.reference
        else:
            newreference = copy.deepcopy(residue.reference)

        patch = self.protein.patchmap[patchname]

        # Add atoms from patch
        for atomname in patch.map:
            newreference.map[atomname] = patch.map[atomname]
            for bond in patch.map[atomname].bonds:
                if bond not in newreference.map:
                    continue
                if atomname not in newreference.map[bond].bonds:
                    newreference.map[bond].bonds.append(atomname)

        # Remove atoms as directed by patch
        for remove in patch.remove:
            if remove in residue.map:
                residue.removeAtom(remove)
            if remove not in newreference.map:
                continue
            removebonds = newreference.map[remove].bonds
            del newreference.map[remove]
            for bond in removebonds:
                index = newreference.map[bond].bonds.index(remove)
                del newreference.map[bond].bonds[index]

        # Add the new dihedrals
        for dihedral in patch.dihedrals:
            newreference.dihedrals.append(dihedral)

        # Point at the new reference
        residue.reference = newreference
        residue.patches.append(patchname)

        # Rename atoms as directed by patch
        for atom in residue.get_atoms():
            if atom.name in patch.altnames:
                residue.renameAtom(atom.name, patch.altnames[atom.name])

        # Replace each atom's reference with the new one
        for atomname in residue.map:
            if newreference.has_atom(atomname):
                atom = residue.get_atom(atomname)
                atom.reference = newreference.map[atomname]

    def set_states(self):
        """Set the state of each residue.

        This is the last step before assigning the forcefield, but is necessary
        so as to distinguish between various protonation states.

        See aa.py for residue-specific functions.
        """
        for residue in self.protein.get_residues():
            if isinstance(residue, (Amino, Nucleic)):
                residue.set_state()

    def assign_termini(self, chain, neutraln=False, neutralc=False):
        """Assign the termini for the given chain by looking at the start and
        end residues.
        """

        if len(chain.residues) == 0:
            text = "Error: chain \"%s\" has 0 residues!" % chain.chain_id
            raise IndexError(text)

        # Set the N-Terminus/ 5' Terminus
        res0 = chain.residues[0]
        if isinstance(res0, Amino):
            res0.set("is_n_term", 1)
            if isinstance(res0, PRO):
                self.apply_patch("NEUTRAL-NTERM", res0)
            elif neutraln:
                self.apply_patch("NEUTRAL-NTERM", res0)
            else:
                self.apply_patch("NTERM", res0)
        elif isinstance(res0, Nucleic):
            res0.set("is5term", 1)
            self.apply_patch("5TERM", res0)

        # Set the C-Terminus/ 3' Terminus
        reslast = chain.residues[-1]
        if isinstance(reslast, Amino):
            reslast.set("is_c_term", 1)
            if neutralc:
                self.apply_patch("NEUTRAL-CTERM", reslast)
            else:
                self.apply_patch("CTERM", reslast)
        elif isinstance(reslast, Nucleic):
            reslast.set("is3term", 1)
            self.apply_patch("3TERM", reslast)
        else:
            for i in range(len(chain.residues)):
                resthis = chain.residues[-1 - i]
                if isinstance(resthis, Amino):
                    resthis.set("is_c_term", 1)
                    if neutralc:
                        self.apply_patch("NEUTRAL-CTERM", resthis)
                    else:
                        self.apply_patch("CTERM", resthis)
                    break
                elif resthis.name in ["NH2", "NME"]:
                    break
                elif isinstance(resthis, Nucleic):
                    resthis.set("is3term", 1)
                    self.apply_patch("3TERM", resthis)
                    break

    def set_termini(self, neutraln=False, neutralc=False):
        """Set the termini for the protein.

        First set all known termini by looking at the ends of the chain. Then
        examine each residue, looking for internal chain breaks.
        """
        # TODO - this function does a lot more than just set termini...
        _LOGGER.info("Setting the termini...")

        # First assign the known termini
        for chain in self.protein.get_chains():
            self.assign_termini(chain, neutraln, neutralc)

        # Now determine if there are any hidden chains
        letters = string.ascii_uppercase + string.ascii_lowercase
        ch_num = 0

        while ch_num < len(self.protein.get_chains()):
            chain = self.protein.chains[ch_num]
            reslist = []
            origlist = []

            # origlist holds the original residue list for the chain
            for residue in chain.get_residues():
                origlist.append(residue)

            for residue in origlist:
                reslist.append(residue)

                # Look for ending termini
                fixflag = 0
                if isinstance(residue, Amino):
                    if (residue.has_atom("OXT") and not residue.is_c_term):
                        fixflag = 1

                elif isinstance(residue, Nucleic):
                    if ((residue.has_atom("H3T") or residue.name.endswith("3"))\
                      and not residue.is3term):
                        fixflag = 1

                if fixflag:
                    # Get an available chain ID
                    chainid = letters[0]
                    id_ = 0
                    id_length = 1
                    while chainid in self.protein.chainmap:
                        id_ += 1
                        if id_ >= len(letters):
                            id_length += 1
                            id_ = 0
                        chainid = letters[id_] * id_length

                    if id_length > 1:
                        message = 'Warning: Reusing chain id: ' + chainid[0] + ''
                        _LOGGER.warn(message)

                    # Make a new chain with these residues
                    newchain = Chain(chainid[0])

                    self.protein.chainmap[chainid] = newchain
                    self.protein.chains.insert(ch_num, newchain)

                    for res in reslist:
                        newchain.addResidue(res)
                        chain.residues.remove(res)
                        res.setChainID(chainid[0])

                    self.assign_termini(chain, neutraln, neutralc)
                    self.assign_termini(newchain, neutraln, neutralc)

                    reslist = []
                    ch_num += 1
            ch_num += 1

        # Update the final chain's chain_id if it is "" unless it's all water
        if "" in self.protein.chainmap:
            notwat = 0
            for res in chain.residues:
                if not isinstance(res, WAT):
                    notwat = 1
                    break

            if notwat == 0:
                _LOGGER.debug("Done setting termini.")
                return

            chain = self.protein.chainmap[""]
            chainid = letters[0]
            id_ = 0
            id_length = 1
            while chainid in self.protein.chainmap:
                id_ += 1
                if id_ >= len(letters):
                    id_length += 1
                    id_ = 0
                chainid = letters[id_] * id_length

            if id_length > 1:
                message = 'Warning: Reusing chain id: ' + chainid[0]
                _LOGGER.warn(message)

            # Use the new chain_id
            self.protein.chainmap[chainid] = chain
            del self.protein.chainmap[""]

            for res in chain.residues:
                res.setChainID(chainid[0])
        _LOGGER.debug("Done setting termini.")

    def find_missing_heavy(self):
        """Repair residues that contain missing heavy (non-Hydrogen) atoms"""
        _LOGGER.info("Checking for missing heavy atoms...")
        misscount = 0
        heavycount = 0
        for residue in self.protein.get_residues():
            if not isinstance(residue, (Amino, Nucleic)):
                continue

            # Check for Missing Heavy Atoms
            for refatomname in residue.reference.map:
                if refatomname.startswith("H"):
                    continue
                if refatomname in ["N+1", "C-1"]:
                    continue
                if refatomname in ["O1P", "O2P"]:
                    if residue.has_atom("OP1") and residue.has_atom("OP2"):
                        continue
                heavycount += 1
                if not residue.has_atom(refatomname):
                    _LOGGER.debug("Missing %s in %s", refatomname, residue)
                    misscount += 1
                    residue.addMissing(refatomname)

            # Check for Extra Atoms
            atomlist = []
            for atom in residue.get("atoms"):
                atomlist.append(atom)

            for atom in atomlist:
                atomname = atom.get("name")
                if atomname in ["OP1", "OP2"] and residue.reference.has_atom("O1P") \
                    and residue.reference.has_atom("O2P"):
                    continue
                if not residue.reference.has_atom(atomname):
                    _LOGGER.debug("Extra atom %s in %s! - ", atomname, residue)
                    residue.removeAtom(atomname)
                    _LOGGER.debug("Deleted this atom.")

        if heavycount == 0:
            raise ValueError(("No heavy atoms found. "
                              "You may also see this message if PDB2PQR does "
                              "not have parameters for any residue in your protein."))

        misspct = 100.0 * float(misscount) / heavycount
        if misspct > REPAIR_LIMIT:
            error = "This PDB file is missing too many (%i out of " % misscount
            error += "%i, %.2f%%) heavy atoms to accurately repair the file.  " % \
                     (heavycount, misspct)
            error += "The current repair limit is set at %i%%. " % REPAIR_LIMIT
            error += "You may also see this message if PDB2PQR does not have "
            error += "parameters for enough residues in your protein."
            raise ValueError(error)
        elif misscount > 0:
            _LOGGER.debug("Missing %i out of %i heavy atoms (%.2f percent) - ",
                          misscount, heavycount, misspct)
            _LOGGER.debug("Will attempt to repair.")
            self.repair_heavy()
        else:
            _LOGGER.debug("No heavy atoms found missing.")
            _LOGGER.debug("Done checking for missing heavy atoms.")

    @staticmethod
    def rebuild_tetrahedral(residue, atomname):
        """Rebuild a tetrahedral hydrogen group.

        This is necessary due to the shortcomings of the quatfit routine - given
        a tetrahedral geometry and two existing hydrogens, the quatfit routines
        have two potential solutions.  This function uses basic tetrahedral
        geometry to fix this issue.

        Parameters
            residue:  The residue in question (residue)
            atomname: The atomname to add (string)
        Returns
            True if successful, False otherwise
        """
        hcount = 0
        nextatomname = None

        atomref = residue.reference.map.get(atomname)
        if atomref is None:
            return False
        bondname = atomref.bonds[0]

        # Return if the bonded atom does not exist
        if not residue.has_atom(bondname):
            return False

        # This group is tetrahedral if bondatom has 4 bonds,
        #  3 of which are hydrogens
        for bond in residue.reference.map[bondname].bonds:
            if bond.startswith("H"):
                hcount += 1
            elif bond != 'C-1' and bond != 'N+1':
                nextatomname = bond

        # Check if this is a tetrahedral group
        if hcount != 3 or nextatomname is None:
            return False

        # Now rebuild according to the tetrahedral geometry
        bondatom = residue.get_atom(bondname)
        nextatom = residue.get_atom(nextatomname)
        numbonds = len(bondatom.bonds)

        if numbonds == 1:

            # Place according to two atoms
            coords = [bondatom.getCoords(), nextatom.getCoords()]
            refcoords = [residue.reference.map[bondname].getCoords(), \
                         residue.reference.map[nextatomname].getCoords()]
            refatomcoords = atomref.getCoords()
            newcoords = find_coordinates(2, coords, refcoords, refatomcoords)
            residue.create_atom(atomname, newcoords)

            # For LEU and ILE residues only: make sure the Hydrogens are in
            # staggered conformation instead of eclipsed.
            if isinstance(residue, LEU):
                hcoords = newcoords
                cbatom = residue.get_atom('CB')
                ang = getDihedral(cbatom.getCoords(), nextatom.getCoords(),
                                  bondatom.getCoords(), hcoords)
                diffangle = 60 - ang
                residue.rotateTetrahedral(nextatom, bondatom, diffangle)

            elif isinstance(residue, ILE):
                hcoords = newcoords
                cg1atom = residue.get_atom('CG1')
                cbatom = residue.get_atom('CB')
                if bondatom.name == 'CD1':
                    ang = getDihedral(cbatom.getCoords(), nextatom.getCoords(),
                                      bondatom.getCoords(), hcoords)
                elif bondatom.name == 'CG2':
                    ang = getDihedral(cg1atom.getCoords(), nextatom.getCoords(),
                                      bondatom.getCoords(), hcoords)
                else:
                    ang = getDihedral(cbatom.getCoords(), nextatom.getCoords(),
                                      bondatom.getCoords(), hcoords)

                diffangle = 60 - ang
                residue.rotateTetrahedral(nextatom, bondatom, diffangle)
            return True

        elif numbonds == 2:

            # Get the single hydrogen coordinates
            hatom = None
            for bond in bondatom.reference.bonds:
                if residue.has_atom(bond) and bond.startswith("H"):
                    hatom = residue.get_atom(bond)
                    break

            # Use the existing hydrogen and rotate about the bond
            residue.rotateTetrahedral(nextatom, bondatom, 120)
            newcoords = hatom.getCoords()
            residue.rotateTetrahedral(nextatom, bondatom, -120)
            residue.create_atom(atomname, newcoords)

            return True

        elif numbonds == 3:

            # Find the one spot the atom can be
            hatoms = []
            for bond in bondatom.reference.bonds:
                if residue.has_atom(bond) and bond.startswith("H"):
                    hatoms.append(residue.get_atom(bond))

            # If this is more than two something is wrong
            if len(hatoms) != 2:
                return False

            # Use the existing hydrogen and rotate about the bond
            residue.rotateTetrahedral(nextatom, bondatom, 120)
            newcoords1 = hatoms[0].getCoords()
            residue.rotateTetrahedral(nextatom, bondatom, 120)
            newcoords2 = hatoms[0].getCoords()
            residue.rotateTetrahedral(nextatom, bondatom, 120)

            # Determine which one hatoms[1] is not in
            if distance(hatoms[1].getCoords(), newcoords1) > 0.1:
                residue.create_atom(atomname, newcoords1)
            else:
                residue.create_atom(atomname, newcoords2)

            return True
        return False

    def add_hydrogens(self):
        """Add the hydrogens to the protein.

        This requires either the rebuild_tetrahedral function for tetrahedral
        geometries or the standard quatfit methods.  These methods use three
        nearby bonds to rebuild the atom; the closer the bonds, the more
        accurate the results.  As such the peptide bonds are used when available.
        """
        count = 0
        _LOGGER.info("Adding hydrogens to the protein...")
        for residue in self.protein.get_residues():
            if not isinstance(residue, (Amino, Nucleic)):
                continue
            for atomname in residue.reference.map:
                if not atomname.startswith("H"):
                    continue
                if residue.has_atom(atomname):
                    continue
                if isinstance(residue, CYS) and residue.ss_bonded and atomname == "HG":
                    continue

                # If this hydrogen is part of a tetrahedral group,
                # follow a different codepath
                if Routines.rebuild_tetrahedral(residue, atomname):
                    count += 1
                    continue

                # Otherwise use the standard quatfit methods
                coords = []
                refcoords = []

                refatomcoords = residue.reference.map[atomname].getCoords()
                bondlist = residue.reference.get_nearest_bonds(atomname)

                for bond in bondlist:
                    if bond == "N+1":
                        atom = residue.peptide_n
                    elif bond == "C-1":
                        atom = residue.peptide_c
                    else:
                        atom = residue.get_atom(bond)

                    if atom is None:
                        continue

                    # Get coordinates, reference coordinates
                    coords.append(atom.getCoords())
                    refcoords.append(residue.reference.map[bond].getCoords())

                    # Exit if we have enough atoms
                    if len(coords) == 3:
                        break

                if len(coords) == 3:
                    newcoords = find_coordinates(3, coords, refcoords, refatomcoords)
                    residue.create_atom(atomname, newcoords)
                    count += 1
                else:
                    _LOGGER.warn("Couldn't rebuild %s in %s!", atomname, residue)
        _LOGGER.debug(" Added %i hydrogen atoms.", count)

    def remove_hydrogens(self):
        """Remove hydrogens from the protein."""
        _LOGGER.info("Stripping hydrogens from the protein...")
        for residue in self.protein.get_residues():
            if not isinstance(residue, (Amino, Nucleic)):
                continue
            for atom in residue.atoms[:]:
                if atom.isHydrogen():
                    residue.removeAtom(atom.name)

    def repair_heavy(self):
        """Repair all heavy atoms.

        Unfortunately the first time we get to an atom we might not be able to
        rebuild it - it might depend on other atoms to be rebuild first (think
        side chains).  As such a 'seenmap' is used to keep track of what we've
        already seen and subsequent attempts to rebuild the atom.
        """
        _LOGGER.info("Rebuilding missing heavy atoms...")
        for residue in self.protein.get_residues():
            if not isinstance(residue, (Amino, Nucleic)):
                continue
            missing = residue.get("missing")
            if missing == []:
                continue

            # Initialize some variables
            seenmap = {}
            nummissing = len(missing)

            while len(missing) > 0:
                coords = []
                refcoords = []

                atomname = missing.pop(0)
                refatomcoords = residue.reference.map[atomname].getCoords()
                bondlist = residue.reference.get_nearest_bonds(atomname)

                for bond in bondlist:
                    if bond == "N+1":
                        atom = residue.peptide_n
                    elif bond == "C-1":
                        atom = residue.peptide_c
                    else: atom = residue.get_atom(bond)

                    if atom is None:
                        continue

                    # Get coordinates, reference coordinates
                    coords.append(atom.getCoords())
                    refcoords.append(residue.reference.map[bond].getCoords())

                    # Exit if we have enough atoms
                    if len(coords) == 3:
                        break

                # We might need other atoms to be rebuilt first
                if len(coords) < 3:
                    try:
                        seenmap[atomname] += 1
                    except KeyError:
                        seenmap[atomname] = 1

                    missing.append(atomname)
                    if seenmap[atomname] > nummissing:
                        text = "Too few atoms present to reconstruct or cap "
                        text += "residue %s in structure! " % (residue)
                        text += "This error is generally caused by missing backbone "
                        text += "atoms in this protein; "
                        text += "you must use an external program to complete gaps "
                        text += "in the protein backbone. "
                        text += "Heavy atoms missing from %s: " % (residue)
                        text += ' '.join(missing)
                        raise ValueError(text)

                else: # Rebuild the atom
                    newcoords = find_coordinates(3, coords, refcoords, refatomcoords)
                    residue.create_atom(atomname, newcoords)
                    _LOGGER.debug("Added %s to %s at coordinates", atomname, residue)
                    _LOGGER.debug(" %.3f %.3f %.3f", newcoords[0], newcoords[1], newcoords[2])
        _LOGGER.debug("Done rebuilding missing atoms.")

    def set_reference_distance(self):
        """Set the distance to the CA atom in the residue.

        This is necessary for determining which atoms are allowed to move during
        rotations.  Uses the shortestPath algorithm found in utilities.py.
        """
        for residue in self.protein.get_residues():
            if not isinstance(residue, Amino):
                continue

            # Initialize some variables
            map_ = {}
            caatom = residue.get_atom("CA")

            if caatom is None:
                text = "Cannot set references to %s without CA atom!"
                raise ValueError(text)

            # Set up the linked map
            for atom in residue.get_atoms():
                map_[atom] = atom.bonds

            # Run the algorithm
            for atom in residue.get_atoms():
                if atom.isBackbone():
                    atom.refdistance = -1
                elif residue.is_c_term and atom.name == "HO":
                    atom.refdistance = 3
                elif residue.is_n_term and (atom.name == "H3" or atom.name == "H2"):
                    atom.refdistance = 2
                else:
                    atom.refdistance = len(shortestPath(map_, atom, caatom)) - 1

    def get_bump_score(self, residue):
        """Get an bump score for the current structure"""

        # Do some setup
        self.cells = Cells(CELL_SIZE)
        self.cells.assign_cells(self.protein)

        self.calculate_dihedral_angles()
        self.set_donors_acceptors()
        self.update_internal_bonds()
        self.set_reference_distance()
        bumpscore = 0.0
        if not isinstance(residue, Amino):
            return 0.0

        # Initialize variables
        for atom in residue.get_atoms():
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
        atom_size = BUMP_HYDROGEN_SIZE if atom.isHydrogen() else BUMP_HEAVY_SIZE

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
            if atom.isHydrogen() and len(atom.bonds) != 0 and atom.bonds[0].hdonor \
               and closeatom.hacceptor:
                continue

            if closeatom.isHydrogen() and len(closeatom.bonds) != 0 and closeatom.bonds[0].hdonor \
                   and atom.hacceptor:
                continue

            dist = distance(atom.getCoords(), closeatom.getCoords())
            other_size = BUMP_HYDROGEN_SIZE if closeatom.isHydrogen() else BUMP_HEAVY_SIZE
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

        self.calculate_dihedral_angles()
        self.set_donors_acceptors()
        self.update_internal_bonds()
        self.set_reference_distance()

        # Determine which residues to debump
        for residue in self.protein.get_residues():
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
                _LOGGER.warn(text)

        _LOGGER.debug("Done checking if we must debump any residues.")

    def find_residue_conflicts(self, residue, write_conflict_info=False):
        """Find conflicts between residues."""
        conflict_names = []
        for atom in residue.get_atoms():
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
        moveablenames = self.get_moveable_names(residue, pivot)
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
        for _ in range(ANGLE_TEST_COUNT):
            anglenum = self.pick_dihedral_angle(residue, curr_conflict_names, anglenum)
            if anglenum == -1:
                return False

            _LOGGER.debug("Using dihedral angle number %i to debump the residue.", anglenum)
            bestscore = self.score_dihedral_angle(residue, anglenum)
            found_improved = False
            bestangle = orig_angle = residue.dihedrals[anglenum]

            # Skip the first angle as it's already known.
            for i in range(1, ANGLE_STEPS):
                newangle = orig_angle + (ANGLE_STEP_SIZE * i)
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
                    if diff > EPSILON:
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

    def calculate_dihedral_angles(self):
        """Calculate the dihedral angle for every residue within the protein"""
        for residue in self.protein.get_residues():
            if not isinstance(residue, Amino):
                continue
            residue.dihedrals = []

            refangles = residue.reference.dihedrals
            for dihed in refangles:
                coords = []
                atoms = dihed.split()
                for i in range(4):
                    atomname = atoms[i]
                    if residue.has_atom(atomname):
                        coords.append(residue.get_atom(atomname).getCoords())

                if len(coords) == 4:
                    angle = getDihedral(coords[0], coords[1], coords[2], coords[3])
                else:
                    angle = None

                residue.add_dihedral_angle(angle)

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
            if atom.isHydrogen() and atom.bonds[0].hdonor and closeatom.hacceptor:
                continue
            if closeatom.isHydrogen() and closeatom.bonds[0].hdonor and atom.hacceptor:
                continue

            dist = distance(atom.getCoords(), closeatom.getCoords())

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
        atom_size = BUMP_HYDROGEN_SIZE if atom.isHydrogen() else BUMP_HEAVY_SIZE

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
            if (atom.isHydrogen() and len(atom.bonds) != 0 and \
                atom.bonds[0].hdonor and closeatom.hacceptor):
                continue

            if (closeatom.isHydrogen() and len(closeatom.bonds) != 0 and \
                closeatom.bonds[0].hdonor and atom.hacceptor):
                continue

            dist = distance(atom.getCoords(), closeatom.getCoords())
            other_size = BUMP_HYDROGEN_SIZE if closeatom.isHydrogen() else BUMP_HEAVY_SIZE
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

            moveablenames = self.get_moveable_names(residue, pivot)

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
                coordlist.append(residue.get_atom(atomname).getCoords())
            else:
                raise ValueError("Error occurred while trying to debump!")

        initcoords = subtract(coordlist[2], coordlist[1])

        moveablenames = self.get_moveable_names(residue, pivot)

        for name in moveablenames:
            atom = residue.get_atom(name)
            movecoords.append(subtract(atom.getCoords(), coordlist[1]))

        newcoords = qchichange(initcoords, movecoords, diff)

        for iatom, atom_name in enumerate(moveablenames):
            atom = residue.get_atom(atom_name)
            self.cells.remove_cell(atom)
            x = (newcoords[iatom][0] + coordlist[1][0])
            y = (newcoords[iatom][1] + coordlist[1][1])
            z = (newcoords[iatom][2] + coordlist[1][2])
            atom.set("x", x)
            atom.set("y", y)
            atom.set("z", z)
            self.cells.add_cell(atom)

        # Set the new angle
        coordlist = []
        for atomname in atomnames:
            if residue.has_atom(atomname):
                coordlist.append(residue.get_atom(atomname).getCoords())
            else:
                raise ValueError("Error occurred while trying to debump!")

        dihed = getDihedral(coordlist[0], coordlist[1], coordlist[2], coordlist[3])
        residue.dihedrals[anglenum] = dihed

    @classmethod
    def get_moveable_names(cls, residue, pivot):
        """Return all atomnames that are further away than the pivot atom.

        Parameters
            residue:  The residue to use
            pivot:    The pivot atomname
        """
        movenames = []
        refdist = residue.get_atom(pivot).refdistance
        for atom in residue.get_atoms():
            if atom.refdistance > refdist:
                movenames.append(atom.name)
        return movenames

    def set_donors_acceptors(self):
        """Set the donors and acceptors within the protein"""
        for residue in self.protein.get_residues():
            residue.set_donors_acceptors()

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
        # _LOGGER.warn(pdb2pka_warnings)

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

        ph = pka_options.pH
        _LOGGER.info("Running propka 3.1 at pH %.2f... ", ph)

        # Initialize some variables
        pkadic = {}

        # Reorder the atoms in each residue to start with N - TONI is this necessary?
        for residue in self.protein.get_residues():
            residue.reorder()

        # TONI Make a string with all non-hydrogen atoms. Previously it was removing the "element"
        # column and hydrogens. This does not seem to be necessary in propKa 3.1 .
        with tempfile.NamedTemporaryFile(mode="w+", suffix=".pdb") as h_free_file:
            for atom in self.protein.get_atoms():
                if not atom.isHydrogen():
                    atomtxt = atom.getPDBString()
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

        for residue in self.protein.get_residues():
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
                            _LOGGER.warn(warn)
                        else:
                            self.apply_patch("NEUTRAL-NTERM", residue)

            if residue.is_c_term:
                key = "C- %i %s" % (resnum, chain_id)
                key = key.strip()
                if key in pkadic:
                    value = pkadic[key]
                    del pkadic[key]
                    if ph < value:
                        if force_field in ["amber", "charmm", "tyl06", "peoepb", "swanson"]:
                            warn = ("C-terminal %s" % key, "neutral")
                            _LOGGER.warn(warn)
                        else:
                            self.apply_patch("NEUTRAL-CTERM", residue)

            key = "%s %i %s" % (resname, resnum, chain_id)
            key = key.strip()
            if key in pkadic:
                value = pkadic[key]
                del pkadic[key]
                if resname == "ARG" and ph >= value:
                    if force_field == "parse":
                        self.apply_patch("AR0", residue)
                        _LOGGER.warn(("Neutral arginines are very rare. Please "
                                      "double-check your setup."))
                    else:
                        warn = (key, "neutral")
                        _LOGGER.warn(warn)
                elif resname == "ASP" and ph < value:
                    if residue.is_c_term and force_field in ["amber", "tyl06", "swanson"]:
                        warn = (key, "Protonated at C-Terminal")
                        _LOGGER.warn(warn)
                    elif residue.is_n_term and force_field in ["amber", "tyl06", "swanson"]:
                        warn = (key, "Protonated at N-Terminal")
                        _LOGGER.warn(warn)
                    else:
                        self.apply_patch("ASH", residue)
                elif resname == "CYS" and ph >= value:
                    if force_field == "charmm":
                        warn = (key, "negative")
                        _LOGGER.warn(warn)
                    else:
                        self.apply_patch("CYM", residue)
                elif resname == "GLU" and ph < value:
                    if residue.is_c_term and force_field in ["amber", "tyl06", "swanson"]:
                        warn = (key, "Protonated at C-Terminal")
                        _LOGGER.warn(warn)
                    elif residue.is_n_term and force_field in ["amber", "tyl06", "swanson"]:
                        warn = (key, "Protonated at N-Terminal")
                        _LOGGER.warn(warn)
                    else:
                        self.apply_patch("GLH", residue)
                elif resname == "HIS" and ph < value:
                    self.apply_patch("HIP", residue)
                elif resname == "LYS" and ph >= value:
                    if force_field == "charmm":
                        warn = (key, "neutral")
                        _LOGGER.warn(warn)
                    elif force_field in ["amber", "tyl06", "swanson"] and residue.get("is_c_term"):
                        warn = (key, "neutral at C-Terminal")
                        _LOGGER.warn(warn)
                    elif force_field == "tyl06" and residue.get("is_n_term"):
                        warn = (key, "neutral at N-Terminal")
                        _LOGGER.warn(warn)
                    else:
                        self.apply_patch("LYN", residue)
                elif resname == "TYR" and ph >= value:
                    if force_field in ["charmm", "amber", "tyl06", "peoepb", "swanson"]:
                        warn = (key, "negative")
                        _LOGGER.warn(warn)
                    else:
                        self.apply_patch("TYM", residue)

        if len(pkadic) > 0:
            warn = ("PDB2PQR could not identify the following residues and residue "
                    "numbers as returned by PROPKA or PDB2PKA")
            _LOGGER.warn(warn)
            for item in pkadic:
                text = "             %s" % item
                _LOGGER.warn(text)

    def hold_residues(self, hlist):
        """Set the stateboolean dictionary to residues in hlist."""
        if not hlist:
            return

        for residue in self.protein.get_residues():
            reskey = (residue.res_seq, residue.chain_id, residue.ins_code)
            if reskey in hlist:
                hlist.remove(reskey)
                if isinstance(residue, Amino):
                    residue.stateboolean = {'FIXEDSTATE': False}
                    _LOGGER.debug("Setting residue {:s} as fixed.".format(str(residue)))
                else:
                    err = "Matched residue {:s} but not subclass of Amino."
                    _LOGGER.warn(err.format(str(residue)))

        if len(hlist) > 0:
            err = "The following fixed residues were not matched (possible internal error): {:s}."
            _LOGGER.warn(err.format(str(hlist)))


class Cells(object):
    """The cells object provides a better way to search for nearby atoms.

    A pure all versus all search is O(n^2) - for every atom, every other atom
    must be searched.  This is rather inefficient, especially for large proteins
    where cells may be tens of angstroms apart.  The cell class breaks down the
    xyz protein space into several 3-D cells of desired size - then by simply
    examining atoms that fall into the adjacent cells one can quickly find nearby
    cells.

    NOTE:  Ideally this should be somehow separated from the routines
            object...
    """

    def __init__(self, cellsize):
        """Initialize the cells.

        Parameters
            cellsize:  The size of each cell (int)
        """
        self.cellmap = {}
        self.cellsize = cellsize

    def assign_cells(self, protein):
        """Place each atom in a virtual cell for easy neighbor comparison."""
        for atom in protein.get_atoms():
            atom.cell = None
            self.add_cell(atom)

    def add_cell(self, atom):
        """Add an atom to the cell

        Parameters
            atom:  The atom to add (atom)
        """
        size = self.cellsize

        x = atom.get("x")
        if x < 0:
            x = (int(x) - 1) // size * size
        else:
            x = int(x) // size * size

        y = atom.get("y")
        if y < 0:
            y = (int(y) - 1) // size * size
        else:
            y = int(y) // size * size

        z = atom.get("z")
        if z < 0:
            z = (int(z) - 1) // size * size
        else:
            z = int(z) // size * size

        key = (x, y, z)
        try:
            self.cellmap[key].append(atom)
        except KeyError:
            self.cellmap[key] = [atom]
        atom.set("cell", key)

    def remove_cell(self, atom):
        """Remove the atom from a cell

        Parameters
            atom:   The atom to add (atom)
        """
        oldcell = atom.get("cell")
        if oldcell is None:
            return
        atom.set("cell", None)
        self.cellmap[oldcell].remove(atom)

    def get_near_cells(self, atom):
        """Find all atoms in bordering cells to an atom

        Parameters
            atom:  The atom to use (atom)
        Returns
            closeatoms:  A list of nearby atoms (list)
        """
        size = self.cellsize
        closeatoms = []
        cell = atom.get("cell")
        if cell is None:
            return closeatoms
        else:
            x = cell[0]
            y = cell[1]
            z = cell[2]
            for i in range(-1 * size, 2 * size, size):
                for j in range(-1 * size, 2 * size, size):
                    for k in range(-1 * size, 2 * size, size):
                        newkey = (x + i, y + j, z + k)
                        try:
                            newatoms = self.cellmap[newkey]
                            for atom2 in newatoms:
                                if atom == atom2:
                                    continue
                                closeatoms.append(atom2)
                        except KeyError:
                            pass

            return closeatoms
