"""Routines for PDB2PQR

This module contains the protein object used in PDB2PQR and methods used to
correct, analyze, and optimize that protein.

Authors:  Jens Erik Nielsen, Todd Dolinsky, Yong Huang
"""

from __future__ import division
import math
import copy
# TODO - need to fix all of the import * statements
from .pdb import *
from .utilities import *
from .quatfit import *
from .forcefield import *
from .structures import *
from .protein import *
from .definitions import *
from io import StringIO
from .errors import PDBInputError, PDBInternalError, PDB2PKAError
from pprint import pformat
import logging


_LOGGER = logging.getLogger(__name__)


CELL_SIZE = 2
BUMP_DIST = 2.0
BUMP_HDIST = 1.5
BUMP_HYDROGEN_SIZE = 0.5
BUMP_HEAVY_SIZE = 1.0
BONDED_SS_LIMIT = 2.5
PEPTIDE_DIST = 1.7
REPAIR_LIMIT = 10
AAS = ["ALA", "ARG", "ASH", "ASN", "ASP", "CYS", "CYM", "GLN", "GLU", "GLH", "GLY", \
       "HIS", "HID", "HIE", "HIP", "HSD", "HSE", "HSP", "ILE", "LEU", "LYS", "LYN", \
       "MET", "PHE", "PRO", "SER", "THR", "TRP", "TYR", "TYM", "VAL"]
NAS = ["A", "A5", "A3", "C", "C5", "C3", "G", "G5", "G3", "T", "T5", "T3", "U", \
       "U5", "U3", "RA", "RG", "RC", "RU", "DA", "DG", "DC", "DT"]


class Routines:
    def __init__(self, protein, definition=None):
        """
            Initialize the Routines class.  The class contains most
            of the main routines that run PDB2PQR

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

    def applyNameScheme(self, forcefield):
        """
            Apply the naming scheme of the given forcefield to the atoms
            within the protein

            Parameters
                forcefield: The forcefield object (forcefield)

        """
        _LOGGER.info("Applying the naming scheme to the protein...")
        for residue in self.protein.getResidues():
            if isinstance(residue, (Amino, WAT, Nucleic)):
                resname = residue.ffname
            else:
                resname = residue.name

            for atom in residue.getAtoms():
                rname, aname = forcefield.getNames(resname, atom.name)
                if resname not in ['LIG', 'WAT', 'ACE', 'NME'] and rname != None:
                    try:
                        if (residue.isNterm or residue.isCterm) and rname != residue.name:
                            rname = residue.name
                    except AttributeError:
                        pass
                if aname != None and rname != None:
                    atom.resName = rname
                    atom.name = aname

        _LOGGER.debug("Done applying naming scheme.")

    def applyForcefield(self, forcefield):
        """
            Apply the forcefield to the atoms within the protein

            Parameters
                forcefield: The forcefield object (forcefield)
            Returns
                hitlist:    A list of atoms that were found in
                            the forcefield (list)
                misslist:   A list of atoms that were not found in
                            the forcefield (list)
        """
        _LOGGER.info("Applying the forcefield to the protein...")
        misslist = []
        hitlist = []
        for residue in self.protein.getResidues():
            if isinstance(residue, (Amino, WAT, Nucleic)):
                resname = residue.ffname
            else:
                resname = residue.name

            # Apply the parameters

            for atom in residue.getAtoms():
                atomname = atom.get("name")
                charge, radius = forcefield.getParams(resname, atomname)
                if charge != None and radius != None:
                    atom.set("ffcharge", charge)
                    atom.set("radius", radius)
                    hitlist.append(atom)
                else:
                    misslist.append(atom)

        _LOGGER.debug("Done applying forcefield.")
        return hitlist, misslist

    def updateResidueTypes(self):
        """
            Find the type of residue as notated in the Amino Acid definition
        """
        _LOGGER.info("Updating residue types... ")
        for chain in self.protein.getChains():
            for residue in chain.get("residues"):
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

    def updateSSbridges(self):
        """
            Check for SS-bridge partners, and if present, set appropriate
            partners
        """
        _LOGGER.info("Updating SS bridges...")
        SGpartners = {}
        for residue in self.protein.getResidues():
            if isinstance(residue, CYS):
                atom = residue.getAtom("SG")
                if atom != None:
                    SGpartners[atom] = []

        for atom in SGpartners:
            for partner in SGpartners:
                if atom == partner or SGpartners[atom] != []: continue
                dist = distance(atom.getCoords(), partner.getCoords())
                if dist < BONDED_SS_LIMIT:
                    SGpartners[atom].append(partner)
                    SGpartners[partner].append(atom)

        for atom in SGpartners:
            res1 = atom.get("residue")
            numpartners = len(SGpartners[atom])
            if numpartners == 1:
                partner = SGpartners[atom][0]
                res2 = partner.get("residue")
                res1.set("SSbonded", 1)
                res1.set("SSbondedpartner", partner)
                self.applyPatch("CYX", res1)
                _LOGGER.debug("%s - %s", res1, res2)
            elif numpartners > 1:
                error = "WARNING: %s has multiple potential " % res1
                error += "SS-bridge partners"
                _LOGGER.warn(error)
            elif numpartners == 0:
                _LOGGER.debug("%s is a free cysteine" % res1, 1)
        _LOGGER.debug("Done updating SS bridges.")

    def updateInternalBonds(self):
        """
            Update the internal bonding network using the reference
            objects in each atom.
        """
        for residue in self.protein.getResidues():
            if isinstance(residue, (Amino, WAT, Nucleic)):
                for atom in residue.getAtoms():
                    if not atom.hasReference(): continue
                    for bond in atom.reference.bonds:
                        if not residue.hasAtom(bond): continue
                        bondatom = residue.getAtom(bond)
                        if bondatom not in atom.bonds:
                            atom.addBond(bondatom)

    def updateBonds(self):
        """
            Update the bonding network of the protein.  This happens
            in 3 steps:
              1.  Applying the PEPTIDE patch to all Amino residues
                  so as to add reference for the N(i+1) and C(i-1)
                  atoms
              2.  UpdateInternalBonds for inter-residue linking
              3.  Set the links to the N(i+1) and C(i-1) atoms
        """

        # Apply the peptide patch

        for residue in self.protein.getResidues():
            if isinstance(residue, Amino):
                if residue.isNterm or residue.isCterm:
                    continue
                else:
                    self.applyPatch("PEPTIDE", residue)

        # Update all internal bonds

        self.updateInternalBonds()

        # Set the peptide bond pointers

        for chain in self.protein.getChains():
            for i in range(chain.numResidues() - 1):
                res1 = chain.residues[i]
                res2 = chain.residues[i + 1]
                if not isinstance(res1, Amino) or not isinstance(res2, Amino):
                    continue
                atom1 = res1.getAtom("C")
                atom2 = res2.getAtom("N")

                if atom1 != None:
                    res2.peptideC = atom1
                if atom2 != None:
                    res1.peptideN = atom2
                if atom1 == None or atom2 == None:
                    continue

                if distance(atom1.getCoords(), atom2.getCoords()) > PEPTIDE_DIST:
                    text = "Gap in backbone detected between %s and %s!" % \
                           (res1, res2)
                    _LOGGER.warn(text)
                    res2.peptideC = None
                    res1.peptideN = None

    def applyPatch(self, patchname, residue):
        """
            Apply a patch to the given residue.  This is one of the key
            functions in PDB2PQR.  A similar function appears in
            definitions.py - that version is needed for residue level
            subtitutions so certain protonation states (i.e. CYM, HSE)
            are detectatble on input.

            This version looks up the particular patch name in the
            patchmap stored in the protein, and then applies the
            various commands to the reference and actual residue
            structures.

            See the inline comments for a more detailed explanation.

            Parameters
                patchname:  The name of the patch (string)
                residue:    The residue to apply the patch to (residue)
        """
        if patchname not in self.protein.patchmap:
            raise PDBInternalError("Unable to find patch %s!" % patchname)

        _LOGGER.debug('PATCH INFO: %s patched with %s', residue,patchname)

        # Make a copy of the reference, i.e. a new reference for
        # this patch.  Two examples:
        #     PEPTIDE is a special case, as it applies to
        #             every residue.
        #     CTERM only applies to one specific residue, so a
        #             deep copy is used.

        if patchname == "PEPTIDE":
            newreference = residue.reference
        else:
            newreference = copy.deepcopy(residue.reference)

        patch = self.protein.patchmap[patchname]

        # Add atoms from patch

        for atomname in patch.map:
            newreference.map[atomname] = patch.map[atomname]
            for bond in patch.map[atomname].bonds:
                if bond not in newreference.map: continue
                if atomname not in newreference.map[bond].bonds:
                    newreference.map[bond].bonds.append(atomname)

        # Remove atoms as directed by patch

        for remove in patch.remove:
            if remove in residue.map: residue.removeAtom(remove)
            if remove not in newreference.map: continue
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

        for atom in residue.getAtoms():
            if atom.name in patch.altnames:
                residue.renameAtom(atom.name, patch.altnames[atom.name])

        # Replace each atom's reference with the new one

        for atomname in residue.map:
            if newreference.hasAtom(atomname):
                atom = residue.getAtom(atomname)
                atom.reference = newreference.map[atomname]

    def setStates(self):
        """
            Set the state of each residue.  This is the last step
            before assigning the forcefield, but is necessary so
            as to distinguish between various protonation states.

            See aa.py for residue-specific functions.
        """
        for residue in self.protein.getResidues():
            if isinstance(residue, (Amino, Nucleic)):
                residue.setState()

    def assignTermini(self, chain, neutraln=False, neutralc=False):
        """
            Assign the termini for the given chain by looking at
            the start and end residues.
        """

        if len(chain.residues) == 0:
            text = "Error: chain \"%s\" has 0 residues!" % chain.chainID
            raise PDBInputError(text)

        # Set the N-Terminus/ 5' Terminus

        res0 = chain.residues[0]
        if isinstance(res0, Amino):
            res0.set("isNterm", 1)
            if isinstance(res0, PRO):
                self.applyPatch("NEUTRAL-NTERM", res0)
            elif neutraln:
                self.applyPatch("NEUTRAL-NTERM", res0)
            else:
                self.applyPatch("NTERM", res0)
        elif isinstance(res0, Nucleic):
            res0.set("is5term", 1)
            self.applyPatch("5TERM", res0)

        # Set the C-Terminus/ 3' Terminus

        reslast = chain.residues[-1]
        if isinstance(reslast, Amino):
            reslast.set("isCterm", 1)
            if neutralc:
                self.applyPatch("NEUTRAL-CTERM", reslast)
            else:
                self.applyPatch("CTERM", reslast)
        elif isinstance(reslast, Nucleic):
            reslast.set("is3term", 1)
            self.applyPatch("3TERM", reslast)
        else:
            for i in range(len(chain.residues)):
                resthis = chain.residues[-1 - i]
                if isinstance(resthis, Amino):
                    resthis.set("isCterm", 1)
                    if neutralc:
                        self.applyPatch("NEUTRAL-CTERM", resthis)
                    else:
                        self.applyPatch("CTERM", resthis)
                    break
                elif resthis.name in ["NH2", "NME"]: break
                elif isinstance(resthis, Nucleic):
                    resthis.set("is3term", 1)
                    self.applyPatch("3TERM", resthis)
                    break

    def setTermini(self, neutraln=False, neutralc=False):
        """
            Set the termini for the protein. First set all known
            termini by looking at the ends of the chain. Then
            examine each residue, looking for internal chain breaks.
        """
        # TODO - this function does a lot more than just set termini...
        _LOGGER.info("Setting the termini...")

        # First assign the known termini

        for chain in self.protein.getChains():
            self.assignTermini(chain, neutraln, neutralc)

        # Now determine if there are any hidden chains

        letters = string.ascii_uppercase + string.ascii_lowercase
        c = 0

        while c < len(self.protein.getChains()):
            chain = self.protein.chains[c]
            reslist = []

            origlist = []

            # origlist holds the original residue list for the chain

            for residue in chain.getResidues():
                origlist.append(residue)

            for residue in origlist:
                reslist.append(residue)
                oldid = residue.chainID

                # Look for ending termini

                fixflag = 0
                if isinstance(residue, Amino):
                    if (residue.hasAtom("OXT") and not residue.isCterm):
                        fixflag = 1

                elif isinstance(residue, Nucleic):
                    if ((residue.hasAtom("H3T") or residue.name.endswith("3"))\
                      and not residue.is3term):
                        fixflag = 1

                if fixflag:

                    # Get an available chain ID
                    chainid = letters[0]
                    id = 0
                    idLength = 1
                    while chainid in self.protein.chainmap:
                        id += 1
                        if id >= len(letters):
                            idLength += 1
                            id = 0
                        chainid = letters[id] * idLength

                    if(idLength > 1):
                        message = 'Warning: Reusing chain id: ' + chainid[0] + ''
                        _LOGGER.warn(message)

                    # Make a new chain with these residues
                    newchain = Chain(chainid[0])

                    self.protein.chainmap[chainid] = newchain
                    self.protein.chains.insert(c, newchain)

                    for res in reslist:
                        newchain.addResidue(res)
                        chain.residues.remove(res)
                        res.setChainID(chainid[0])

                    self.assignTermini(chain, neutraln, neutralc)
                    self.assignTermini(newchain, neutraln, neutralc)

                    reslist = []
                    c += 1

            c += 1

        # Update the final chain's chainID if it is "" unless it's all water

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
            id = 0
            idLength = 1
            while chainid in self.protein.chainmap:
                id += 1
                if id >= len(letters):
                    idLength += 1
                    id = 0
                chainid = letters[id] * idLength

            if(idLength > 1):
                message = 'Warning: Reusing chain id: ' + chainid[0]
                _LOGGER.warn(message)

            # Use the new chainID

            self.protein.chainmap[chainid] = chain
            del self.protein.chainmap[""]

            for res in chain.residues:
                res.setChainID(chainid[0])

        _LOGGER.debug("Done setting termini.")


    def findMissingHeavy(self):
        """
            Repair residues that contain missing heavy (non-Hydrogen) atoms
        """
        _LOGGER.info("Checking for missing heavy atoms...")
        misscount = 0
        heavycount = 0
        for residue in self.protein.getResidues():
            if not isinstance(residue, (Amino, Nucleic)): continue

            # Check for Missing Heavy Atoms

            for refatomname in residue.reference.map:
                if refatomname.startswith("H"): continue
                if refatomname in ["N+1", "C-1"]: continue
                if refatomname in ["O1P", "O2P"]:
                    if residue.hasAtom("OP1") and residue.hasAtom("OP2"): continue
                heavycount += 1
                if not residue.hasAtom(refatomname):
                    _LOGGER.debug("Missing %s in %s", refatomname, residue)
                    misscount += 1
                    residue.addMissing(refatomname)

            # Check for Extra Atoms

            atomlist = []
            for atom in residue.get("atoms"):
                atomlist.append(atom)

            for atom in atomlist:
                atomname = atom.get("name")
                if atomname in ["OP1", "OP2"] and residue.reference.hasAtom("O1P") \
                    and residue.reference.hasAtom("O2P"): continue
                if not residue.reference.hasAtom(atomname):
                    _LOGGER.debug("Extra atom %s in %s! - ", atomname, residue)
                    residue.removeAtom(atomname)
                    _LOGGER.debug("Deleted this atom.")

        if heavycount == 0:
            raise PDBInputError("No heavy atoms found. " +
                                "You may also see this message if PDB2PQR does not have parameters for any residue in your protein.")

        misspct = 100.0 * float(misscount) / heavycount
        if misspct > REPAIR_LIMIT:
            error = "This PDB file is missing too many (%i out of " % misscount
            error += "%i, %.2f%%) heavy atoms to accurately repair the file.  " % \
                     (heavycount, misspct)
            error += "The current repair limit is set at %i%%. " % REPAIR_LIMIT
            error += "You may also see this message if PDB2PQR does not have parameters for enough residues in your protein."
            raise PDBInputError(error)
        elif misscount > 0:
            _LOGGER.debug("Missing %i out of %i heavy atoms (%.2f percent) - ",
                          misscount, heavycount, misspct)
            _LOGGER.debug("Will attempt to repair.")
            self.repairHeavy()
        else:
            _LOGGER.debug("No heavy atoms found missing.")
            _LOGGER.debug("Done checking for missing heavy atoms.")

    @staticmethod
    def rebuildTetrahedral(residue, atomname):
        """
            Rebuild a tetrahedral hydrogen group.  This is necessary
            due to the shortcomings of the quatfit routine - given a
            tetrahedral geometry and two existing hydrogens, the
            quatfit routines have two potential solutions.  This function
            uses basic tetrahedral geometry to fix this issue.

            Parameters
                residue:  The residue in question (residue)
                atomname: The atomname to add (string)
            Returns
                1 if successful, 0 otherwise
        """

        hcount = 0
        nextatomname = None

        atomref = residue.reference.map.get(atomname)
        if atomref is None:
            return False
        bondname = atomref.bonds[0]

        # Return if the bonded atom does not exist

        if not residue.hasAtom(bondname):
            return False

        # This group is tetrahedral if bondatom has 4 bonds,
        #  3 of which are hydrogens

        for bond in residue.reference.map[bondname].bonds:
            if bond.startswith("H"):
                hcount += 1
            elif bond != 'C-1' and bond != 'N+1':
                nextatomname = bond

        # Check if this is a tetrahedral group

        if hcount != 3 or nextatomname == None:
            return False

        # Now rebuild according to the tetrahedral geometry

        bondatom = residue.getAtom(bondname)
        nextatom = residue.getAtom(nextatomname)
        numbonds = len(bondatom.bonds)

        if numbonds == 1:

            # Place according to two atoms

            coords = [bondatom.getCoords(), nextatom.getCoords()]
            refcoords = [residue.reference.map[bondname].getCoords(), \
                         residue.reference.map[nextatomname].getCoords()]
            refatomcoords = atomref.getCoords()
            newcoords = findCoordinates(2, coords, refcoords, refatomcoords)
            residue.createAtom(atomname, newcoords)

            # For LEU and ILE residues only: make sure the Hydrogens are in staggered conformation instead of eclipsed.

            if isinstance(residue, LEU):
                hcoords = newcoords
                cbatom = residue.getAtom('CB')
                ang = getDihedral(cbatom.getCoords(), nextatom.getCoords(), bondatom.getCoords(), hcoords)
                diffangle = 60 - ang
                residue.rotateTetrahedral(nextatom, bondatom, diffangle)

            elif isinstance(residue, ILE):
                hcoords = newcoords
                cg1atom = residue.getAtom('CG1')
                cbatom = residue.getAtom('CB')
                if bondatom.name == 'CD1':
                    ang = getDihedral(cbatom.getCoords(), nextatom.getCoords(), bondatom.getCoords(), hcoords)
                elif bondatom.name == 'CG2':
                    ang = getDihedral(cg1atom.getCoords(), nextatom.getCoords(), bondatom.getCoords(), hcoords)
                else:
                    ang = getDihedral(cbatom.getCoords(), nextatom.getCoords(), bondatom.getCoords(), hcoords)

                diffangle = 60 - ang
                residue.rotateTetrahedral(nextatom, bondatom, diffangle)

            return 1

        elif numbonds == 2:

            # Get the single hydrogen coordinates

            hatom = None
            for bond in bondatom.reference.bonds:
                if residue.hasAtom(bond) and bond.startswith("H"):
                    hatom = residue.getAtom(bond)
                    break

            # Use the existing hydrogen and rotate about the bond

            residue.rotateTetrahedral(nextatom, bondatom, 120)
            newcoords = hatom.getCoords()
            residue.rotateTetrahedral(nextatom, bondatom, -120)
            residue.createAtom(atomname, newcoords)

            return 1

        elif numbonds == 3:

            # Find the one spot the atom can be

            hatoms = []
            for bond in bondatom.reference.bonds:
                if residue.hasAtom(bond) and bond.startswith("H"):
                    hatoms.append(residue.getAtom(bond))

            # If this is more than two something is wrong

            if len(hatoms) != 2: return 0

            # Use the existing hydrogen and rotate about the bond

            residue.rotateTetrahedral(nextatom, bondatom, 120)
            newcoords1 = hatoms[0].getCoords()
            residue.rotateTetrahedral(nextatom, bondatom, 120)
            newcoords2 = hatoms[0].getCoords()
            residue.rotateTetrahedral(nextatom, bondatom, 120)

            # Determine which one hatoms[1] is not in

            if distance(hatoms[1].getCoords(), newcoords1) > 0.1:
                residue.createAtom(atomname, newcoords1)
            else:
                residue.createAtom(atomname, newcoords2)

            return 1

    def addHydrogens(self):
        """
            Add the hydrogens to the protein.  This requires either
            the rebuildTetrahedral function for tetrahedral geometries
            or the standard quatfit methods.  These methods use three
            nearby bonds to rebuild the atom; the closer the bonds, the
            more accurate the results.  As such the peptide bonds are
            used when available.
        """
        count = 0
        _LOGGER.info("Adding hydrogens to the protein...")
        for residue in self.protein.getResidues():
            if not isinstance(residue, (Amino, Nucleic)):
                continue
            for atomname in residue.reference.map:
                if not atomname.startswith("H"):
                    continue
                if residue.hasAtom(atomname):
                    continue
                if isinstance(residue, CYS) and residue.SSbonded and atomname == "HG":
                    continue

                # If this hydrogen is part of a tetrahedral group,
                #  follow a different codepath

                if Routines.rebuildTetrahedral(residue, atomname):
                    count += 1
                    continue

                # Otherwise use the standard quatfit methods

                coords = []
                refcoords = []

                refatomcoords = residue.reference.map[atomname].getCoords()
                bondlist = residue.reference.getNearestBonds(atomname)

                for bond in bondlist:
                    if bond == "N+1":
                        atom = residue.peptideN
                    elif bond == "C-1":
                        atom = residue.peptideC
                    else:
                        atom = residue.getAtom(bond)

                    if atom == None:
                        continue

                    # Get coordinates, reference coordinates

                    coords.append(atom.getCoords())
                    refcoords.append(residue.reference.map[bond].getCoords())

                    # Exit if we have enough atoms

                    if len(coords) == 3:
                        break

                if len(coords) == 3:
                    newcoords = findCoordinates(3, coords, refcoords, refatomcoords)
                    residue.createAtom(atomname, newcoords)
                    count += 1
                else:
                    _LOGGER.warn("Couldn't rebuild %s in %s!", atomname, residue)

        _LOGGER.debug(" Added %i hydrogen atoms.", count)

    def removeHydrogens(self):
        _LOGGER.info("Stripping hydrogens from the protein...")

        for residue in self.protein.getResidues():
            if not isinstance(residue, (Amino, Nucleic)):
                continue
            for atom in residue.atoms[:]:
                if atom.isHydrogen():
                    residue.removeAtom(atom.name)

    def repairHeavy(self):
        """
            Repair all heavy atoms.  Unfortunately the first time we
            get to an atom we might not be able to rebuild it - it
            might depend on other atoms to be rebuild first (think side
            chains).  As such a 'seenmap' is used to keep track of what
            we've already seen and subsequent attempts to rebuild the
            atom.
        """
        _LOGGER.info("Rebuilding missing heavy atoms...")
        for residue in self.protein.getResidues():
            if not isinstance(residue, (Amino, Nucleic)):
                continue
            missing = residue.get("missing")
            if missing == []: continue

            # Initialize some variables

            seenmap = {}
            nummissing = len(missing)

            while len(missing) > 0:
                coords = []
                refcoords = []

                atomname = missing.pop(0)
                refatomcoords = residue.reference.map[atomname].getCoords()
                bondlist = residue.reference.getNearestBonds(atomname)

                for bond in bondlist:
                    if bond == "N+1": atom = residue.peptideN
                    elif bond == "C-1": atom = residue.peptideC
                    else: atom = residue.getAtom(bond)

                    if atom == None: continue

                    # Get coordinates, reference coordinates

                    coords.append(atom.getCoords())
                    refcoords.append(residue.reference.map[bond].getCoords())

                    # Exit if we have enough atoms

                    if len(coords) == 3: break

                # We might need other atoms to be rebuilt first

                if len(coords) < 3:
                    try:
                        seenmap[atomname] += 1
                    except KeyError:
                        seenmap[atomname] = 1

                    missing.append(atomname)
                    if seenmap[atomname] > nummissing:
                        text = "Too few atoms present to reconstruct or cap residue %s in structure! " % (residue)
                        text += "This error is generally caused by missing backbone atoms in this protein; "
                        text += "you must use an external program to complete gaps in the protein backbone. "
                        text += "Heavy atoms missing from %s: " % (residue)
                        text += ' '.join(missing)
                        raise PDBInputError(text)

                else: # Rebuild the atom
                    newcoords = findCoordinates(3, coords, refcoords, refatomcoords)
                    residue.createAtom(atomname, newcoords)
                    _LOGGER.debug("Added %s to %s at coordinates", atomname, residue)
                    _LOGGER.debug(" %.3f %.3f %.3f", newcoords[0], newcoords[1], newcoords[2])

        _LOGGER.debug("Done rebuilding missing atoms.")

    def setReferenceDistance(self):
        """
            Set the distance to the CA atom in the residue.
            This is necessary for determining which atoms are
            allowed to move during rotations.  Uses the
            shortestPath algorithm found in utilities.py.
        """
        for residue in self.protein.getResidues():
            if not isinstance(residue, Amino): continue

            # Initialize some variables

            map = {}
            caatom = residue.getAtom("CA")

            if caatom == None:
                text = "Cannot set references to %s without CA atom!"
                raise PDBInputError(text)

            # Set up the linked map

            for atom in residue.getAtoms():
                map[atom] = atom.bonds

            # Run the algorithm

            for atom in residue.getAtoms():
                if atom.isBackbone():
                    atom.refdistance = -1
                elif residue.isCterm and atom.name == "HO":   # special case for HO in Cterm
                    atom.refdistance = 3
                elif residue.isNterm and (atom.name == "H3" or atom.name == "H2"):  # special case for H2 or H3 in Nterm
                    atom.refdistance = 2
                else:
                    atom.refdistance = len(shortestPath(map, atom, caatom)) - 1

    def getbumpscore(self, residue):
        """Get an bump score for the current structure"""

        # Do some setup

        self.cells = Cells(CELL_SIZE)
        self.cells.assignCells(self.protein)

        self.calculateDihedralAngles()
        self.setDonorsAndAcceptors()
        self.updateInternalBonds()
        self.setReferenceDistance()
        bumpscore = 0.0
        #for residue in self.protein.getResidues():
        if not isinstance(residue, Amino): return 0.0
        # Initialize variables

        conflictnames = []

        for atom in residue.getAtoms():
            atomname = atom.name
            #if not atom.added: continue
            if atomname[0] != "H": continue
            #if atom.optimizeable: continue
            #print atomname,atom.optimizeable,atom.added
            bumpscore = bumpscore + self.getbumpscore_atom(atom)
        return bumpscore


    def getbumpscore_atom(self, atom):
        """
            Find nearby atoms for conflict-checking.  Uses
            neighboring cells to compare atoms rather than an all
            versus all O(n^2) algorithm, which saves a great deal
            of time.  There are several instances where we ignore
            potential conflicts; these include donor/acceptor pairs,
            atoms in the same residue, and bonded CYS bridges.

            Parameters
                atom:  Find nearby atoms to this atom (Atom)
            Returns
                bumpscore: a bump score sum((dist-cutoff)**20 for all near atoms

            Jens rewrote this function from findNearbyAtoms to
            be usable for detecting bumps for optimzable hydrogens
        """
        # Initialize some variables
        nearatoms = {}
        residue = atom.residue
        #print 'Cutoff distance',cutoff
        atom_size = BUMP_HYDROGEN_SIZE if atom.isHydrogen() else BUMP_HEAVY_SIZE

        # Get atoms from nearby cells

        closeatoms = self.cells.getNearCells(atom)

        # Loop through and see if any are within the cutoff

        bumpscore = 0.0
        for closeatom in closeatoms:
            closeresidue = closeatom.residue
            if closeresidue == residue and (closeatom in atom.bonds or atom in closeatom.bonds):
                continue

            if not isinstance(closeresidue, Amino):
                continue
            if isinstance(residue, CYS):
                if residue.SSbondedpartner == closeatom: continue

            # Also ignore if this is a donor/acceptor pair
            pair_ignored = False
            if atom.isHydrogen() and len(atom.bonds) != 0 and atom.bonds[0].hdonor \
               and closeatom.hacceptor:
                #pair_ignored = True
                continue
            if closeatom.isHydrogen() and len(closeatom.bonds) != 0 and closeatom.bonds[0].hdonor \
                   and atom.hacceptor:
                #pair_ignored = True
                continue

            dist = distance(atom.getCoords(), closeatom.getCoords())
            other_size = BUMP_HYDROGEN_SIZE if closeatom.isHydrogen() else BUMP_HEAVY_SIZE
            cutoff = atom_size + other_size
            if dist < cutoff:
                bumpscore = bumpscore + 1000.0
                if pair_ignored:
                    _LOGGER.debug('This bump is a donor/acceptor pair.')
#                if heavy_not_ignored:
#                    print 'Bumped {0} against {1} within residue'.format(atom.name, closeatom.name)
                #nearatoms[closeatom] = (dist-cutoff)**2
        _LOGGER.debug('BUMPSCORE ' + str(bumpscore))
        return bumpscore


    def debumpProtein(self):
        """
            Make sure that none of the added atoms were rebuilt
            on top of existing atoms.  See each called function
            for more information.
        """

        _LOGGER.info("Checking if we must debump any residues... ")

        # Do some setup

        self.cells = Cells(CELL_SIZE)
        self.cells.assignCells(self.protein)

        self.calculateDihedralAngles()
        self.setDonorsAndAcceptors()
        self.updateInternalBonds()
        self.setReferenceDistance()

        # Determine which residues to debump

        for residue in self.protein.getResidues():
            if not isinstance(residue, Amino): continue

            # Initialize variables

            conflictnames = self.findResidueConflicts(residue, True)

            if not conflictnames:
                continue

            # Otherwise debump the residue

            _LOGGER.debug("Starting to debump %s...", residue)
            _LOGGER.debug("Debumping cutoffs: %2.1f for heavy-heavy, %2.1f for hydrogen-heavy, and %2.1f for hydrogen-hydrogen." %
                       (BUMP_HEAVY_SIZE*2,
                        BUMP_HYDROGEN_SIZE+BUMP_HEAVY_SIZE,
                        BUMP_HYDROGEN_SIZE*2), 1)
            if self.debumpResidue(residue, conflictnames):
                _LOGGER.debug("Debumping Successful!", 1)
            else:
                text = "WARNING: Unable to debump %s" % residue
                _LOGGER.warn(text)

        _LOGGER.debug("Done checking if we must debump any residues.")

    def findResidueConflicts(self, residue, writeConflictInfo=False):
        conflictnames = []
        for atom in residue.getAtoms():
            atomname = atom.name
            if not atom.added: continue
            if atomname == "H": continue
            if atom.optimizeable: continue

            nearatoms = self.findNearbyAtoms(atom)

            # If something is too close, we must debump the residue

            if nearatoms != {}:
                conflictnames.append(atomname)
                if writeConflictInfo:
                    for repatom in nearatoms:
                        _LOGGER.debug("%s %s is too close to %s %s", residue, atomname, repatom.residue, repatom.name)

        return conflictnames

    def scoreDihedralAngle(self, residue, anglenum):
        score = 0
        atomnames = residue.reference.dihedrals[anglenum].split()
        pivot = atomnames[2]
        moveablenames = self.getMoveableNames(residue, pivot)
        for name in moveablenames:
            nearatoms = self.findNearbyAtoms(residue.getAtom(name))
            for v in nearatoms.values():
                score += v

        return score


    def debumpResidue(self, residue, conflictnames):
        """
            Debump a specific residue.  Only should be called
            if the residue has been detected to have a conflict.
            If called, try to rotate about dihedral angles to
            resolve the conflict.

            Parameters
                residue:  The residue in question
                conflictnames:  A list of atomnames that were
                                rebuilt too close to other atoms
            Returns
                True if successful, False otherwise
        """

        # Initialize some variables

        ANGLE_STEPS = 72
        ANGLE_STEP_SIZE = float(360 // ANGLE_STEPS)

        ANGLE_TEST_COUNT = 10

        anglenum = -1
        currentConflictNames = conflictnames

        EPSILON = 0.0000001

        # Try (up to 10 times) to find a workable solution

        for _ in range(ANGLE_TEST_COUNT):

            anglenum = self.pickDihedralAngle(residue, currentConflictNames, anglenum)

            if anglenum == -1: return False

            _LOGGER.debug("Using dihedral angle number %i to debump the residue.", anglenum)

            bestscore = self.scoreDihedralAngle(residue, anglenum)
            foundImprovement = False
            bestangle = originalAngle = residue.dihedrals[anglenum]

            #Skip the first angle as it's already known.
            for i in range(1, ANGLE_STEPS):
                newangle = originalAngle + (ANGLE_STEP_SIZE * i)
                self.setDihedralAngle(residue, anglenum, newangle)

                # Check for conflicts

                score = self.scoreDihedralAngle(residue, anglenum)

                if score == 0:
                    if not self.findResidueConflicts(residue):
                        _LOGGER.debug("No conflicts found at angle " + repr(newangle))
                        return True
                    else:
                        bestangle = newangle
                        foundImprovement = True
                        break

                # Set the best angle
                elif score < bestscore:
                    diff = abs(bestscore - score)
                    #Don't update if it's effectively a tie
                    if diff > EPSILON:
                        bestscore = score
                        bestangle = newangle
                        foundImprovement = True

            self.setDihedralAngle(residue, anglenum, bestangle)
            currentConflictNames = self.findResidueConflicts(residue)

            if foundImprovement:
                _LOGGER.debug("Best score of "+repr(bestscore)+" at angle "+repr(bestangle)+". New conflict set: ")
                _LOGGER.debug(str(currentConflictNames))
            else:
                _LOGGER.debug("No improvement found for this dihedral angle.")

        # If we're here, debumping was unsuccessful

        return False

    def calculateDihedralAngles(self):
        """
            Calculate the dihedral angle for every residue within the protein
        """
        for residue in self.protein.getResidues():
            if not isinstance(residue, Amino): continue
            residue.dihedrals = []

            refangles = residue.reference.dihedrals
            for di in refangles:
                coords = []
                atoms = di.split()
                for i in range(4):
                    atomname = atoms[i]
                    if residue.hasAtom(atomname):
                        coords.append(residue.getAtom(atomname).getCoords())

                if len(coords) == 4: angle = getDihedral(coords[0], coords[1], coords[2], coords[3])
                else: angle = None

                residue.addDihedralAngle(angle)

    def getClosestAtom(self, atom):
        """
            Get the closest atom that does not form a donor/acceptor pair.
            Used to detect potential conflicts.

            NOTE:  Cells must be set before using this function.

            Parameters
                atom:  The atom in question (Atom)
            Returns
                bestatom:  The closest atom to the input atom that does not
                           satisfy a donor/acceptor pair.
        """
        # Initialize some variables

        bestdist = 999.99
        bestwatdist = 999.99
        bestatom = None
        bestwatatom = None
        residue = atom.residue

        # Get atoms from nearby cells

        closeatoms = self.cells.getNearCells(atom)

        # Loop through and see which is the closest

        for closeatom in closeatoms:
            closeresidue = closeatom.residue
            if closeresidue == residue: continue
            if not isinstance(closeresidue, (Amino,WAT)): continue
            if isinstance(residue, CYS):
                if residue.SSbondedpartner == closeatom: continue

            # Also ignore if this is a donor/acceptor pair

            if atom.isHydrogen() and atom.bonds[0].hdonor \
               and closeatom.hacceptor: continue
            if closeatom.isHydrogen() and closeatom.bonds[0].hdonor \
                   and atom.hacceptor:
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
            txt = "Warning: %s in %s skipped when optimizing %s in %s" % (bestwatatom.name,
                                                                           bestwatatom.residue,
                                                                           atom.name, residue)
            _LOGGER.warn(txt)

        return bestatom

    def findNearbyAtoms(self, atom):
        """
            Find nearby atoms for conflict-checking.  Uses
            neighboring cells to compare atoms rather than an all
            versus all O(n^2) algorithm, which saves a great deal
            of time.  There are several instances where we ignore
            potential conflicts; these include donor/acceptor pairs,
            atoms in the same residue, and bonded CYS bridges.

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

        closeatoms = self.cells.getNearCells(atom)

        # Loop through and see if any are within the cutoff

        for closeatom in closeatoms:
            closeresidue = closeatom.residue
            if closeresidue == residue and (closeatom in atom.bonds or atom in closeatom.bonds):
                continue

            if not isinstance(closeresidue, (Amino, WAT)):
                continue
            if isinstance(residue, CYS) and residue.SSbondedpartner == closeatom:
                continue

            # Also ignore if this is a donor/acceptor pair

            if (atom.isHydrogen() and len(atom.bonds) != 0 and
                atom.bonds[0].hdonor and closeatom.hacceptor):
                continue

            if (closeatom.isHydrogen() and len(closeatom.bonds) != 0 and
                closeatom.bonds[0].hdonor and atom.hacceptor):
                continue

            dist = distance(atom.getCoords(), closeatom.getCoords())
            other_size = BUMP_HYDROGEN_SIZE if closeatom.isHydrogen() else BUMP_HEAVY_SIZE
            cutoff = atom_size + other_size
            if dist < cutoff:
                nearatoms[closeatom] = cutoff - dist

        return nearatoms


    def pickDihedralAngle(self, residue, conflictnames, oldnum=None):
        """
            Choose an angle number to use in debumping

            Algorithm
                Instead of simply picking a random chiangle, this function
                uses a more intelligent method to improve efficiency.
                The algorithm uses the names of the conflicting atoms
                within the residue to determine which angle number
                has the best chance of fixing the problem(s). The method
                also insures that the same chiangle will not be run twice
                in a row.
            Parameters
                residue:    The residue that is being debumped (Residue)
                conflictnames: A list of atom names that are currently
                               conflicts (list)
                oldnum    : The old dihedral angle number (int)
            Returns
                bestnum    : The new dihedral angle number (int)
        """
        bestnum = -1
        best = 0

        iList = list(range(len(residue.dihedrals)))
        #Make sure our testing is done round robin.
        if oldnum is not None and oldnum >= 0 and len(iList) > 0:
            del iList[oldnum]
            testDihedralIndecies = iList[oldnum:] + iList[:oldnum]
        else:
            testDihedralIndecies = iList

        for i in testDihedralIndecies:
            if i == oldnum: continue
            if residue.dihedrals[i] is None: continue

            score = 0
            atomnames = residue.reference.dihedrals[i].split()
            pivot = atomnames[2]

            moveablenames = self.getMoveableNames(residue, pivot)

            # If this pivot only moves the conflict atoms, pick it

            if conflictnames == moveablenames: return i

            # Otherwise find the pivot with the most matches

            for name in conflictnames:
                if name in moveablenames:
                    score += 1
                    if score > best:
                        best = score
                        bestnum = i

        # Return the best angle.  If none were found, return -1.

        return bestnum

    def setDihedralAngle(self, residue, anglenum, angle):
        """
            Rotate a residue about a given angle. Uses the quatfit
            methods to perform the matrix mathematics.

            Parameters
                residue:   The residue to rotate
                anglenum:  The number of the angle to rotate as
                           listed in residue.dihedrals
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
            if residue.hasAtom(atomname):
                coordlist.append(residue.getAtom(atomname).getCoords())
            else:
                raise PDBInputError("Error occurred while trying to debump!")

        initcoords = subtract(coordlist[2], coordlist[1])

        moveablenames = self.getMoveableNames(residue, pivot)

        for name in moveablenames:
            atom = residue.getAtom(name)
            movecoords.append(subtract(atom.getCoords(), coordlist[1]))

        newcoords = qchichange(initcoords, movecoords, diff)

        for i in range(len(moveablenames)):
            atom = residue.getAtom(moveablenames[i])
            self.cells.removeCell(atom)
            x = (newcoords[i][0] + coordlist[1][0])
            y = (newcoords[i][1] + coordlist[1][1])
            z = (newcoords[i][2] + coordlist[1][2])
            atom.set("x", x)
            atom.set("y", y)
            atom.set("z", z)
            self.cells.addCell(atom)


        # Set the new angle

        coordlist = []
        for atomname in atomnames:
            if residue.hasAtom(atomname):
                coordlist.append(residue.getAtom(atomname).getCoords())
            else:
                raise PDBInputError("Error occurred while trying to debump!")

        di = getDihedral(coordlist[0], coordlist[1], coordlist[2], coordlist[3])
        residue.dihedrals[anglenum] = di

    def getMoveableNames(self, residue, pivot):
        """
            Return all atomnames that are further away than the
            pivot atom.

            Parameters
                residue:  The residue to use
                pivot:    The pivot atomname
        """
        movenames = []
        refdist = residue.getAtom(pivot).refdistance
        for atom in residue.getAtoms():
            if atom.refdistance > refdist:
                movenames.append(atom.name)
        return movenames

    def setDonorsAndAcceptors(self):
        """
            Set the donors and acceptors within the protein
        """
        for residue in self.protein.getResidues():
            residue.setDonorsAndAcceptors()

    def runPDB2PKA(self, ph, ff, pdblist, ligand, pdb2pka_params):
        if ff.lower() != 'parse':
            PDB2PKAError('PDB2PKA can only be run with the PARSE force field.')

        _LOGGER.info("Running PDB2PKA and applying at pH %.2f... ", ph)

        import pka
        from pdb2pka import pka_routines
        
        init_params = pdb2pka_params.copy()
        init_params.pop('pairene')
        init_params.pop('clean_output')
        
        results = pka.pre_init(original_pdb_list=pdblist, ff=ff, ligand=ligand,
                               **init_params)
        output_dir, protein, routines, forcefield,apbs_setup, ligand_titratable_groups, maps, sd = results

        mypkaRoutines = pka_routines.pKaRoutines(protein, routines, forcefield, apbs_setup, output_dir, maps, sd,
                                                 restart=pdb2pka_params.get('clean_output'),
                                                 pairene=pdb2pka_params.get('pairene'))

        _LOGGER.info('Doing full pKa calculation')
        mypkaRoutines.runpKa()

        pdb2pka_warnings = mypkaRoutines.warnings[:]

        _LOGGER.warn(pdb2pka_warnings)

        residue_ph = {}
        for pka_residue_tuple, calc_ph in mypkaRoutines.ph_at_0_5.items():
            tit_type, chain_id, number_str = pka_residue_tuple
            if tit_type == 'NTR':
                tit_type = 'N+'
            elif tit_type == 'CTR':
                tit_type = 'C-'

            key = ' '.join([tit_type, number_str, chain_id])
            residue_ph[key] = calc_ph

        pformat(residue_ph)

        self.apply_pka_values(ff, ph, residue_ph)

        _LOGGER.debug('Finished running PDB2PKA.')


    def runPROPKA31(self, pka_options):
        """
            Run PROPKA 3.1 on the current protein, setting protonation states to
            the correct values. pH is set in pka_options

            Parameters
               pka_options: Options for propKa 3.1, including pH

            Returns
               pka_molecule: pKa's internal molecule object (including pKa's, etc)
               not_found:    dict of residues found in pka_molecule but not in PDB2PQR (with pKa)
        """

        # See https://github.com/jensengroup/propka-3.1/blob/master/scripts/propka31.py
        import propka.molecular_container
        import tempfile

        def delete_propka_input(fn):
            import os
            p, f = os.path.split(fn)
            f = f.replace('.pdb', '.propka_input')
            os.remove(f)

        ph = pka_options.pH
        _LOGGER.info("Running propka 3.1 at pH %.2f... ", ph)

        # Initialize some variables
        pkadic = {}

        # Reorder the atoms in each residue to start with N - TONI is this necessary?
        for residue in self.protein.getResidues():
            residue.reorder()

        # TONI Make a string with all non-hydrogen atoms. Previously it was removing the "element"
        # column and hydrogens. This does not seem to be necessary in propKa 3.1 .
        HFreeProteinFile = tempfile.NamedTemporaryFile(mode="w+", suffix=".pdb", dir=os.getcwd(), delete=False)
        for atom in self.protein.getAtoms():
            if not atom.isHydrogen():
                atomtxt = atom.getPDBString()
                HFreeProteinFile.write(atomtxt + '\n')
        #HFreeProteinFile.seek(0)
        HFreeProteinFile.close();

        # Run PropKa 3.1 -------------

        # Creating protein object. Annoyingly, at this stage propka generates a
        # *.propka_input file in PWD and does not delete it (irregardless of the original .pdb location)
        pka_molecule = propka.molecular_container.Molecular_container(HFreeProteinFile.name, pka_options)
        #HFreeProteinFile.open()
        delete_propka_input(HFreeProteinFile.name)
        HFreeProteinFile.close()

        # calculating pKa values for ionizable residues -
        pka_molecule.calculate_pka()

        ##  pka_molecule.write_pka()

        for grp in pka_molecule.conformations['AVR'].groups:
            key = str.strip('%s %s %s' % (grp.residue_type, grp.atom.resNumb, grp.atom.chainID))
            pkadic[key] = grp.pka_value

        self.protein.pka_protein = pka_molecule

        return pkadic



    def runPROPKA(self, ph, ff, rootname, outname, options, version=30):
        """
            Run PROPKA on the current protein, setting protonation states to
            the correct values

            Parameters
               ph:  The desired pH of the system
               ff:  The forcefield name to be used
               outname: The name of the PQR outfile
               options: Options to propka
               version: may be 30 or 31 (uses external propka 3.1)
        """
        _LOGGER.info("Running PROPKA v%d and applying at pH %.2f... ", version, ph)

        pkadic = self.runPROPKA31(options)

        if len(pkadic) == 0:
            raise ValueError("PROPKA returned empty results!")
            return

        # Now apply each pka to the appropriate residue
        self.apply_pka_values(ff, ph, pkadic)

        _LOGGER.debug("Done running PROPKA")


    def apply_pka_values(self, ff, ph, pkadic):
        _LOGGER.info('Applying pKa values at a pH of %.2f:', ph)
        formatted_pkadict = pformat(pkadic)
        _LOGGER.debug(formatted_pkadict+'')

        for residue in self.protein.getResidues():
            if not isinstance(residue, Amino):
                continue
            resname = residue.name
            resnum = residue.resSeq
            chainID = residue.chainID

            if residue.isNterm:
                key = "N+ %i %s" % (resnum, chainID)
                key = key.strip()
                if key in pkadic:
                    value = pkadic[key]
                    del pkadic[key]
                    if ph >= value:
                        if ff in ["amber", "charmm", "tyl06", "peoepb", "swanson"]:
                            warn = ("N-terminal %s" % key, "neutral")
                            _LOGGER.warn(warn)
                        else:
                            self.applyPatch("NEUTRAL-NTERM", residue)

            if residue.isCterm:
                key = "C- %i %s" % (resnum, chainID)
                key = key.strip()
                if key in pkadic:
                    value = pkadic[key]
                    del pkadic[key]
                    if ph < value:
                        if ff in ["amber", "charmm", "tyl06", "peoepb", "swanson"]:
                            warn = ("C-terminal %s" % key, "neutral")
                            _LOGGER.warn(warn)
                        else:
                            self.applyPatch("NEUTRAL-CTERM", residue)

            key = "%s %i %s" % (resname, resnum, chainID)
            key = key.strip()
            if key in pkadic:
                value = pkadic[key]
                del pkadic[key]
                if resname == "ARG" and ph >= value:
                    if ff == "parse":
                        self.applyPatch("AR0", residue)
                        _LOGGER.warn("Neutral arginines are very rare. Please double-check your setup.")
                    else:
                        warn = (key, "neutral")
                        _LOGGER.warn(warn)
                elif resname == "ASP" and ph < value:
                    if residue.isCterm and ff in ["amber", "tyl06", "swanson"]:
                        warn = (key, "Protonated at C-Terminal")
                        _LOGGER.warn(warn)
                    elif residue.isNterm and ff in ["amber", "tyl06", "swanson"]:
                        warn = (key, "Protonated at N-Terminal")
                        _LOGGER.warn(warn)
                    else:
                        self.applyPatch("ASH", residue)
                elif resname == "CYS" and ph >= value:
                    if ff == "charmm":
                        warn = (key, "negative")
                        _LOGGER.warn(warn)
                    else:
                        self.applyPatch("CYM", residue)
                elif resname == "GLU" and ph < value:
                    if residue.isCterm and ff in ["amber", "tyl06", "swanson"]:
                        warn = (key, "Protonated at C-Terminal")
                        _LOGGER.warn(warn)
                    elif residue.isNterm and ff in ["amber", "tyl06", "swanson"]:
                        warn = (key, "Protonated at N-Terminal")
                        _LOGGER.warn(warn)
                    else:
                        self.applyPatch("GLH", residue)
                elif resname == "HIS" and ph < value:
                    self.applyPatch("HIP", residue)
                elif resname == "LYS" and ph >= value:
                    if ff == "charmm":
                        warn = (key, "neutral")
                        _LOGGER.warn(warn)
                    elif ff in ["amber", "tyl06", "swanson"] and residue.get("isCterm"):
                        warn = (key, "neutral at C-Terminal")
                        _LOGGER.warn(warn)
                    elif ff == "tyl06" and residue.get("isNterm"):
                        warn = (key, "neutral at N-Terminal")
                        _LOGGER.warn(warn)
                    else:
                        self.applyPatch("LYN", residue)
                elif resname == "TYR" and ph >= value:
                    if ff in ["charmm", "amber", "tyl06", "peoepb", "swanson"]:
                        warn = (key, "negative")
                        _LOGGER.warn(warn)
                    else:
                        self.applyPatch("TYM", residue)

        if len(warnings) > 0:
            init = "PDB2PKA determined the following residues to be in a protonation state not supported by the {ff} forcefield.  All were reset to their pH 7.0 model titration state."
            _LOGGER.warn(init.format(ff=options.ff))
            _LOGGER.warn(init)
            for warn in warnings:
                text = "             %s (%s)" % (warn[0], warn[1])
                _LOGGER.warn(text)

        if len(pkadic) > 0:
            warn = "PDB2PQR could not identify the following residues and residue numbers as returned by PROPKA or PDB2PKA:"
            _LOGGER.warn(warn)
            for item in pkadic:
                text = "             %s" % item
                _LOGGER.warn(text)

    def holdResidues(self, hlist):
        """Set the stateboolean dictionary to residues in hlist."""

        if not hlist:
            return

        hlist_copy = hlist.copy()
        for residue in self.protein.getResidues():
            reskey = (residue.resSeq, residue.chainID, residue.iCode)
            if reskey in hlist:
                hlist.remove(reskey)
                if isinstance(residue, Amino):
                    residue.stateboolean = {'FIXEDSTATE': False}
                    _LOGGER.debug("Setting residue {:s} as fixed.".format(str(residue)))
                else:
                    _LOGGER.warn("Matched residue {:s} but not subclass of Amino.".format(str(residue)))

        if len(hlist) > 0:
            _LOGGER.warn("The following fixed residues were not matched (possible internal error): {:s}."
                       .format(str(hlist)))


class Cells:
    """
        The cells object provides a better way to search for nearby atoms. A
        pure all versus all search is O(n^2) - for every atom, every other atom
        must be searched.  This is rather inefficient, especially for large
        proteins where cells may be tens of angstroms apart.  The cell class
        breaks down the xyz protein space into several 3-D cells of desired
        size - then by simply examining atoms that fall into the adjacent
        cells one can quickly find nearby cells.

        NOTE:  Ideally this should be somehow separated from the routines
               object...
        """
    def __init__(self, cellsize):
        """
            Initialize the cells.

            Parameters
                cellsize:  The size of each cell (int)
        """
        self.cellmap = {}
        self.cellsize = cellsize

    def assignCells(self, protein):
        """
            Place each atom in a virtual cell for easy neighbor comparison
        """
        for atom in protein.getAtoms():
            atom.cell = None
            self.addCell(atom)

    def addCell(self, atom):
        """
            Add an atom to the cell

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

    def removeCell(self, atom):
        """
             Remove the atom from a cell

             Parameters
                 atom:   The atom to add (atom)
        """
        oldcell = atom.get("cell")
        if oldcell == None: return
        atom.set("cell", None)
        self.cellmap[oldcell].remove(atom)

    def getNearCells(self, atom):
        """
            Find all atoms in bordering cells to an atom

            Parameters
                atom:  The atom to use (atom)
            Returns
                closeatoms:  A list of nearby atoms (list)
        """
        size = self.cellsize
        closeatoms = []
        cell = atom.get("cell")
        if cell == None:
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
                                if atom == atom2: continue
                                closeatoms.append(atom2)
                        except KeyError: pass

            return closeatoms
