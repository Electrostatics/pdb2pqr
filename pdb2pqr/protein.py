"""Routines for PDB2PQR

This module contains the protein object used in PDB2PQR and associated methods

TODO - this module should be broken into separate files.

Authors:  Todd Dolinsky, Yong Huang
"""
import string
import logging
import copy
import pprint
from . import residue as residue_
from . import aa
from . import na
from . import utilities as util
from . import structures as struct
from . import pdb
from . import forcefield
from . import quatfit as quat
from .config import AA_NAMES, NA_NAMES, BONDED_SS_LIMIT, PEPTIDE_DIST


_LOGGER = logging.getLogger(__name__)


class Protein(object):
    """Protein class

    The protein class represents the parsed PDB, and provides a hierarchy of
    information - each Protein contains a list of Chain objects as provided in
    the PDB file.  Each Chain then contains its associated list of Residue
    objects, and each Residue contains a list of Atom objects, completing the
    hierarchy.
    """

    def __init__(self, pdblist, definition):
        """Initialize using parsed PDB file

        Args:
            pdblist: List of Classes of PDB lines as created
        """
        self.chainmap = {}
        self.chains = []
        self.residues = []
        self.definition = definition

        chain_dict = {}
        previous_atom = None
        residue = []
        num_models = 0
        num_chains = 1
        count = 0

        for record in pdblist: # Find number of chains
            if isinstance(record, pdb.TER):
                num_chains += 1

        for record in pdblist:
            if isinstance(record, (pdb.ATOM, pdb.HETATM)):
                if record.chain_id == "":
                    if num_chains > 1 and record.res_name not in ["WAT", "HOH"]:
                        # Assign a chain ID
                        record.chain_id = string.ascii_uppercase[count]

                chain_id = record.chain_id
                res_seq = record.res_seq
                ins_code = record.ins_code

                if previous_atom is None:
                    previous_atom = record

                if chain_id not in chain_dict:
                    my_chain = struct.Chain(chain_id)
                    chain_dict[chain_id] = my_chain

                if res_seq != previous_atom.res_seq or \
                      ins_code != previous_atom.ins_code or \
                      chain_id != previous_atom.chain_id:
                    my_residue = self.create_residue(residue, previous_atom.res_name)
                    chain_dict[previous_atom.chain_id].add_residue(my_residue)
                    residue = []

                residue.append(record)
                previous_atom = record

            elif isinstance(record, pdb.END):
                my_residue = self.create_residue(residue, previous_atom.res_name)
                chain_dict[previous_atom.chain_id].add_residue(my_residue)
                residue = []

            elif isinstance(record, pdb.MODEL):
                num_models += 1
                if residue == []:
                    continue
                if num_models > 1:
                    my_residue = self.create_residue(residue,
                                                     previous_atom.res_name)
                    chain_dict[previous_atom.chain_id].add_residue(my_residue)
                    break

            elif isinstance(record, pdb.TER):
                count += 1

        if residue != [] and num_models <= 1:
            my_residue = self.create_residue(residue, previous_atom.res_name)
            chain_dict[previous_atom.chain_id].add_residue(my_residue)

        # Keep a map for accessing chains via chain_id
        self.chainmap = chain_dict.copy()

        # Make a list for sequential ordering of chains
        if "" in chain_dict:
            chain_dict["ZZ"] = chain_dict[""]
            del chain_dict[""]

        keys = list(chain_dict.keys())
        keys.sort()

        for key in keys:
            self.chains.append(chain_dict[key])

        for chain in self.chains:
            for residue in chain.residues:
                self.residues.append(residue)

    @property
    def num_missing_heavy(self):
        """Return number of missing biomolecular heavy atoms in structure."""
        natom = 0
        for residue in self.residues:
            if not isinstance(residue, (aa.Amino, na.Nucleic)):
                continue
            residue.missing = []
            for refatomname in residue.reference.map:
                if refatomname.startswith("H"):
                    continue
                if refatomname in ["N+1", "C-1"]:
                    continue
                if refatomname in ["O1P", "O2P"]:
                    if residue.has_atom("OP1") and residue.has_atom("OP2"):
                        continue
                if not residue.has_atom(refatomname):
                    _LOGGER.debug("Missing %s in %s", refatomname, residue)
                    natom += 1
                    residue.missing.append(refatomname)
        return natom

    @property
    def num_heavy(self):
        """Return number of biomolecular heavy atoms in structure.

        TODO - figure out if this is redundant with self.num_bio_atoms
        num_bio_atoms include hydrogens but those are stripped off eventually
        """
        natom = 0
        for residue in self.residues:
            if not isinstance(residue, (aa.Amino, na.Nucleic)):
                continue
            for refatomname in residue.reference.map:
                if refatomname.startswith("H"):
                    continue
                if refatomname in ["N+1", "C-1"]:
                    continue
                if refatomname in ["O1P", "O2P"]:
                    if residue.has_atom("OP1") and residue.has_atom("OP2"):
                        continue
                natom += 1
        return natom

    @property
    def reference_map(self):
        """Return definition reference map."""
        return self.definition.map

    @property
    def patch_map(self):
        """Return definition patch map."""
        return self.definition.patches

    @property
    def num_bio_atoms(self):
        """Return the number of ATOM (not HETATM) records protein."""
        natom = 0
        for atom in self.atoms:
            if atom.type == "ATOM":
                natom += 1
        return natom

    def set_hip(self):
        """Set all HIS states to HIP."""
        for residue in self.residues:
            if isinstance(residue, aa.HIS):
                self.apply_patch("HIP", residue)

    def set_termini(self, neutraln=False, neutralc=False):
        """Set the termini for the protein.

        First set all known termini by looking at the ends of the chain. Then
        examine each residue, looking for internal chain breaks.

        TODO:  this function needs to be cleaned and simplified

        Args:
            neutraln:  whether N-terminus is neutral
            neutralc:  whether C-terminus is neutral
        """
        # First assign the known termini
        for chain in self.chains:
            self.assign_termini(chain, neutraln, neutralc)

        # Now determine if there are any hidden chains
        letters = string.ascii_uppercase + string.ascii_lowercase
        ch_num = 0

        while ch_num < len(self.chains):
            chain = self.chains[ch_num]
            reslist = []
            origlist = []

            # origlist holds the original residue list for the chain
            for residue in chain.residues:
                origlist.append(residue)

            for residue in origlist:
                reslist.append(residue)

                # Look for ending termini
                fixflag = 0
                if isinstance(residue, aa.Amino):
                    if (residue.has_atom("OXT") and not residue.is_c_term):
                        fixflag = 1

                elif isinstance(residue, na.Nucleic):
                    if ((residue.has_atom("H3T") or residue.name.endswith("3"))\
                      and not residue.is3term):
                        fixflag = 1

                if fixflag:
                    # Get an available chain ID
                    chainid = letters[0]
                    id_ = 0
                    id_length = 1
                    while chainid in self.chainmap:
                        id_ += 1
                        if id_ >= len(letters):
                            id_length += 1
                            id_ = 0
                        chainid = letters[id_] * id_length

                    if id_length > 1:
                        message = 'Warning: Reusing chain id: ' + chainid[0] + ''
                        _LOGGER.warning(message)

                    # Make a new chain with these residues
                    newchain = struct.Chain(chainid[0])

                    self.chainmap[chainid] = newchain
                    self.chains.insert(ch_num, newchain)

                    for res in reslist:
                        newchain.add_residue(res)
                        chain.residues.remove(res)
                        res.set_chain_id(chainid[0])

                    self.assign_termini(chain, neutraln, neutralc)
                    self.assign_termini(newchain, neutraln, neutralc)

                    reslist = []
                    ch_num += 1
            ch_num += 1

        # Update the final chain's chain_id if it is "" unless it's all water
        if "" in self.chainmap:
            notwat = 0
            for res in chain.residues:
                if not isinstance(res, aa.WAT):
                    notwat = 1
                    break

            if notwat == 0:
                _LOGGER.debug("Done setting termini.")
                return

            chain = self.chainmap[""]
            chainid = letters[0]
            id_ = 0
            id_length = 1
            while chainid in self.chainmap:
                id_ += 1
                if id_ >= len(letters):
                    id_length += 1
                    id_ = 0
                chainid = letters[id_] * id_length

            if id_length > 1:
                message = 'Warning: Reusing chain id: ' + chainid[0]
                _LOGGER.warning(message)

            # Use the new chain_id
            self.chainmap[chainid] = chain
            del self.chainmap[""]

            for res in chain.residues:
                res.set_chain_id(chainid[0])
        _LOGGER.debug("Done setting termini.")

    def set_states(self):
        """Set the state of each residue.

        This is the last step before assigning the forcefield, but is necessary
        so as to distinguish between various protonation states.

        See aa.py for residue-specific functions.
        """
        for residue in self.residues:
            if isinstance(residue, (aa.Amino, na.Nucleic)):
                residue.set_state()

    def add_hydrogens(self):
        """Add the hydrogens to the protein.

        This requires either the rebuild_tetrahedral function for tetrahedral
        geometries or the standard quatfit methods.  These methods use three
        nearby bonds to rebuild the atom; the closer the bonds, the more
        accurate the results.  As such the peptide bonds are used when available.
        """
        count = 0
        for residue in self.residues:
            if not isinstance(residue, (aa.Amino, na.Nucleic)):
                continue
            for atomname in residue.reference.map:
                if not atomname.startswith("H"):
                    continue
                if residue.has_atom(atomname):
                    continue
                if isinstance(residue, aa.CYS) and residue.ss_bonded and atomname == "HG":
                    continue

                # If this hydrogen is part of a tetrahedral group,
                # follow a different codepath
                if residue.rebuild_tetrahedral(atomname):
                    count += 1
                    continue

                # Otherwise use the standard quatfit methods
                coords = []
                refcoords = []

                refatomcoords = residue.reference.map[atomname].coords
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
                    coords.append(atom.coords)
                    refcoords.append(residue.reference.map[bond].coords)

                    # Exit if we have enough atoms
                    if len(coords) == 3:
                        break

                if len(coords) == 3:
                    newcoords = quat.find_coordinates(3, coords, refcoords, refatomcoords)
                    residue.create_atom(atomname, newcoords)
                    count += 1
                else:
                    _LOGGER.warning("Couldn't rebuild %s in %s!", atomname, residue)
        _LOGGER.debug(" Added %i hydrogen atoms.", count)

    def set_donors_acceptors(self):
        """Set the donors and acceptors within the protein"""
        for residue in self.residues:
            residue.set_donors_acceptors()

    def calculate_dihedral_angles(self):
        """Calculate the dihedral angle for every residue within the protein"""
        for residue in self.residues:
            if not isinstance(residue, aa.Amino):
                continue
            residue.dihedrals = []

            refangles = residue.reference.dihedrals
            for dihed in refangles:
                coords = []
                atoms = dihed.split()
                for i in range(4):
                    atomname = atoms[i]
                    if residue.has_atom(atomname):
                        coords.append(residue.get_atom(atomname).coords)

                if len(coords) == 4:
                    angle = util.dihedral(coords[0], coords[1], coords[2], coords[3])
                else:
                    angle = None

                residue.add_dihedral_angle(angle)

    def set_reference_distance(self):
        """Set the distance to the CA atom in the residue.

        This is necessary for determining which atoms are allowed to move during
        rotations.  Uses the shortest_path algorithm found in utilities.py.
        """
        for residue in self.residues:
            if not isinstance(residue, aa.Amino):
                continue

            # Initialize some variables
            map_ = {}
            caatom = residue.get_atom("CA")

            if caatom is None:
                text = "Cannot set references to %s without CA atom!"
                raise ValueError(text)

            # Set up the linked map
            for atom in residue.atoms:
                map_[atom] = atom.bonds

            # Run the algorithm
            for atom in residue.atoms:
                if atom.is_backbone:
                    atom.refdistance = -1
                elif residue.is_c_term and atom.name == "HO":
                    atom.refdistance = 3
                elif residue.is_n_term and (atom.name == "H3" or atom.name == "H2"):
                    atom.refdistance = 2
                else:
                    atom.refdistance = len(util.shortest_path(map_, atom, caatom)) - 1

    def remove_hydrogens(self):
        """Remove hydrogens from the protein."""
        for residue in self.residues:
            if not isinstance(residue, (aa.Amino, na.Nucleic)):
                continue
            for atom in residue.atoms[:]:
                if atom.is_hydrogen:
                    residue.remove_atom(atom.name)

    def assign_termini(self, chain, neutraln=False, neutralc=False):
        """Assign the termini for the given chain by looking at the start and
        end residues.

        Args:
            chain:  chain of protein
            neutraln:  whether to neutralize N-terminus
            neutralc:  whether to neutralize C-terminus
        """

        if len(chain.residues) == 0:
            text = "Error: chain \"%s\" has 0 residues!" % chain.chain_id
            raise IndexError(text)

        # Set the N-Terminus/ 5' Terminus
        res0 = chain.residues[0]
        if isinstance(res0, aa.Amino):
            res0.is_n_term = True
            if isinstance(res0, aa.PRO):
                self.apply_patch("NEUTRAL-NTERM", res0)
            elif neutraln:
                self.apply_patch("NEUTRAL-NTERM", res0)
            else:
                self.apply_patch("NTERM", res0)
        elif isinstance(res0, na.Nucleic):
            res0.is5term = True
            self.apply_patch("5TERM", res0)

        # Set the C-Terminus/ 3' Terminus
        reslast = chain.residues[-1]
        if isinstance(reslast, aa.Amino):
            reslast.is_c_term = True
            if neutralc:
                self.apply_patch("NEUTRAL-CTERM", reslast)
            else:
                self.apply_patch("CTERM", reslast)
        elif isinstance(reslast, na.Nucleic):
            reslast.is3term = True
            self.apply_patch("3TERM", reslast)
        else:
            for i in range(len(chain.residues)):
                resthis = chain.residues[-1 - i]
                if isinstance(resthis, aa.Amino):
                    resthis.is_c_term = True
                    if neutralc:
                        self.apply_patch("NEUTRAL-CTERM", resthis)
                    else:
                        self.apply_patch("CTERM", resthis)
                    break
                elif resthis.name in ["NH2", "NME"]:
                    break
                elif isinstance(resthis, na.Nucleic):
                    resthis.is3term = True
                    self.apply_patch("3TERM", resthis)
                    break

    def update_internal_bonds(self):
        """Update the internal bonding network using the reference objects in each atom."""
        for residue in self.residues:
            if isinstance(residue, (aa.Amino, aa.WAT, na.Nucleic)):
                for atom in residue.atoms:
                    if not atom.has_reference:
                        continue
                    for bond in atom.reference.bonds:
                        if not residue.has_atom(bond):
                            continue
                        bondatom = residue.get_atom(bond)
                        if bondatom not in atom.bonds:
                            atom.add_bond(bondatom)

    def update_bonds(self):
        """Update the bonding network of the protein.

        This happens in 3 steps:
            1.  Applying the PEPTIDE patch to all Amino residues so as to add
                reference for the N(i+1) and C(i-1) atoms
            2.  UpdateInternal_bonds for inter-residue linking
            3.  Set the links to the N(i+1) and C(i-1) atoms
        """
        # Apply the peptide patch
        for residue in self.residues:
            if isinstance(residue, aa.Amino):
                if residue.is_n_term or residue.is_c_term:
                    continue
                else:
                    self.apply_patch("PEPTIDE", residue)

        # Update all internal bonds
        self.update_internal_bonds()

        # Set the peptide bond pointers
        for chain in self.chains:
            for i in range(len(chain.residues) - 1):
                res1 = chain.residues[i]
                res2 = chain.residues[i + 1]
                if not isinstance(res1, aa.Amino) or not isinstance(res2, aa.Amino):
                    continue
                atom1 = res1.get_atom("C")
                atom2 = res2.get_atom("N")

                if atom1 is not None:
                    res2.peptide_c = atom1
                if atom2 is not None:
                    res1.peptide_n = atom2
                if atom1 is None or atom2 is None:
                    continue

                if util.distance(atom1.coords, atom2.coords) > PEPTIDE_DIST:
                    text = "Gap in backbone detected between %s and %s!" % \
                           (res1, res2)
                    _LOGGER.warning(text)
                    res2.peptide_c = None
                    res1.peptide_n = None

    def apply_patch(self, patchname, residue):
        """Apply a patch to the given residue.

        This is one of the key functions in PDB2PQR.  A similar function appears
        in definitions.py - that version is needed for residue level subtitutions
        so certain protonation states (i.e. CYM, HSE) are detectatble on input.

        This version looks up the particular patch name in the patch_map stored
        in the protein, and then applies the various commands to the reference
        and actual residue structures.

        See the inline comments for a more detailed explanation.

        Parameters
            patchname:  The name of the patch (string)
            residue:    The residue to apply the patch to (residue)
        """
        if patchname not in self.patch_map:
            raise KeyError("Unable to find patch %s!" % patchname)

        _LOGGER.debug('PATCH INFO: %s patched with %s', residue, patchname)

        if patchname == "PEPTIDE":
            newreference = residue.reference
        else:
            newreference = copy.deepcopy(residue.reference)

        patch = self.patch_map[patchname]

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
                residue.remove_atom(remove)
            if remove not in newreference.map:
                continue
            removebonds = newreference.map[remove].bonds
            del newreference.map[remove]
            for bond in removebonds:
                index = newreference.map[bond].bonds.index(remove)
                del newreference.map[bond].bonds[index]

        # Add the new dihedrals
        for dihedral_ in patch.dihedrals:
            newreference.dihedrals.append(dihedral_)

        # Point at the new reference
        residue.reference = newreference
        residue.patches.append(patchname)

        # Rename atoms as directed by patch
        for atom in residue.atoms:
            if atom.name in patch.altnames:
                residue.rename_atom(atom.name, patch.altnames[atom.name])

        # Replace each atom's reference with the new one
        for atomname in residue.map:
            if newreference.has_atom(atomname):
                atom = residue.get_atom(atomname)
                atom.reference = newreference.map[atomname]

    def update_ss_bridges(self):
        """Check for SS-bridge partners, and if present, set appropriate partners."""
        sg_partners = {}
        for residue in self.residues:
            if isinstance(residue, aa.CYS):
                atom = residue.get_atom("SG")
                if atom != None:
                    sg_partners[atom] = []

        for atom in sg_partners:
            for partner in sg_partners:
                if atom == partner or sg_partners[atom] != []:
                    continue
                dist = util.distance(atom.coords, partner.coords)
                if dist < BONDED_SS_LIMIT:
                    sg_partners[atom].append(partner)
                    sg_partners[partner].append(atom)

        for atom in sg_partners:
            res1 = atom.residue
            numpartners = len(sg_partners[atom])
            if numpartners == 1:
                partner = sg_partners[atom][0]
                res2 = partner.residue
                res1.ss_bonded = True
                res1.ss_bonded_partner = partner
                self.apply_patch("CYX", res1)
                _LOGGER.debug("%s - %s", res1, res2)
            elif numpartners > 1:
                error = "WARNING: %s has multiple potential " % res1
                error += "SS-bridge partners"
                _LOGGER.warning(error)
            elif numpartners == 0:
                _LOGGER.debug("%s is a free cysteine", res1)

    def update_residue_types(self):
        """Find the type of residue as notated in the Amino Acid definition"""
        for chain in self.chains:
            for residue in chain.residues:
                # TODO - why are we setting residue types to numeric values?
                name = residue.name
                if name in AA_NAMES:
                    residue.type = 1
                elif name == "WAT":
                    residue.type = 3
                elif name in NA_NAMES:
                    residue.type = 4
                else: # Residue is a ligand or unknown
                    residue.type = 2

    def apply_force_field(self, forcefield_):
        """Apply the forcefield to the atoms within the protein

        Parameters
            forcefield: forcefield object (forcefield)
        Returns
            hitlist:  list of atoms that were found in the forcefield (list)
            misslist:  list of atoms that were not found in the forcefield (list)
        """
        misslist = []
        hitlist = []
        for residue in self.residues:
            if isinstance(residue, (aa.Amino, aa.WAT, na.Nucleic)):
                resname = residue.ffname
            else:
                resname = residue.name

            for atom in residue.atoms:
                atomname = atom.name
                charge, radius = forcefield_.get_params(resname, atomname)
                if charge is not None and radius is not None:
                    atom.ffcharge = charge
                    atom.radius = radius
                    hitlist.append(atom)
                else:
                    misslist.append(atom)
        return hitlist, misslist

    def apply_name_scheme(self, forcefield_):
        """Apply the naming scheme of the given forcefield to the atoms within the protein

        Args:
            forcefield: forcefield object (forcefield)
        """
        for residue in self.residues:
            if isinstance(residue, (aa.Amino, aa.WAT, na.Nucleic)):
                resname = residue.ffname
            else:
                resname = residue.name

            for atom in residue.atoms:
                rname, aname = forcefield_.get_names(resname, atom.name)
                if resname not in ['LIG', 'WAT', 'ACE', 'NME'] and rname != None:
                    try:
                        if (residue.is_n_term or residue.is_c_term) and rname != residue.name:
                            rname = residue.name
                    except AttributeError:
                        pass
                if aname != None and rname != None:
                    atom.res_name = rname
                    atom.name = aname

    def apply_pka_values(self, force_field, ph, pkadic):
        """Apply calculated pKa values to assign titration states."""
        _LOGGER.info('Applying pKa values at a pH of %.2f:', ph)
        formatted_pkadict = pprint.pformat(pkadic)
        _LOGGER.debug("%s", formatted_pkadict)

        for residue in self.residues:
            if not isinstance(residue, aa.Amino):
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
                            _LOGGER.warning(warn)
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
                        self.apply_patch("ASH", residue)
                elif resname == "CYS" and ph >= value:
                    if force_field == "charmm":
                        warn = (key, "negative")
                        _LOGGER.warning(warn)
                    else:
                        self.apply_patch("CYM", residue)
                elif resname == "GLU" and ph < value:
                    if residue.is_c_term and force_field in ["amber", "tyl06", "swanson"]:
                        warn = (key, "Protonated at C-Terminal")
                        _LOGGER.warning(warn)
                    elif residue.is_n_term and force_field in ["amber", "tyl06", "swanson"]:
                        warn = (key, "Protonated at N-Terminal")
                        _LOGGER.warning(warn)
                    else:
                        self.apply_patch("GLH", residue)
                elif resname == "HIS" and ph < value:
                    self.apply_patch("HIP", residue)
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
                        self.apply_patch("LYN", residue)
                elif resname == "TYR" and ph >= value:
                    if force_field in ["charmm", "amber", "tyl06", "peoepb", "swanson"]:
                        warn = (key, "negative")
                        _LOGGER.warning(warn)
                    else:
                        self.apply_patch("TYM", residue)

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
        for residue in self.residues:
            reskey = (residue.res_seq, residue.chain_id, residue.ins_code)
            if reskey in hlist:
                hlist.remove(reskey)
                if isinstance(residue, aa.Amino):
                    residue.stateboolean = {'FIXEDSTATE': False}
                    _LOGGER.debug("Setting residue {:s} as fixed.".format(str(residue)))
                else:
                    err = "Matched residue {:s} but not subclass of Amino."
                    _LOGGER.warning(err.format(str(residue)))
        if len(hlist) > 0:
            err = "The following fixed residues were not matched (possible internal error): {:s}."
            _LOGGER.warning(err.format(str(hlist)))

    def create_residue(self, residue, resname):
        """Create a residue object.

        If the resname is a known residue type, try to make that specific
        object, otherwise just make a standard residue object.

        Args:
            residue:  A list of atoms (list)
            resname:  The name of the residue (string)
        Returns:
            residue:  The residue object (Residue)
        """
        try:
            refobj = self.definition.map[resname]
            if refobj.name != resname:
                klass = getattr(aa, refobj.name)
                residue = klass(residue, refobj)
                residue.reference = refobj
            else:
                klass = getattr(aa, resname)
                residue = klass(residue, refobj)
        except (KeyError, NameError):
            _LOGGER.debug("Parsing %s as new residue", resname)
            residue = residue_.Residue(residue)
        return residue

    def repair_heavy(self):
        """Repair all heavy atoms.

        Unfortunately the first time we get to an atom we might not be able to
        rebuild it - it might depend on other atoms to be rebuild first (think
        side chains).  As such a 'seenmap' is used to keep track of what we've
        already seen and subsequent attempts to rebuild the atom.
        """
        total_missing = self.num_missing_heavy
        if total_missing == 0:
            _LOGGER.warning("No heavy atoms need to be repaired.")
            return
        for residue in self.residues:
            if not isinstance(residue, (aa.Amino, na.Nucleic)):
                continue
            atomlist = list(residue.atoms)
            for atom in atomlist:
                atomname = atom.name
                if atomname in ["OP1", "OP2"] and residue.reference.has_atom("O1P") \
                    and residue.reference.has_atom("O2P"):
                    continue
                if not residue.reference.has_atom(atomname):
                    _LOGGER.warning("Extra atom %s in %s! - ", atomname, residue)
                    residue.remove_atom(atomname)
                    _LOGGER.warning("Deleted this atom.")

            missing = residue.missing
            if missing == []:
                continue
            seenmap = {}
            nummissing = len(missing)
            while len(missing) > 0:
                coords = []
                refcoords = []
                atomname = missing.pop(0)
                refatomcoords = residue.reference.map[atomname].coords
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
                    coords.append(atom.coords)
                    refcoords.append(residue.reference.map[bond].coords)
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
                    newcoords = quat.find_coordinates(3, coords, refcoords, refatomcoords)
                    residue.create_atom(atomname, newcoords)
                    _LOGGER.debug("Added %s to %s at coordinates", atomname, residue)
                    _LOGGER.debug(" %.3f %.3f %.3f", newcoords[0], newcoords[1], newcoords[2])

    def create_html_typemap(self, definition, outfilename):
        """Create an HTML typemap file at the desired location.

        If a type cannot be found for an atom a blank is listed.

        Args:
            definition: The definition objects.
            outfilename:  The name of the file to write (string)
        """
        # Cache the initial atom numbers
        numcache = {}
        for atom in self.atoms:
            numcache[atom] = atom.serial
        self.reserialize()

        amberff = forcefield.Forcefield("amber", definition, None)
        charmmff = forcefield.Forcefield("charmm", definition, None)

        with open(outfilename, "w") as file_:
            file_.write("<HTML>\n")
            file_.write("<HEAD>\n")
            file_.write("<TITLE>PQR Typemap (beta)</TITLE>\n")
            file_.write("</HEAD>\n")
            file_.write("<BODY>\n")
            file_.write(("<H3>This is a developmental page including the atom "
                         "type for the atoms in the PQR file.</H3><P>\n"))
            file_.write("<TABLE CELLSPACING=2 CELLPADDING=2 BORDER=1>\n")
            file_.write(("<tr><th>Atom Number</th><th>Atom Name</th><th>Residue "
                         "Name</th><th>Chain ID</th><th>AMBER Atom Type</th><th>"
                         "CHARMM Atom Type</th></tr>\n"))

            for atom in self.atoms:
                if isinstance(atom.residue, (aa.Amino, aa.WAT, na.Nucleic)):
                    resname = atom.residue.ffname
                else:
                    resname = atom.residue.name
                ambergroup = amberff.get_group(resname, atom.name)
                charmmgroup = charmmff.get_group(resname, atom.name)
                file_.write(("<tr><td>%s</td><td>%s</td><td>%s</td><td>%s</td>"
                             "<td>%s</td><td>%s</td></tr>\n") % (atom.serial, atom.name,
                                                                 resname, atom.chain_id,
                                                                 ambergroup, charmmgroup))

            file_.write("</table>\n")
            file_.write("</BODY></HTML>\n")

        # Return the original numbers back
        for atom in self.atoms:
            atom.serial = numcache[atom]

    def reserialize(self):
        """Generate new serial numbers for atoms in the protein"""
        count = 1
        for atom in self.atoms:
            atom.serial = count
            count += 1

    @property
    def atoms(self):
        """Return all Atom objects in list format

        Returns:
            atomlist:  List of Atom objects in the protein (list)
        """
        atomlist = []
        for chain in self.chains:
            for atom in chain.atoms:
                atomlist.append(atom)
        return atomlist

    @property
    def charge(self):
        """Get the total charge on the protein

        NOTE:  Since the misslist is used to identify incorrect charge
        assignments, this routine does not list the 3 and 5 termini of nucleic
        acid chains as having non-integer charge even though they are
        (correctly) non-integer.

        Returns:
            misslist:  List of residues with non-integer charges (list)
            charge:  The total charge on the protein (float)
        """
        charge = 0.0
        misslist = []
        for chain in self.chains:
            for residue in chain.residues:
                rescharge = residue.charge
                charge += rescharge
                if isinstance(residue, na.Nucleic):
                    if residue.is3term or residue.is5term:
                        continue
                if float("%i" % rescharge) != rescharge:
                    misslist.append(residue)
        return misslist, charge

    def __str__(self):
        output = []
        for chain in self.chains:
            output.append(chain.get_summary())
        return ' '.join(output)
