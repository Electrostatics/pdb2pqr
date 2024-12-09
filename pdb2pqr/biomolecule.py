"""The biomolecule object used in PDB2PQR and associated methods.

.. todo:: This module should be broken into separate files.

Authors:  Todd Dolinsky, Yong Huang
"""

import copy
import logging
import pprint
import string

from . import aa, forcefield, io, na, pdb
from . import quatfit as quat
from . import residue as residue_
from . import structures as struct
from . import utilities as util
from .config import BONDED_SS_LIMIT, PEPTIDE_DIST, RNA_MAPPING

_LOGGER = logging.getLogger(__name__)
_LOGGER.addFilter(io.DuplicateFilter())


class Biomolecule:
    """Biomolecule class.

    This class represents the parsed PDB, and provides a hierarchy of
    information - each Biomolecule object contains a list of Chain objects as
    provided in the PDB file.  Each Chain then contains its associated list of
    Residue objects, and each Residue contains a list of Atom objects,
    completing the hierarchy.
    """

    def __init__(self, pdblist, definition):
        """Initialize using parsed PDB file

        :param pdblist:  list of objects from :mod:`pdb` from lines of PDB file
        :type pdblist:  list
        :param definition:  topology definition object
        :type definition:  Definition
        """
        self.chainmap = {}
        self.chains = []
        self.residues = []
        self.definition = definition
        self.pdblist = pdblist
        chain_dict = {}
        previous_atom = None
        residue = []
        num_models = 0
        num_chains = 1
        count = 0
        for record in pdblist:  # Find number of chains
            if isinstance(record, pdb.TER):
                num_chains += 1
        for record in pdblist:
            if isinstance(record, (pdb.ATOM, pdb.HETATM)):
                if (
                    record.chain_id == ""
                    and num_chains > 1
                    and record.res_name
                    not in [
                        "WAT",
                        "HOH",
                    ]
                ):
                    # Assign a chain ID
                    try:
                        record.chain_id = (
                            string.ascii_uppercase
                            + string.ascii_lowercase
                            + string.digits
                        )[count]
                    except IndexError:
                        raise Exception(
                            "Too many chains exist in biomolecule. "
                            "Consider preparing subsets."
                        )

                chain_id = record.chain_id
                res_seq = record.res_seq
                ins_code = record.ins_code
                if previous_atom is None:
                    previous_atom = record
                if chain_id not in chain_dict:
                    my_chain = struct.Chain(chain_id)
                    chain_dict[chain_id] = my_chain
                if (
                    res_seq != previous_atom.res_seq
                    or ins_code != previous_atom.ins_code
                    or chain_id != previous_atom.chain_id
                ):
                    my_residue = self.create_residue(
                        residue, previous_atom.res_name
                    )
                    chain_dict[previous_atom.chain_id].add_residue(my_residue)
                    residue = []
                residue.append(record)
                previous_atom = record
            elif isinstance(record, pdb.END):
                my_residue = self.create_residue(
                    residue, previous_atom.res_name
                )
                chain_dict[previous_atom.chain_id].add_residue(my_residue)
                residue = []
            elif isinstance(record, pdb.MODEL):
                num_models += 1
                if residue == []:
                    continue
                if num_models > 1:
                    my_residue = self.create_residue(
                        residue, previous_atom.res_name
                    )
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
        """Return number of missing biomolecular heavy atoms in structure.

        :return:  number of missing heavy atoms in structure
        :rtype:  int
        """
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
                if refatomname == "O1P" and residue.has_atom("OP1"):
                    continue
                if refatomname == "O2P" and residue.has_atom("OP2"):
                    continue
                if not residue.has_atom(refatomname):
                    _LOGGER.warning(
                        f"Missing atom {refatomname} in residue {residue}"
                    )
                    natom += 1
                    residue.missing.append(refatomname)
        return natom

    @property
    def num_heavy(self):
        """Return number of biomolecular heavy atoms in structure.

        .. todo::
           Figure out if this is redundant with
           :func:`Biomolecule.num_bio_atoms`

        .. note::
           Includes hydrogens (but those are stripped off eventually)

        :return:  number of heavy atoms
        :rtype:  int
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
                if refatomname == "O1P" and residue.has_atom("OP1"):
                    continue
                if refatomname == "O2P" and residue.has_atom("OP2"):
                    continue
                natom += 1
        return natom

    @property
    def reference_map(self):
        """Return definition reference map.

        :return:  definition reference map
        :rtype:  dict
        """
        _LOGGER.error(type(self.definition.map))
        return self.definition.map

    @property
    def patch_map(self):
        """Return definition patch maps.

        :return:  definition patch maps
        :rtype:  list
        """
        return self.definition.patches

    @property
    def num_bio_atoms(self):
        """Return the number of ATOM (not HETATM) records in the biomolecule.

        :return:  number of ATOM records
        :rtype:  int
        """
        return sum(1 for atom in self.atoms if atom.type == "ATOM")

    def set_hip(self):
        """Set all HIS states to HIP."""
        for residue in self.residues:
            if isinstance(residue, aa.HIS):
                self.apply_patch("HIP", residue)

    def set_termini(self, *, neutraln=False, neutralc=False):
        """Set the termini for a protein.

        First set all known termini by looking at the ends of the chain. Then
        examine each residue, looking for internal chain breaks.

        .. todo::  This function needs to be cleaned and simplified

        :param neutraln:  indicate whether N-terminus is neutral
        :type neutraln:  bool
        :param neutralc:  indicate whether C-terminus is neutral
        :type neutralc:  bool
        """
        # First assign the known termini
        chain = None
        for chain in self.chains:
            self.assign_termini(chain, neutraln=neutraln, neutralc=neutralc)
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
                    if residue.has_atom("OXT") and not residue.is_c_term:
                        fixflag = 1
                elif (
                    isinstance(residue, na.Nucleic)
                    and (residue.has_atom("H3T") or residue.name.endswith("3"))
                    and not residue.is3term
                ):
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
                        message = (
                            "Warning: Reusing chain id: " + chainid[0] + ""
                        )
                        _LOGGER.warning(message)
                    # Make a new chain with these residues
                    newchain = struct.Chain(chainid[0])
                    self.chainmap[chainid] = newchain
                    self.chains.insert(ch_num, newchain)
                    for res in reslist:
                        newchain.add_residue(res)
                        chain.residues.remove(res)
                        res.set_chain_id(chainid[0])
                    self.assign_termini(
                        chain, neutraln=neutraln, neutralc=neutralc
                    )
                    self.assign_termini(
                        newchain, neutraln=neutraln, neutralc=neutralc
                    )
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
                message = "Warning: Reusing chain id: " + chainid[0]
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

        See :mod:`aa` for residue-specific functions.
        """
        for residue in self.residues:
            if isinstance(residue, (aa.Amino, na.Nucleic)):
                residue.set_state()

    def add_hydrogens(self, hlist=None):
        """Add the hydrogens to the biomolecule.

        This requires either the rebuild_tetrahedral function for tetrahedral
        geometries or the standard quatfit methods.  These methods use three
        nearby bonds to rebuild the atom; the closer the bonds, the more
        accurate the results.  As such the peptide bonds are used when
        available.
        """
        count = 0
        for residue in self.residues:
            if not isinstance(residue, (aa.Amino, na.Nucleic)):
                continue

            reskey = (residue.res_seq, residue.chain_id, residue.ins_code)
            if hlist is not None and reskey in hlist:
                continue

            for atomname in residue.reference.map:
                if not atomname.startswith("H"):
                    continue
                if residue.has_atom(atomname):
                    continue
                if isinstance(residue, aa.CYS) and (
                    residue.ss_bonded and atomname == "HG"
                ):
                    continue
                if hasattr(residue, "rebuild_tetrahedral"):
                    # If this hydrogen is part of a tetrahedral group,
                    # follow a different codepath
                    if residue.rebuild_tetrahedral(atomname):
                        count += 1
                        continue
                else:
                    _LOGGER.warning(
                        "Tetrahedral hydrogen reconstruction not available "
                        "for nucleic acids. Some hydrogens may be missing "
                        "(if so, this is a bug)."
                    )
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
                    newcoords = quat.find_coordinates(
                        3, coords, refcoords, refatomcoords
                    )
                    residue.create_atom(atomname, newcoords)
                    count += 1
                else:
                    _LOGGER.warning(
                        f"Couldn't rebuild {atomname} in {residue}!"
                    )
        _LOGGER.debug(f" Added {count} hydrogen atoms.")

    def set_donors_acceptors(self):
        """Set the donors and acceptors within the biomolecule."""
        for residue in self.residues:
            residue.set_donors_acceptors()

    def calculate_dihedral_angles(self):
        """Calculate dihedral angles for every residue in the biomolecule."""
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
                    angle = util.dihedral(
                        coords[0], coords[1], coords[2], coords[3]
                    )
                else:
                    angle = None
                residue.add_dihedral_angle(angle)

    def set_reference_distance(self):
        """Set the distance to the CA atom in the residue.

        This is necessary for determining which atoms are allowed to move
        during rotations.  Uses the :func:`shortest_path` algorithm found in
        :mod:`utilities`.

        :raises ValueError:  if shortest path cannot be found (e.g., if the
            atoms are not connected)
        """
        for residue in self.residues:
            if not isinstance(residue, aa.Amino):
                continue
            # Initialize some variables
            map_ = {}
            caatom = residue.get_atom("CA")
            if caatom is None:
                # TODO: What does the %s mean? Is it the residue name?
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
                elif residue.is_n_term and (atom.name in ("H3", "H2")):
                    atom.refdistance = 2
                else:
                    path = util.shortest_path(map_, atom, caatom)
                    if path is not None:
                        atom.refdistance = len(path) - 1
                    else:
                        raise ValueError(
                            "Found gap in biomolecule structure for atom "
                            f"{atom}"
                        )

    def remove_hydrogens(self):
        """Remove hydrogens from the biomolecule."""
        for residue in self.residues:
            if not isinstance(residue, (aa.Amino, na.Nucleic)):
                continue
            for atom in residue.atoms[:]:
                if atom.is_hydrogen:
                    residue.remove_atom(atom.name)

    def assign_termini(self, chain, *, neutraln=False, neutralc=False):
        """Assign the termini for the given chain.

        Assignment made by looking at the start and end residues.

        :param chain:  chain of biomolecule
        :type chain:  Chain
        :param neutraln:  indicate whether to neutralize N-terminus
        :type neutraln:  bool
        :param neutralc:  indicate whether to neutralize C-terminus
        :type neutralc:  bool
        """
        if len(chain.residues) == 0:
            text = f'Error: chain "{chain.chain_id}" has 0 residues!'
            raise IndexError(text)

        res0 = chain.residues[0]
        reslast = chain.residues[-1]
        # Check if chain is cyclic. Amide distance ranges between 1.325 - 1.346
        if "N" in res0.map and "C" in reslast.map:
            dist = util.distance(res0.map["N"].coords, reslast.map["C"].coords)
            if dist < 1.35:
                # If the chain is cyclic, don't apply termini.
                return

        # Set the N-Terminus/ 5' Terminus
        if isinstance(res0, aa.Amino):
            res0.is_n_term = True
            # If N is bonded to more than one heavy atom switch to neutral-nterm
            heavy_n_bonds = []
            if "N" in res0.map:
                heavy_n_bonds = [
                    a for a in res0.map["N"].bonds if a.name[0] != "H"
                ]

            if neutraln or len(heavy_n_bonds) > 1:
                self.apply_patch("NEUTRAL-NTERM", res0)
            else:
                self.apply_patch("NTERM", res0)
        elif isinstance(res0, na.Nucleic):
            res0.is5term = True
            self.apply_patch("5TERM", res0)
        # Set the C-Terminus/ 3' Terminus
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
        """Update the internal bonding network.

        Update using the reference objects in each atom.
        """
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
        """Update the bonding network of the biomolecule.

        This happens in 3 steps:

        1. Apply the PEPTIDE patch to all Amino residues to add
           reference for the N(i+1) and C(i-1) atoms
        2. UpdateInternal_bonds for inter-residue linking
        3. Set the links to the N(i+1) and C(i-1) atoms
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
                if not isinstance(res1, aa.Amino) or not isinstance(
                    res2, aa.Amino
                ):
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
                    text = (
                        f"Gap in backbone detected between {res1} and {res2}!"
                    )
                    _LOGGER.warning(text)
                    res2.peptide_c = None
                    res1.peptide_n = None

    def apply_patch(self, patchname: str, residue: residue_.Residue):
        """Apply a patch to the given residue.

        This is one of the key functions in PDB2PQR.  A similar function
        appears in :mod:`definitions` - that version is needed for residue
        level subtitutions so certain protonation states (i.e. CYM, HSE) are
        detectatble on input.

        This version looks up the particular patch name in the patch_map
        stored in the biomolecule, and then applies the various commands to the
        reference and actual residue structures.

        See the inline comments for a more detailed explanation.

        :param patchname:  the name of the patch
        :type patchname:  str
        :param residue:  the residue to apply the patch to
        :type residue:  Residue
        """
        if patchname not in self.patch_map:
            raise KeyError(f"Unable to find patch {patchname}!")
        _LOGGER.debug(f"PATCH INFO: {residue} patched with {patchname}")
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
        """Check and set SS-bridge partners."""
        sg_partners = {}
        for residue in self.residues:
            if isinstance(residue, aa.CYS):
                atom = residue.get_atom("SG")
                if atom is not None:
                    sg_partners[atom] = []
        for atom in sg_partners:
            for partner, value in sg_partners.items():
                if atom == partner or sg_partners[atom] != []:
                    continue
                dist = util.distance(atom.coords, partner.coords)
                if dist < BONDED_SS_LIMIT:
                    sg_partners[atom].append(partner)
                    value.append(atom)
        for atom in sg_partners:
            res1 = atom.residue
            numpartners = len(sg_partners[atom])
            if numpartners == 1:
                partner = sg_partners[atom][0]
                res2 = partner.residue
                res1.ss_bonded = True
                res1.ss_bonded_partner = partner
                self.apply_patch("CYX", res1)
                _LOGGER.debug(f"{res1} - {res2}")
            elif numpartners > 1:
                error = f"WARNING: {res1} has multiple potential "
                error += "SS-bridge partners"
                _LOGGER.warning(error)
            elif numpartners == 0:
                _LOGGER.debug(f"{res1} is a free cysteine")

    def apply_force_field(self, forcefield_):
        """Apply the forcefield to the atoms within the biomolecule.

        :param forcefield_:  forcefield object
        :type forcefield_:  Forcefield
        :return:  (list of atoms that were found in the forcefield,
            list of atoms that were not found in the forcefield)
        :rtype:  (list, list)
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
            charge_err = util.noninteger_charge(residue.charge)
            if charge_err:
                _LOGGER.warning(
                    f"Residue {residue} has non-integer charge: {charge_err}. "
                )
        return hitlist, misslist

    def apply_name_scheme(self, forcefield_):
        """Apply the naming scheme of the given forcefield.

        :param forcefield_:  forcefield object
        :type forcefield_:  Forcefield
        """
        for residue in self.residues:
            if isinstance(residue, (aa.Amino, aa.WAT, na.Nucleic)):
                resname = residue.ffname
            else:
                resname = residue.name
            for atom in residue.atoms:
                rname, aname = forcefield_.get_names(resname, atom.name)
                if resname not in ["LIG", "WAT", "ACE", "NME"] and (
                    rname is not None
                ):
                    try:
                        if (residue.is_n_term or residue.is_c_term) and (
                            rname != residue.name
                        ):
                            if len(rname) == 4 and rname[0] in ("C", "N"):
                                # Remove the C/N prefix to keep the protonation state of the residue
                                # in the terminal residues
                                rname = rname[1:]
                            elif rname.startswith("NEUTRAL-"):
                                # Remove the NEUTRAL-C and NEUTRAL-N prefixes to keep protonation state
                                rname = rname[9:]
                            else:
                                # This is the old code which will overwrite the protonation state
                                # of the residue but have wrong hydrogens. Not sure if needed but
                                # nobody seems to know.
                                rname = residue.name
                    except AttributeError:
                        pass
                if aname is not None and rname is not None:
                    atom.res_name = rname
                    atom.name = aname

    def apply_pka_values(self, force_field, ph, pkadic):
        """Apply calculated pKa values to assign titration states.

        :param force_field:  force field name (determines naming of residues)
        :type force_field:  str
        :param ph:  pH value
        :type ph:  float
        :param pkadic:  dictionary of pKa values for residues
        :type pkadic:  dict
        """
        _LOGGER.info(f"Applying pKa values at a pH of {ph:.2f}:")
        formatted_pkadict = pprint.pformat(pkadic)
        _LOGGER.debug(f"{formatted_pkadict}")
        for residue in self.residues:
            if not isinstance(residue, aa.Amino):
                continue
            resname = residue.name
            resnum = residue.res_seq
            chain_id = residue.chain_id
            if residue.is_n_term:
                key = f"N+  {resnum:>3} {chain_id}"
                key = key.strip()
                if key in pkadic:
                    value = pkadic[key]
                    del pkadic[key]
                    if ph >= value:
                        if force_field in [
                            "amber",
                            "charmm",
                            "tyl06",
                            "peoepb",
                            "swanson",
                        ]:
                            warn = f"N-terminal {key} neutral"
                            _LOGGER.warning(warn)
                        else:
                            self.apply_patch("NEUTRAL-NTERM", residue)
            if residue.is_c_term:
                key = f"C-  {resnum:>3} {chain_id}"
                key = key.strip()
                if key in pkadic:
                    value = pkadic[key]
                    del pkadic[key]
                    if ph < value:
                        if force_field in [
                            "amber",
                            "charmm",
                            "tyl06",
                            "peoepb",
                            "swanson",
                        ]:
                            warn = f"C-terminal {key} neutral"
                            _LOGGER.warning(warn)
                        else:
                            self.apply_patch("NEUTRAL-CTERM", residue)
            key = f"{resname} {resnum} {chain_id}"
            key = key.strip()
            if key in pkadic:
                value = pkadic[key]
                del pkadic[key]
                if resname == "ARG" and ph >= value:
                    if force_field == "parse":
                        self.apply_patch("AR0", residue)
                        _LOGGER.warning(
                            "Neutral arginines are very rare. Please "
                            "double-check your setup."
                        )
                    else:
                        warn = (key, "neutral")
                        _LOGGER.warning(warn)
                elif resname == "ASP" and ph < value:
                    if residue.is_c_term and force_field in [
                        "amber",
                        "tyl06",
                        "swanson",
                    ]:
                        warn = (key, "Protonated at C-Terminal")
                        _LOGGER.warning(warn)
                    elif residue.is_n_term and force_field in [
                        "amber",
                        "tyl06",
                        "swanson",
                    ]:
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
                    if residue.is_c_term and force_field in [
                        "amber",
                        "tyl06",
                        "swanson",
                    ]:
                        warn = (key, "Protonated at C-Terminal")
                        _LOGGER.warning(warn)
                    elif residue.is_n_term and force_field in [
                        "amber",
                        "tyl06",
                        "swanson",
                    ]:
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
                    elif (
                        force_field in ["amber", "tyl06", "swanson"]
                        and residue.is_c_term
                    ):
                        warn = (key, "neutral at C-Terminal")
                        _LOGGER.warning(warn)
                    elif force_field == "tyl06" and residue.is_n_term:
                        warn = (key, "neutral at N-Terminal")
                        _LOGGER.warning(warn)
                    else:
                        self.apply_patch("LYN", residue)
                elif resname == "TYR" and ph >= value:
                    if force_field in [
                        "charmm",
                        "amber",
                        "tyl06",
                        "peoepb",
                        "swanson",
                    ]:
                        warn = (key, "negative")
                        _LOGGER.warning(warn)
                    else:
                        self.apply_patch("TYM", residue)
        if len(pkadic) > 0:
            warn = (
                "PDB2PQR could not identify the following residues and residue"
                " numbers as returned by PROPKA or PDB2PKA"
            )
            _LOGGER.warning(warn)
            for item in pkadic:
                text = f"             {item}"
                _LOGGER.warning(text)

    def hold_residues(self, hlist):
        """Set fixed state of specified residues.

        :param hlist:  list of (res_seq, chainid, ins_code) specifying the
            residues for altering fixed state status.
        :type hlist:  [(str, str, str)]
        """
        if not hlist:
            return
        for residue in self.residues:
            reskey = (residue.res_seq, residue.chain_id, residue.ins_code)
            if reskey in hlist:
                hlist.remove(reskey)
                if isinstance(residue, aa.Amino):
                    residue.stateboolean = {"FIXEDSTATE": False}
                    _LOGGER.debug(f"Setting residue {residue} as fixed.")
                else:
                    err = f"Matched residue {residue} but not subclass of "
                    err += "Amino."
                    _LOGGER.warning(err)
        if len(hlist) > 0:
            err = (
                "The following fixed residues were not matched (possible "
                f"internal error): {hlist}."
            )
            _LOGGER.warning(err)

    def create_residue(self, residue, resname):
        """Create a residue object.

        If the resname is a known residue type, try to make that specific
        object, otherwise just make a standard residue object.

        :param residue:  a list of atoms
        :type residue:  list
        :param resname:  the name of the residue
        :type resname:  str
        :return:  the residue object
        :rtype:  Residue
        """
        if (resname not in self.definition.map) and (resname in RNA_MAPPING):
            resname = RNA_MAPPING[resname]
        try:
            refobj = self.definition.map[resname]
            if refobj.name != resname:
                try:
                    klass = getattr(aa, refobj.name)
                except AttributeError:
                    klass = getattr(na, refobj.name)
                residue = klass(residue, refobj)
                residue.reference = refobj
            else:
                try:
                    klass = getattr(aa, resname)
                except AttributeError:
                    klass = getattr(na, resname)
                residue = klass(residue, refobj)
            residue.rename_residue(resname)
        except (KeyError, NameError):
            _LOGGER.warning(
                f"Unable to find amino or nucleic acid definition "
                f"for {resname}.  Parsing as new residue."
            )
            residue = residue_.Residue(residue)
        return residue

    def repair_heavy(self):
        """Repair all heavy atoms.

        Unfortunately the first time we get to an atom we might not be able to
        rebuild it - it might depend on other atoms to be rebuild first (think
        side chains).  As such a 'seenmap' is used to keep track of what we've
        already seen and subsequent attempts to rebuild the atom.

        :raise ValueError:  missing atoms prevent reconstruction
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
                if atomname in ["O1P", "OP1"] and residue.has_atom("OP1"):
                    continue
                if atomname in ["O2P", "OP2"] and residue.has_atom("OP2"):
                    continue
                if not residue.reference.has_atom(atomname):
                    _LOGGER.warning(f"Extra atom {atomname} in {residue}! - ")
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
                        missing_str = " ".join(missing)
                        err = (
                            "Too few atoms present to reconstruct or cap "
                            f"residue {residue} in structure! This error is "
                            "generally caused by missing backbone atoms in "
                            "this biomolecule; you must use an external "
                            "program to complete gaps in the biomolecule "
                            f"backbone. Heavy atoms missing from {residue}:  "
                            f"{missing_str}"
                        )
                        raise ValueError(err)
                else:  # Rebuild the atom
                    newcoords = quat.find_coordinates(
                        3, coords, refcoords, refatomcoords
                    )
                    residue.create_atom(atomname, newcoords)
                    _LOGGER.info(
                        f"Added atom {atomname} to residue {residue} at "
                        f"coordinates {newcoords[0]:.3f}, {newcoords[1]:.3f}, "
                        f"{newcoords[2]:.3f}"
                    )

    def create_html_typemap(self, definition, outfilename):
        """Create an HTML typemap file at the desired location.

        If a type cannot be found for an atom a blank is listed.

        :param definition:  the definition objects.
        :type definition:  Definition
        :param outfilename:  the name of the file to write
        :type outfilename:  str
        """
        # Cache the initial atom numbers
        numcache = {atom: atom.serial for atom in self.atoms}
        self.reserialize()
        amberff = forcefield.Forcefield("amber", definition, None)
        charmmff = forcefield.Forcefield("charmm", definition, None)
        with open(outfilename, "w") as file_:
            file_.write("<HTML>\n")
            file_.write("<HEAD>\n")
            file_.write("<TITLE>PQR Typemap (beta)</TITLE>\n")
            file_.write("</HEAD>\n")
            file_.write("<BODY>\n")
            file_.write(
                "<H3>This is a developmental page including the atom "
                "type for the atoms in the PQR file.</H3><P>\n"
            )
            file_.write("<TABLE CELLSPACING=2 CELLPADDING=2 BORDER=1>\n")
            file_.write(
                "<tr><th>Atom Number</th><th>Atom Name</th><th>Residue "
                "Name</th><th>Chain ID</th><th>AMBER Atom Type</th><th>"
                "CHARMM Atom Type</th></tr>\n"
            )
            for atom in self.atoms:
                if isinstance(atom.residue, (aa.Amino, aa.WAT, na.Nucleic)):
                    resname = atom.residue.ffname
                else:
                    resname = atom.residue.name
                ambergroup = amberff.get_group(resname, atom.name)
                charmmgroup = charmmff.get_group(resname, atom.name)
                file_.write(
                    f"<tr><td>{atom.serial}</td>"
                    f"<td>{atom.name}</td>"
                    f"<td>{resname}</td>"
                    f"<td>{atom.chain_id}</td>"
                    f"<td>{ambergroup}</td>"
                    f"<td>{charmmgroup}</td></tr>\n"
                )
            file_.write("</table>\n")
            file_.write("</BODY></HTML>\n")
        # Return the original numbers back
        for atom in self.atoms:
            atom.serial = numcache[atom]

    def reserialize(self):
        """Generate new serial numbers for atoms in the biomolecule."""
        for count, atom in enumerate(self.atoms, start=1):
            atom.serial = count

    @property
    def atoms(self):
        """Return all Atom objects in list format.

        :return:  all atom objects
        :rtype:  [Atom]
        """
        atomlist = []
        for chain in self.chains:
            for atom in chain.atoms:
                atomlist.append(atom)
        return atomlist

    @property
    def charge(self):
        """Get the total charge on the biomolecule

        .. todo::
           Since the misslist is used to identify incorrect charge
           assignments, this routine does not list the 3 and 5 termini of
           nucleic acid chains as having non-integer charge even though they
           are (correctly) non-integer.

        :return:  (list of residues with non-integer charges,
            the total charge on the biomolecule)
        :rtype:  (list, float)
        """
        charge = 0.0
        misslist = []
        for chain in self.chains:
            for residue in chain.residues:
                rescharge = residue.charge
                charge += rescharge
                if isinstance(residue, na.Nucleic) and (
                    residue.is3term or residue.is5term
                ):
                    continue
                if util.noninteger_charge(rescharge):
                    misslist.append(residue)
        return misslist, charge

    def __str__(self):
        output = [chain.get_summary() for chain in self.chains]
        return " ".join(output)
