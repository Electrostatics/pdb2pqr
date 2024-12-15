"""Hydrogen optimization module for PDB2PQR.

This is an module for hydrogen optimization routines.

.. todo::
    This module has too many lines and should be simplified.

.. codeauthor:: Todd Dolinsky
.. codeauthor:: Jens Erik Nielsen
.. codeauthor:: Yong Huang
.. codeauthor:: Nathan Baker
"""

__author__ = "Todd Dolinsky, Jens Erik Nielsen, Yong Huang, Nathan Baker"
import logging
from xml import sax

from .. import aa, cells, io
from .. import definitions as defns
from .. import quatfit as quat
from .. import utilities as util
from ..config import HYD_DEF_PATH
from . import structures
from .structures import (
    HydrogenConformation,
    HydrogenDefinition,
    HydrogenHandler,
    PotentialBond,
)

_LOGGER = logging.getLogger(__name__)
_LOGGER.addFilter(io.DuplicateFilter())


#: A dictionary of the titration states for residues
TITRATION_DICT = {
    "ASH1c": "1",
    "ASH1t": "2",
    "ASH2c": "3",
    "ASH2t": "4",
    "ASP": "0",
    "GLH1c": "1",
    "GLH1t": "2",
    "GLH2c": "3",
    "GLH2t": "4",
    "GLU": "0",
    "ARG0": "1+2+3+4",
    "ARG": "1+2+3+4+5",
    "LYS": "1",
    "LYS0": "0",
    "TYR": "1",
    "TYR-": "0",
    "HSD": "1",
    "HSE": "2",
    "HSP": "1+2",
    "H3": "1",
    "H2": "2",
    "H3+H2": "1+2",
    "CTR01c": "1",
    "CTR01t": "2",
    "CTR02c": "3",
    "CTR02t": "4",
    "CTR-": "0",
}


def create_handler(hyd_path=HYD_DEF_PATH):
    """Create and populate a hydrogen handler.

    :param hyd_def_file:  path to hydrogen definition file
    :type hyd_def_file:  string or pathlib.Path object
    :return: HydrogenHandler object
    :rtype: HydrogenHandler
    """
    handler = HydrogenHandler()
    hyd_path = io.test_dat_file(hyd_path)
    with open(hyd_path) as hyd_file:
        sax.make_parser()
        sax.parseString(hyd_file.read(), handler)
    return handler


class HydrogenRoutines:
    """The main routines for hydrogen optimization.

    .. todo::
       This class really needs to be refactored.
    """

    def __init__(self, debumper, handler):
        """Initialize object.

        :param debumper:  Debump object
        :type debumper:  debump.Debump
        :param handler:  HydrogenHandler object
        :type handler:  HydrogenHandler
        """
        self.debumper = debumper
        self.biomolecule = debumper.biomolecule
        self.optlist = []
        self.atomlist = []
        self.resmap = {}
        self.hydrodefs = []
        self.map = handler.map

    def switchstate(self, states, amb, state_id):
        """Switch a residue to a new state by first removing all hydrogens.

        :param states:  the list of states
        :type states:  list
        :param amb:  the amibiguity to switch
        :type amb:  tup
        :param state_id:  the state id to switch to
        :type state_id:  int
        """
        if states == "pKa":
            return self.pka_switchstate(amb, state_id)

        if state_id > len(states):
            raise IndexError("Invalid State ID!")

        # First Remove all Hs
        residue = amb.residue
        hdef = amb.hdef
        for conf in hdef.conformations:
            hname = conf.hname
            boundname = conf.boundatom
            if residue.get_atom(hname) is not None:
                _LOGGER.debug(
                    f"Removing {residue.name} {residue.res_seq} {hname}"
                )
                residue.remove_atom(hname)
            residue.get_atom(boundname).hacceptor = 1
            residue.get_atom(boundname).hdonor = 0

        # Update the IntraBonds
        name = residue.name
        defresidue = self.debumper.aadef.get_residue(name)
        residue.updateIntraBonds(defresidue)

        # Now build appropriate atoms
        state = states[state_id]
        for conf in state:
            _LOGGER.debug(conf)
            refcoords = []
            defcoords = []
            defatomcoords = []
            if conf == ():
                continue  # Nothing to add
            hname = conf.hname
            for atom in conf.atoms:
                # print confatoms
                atomname = atom.name
                resatom = residue.get_atom(atomname)
                if atomname == hname:
                    defatomcoords = atom.coords
                elif resatom is not None:
                    refcoords.append(resatom.coords)
                    defcoords.append(atom.coords)
                else:
                    raise KeyError("Could not find necessary atom!")

            newcoords = quat.find_coordinates(
                3, refcoords, defcoords, defatomcoords
            )
            boundname = conf.boundatom
            residue.create_atom(hname, newcoords, "ATOM")
            residue.addDebumpAtom(residue.get_atom(hname))
            residue.get_atom(boundname).addIntraBond(hname)
            residue.get_atom(boundname).hacceptor = 0
            residue.get_atom(boundname).hdonor = 1
            # Setting the SybylType for the newly built H
            residue.get_atom(hname).sybyl_type = "H"
            # formal charge for PEOE_PB
            residue.get_atom(hname).formalcharge = 0.0
            # flag the added hydrogen
            residue.get_atom(hname).titratableH = True
            residue.get_atom(hname).addIntraBond(boundname)

    @classmethod
    def pka_switchstate(cls, amb, state_id_):
        """Switch a residue to a new state by first removing all hydrogens.

        This routine is used in pKa calculations only!

        :param amb:  the amibiguity to switch
        :type amb:  tup
        :param state_id_:  the state id to switch to
        :type state_id_:  int
        """
        titrationdict = TITRATION_DICT
        state_id = titrationdict[state_id_]
        state_id = state_id.split("+")
        new_state_id = [int(i) for i in state_id]
        residue = amb.residue
        hdef = amb.hdef
        for conf in hdef.conformations:
            hname = conf.hname
            boundname = conf.boundatom
            if residue.get_atom(hname) is not None:
                residue.remove_atom(hname)
            residue.get_atom(boundname).hacceptor = 1
            residue.get_atom(boundname).hdonor = 0

        # Update the IntraBonds
        for state_id in new_state_id:
            if state_id == 0:
                continue
            conf = hdef.conformations[state_id - 1]
            refcoords = []
            defcoords = []
            defatomcoords = []
            if conf == ():
                continue
            hname = conf.hname
            for atom in conf.atoms:
                if (
                    residue.is_n_term
                    and residue.name == "PRO"
                    and atom.name == "H"
                ):
                    atom.name = "CD"
                    atom.x = 1.874
                    atom.y = 0.862
                    atom.z = 1.306

            if not residue.rebuild_tetrahedral(hname):
                for atom in conf.atoms:
                    atomname = atom.name
                    resatom = residue.get_atom(atomname)
                    if atomname == hname:
                        defatomcoords = atom.coords
                    elif resatom is not None:
                        refcoords.append(resatom.coords)
                        defcoords.append(atom.coords)
                    else:
                        raise KeyError("Could not find necessary atom!")

                newcoords = quat.find_coordinates(
                    3, refcoords, defcoords, defatomcoords
                )
                residue.create_atom(hname, newcoords)

            boundname = conf.boundatom
            residue.get_atom(boundname).hacceptor = 0
            residue.get_atom(boundname).hdonor = 1

            # Setting the SybylType for the newly built H
            residue.get_atom(hname).sybyl_type = "H"

            # formal charge for PEOE_PB
            residue.get_atom(hname).formalcharge = 0.0

            # flag the added hydrogen
            residue.get_atom(hname).titratableH = True

        # Update intrabonds again
        if residue.is_n_term and residue.name == "PRO":
            for atom in residue.atoms:
                if atom.name == "H":
                    residue.remove_atom("H")
        residue.update_terminus_status()

    def cleanup(self):
        """Delete extra carboxylic atoms.

        If there are any extra carboxlyic ``*1`` atoms, delete them.
        This may occur when no optimization is chosen.
        """
        for residue in self.debumper.biomolecule.residues:
            if not isinstance(residue, aa.Amino):
                continue
            if residue.name == "GLH" or "GLH" in residue.patches:
                if residue.has_atom("HE1") and residue.has_atom("HE2"):
                    residue.remove_atom("HE1")
            elif residue.name == "ASH" or "ASH" in residue.patches:
                if residue.has_atom("HD1") and residue.has_atom("HD2"):
                    residue.remove_atom("HD1")

    def is_optimizeable(self, residue):
        """Check to see if the given residue is optimizeable.

        There are three ways to identify a residue:

        1. By name (i.e., HIS)
        2. By reference name - a PDB file HSP has a HIS reference name
        3. By patch - applied by :program:`propka`, terminal selection

        :param residue:  the residue in question
        :type residue:  Residue
        :return: None if not optimizeable, otherwise the OptimizationHolder
            instance that corresponds to the residue.
        :rtype: None or OptimizationHolder
        """
        optinstance = None
        if not isinstance(residue, (aa.Amino, aa.WAT)):
            return optinstance

        if residue.name in self.map:
            optinstance = self.map[residue.name]
        elif residue.reference.name in self.map:
            optinstance = self.map[residue.reference.name]
        else:
            for patch in residue.patches:
                if patch in self.map:
                    optinstance = self.map[patch]
                    break

        # If alcoholic, make sure the hydrogen is present
        if optinstance is not None and optinstance.opttype == "Alcoholic":
            atomname = next(iter(optinstance.map.keys()))
            if not residue.reference.has_atom(atomname):
                optinstance = None
        return optinstance

    def set_optimizeable_hydrogens(self):
        """Set any hydrogen listed in HYDROGENS.xml that is optimizeable.

        Used BEFORE hydrogen optimization to label atoms so that they won't be
        debumped - i.e. if SER HG is too close to another atom, don't debump
        but wait for optimization.

        .. note::
           This function should not be used if full optimization is not taking
           place.
        """
        for residue in self.biomolecule.residues:
            optinstance = self.is_optimizeable(residue)
            if optinstance is None:
                continue
            for atom in residue.atoms:
                if atom.name in optinstance.map:
                    atom.optimizeable = 1

    def initialize_full_optimization(self):
        """Initialize the full optimization.

        Detects all optimizeable donors and acceptors and sets the internal
        optlist.
        """
        # Do some setup
        self.debumper.cells = cells.Cells(5)
        self.debumper.cells.assign_cells(self.biomolecule)
        self.biomolecule.calculate_dihedral_angles()
        self.biomolecule.set_donors_acceptors()
        self.biomolecule.update_internal_bonds()
        self.biomolecule.set_reference_distance()
        self.optlist = []
        self.atomlist = []
        # First initialize the various types
        for residue in self.biomolecule.residues:
            optinstance = self.is_optimizeable(residue)
            if isinstance(residue, aa.Amino):
                residue.fixed = (
                    1 if False in residue.stateboolean.values() else 0
                )
            if optinstance is None:
                continue

            type_ = optinstance.opttype
            if residue.fixed != 1:
                klass = getattr(structures, type_)
                myobj = klass(residue, optinstance, self.debumper)
                self.atomlist += myobj.atomlist
                self.optlist.append(myobj)
                self.resmap[residue] = myobj
        _LOGGER.debug("Done.")

    def initialize_wat_optimization(self):
        """Initialize optimization for waters only.

        Detects all optimizeable donors and acceptors and sets the internal
        optlist.
        """
        _LOGGER.info("Initializing water bonding optimization...")
        # Do some setup
        self.debumper.cells = cells.Cells(5)
        self.debumper.cells.assign_cells(self.biomolecule)
        self.biomolecule.calculate_dihedral_angles()
        self.biomolecule.set_donors_acceptors()
        self.biomolecule.update_internal_bonds()
        self.biomolecule.set_reference_distance()
        self.optlist = []
        # First initialize the various types
        for residue in self.biomolecule.residues:
            optinstance = self.is_optimizeable(residue)
            if optinstance is None:
                continue
            type_ = optinstance.opttype
            if type_ == "Water":
                klass = getattr(structures, type_)
                myobj = klass(residue, optinstance, self.debumper)
                self.atomlist += myobj.atomlist
                self.optlist.append(myobj)
                self.resmap[residue] = myobj
        _LOGGER.debug("Done.")

    def optimize_hydrogens(self):
        """The main driver for the optimization.

        .. note::
           Should be called only after the optlist has been initialized.

        .. todo::
           Remove hard-coded :makevar:`progress` threshold and increment
           values.

        .. todo::
           This function needs to be simplified.
        """
        _LOGGER.debug("Optimization progress:")
        optlist = self.optlist
        connectivity = {}
        # Initialize the detection progress
        if len(optlist) == 0:
            return
        _LOGGER.debug("  Detecting potential hydrogen bonds")
        progress = 0.0
        increment = 1.0 / len(optlist)
        for obj in optlist:
            connectivity[obj] = []
            for atom in obj.atomlist:
                closeatoms = self.debumper.cells.get_near_cells(atom)
                for closeatom in closeatoms:
                    # Conditions for continuing
                    if atom.residue == closeatom.residue:
                        continue
                    if not (closeatom.hacceptor or closeatom.hdonor):
                        continue
                    if (
                        atom.hdonor
                        and not atom.hacceptor
                        and not closeatom.hacceptor
                    ):
                        continue
                    if (
                        atom.hacceptor
                        and not atom.hdonor
                        and not closeatom.hdonor
                    ):
                        continue
                    dist = util.distance(atom.coords, closeatom.coords)
                    if dist < 4.3:
                        residue = atom.residue
                        hbond = PotentialBond(atom, closeatom, dist)
                        # Store the potential bond
                        obj.hbonds.append(hbond)
                        # Keep track of connectivity
                        if closeatom in self.atomlist:
                            closeobj = self.resmap[closeatom.residue]
                            if closeobj not in connectivity[obj]:
                                connectivity[obj].append(closeobj)
            progress += increment
            while progress >= 0.0499:
                progress -= 0.05
        # Some residues might have no nearby hbonds - if so, place at
        # default state
        for obj in optlist:
            if len(obj.hbonds) == 0:
                if obj.residue.fixed:
                    continue
                _LOGGER.debug(
                    f"{obj.residue} has no nearby partners - fixing."
                )
                obj.finalize()
        # Determine the distinct networks
        networks = []
        seen = []
        for obj1 in optlist:
            if obj1.residue.fixed:
                continue
            if obj1 in seen:
                continue
            network = util.analyze_connectivity(connectivity, obj1)
            for obj2 in network:
                if obj2 not in seen:
                    seen.append(obj2)
            networks.append(network)
        # Initialize the output progress
        if len(networks) > 0:
            _LOGGER.debug("Optimizing hydrogen bonds")
            progress = 0.0
            increment = 1.0 / len(networks)
        # Work on the networks
        for network in networks:
            txt = ""
            for obj in network:
                txt += f"{obj}, "
            _LOGGER.debug(f"Starting network {txt[:-2]}")
            #  FIRST:  Only optimizeable to backbone atoms
            _LOGGER.debug("* Optimizeable to backbone *")
            hbondmap = {}
            for obj in network:
                for hbond in obj.hbonds:
                    if hbond.atom2 not in self.atomlist:
                        hbondmap[hbond] = hbond.dist
            hbondlist = util.sort_dict_by_value(hbondmap)
            hbondlist.reverse()
            for hbond in hbondlist:
                atom = hbond.atom1
                atom2 = hbond.atom2
                obj = self.resmap[atom.residue]

                if atom.residue.fixed:
                    continue
                if atom.hdonor:
                    obj.try_donor(atom, atom2)
                if atom.hacceptor:
                    obj.try_acceptor(atom, atom2)
            # SECOND:  Non-dual water Optimizeable to Optimizeable
            _LOGGER.debug("* Optimizeable to optimizeable *")
            hbondmap = {}
            seenlist = []
            for obj in network:
                for hbond in obj.hbonds:
                    if hbond.atom2 in self.atomlist and not isinstance(
                        hbond.atom1.residue, aa.WAT
                    ):
                        if not isinstance(hbond.atom2.residue, aa.WAT):
                            # Only get one hbond pair
                            if (hbond.atom2, hbond.atom1) not in seenlist:
                                hbondmap[hbond] = hbond.dist
                                seenlist.append((hbond.atom1, hbond.atom2))
            hbondlist = util.sort_dict_by_value(hbondmap)
            hbondlist.reverse()
            for hbond in hbondlist:
                atom = hbond.atom1
                atom2 = hbond.atom2
                obj1 = self.resmap[atom.residue]
                obj2 = self.resmap[atom2.residue]
                # Atoms may no longer exist if already optimized
                if not atom.residue.has_atom(atom.name):
                    continue
                if not atom2.residue.has_atom(atom2.name):
                    continue
                res = 0
                if atom.hdonor and atom2.hacceptor:
                    res = obj1.try_both(atom, atom2, obj2)
                if atom.hacceptor and atom2.hdonor and res == 0:
                    obj2.try_both(atom2, atom, obj1)
            # THIRD:  All water-water residues
            _LOGGER.debug("* Water to Water *")
            hbondmap = {}
            seenlist = []
            for obj in network:
                for hbond in obj.hbonds:
                    residue = hbond.atom1.residue
                    if isinstance(residue, aa.WAT) and isinstance(
                        hbond.atom2.residue, aa.WAT
                    ):
                        if (hbond.atom2, hbond.atom1) not in seenlist:
                            hbondmap[hbond] = hbond.dist
                            seenlist.append((hbond.atom1, hbond.atom2))
            hbondlist = util.sort_dict_by_value(hbondmap)
            hbondlist.reverse()
            for hbond in hbondlist:
                atom = hbond.atom1
                atom2 = hbond.atom2
                obj1 = self.resmap[atom.residue]
                obj2 = self.resmap[atom2.residue]
                res = 0
                if atom.hdonor and atom2.hacceptor:
                    res = obj1.try_both(atom, atom2, obj2)
                if atom.hacceptor and atom2.hdonor and res == 0:
                    obj2.try_both(atom2, atom, obj1)
            # FOURTH: Complete all residues
            for obj in network:
                obj.complete()
            # STEP 5:  Update progress meter
            progress += 100.0 * increment
            while progress >= 5.0:
                progress -= 5.0

    def parse_hydrogen(self, res, topo):
        """Parse a list of lines in order to make a hydrogen definition.

        This is the current definition:
        ``Name Ttyp  A R # Stdconf HT Chi OPTm``

        .. todo::
           The type of the :makevar:`res` appears to be incorrect.

        .. todo::
           This function is too long and needs to be simplified.

        :param res:  the lines to parse (list)
        :type res:  unknown
        :param topo:  Topology object
        :type topo:  pdb2pqr.topology.Topology
        :return:  the hydrogen definition object
        :rtype: HydrogenDefinition

        """
        name = self.map[res].name
        opttype = self.map[res].opttype
        optangle = self.map[res].optangle
        map_ = self.map[res].map
        mydef = HydrogenDefinition(name, opttype, optangle, map_)
        patch_map = []
        atoms = []
        refatoms = []
        conformernames = []
        refmap = {}
        titrationstatemap = {}
        tautomermap = {}
        conformermap = {}
        atommap = {}
        ntrmap = {}
        hmap = {}
        nonhmap = {}
        # reference map from TOPOLOGY.xml
        for res_ in topo.residues:
            refmap[res_.name] = res_.reference
            for atom in refmap[res_.name].atoms:
                atommap[res_.name, atom.name] = atom
            for titrationstate in res_.titration_states:
                titrationstatemap[titrationstate.name] = titrationstate
                for tautomer in titrationstate.tautomers:
                    tautomermap[tautomer.name] = tautomer
                    for conformer in tautomer.conformers:
                        conformermap[conformer.name] = conformer
        if name == "CYS":
            _ = refmap["CYS"]
            atoms = ["HG"]
            refatoms = ["SG", "CB"]
        elif name == "HIS":
            _ = refmap["HIS"]
            atoms = ["HD1", "HE2"]
            for atom in atoms:
                refatoms = ["ND1", "CG", "CE1"]
        elif name == "LYS":
            _ = self.debumper.biomolecule.reference_map[name]
            patch_map = self.debumper.biomolecule.patch_map["LYN"]
            atoms = patch_map.remove
            refatoms = ["HZ1", "HZ2", "NZ"]
        elif name == "TYR":
            _ = self.debumper.biomolecule.reference_map[name]
            patch_map = self.debumper.biomolecule.patch_map["TYM"]
            atoms = patch_map.remove
            refatoms = ["OH", "CZ", "CE2"]
        elif name == "WAT":
            _ = self.debumper.biomolecule.reference_map[name]
            patch_map = self.debumper.biomolecule.patch_map["HOH"]
            atoms = ["H1", "H2"]
            refatoms = None
        elif name == "NTR":
            ntrmap = {}  # map for N-TERM
            for tautomer in titrationstatemap["NTER"].tautomers:
                for conformer in tautomermap[tautomer.name].conformers:
                    for conformeradds in conformermap[
                        conformer.name
                    ].conformer_adds:
                        for atom in conformeradds.atoms:
                            ntrmap[atom.name] = atom
            atoms = ["H3", "H2"]
            refatoms = ["CA", "H", "N"]
        elif name == "CTR":
            hmap = {}  # map for h atoms
            nonhmap = {}  # map for refatoms
            conformernames = []
            for tautomer in titrationstatemap["CTER"].tautomers:
                for conformer in tautomermap[tautomer.name].conformers:
                    for conformeradds in conformermap[
                        conformer.name
                    ].conformer_adds:
                        for atom in conformeradds.atoms:
                            nonhmap[atom.name] = atom
            for tautomer in titrationstatemap["CTER0"].tautomers:
                for conformer in tautomermap[tautomer.name].conformers:
                    conformernames.append(conformer.name)
                    for conformeradds in conformermap[
                        conformer.name
                    ].conformer_adds:
                        for atom in conformeradds.atoms:
                            hmap[conformer.name, atom.name] = atom
            atoms = ["HO"]
            refatoms = ["O", "C", "OXT"]
        elif name in ["SER", "GLN", "THR", "ARG", "ASN"]:
            _ = refmap[name]
            if name == "SER":
                atoms = ["HG"]
                refatoms = ["OG", "CB"]
            elif name == "GLN":
                atoms = ["HE21"]
                refatoms = ["NE2"]
            elif name == "THR":
                atoms = ["HG1"]
                refatoms = ["OG1", "CB"]
            elif name == "ARG":
                atoms = ["HH11", "HH12", "HH21", "HH22", "HE"]
                for atom in atoms:
                    refatoms = ["NH1", "NH2", "CZ"]
            elif name == "ASN":
                atoms = ["HD21"]
                refatoms = ["ND2"]
        elif name == "ASH":
            hmap = {}  # map for h atoms
            nonhmap = {}  # map for refatoms
            conformernames = []
            _ = refmap["ASP"]
            for tautomer in titrationstatemap["ASH"].tautomers:
                for conformer in tautomermap[tautomer.name].conformers:
                    for conformeradds in conformermap[
                        conformer.name
                    ].conformer_adds:
                        for atom in conformeradds.atoms:
                            hmap[conformer.name, atom.name] = atom
                            conformernames.append(conformer.name)
            atoms = ["HD1", "HD2"]
            refatoms = ["OD1", "CG", "OD2"]
        elif name == "GLH":
            hmap = {}  # map for h atoms
            nonhmap = {}  # map for refatoms
            conformernames = []
            _ = refmap["GLU"]
            for tautomer in titrationstatemap["GLH"].tautomers:
                for conformer in tautomermap[tautomer.name].conformers:
                    for conformeradds in conformermap[
                        conformer.name
                    ].conformer_adds:
                        for atom in conformeradds.atoms:
                            hmap[conformer.name, atom.name] = atom
                            conformernames.append(conformer.name)
            atoms = ["HE1", "HE2"]
            refatoms = ["OE1", "CD", "OE2"]
        else:
            patch_map = self.debumper.biomolecule.patch_map[name]
            atoms = list(patch_map.map.keys())
            atoms.sort()
        if name in ["NTR"]:
            bondlength = 1.0
            for atom in atoms:
                hname = atom
                x = ntrmap[hname].x
                y = ntrmap[hname].y
                z = ntrmap[hname].z
                bondatom = ntrmap[hname].bonds[0]
                myconf = HydrogenConformation(hname, bondatom, bondlength)
                atom = defns.DefinitionAtom(hname, x, y, z)
                myconf.add_atom(atom)

                # TODO - lots of arbitrary undefined numbers in this section
                for atom_ in refatoms:
                    if atom_ == "N":
                        natom = defns.DefinitionAtom(atom_, 1.201, 0.847, 0.0)
                        myconf.add_atom(natom)
                    elif atom_ == "CA":
                        caatom = defns.DefinitionAtom(atom_, 0.0, 0.0, 0.0)
                        myconf.add_atom(caatom)
                    elif atom_ == "H":
                        caatom = defns.DefinitionAtom(
                            atom_, 1.201, 1.847, 0.000
                        )
                        myconf.add_atom(caatom)
                mydef.add_conf(myconf)
        elif name in ["CTR"]:
            bondlength = 1.0
            for conformer in conformernames:
                for atom in atoms:
                    hname = atom
                    x = hmap[conformer, hname].x
                    y = hmap[conformer, hname].y
                    z = hmap[conformer, hname].z
                    bondatom = hmap[conformer, hname].bonds[0]
                    myconf = HydrogenConformation(hname, bondatom, bondlength)
                    atom = defns.DefinitionAtom(hname, x, y, z)
                    myconf.add_atom(atom)

                    # TODO - the following code is almost nonsensical
                    for atom_ in refatoms:
                        if atom_ == "C":
                            catom = defns.DefinitionAtom(
                                atom_, -1.250, 0.881, 0.000
                            )
                            myconf.add_atom(catom)
                        else:
                            atomname = atom_
                            x = nonhmap[atom_].x
                            y = nonhmap[atom_].y
                            z = nonhmap[atom_].z
                            atom2 = defns.DefinitionAtom(atomname, x, y, z)
                            myconf.add_atom(atom2)
                    mydef.add_conf(myconf)

        elif name in ["ASH", "GLH"]:
            for conformer in conformernames:
                for atom in atoms:
                    hname = atom
                    if ("1" in conformer and "1" in atom) or (
                        "2" in conformer and "2" in atom
                    ):
                        x = hmap[conformer, hname].x
                        y = hmap[conformer, hname].y
                        z = hmap[conformer, hname].z
                        bondatom = hmap[conformer, hname].bonds[0]
                        bondlength = 1.0
                        myconf = HydrogenConformation(
                            hname, bondatom, bondlength
                        )
                        atom = defns.DefinitionAtom(hname, x, y, z)
                        myconf.add_atom(atom)

                        for atom_ in refatoms:
                            atomname = atom_
                            refresname = ""
                            if name == "ASH":
                                refresname = "ASP"
                            elif name == "GLH":
                                refresname = "GLU"
                            x = atommap[refresname, atom_].x
                            y = atommap[refresname, atom_].y
                            z = atommap[refresname, atom_].z
                            atom2 = defns.DefinitionAtom(atomname, x, y, z)
                            myconf.add_atom(atom2)
                        mydef.add_conf(myconf)
        elif name not in ["WAT"]:
            bondlength = 1.0
            for atom in atoms:
                hname = atom
                x = atommap[name, hname].x
                y = atommap[name, hname].y
                z = atommap[name, hname].z
                bondatom = atommap[name, hname].bonds[0]
                myconf = HydrogenConformation(hname, bondatom, bondlength)
                atom = defns.DefinitionAtom(hname, x, y, z)
                myconf.add_atom(atom)
                if refatoms is not None:
                    if name == "HIS" and atom.name == "HE2":
                        refatoms = ["NE2", "CE1", "CD2"]
                    if name == "ARG" and atom.name == "HE":
                        refatoms = ["NE", "CZ", "NH1"]
                    # FIXME: 2020/07/06 intendo -  the "atom" is reused in
                    # the outer for loop and ambiguous
                    for atom in refatoms:
                        atomname = atom
                        x = atommap[name, atomname].x
                        y = atommap[name, atomname].y
                        z = atommap[name, atomname].z
                        atom = defns.DefinitionAtom(atomname, x, y, z)
                        myconf.add_atom(atom)
                    mydef.add_conf(myconf)
        return mydef

    def read_hydrogen_def(self, topo):
        """Read the hydrogen definition file

        :param topo:  Topology object
        :type topo:  Topology object
        """
        self.hydrodefs = []
        for mapping in self.map:
            res = mapping
            mydef = self.parse_hydrogen(res, topo)
            self.hydrodefs.append(mydef)
            res = ""
