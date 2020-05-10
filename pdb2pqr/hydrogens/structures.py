"""Topology-related classes for hydrogen optimization."""
import logging
from xml import sax
from .. import definitions as defns
from .. import utilities as util
from .. import aa
from ..config import ANGLE_CUTOFF, DIST_CUTOFF
from . import optimize


_LOGGER = logging.getLogger(__name__)


class Flip(optimize.Optimize):
    """The holder for optimization of flippable residues."""

    def __init__(self, residue, optinstance, routines):
        """Initialize a potential flip.

        Rather than flipping the given residue back and forth, take each atom
        that would be flipped and pre-flip it, making a new *FLIP atom in its
        place.

        Args:
            residue:      The residue to flip (residue)
            optinstance:  The optimization instance containing information about
                what to optimize
        """
        # Initialize some variables
        self.optinstance = optinstance
        self.residue = residue
        self.routines = routines
        self.atomlist = []
        self.hbonds = []
        map_ = {}

        # Get all moveable names for this angle/residue pair
        dihedral = optinstance.optangle
        pivot = dihedral.split()[2]
        moveablenames = residue.get_moveable_names(pivot)
        # HO in CTERM shouldn't be in the list of flip atoms
        if residue.is_c_term:
            newmoveablenames = []
            for name in moveablenames:
                if name == "HO":
                    pass
                else:
                    newmoveablenames.append(name)
            moveablenames = newmoveablenames

        # Cache current coordinates
        for name in moveablenames:
            atom = residue.get_atom(name)
            map_[name] = atom.coords

        # Flip the residue about the angle
        anglenum = residue.reference.dihedrals.index(dihedral)
        if anglenum == -1:
            raise ValueError("Unable to find Flip dihedral angle!")
        newangle = 180.0 + residue.dihedrals[anglenum]
        self.routines.set_dihedral_angle(residue, anglenum, newangle)

        # Create new atoms at cached positions
        for name in map_:
            newname = "%sFLIP" % name
            residue.create_atom(newname, map_[name])
            newatom = residue.get_atom(newname)
            self.routines.cells.add_cell(newatom)

            # Set the bonds
            newatom.reference = residue.reference.map[name]
            for bond in newatom.reference.bonds:
                newbond = "%sFLIP" % bond
                if residue.has_atom(newbond):
                    bondatom = residue.map[newbond]
                    if bondatom not in newatom.bonds:
                        newatom.bonds.append(bondatom)
                    if newatom not in bondatom.bonds:
                        bondatom.bonds.append(newatom)

                # And connect back to the existing structure
                newbond = bond
                if residue.has_atom(newbond):
                    bondatom = residue.map[newbond]
                    if bondatom not in newatom.bonds:
                        newatom.bonds.append(bondatom)
                    if newatom not in bondatom.bonds:
                        bondatom.bonds.append(newatom)

        residue.set_donors_acceptors()

        # Add to the optimization list
        for name in moveablenames:

            # Get the atom
            atom = residue.get_atom(name)
            if not atom.is_hydrogen:
                if atom.hdonor or atom.hacceptor:
                    self.atomlist.append(atom)

            # And the FLIP
            atom = residue.get_atom("%sFLIP" % name)
            if not atom.is_hydrogen:
                if atom.hdonor or atom.hacceptor:
                    self.atomlist.append(atom)

        # Special case: Neutral unassigned HIS can be acceptors
        if isinstance(residue, aa.HIS):
            if residue.name == "HIS" and len(residue.patches) == 1:
                for atom in self.atomlist:
                    if atom.name.startswith("N"):
                        atom.hacceptor = 1

    def try_both(self, donor, acc, accobj):
        """Called when both the donor and acceptor are optimizeable.

        If one is fixed, we only need to try one side.  Otherwise first try to
        satisfy the donor - if that's succesful, try to satisfy the acceptor.
        An undo may be necessary if the donor is satisfied and the acceptor
        isn't.
        """
        # If one residue if fixed, use other functions
        if donor.residue.fixed:
            if accobj.try_acceptor(acc, donor):
                return 1
            return 0
        if acc.residue.fixed:
            if self.try_donor(donor, acc):
                return 1
            return 0

        _LOGGER.debug("Working on %s %s (donor) to %s %s (acceptor)",
                      donor.residue, donor.name, acc.residue, acc.name)
        if self.is_hbond(donor, acc):
            if accobj.try_acceptor(acc, donor):
                self.fix_flip(donor)
                donor.hacceptor = 0
                _LOGGER.debug("NET BOND SUCCESSFUL!")
                return 1
            return 0
        return 0

    def try_donor(self, donor, acc):
        """The main driver for adding a hydrogen to an optimizeable residue."""
        residue = self.residue

        # Do some error checking
        if not acc.hacceptor:
            return 0

        _LOGGER.debug("Working on %s %s (donor) to %s %s (acceptor)",
                      donor.residue, donor.name, acc.residue, acc.name)

        if self.is_hbond(donor, acc):
            residue.fixed = donor.name
            self.fix_flip(donor)
            donor.hacceptor = 0
            return 1
        return 0

    def try_acceptor(self, acc, donor):
        """The main driver for adding an LP to an optimizeable residue. """
        residue = acc.residue

        # Do some error checking
        if not donor.hdonor:
            return 0

        _LOGGER.debug("Working on %s %s (acceptor) to %s %s (donor)",
                      acc.residue, acc.name, donor.residue, donor.name)
        if self.is_hbond(donor, acc):
            residue.fixed = acc.name
            self.fix_flip(acc)
            acc.hdonor = 0
            return 1
        return 0

    def fix_flip(self, bondatom):
        """Called if a hydrogen bond has been found using the bondatom.
        If bondatom is *FLIP, remove all * atoms, otherwise remove all *FLIP atoms.
        """

        # Initialize some variables
        atomlist = []
        residue = bondatom.residue
        for atom in residue.atoms:
            atomlist.append(atom)

        # Set a flag to see whether to delete the FLIPs or not
        flag = 0
        if bondatom.name.endswith("FLIP"):
            flag = 1

        dstr = "fix_flip called for residue {:s}, bondatom {:s} and flag {:d}"
        _LOGGER.debug(dstr.format(str(residue), str(bondatom), flag))
        residue.wasFlipped = (flag == 0)

        # Delete the appropriate atoms
        for atom in atomlist:
            atomname = atom.name
            if atomname.endswith("FLIP") and flag: # Delete the other list
                if residue.has_atom(atomname[:-4]):
                    self.routines.cells.remove_cell(residue.get_atom(atomname[:-4]))
                    residue.remove_atom(atomname[:-4])
            elif atomname.endswith("FLIP"):  # Delete the flip
                self.routines.cells.remove_cell(atom)
                residue.remove_atom(atomname)
            else: continue

        residue.fixed = 1

    def finalize(self):
        """Finalizes a flippable back to its original state - since the original
        atoms are now *FLIP, it deletes the * atoms and renames the *FLIP atoms back to *.
        """
        residue = self.residue
        if residue.fixed:
            return
        atomlist = []
        for atom in residue.atoms:
            atomlist.append(atom)
        for atom in atomlist:
            if atom.name.endswith("FLIP"):
                self.routines.cells.remove_cell(atom)
                residue.remove_atom(atom.name[:-4])
                residue.rename_atom(atom.name, atom.name[:-4])
        residue.fixed = 1

    def complete(self):
        """Complete the flippable residue optimization.

        Call the finalize function, and then rename all FLIP atoms back to their
        standard names.
        """
        residue = self.residue
        self.finalize()
        # Rename all *FLIP atoms
        for atom in residue.atoms:
            atomname = atom.name
            if atomname.endswith("FLIP"):
                residue.rename_atom(atomname, atomname[:-4])


class Generic(optimize.Optimize):
    """Generic optimization class"""

    def __init__(self, residue, optinstance, routines):
        self.optinstance = optinstance
        self.residue = residue
        self.routines = routines
        self.hbonds = []
        self.atomlist = []

    def finalize(self):
        """Initialize some variable and pass?"""
        pass

    def complete(self):
        """If not already fixed, finalize"""
        if not self.residue.fixed:
            self.finalize()


class Alcoholic(optimize.Optimize):
    """The class for alcoholic residues"""

    def __init__(self, residue, optinstance, routines):
        """Initialize the alcoholic class by removing the alcoholic hydrogen if
        it exists."""

        # Initialize some variables
        self.optinstance = optinstance
        self.residue = residue
        self.routines = routines
        self.atomlist = []
        self.hbonds = []

        name = list(optinstance.map.keys())[0]
        self.hname = name

        bondname = residue.reference.get_atom(name).bonds[0]
        self.atomlist.append(residue.get_atom(bondname))
        if residue.has_atom(name):
            atom = residue.get_atom(name)
            self.routines.cells.remove_cell(atom)
            residue.remove_atom(name)

    def try_both(self, donor, acc, accobj):
        """Called when both the donor and acceptor are optimizeable.

           If one is fixed, we only need to try one side.  Otherwise
           first try to satisfy the donor - if that's succesful,
           try to satisfy the acceptor.  An undo may be necessary
           if the donor is satisfied and the acceptor isn't.
        """
        # If one residue if fixed, use other functions

        residue = donor.residue
        if donor.residue.fixed:
            if accobj.try_acceptor(acc, donor):
                return 1
            return 0
        if acc.residue.fixed:
            if self.try_donor(donor, acc):
                return 1
            return 0

        if self.try_donor(donor, acc):
            if accobj.try_acceptor(acc, donor):
                _LOGGER.debug("NET BOND SUCCESSFUL!")
                return 1
            residue.remove_atom(self.hname)
            _LOGGER.debug("REMOVED NET HBOND")
            return 0
        return 0

    def try_donor(self, donor, acc):
        """The main driver for adding a hydrogen to an optimizeable residue."""
        residue = self.residue

        # Do some error checking
        if not acc.hacceptor:
            return 0

        # Get the name of the atom to add
        newname = self.hname
        if residue.has_atom(newname):
            return 0

        _LOGGER.debug("Working on %s %s (donor) to %s %s (acceptor)",
                      donor.residue, donor.name, acc.residue, acc.name)

        # Act depending on the number of bonds
        if len(donor.bonds) == 1: # No H or LP attached
            self.make_atom_with_one_bond_h(donor, newname)
            newatom = donor.residue.get_atom(newname)
            return self.try_single_alcoholic_h(donor, acc, newatom)
        elif len(donor.bonds) == 2:
            loc1, loc2 = self.get_positions_with_two_bonds(donor)
            return self.try_positions_with_two_bonds_h(donor, acc, newname, loc1, loc2)
        elif len(donor.bonds) == 3:
            loc = self.get_position_with_three_bonds(donor)
            return self.try_positions_three_bonds_h(donor, acc, newname, loc)

        return 0

    def try_acceptor(self, acc, donor):
        """The main driver for adding an LP to an optimizeable residue."""
        residue = acc.residue

        # Do some error checking
        if not donor.hdonor:
            return 0

        # Get the name of the LP to add
        if residue.has_atom("LP2"):
            return 0
        elif residue.has_atom("LP1"):
            newname = "LP2"
        else: newname = "LP1"

        _LOGGER.debug("Working on %s %s (acceptor) to %s %s (donor)",
                      acc.residue, acc.name, donor.residue, donor.name)

        # Act depending on the number of bonds
        if len(acc.bonds) == 1: # No H or LP attached
            self.make_atom_with_one_bond_lp(acc, newname)
            newatom = acc.residue.get_atom(newname)
            return self.try_single_alcoholic_lp(acc, donor, newatom)
        elif len(acc.bonds) == 2:
            loc1, loc2 = self.get_positions_with_two_bonds(acc)
            return self.try_positions_with_two_bonds_lp(acc, donor, newname, loc1, loc2)
        elif len(acc.bonds) == 3:
            loc = self.get_position_with_three_bonds(acc)
            return self.try_positions_three_bonds_lp(acc, donor, newname, loc)

        return 0

    def finalize(self):
        """Finalize an alcoholic residue.
        Try to minimize conflict with nearby atoms by building away from them.
        Called when LPs are still present so as to account for their bonds.
        """

        # Initialize some variables
        residue = self.residue
        atom = self.atomlist[0]

        # Conditions for return
        addname = self.hname
        if residue.fixed:
            return
        if residue.has_atom(addname):
            return

        if len(atom.bonds) == 1:

            # Initialize variables
            pivot = atom.bonds[0]
            #bestdist = 0.0
            bestcoords = []
            bestenergy = 999.99

            # Add atom and debump
            self.make_atom_with_one_bond_h(atom, addname)
            newatom = residue.get_atom(addname)
            self.routines.cells.add_cell(newatom)

            for _ in range(18):
                residue.rotate_tetrahedral(pivot, atom, 20.0)

                closeatoms = self.routines.cells.get_near_cells(atom)
                energy = 0.0
                for catom in closeatoms:
                    energy += self.get_pair_energy(atom, catom) + self.get_pair_energy(catom, atom)
                if energy < bestenergy:
                    bestenergy = energy
                    bestcoords = newatom.coords

            if bestcoords != []:
                newatom.x = bestcoords[0]
                newatom.y = bestcoords[1]
                newatom.z = bestcoords[2]

        elif len(atom.bonds) == 2:
            loc1, loc2 = self.get_positions_with_two_bonds(atom)
            residue.create_atom(addname, loc1)
            newatom = residue.get_atom(addname)
            self.routines.cells.add_cell(newatom)

            # Debump residue if necessary by trying the other location
            closeatoms = self.routines.cells.get_near_cells(atom)
            energy1 = 0.0
            for catom in closeatoms:
                energy1 += self.get_pair_energy(atom, catom) + self.get_pair_energy(catom, atom)

            # Place at other location
            self.routines.cells.remove_cell(newatom)
            newatom.x = loc2[0]
            newatom.y = loc2[1]
            newatom.z = loc2[2]
            self.routines.cells.add_cell(newatom)

            energy2 = 0.0
            for catom in closeatoms:
                energy2 += self.get_pair_energy(atom, catom) + self.get_pair_energy(catom, atom)

            # If this is worse, switch back
            if energy2 > energy1:
                self.routines.cells.remove_cell(newatom)
                newatom.x = loc1[0]
                newatom.y = loc1[1]
                newatom.z = loc1[2]
                self.routines.cells.add_cell(newatom)

        elif len(atom.bonds) == 3:
            loc = self.get_position_with_three_bonds(atom)
            residue.create_atom(addname, loc)
            self.routines.cells.add_cell(residue.get_atom(addname))

    def complete(self):
        """Complete an alcoholic optimization.
        Call finalize(), and then remove all extra LP atoms.
        """
        # Initialize some variables
        residue = self.residue
        self.finalize()
        residue.fixed = 1

        # Remove all LP atoms
        atomlist = []
        for atom in residue.atoms:
            atomlist.append(atom)
        for atom in atomlist:
            if atom.name.startswith("LP"):
                residue.remove_atom(atom.name)


class Water(optimize.Optimize):
    """The class for water residues"""

    def __init__(self, residue, optinstance, routines):
        """Initialize the water optimization class"""
        self.optinstance = optinstance
        self.residue = residue
        self.routines = routines
        self.hbonds = []

        oxatom = residue.get_atom("O")
        if oxatom is None:
            raise KeyError("Unable to find oxygen atom in %s!" % residue)

        oxatom.hdonor = 1
        oxatom.hacceptor = 1
        self.atomlist = [oxatom]

    def try_both(self, donor, acc, accobj):
        """Called when both the donor and acceptor are optimizeable.
        If one is fixed, we only need to try one side.  Otherwise
        first try to satisfy the donor - if that's succesful,
        try to satisfy the acceptor.  An undo may be necessary
        if the donor is satisfied and the acceptor isn't.
        """
        # If one residue if fixed, use other functions
        residue = donor.residue

        if donor.residue.fixed:
            if accobj.try_acceptor(acc, donor):
                return 1
            return 0
        if acc.residue.fixed:
            if self.try_donor(donor, acc):
                return 1
            return 0

        if self.try_donor(donor, acc):
            if accobj.try_acceptor(acc, donor):
                _LOGGER.debug("NET BOND SUCCESSFUL!")
                return 1
            # We need to undo what we did to the donor
            _LOGGER.debug("REMOVED NET HBOND")
            if residue.has_atom("H2"):
                residue.remove_atom("H2")
            elif residue.has_atom("H1"):
                residue.remove_atom("H1")
            return 0
        return 0

    def try_acceptor(self, acc, donor):
        """The main driver for adding an LP to an optimizeable residue."""
        residue = acc.residue

        # Do some error checking
        if not donor.hdonor:
            return 0

        # Get the name of the LP to add
        if residue.has_atom("LP2"):
            return 0
        elif residue.has_atom("LP1"):
            newname = "LP2"
        else:
            newname = "LP1"

        _LOGGER.debug("Working on %s %s (acceptor) to %s %s (donor)",
                      acc.residue, acc.name, donor.residue, donor.name)

        # Act depending on the number of bonds
        if len(acc.bonds) == 0:

            if self.is_hbond(donor, acc):

                # Find the best donor hydrogen and use that
                bestdist = util.distance(acc.coords, donor.coords)
                for donorh in donor.bonds:
                    dist = util.distance(acc.coords, donorh.coords)
                    if dist < bestdist:
                        bestdist = dist

                # Point the LP to the best H
                # TODO - this looks like a bug...
                # donorh ends up being the last item in donor.bonds
                # This may be fixed by setting a best_donorh to go with bestdist
                # and using best_donorh in the function below
                self.make_atom_with_no_bonds(acc, donorh, newname)
                _LOGGER.warning("The best donorH was not picked (BUG?).")
                _LOGGER.debug("Added %s to %s", newname, acc.residue)
                return 1
            return 0

        elif len(acc.bonds) == 1: # No H or LP attached
            _LOGGER.debug("Trying to add %s to %s with one bond", newname, acc.residue)
            self.make_water_with_one_bond(acc, newname)
            newatom = acc.residue.get_atom(newname)
            return self.try_single_alcoholic_lp(acc, donor, newatom)
        elif len(acc.bonds) == 2:
            _LOGGER.debug("Trying to add %s to %s with two bonds", newname, acc.residue)
            loc1, loc2 = self.get_positions_with_two_bonds(acc)
            return self.try_positions_with_two_bonds_lp(acc, donor, newname, loc1, loc2)
        elif len(acc.bonds) == 3:
            _LOGGER.debug("Trying to add %s to %s with three bonds", newname, acc.residue)
            loc = self.get_position_with_three_bonds(acc)
            return self.try_positions_three_bonds_lp(acc, donor, newname, loc)
        return 0

    def try_donor(self, donor, acc):
        """The main driver for adding a hydrogen to an optimizeable residue."""
        residue = self.residue

        # Do some error checking
        if not acc.hacceptor:
            return 0

        # Get the name of the atom to add
        if residue.has_atom("H2"):
            return 0
        elif residue.has_atom("H1"):
            newname = "H2"
        else:
            newname = "H1"

        _LOGGER.debug("Working on %s %s (donor) to %s %s (acceptor)",
                      donor.residue, donor.name, acc.residue, acc.name)

        # Act depending on the number of bonds
        if len(donor.bonds) == 0:
            self.make_atom_with_no_bonds(donor, acc, newname)
            if self.is_hbond(donor, acc):
                return 1
            self.routines.cells.remove_cell(residue.get_atom(newname))
            residue.remove_atom(newname)
            return 0
        if len(donor.bonds) == 1:
            self.make_water_with_one_bond(donor, newname)
            newatom = donor.residue.get_atom(newname)
            return self.try_single_alcoholic_h(donor, acc, newatom)
        elif len(donor.bonds) == 2:
            loc1, loc2 = self.get_positions_with_two_bonds(donor)
            return self.try_positions_with_two_bonds_h(donor, acc, newname, loc1, loc2)
        elif len(donor.bonds) == 3:
            loc = self.get_position_with_three_bonds(donor)
            return self.try_positions_three_bonds_h(donor, acc, newname, loc)
        return 0

    def finalize(self):
        """Finalize a water residue.
        Try to minimize conflict with nearby atoms by building away from them.
        Called when LPs are still present so as to account for their bonds."""

        residue = self.residue

        # Conditions for return
        if residue.fixed:
            _LOGGER.debug("Residue %s already fixed", residue)
            return
        if residue.has_atom("H2"):
            _LOGGER.debug("Residue %s already has H2", residue)
            return

        atom = residue.get_atom("O")
        if not residue.has_atom("H1"):
            addname = "H1"
        else:
            addname = "H2"

        _LOGGER.debug("Finalizing %s by adding %s (%i current O bonds)",
                      residue, addname, len(atom.bonds))

        if len(atom.bonds) == 0:

            newcoords = []

            # Build hydrogen away from closest atom
            closeatom = self.routines.get_closest_atom(atom)
            if closeatom != None:
                vec = util.subtract(atom.coords, closeatom.coords)
                dist = util.distance(atom.coords, closeatom.coords)

                for i in range(3):
                    newcoords.append(vec[i]/dist + atom.coords[i])

            else:
                newcoords = util.add(atom.coords, [1.0, 0.0, 0.0])

            residue.create_atom(addname, newcoords)
            self.routines.cells.add_cell(residue.get_atom(addname))

            self.finalize()

        elif len(atom.bonds) == 1:

            # Initialize variables
            pivot = atom.bonds[0]
            bestdist = 0.0
            bestcoords = []

            # Add atom and debump
            self.make_water_with_one_bond(atom, addname)
            newatom = residue.get_atom(addname)
            self.routines.cells.add_cell(newatom)

            for i in range(18):
                residue.rotate_tetrahedral(pivot, atom, 20.0)
                nearatom = self.routines.get_closest_atom(newatom)

                # If there is no closest conflict, continue
                if nearatom is None:
                    continue

                dist = util.distance(nearatom.coords, newatom.coords)

                if dist > bestdist:
                    bestdist = dist
                    bestcoords = newatom.coords

            if bestcoords != []:
                newatom.x = bestcoords[0]
                newatom.y = bestcoords[1]
                newatom.z = bestcoords[2]

            if addname == "H1":
                self.finalize()

            residue.fixed = 1

        elif len(atom.bonds) == 2:

            loc1, loc2 = self.get_positions_with_two_bonds(atom)
            residue.create_atom(addname, loc1)
            newatom = residue.get_atom(addname)
            self.routines.cells.add_cell(newatom)

            # Debump residue if necessary by trying the other location
            nearatom = self.routines.get_closest_atom(newatom)
            if nearatom != None:
                dist1 = util.distance(newatom.coords, nearatom.coords)

                # Place at other location
                self.routines.cells.remove_cell(atom)
                newatom.x = loc2[0]
                newatom.y = loc2[1]
                newatom.z = loc2[2]
                self.routines.cells.add_cell(atom)

                nearatom = self.routines.get_closest_atom(newatom)
                if nearatom != None:

                    # If this is worse, switch back
                    if util.distance(newatom.coords, nearatom.coords) < dist1:
                        self.routines.cells.remove_cell(atom)
                        newatom.x = loc1[0]
                        newatom.y = loc1[1]
                        newatom.z = loc1[2]
                        self.routines.cells.add_cell(atom)

            if addname == "H1":
                self.finalize()

        elif len(atom.bonds) == 3:

            loc = self.get_position_with_three_bonds(atom)
            residue.create_atom(addname, loc)
            self.routines.cells.add_cell(residue.get_atom(addname))

    def complete(self):
        """Complete the water optimization class"""
        self.finalize()
        residue = self.residue

        atomlist = []
        for atom in residue.atoms:
            atomlist.append(atom)
        for atom in atomlist:
            if atom.name.startswith("LP"):
                residue.remove_atom(atom.name)


class Carboxylic(optimize.Optimize):
    """The class for carboxylic residues"""

    def __init__(self, residue, optinstance, routines):
        """Initialize a case where the lone hydrogen atom can have four
        different orientations.  Works similar to initializeFlip by preadding
        the necessary atoms.

        This also takes into account that the carboxyl group
        has different bond lengths for the two C-O bonds -
        this is probably due to one bond being assigned
        as a C=O.  As a result hydrogens are only added to
        the C-O (longer) bond.

        Args:
            residue:  The residue to flip (residue)
            dihedral: The angle to flip about
            hname:    The name of one of the hydrogens to add

        Returns:
            optlist:  A list of optimizeable donors and
                        acceptors in the residue (list)
        """
        # Initialize some variables
        self.optinstance = optinstance
        self.residue = residue
        self.routines = routines
        self.atomlist = []
        self.hbonds = []
        self.hlist = []

        hname2 = ""
        hname1 = ""
        for name in optinstance.map.keys():
            if name.endswith("2"):
                hname2 = name
            else:
                hname1 = name

        bondatom1 = residue.get_atom(optinstance.map[hname1].bond)
        bondatom2 = residue.get_atom(optinstance.map[hname2].bond)
        longflag = 0

        # If one bond in the group is significantly (0.05 A)
        # longer than the other, use that group only
        for pivotatom in bondatom1.bonds:
            if not pivotatom.is_hydrogen:
                the_pivatom = pivotatom
                break

        dist1 = util.distance(the_pivatom.coords, bondatom1.coords)
        dist2 = util.distance(the_pivatom.coords, bondatom2.coords)

        order = [hname1, hname2]

        if dist2 > dist1 and abs(dist1 - dist2) > 0.05:
            longflag = 1
            order = [hname2, hname1]

        elif dist1 > dist2 and abs(dist1 - dist2) > 0.05:
            longflag = 1
            order = [hname1, hname2]

        for hname in order:
            bondatom = residue.get_atom(optinstance.map[hname].bond)

            # First mirror the hydrogen about the same donor
            for dihedral in residue.reference.dihedrals:
                if dihedral.endswith(hname):
                    the_dihedral = dihedral
                    break

            anglenum = residue.reference.dihedrals.index(the_dihedral)
            if anglenum == -1:
                raise IndexError("Unable to find Carboxylic dihedral angle!")

            if residue.dihedrals[anglenum] is None:
                self.atomlist.append(bondatom)
                continue

            newangle = 180.0 + residue.dihedrals[anglenum]
            self.routines.set_dihedral_angle(residue, anglenum, newangle)

            hatom = residue.get_atom(hname)
            newcoords = hatom.coords

            # Flip back to return original atom
            newangle = 180.0 + residue.dihedrals[anglenum]
            self.routines.set_dihedral_angle(residue, anglenum, newangle)

            # Rename the original atom and rebuild the new atom
            residue.rename_atom(hname, "%s1" % hname)
            newname = "%s2" % hname
            residue.create_atom(newname, newcoords)
            newatom = residue.get_atom(newname)
            self.routines.cells.add_cell(newatom)
            newatom.refdistance = hatom.refdistance

            # Set the bonds for the new atom
            if bondatom not in newatom.bonds:
                newatom.bonds.append(bondatom)
            if newatom not in bondatom.bonds:
                bondatom.bonds.append(newatom)

            # Break if this is the only atom to add
            self.atomlist.append(bondatom)
            self.hlist.append(residue.get_atom("%s1" % hname))
            self.hlist.append(residue.get_atom("%s2" % hname))

            if longflag:
                break

        residue.set_donors_acceptors()

    def try_both(self, donor, acc, accobj):
        """Called when both the donor and acceptor are optimizeable.
           If one is fixed, we only need to try one side.  Otherwise
           first try to satisfy the donor - if that's succesful,
           try to satisfy the acceptor.  An undo may be necessary
           if the donor is satisfied and the acceptor isn't.
        """
        # If one residue if fixed, use other functions
        if donor.residue.fixed:
            if accobj.try_acceptor(acc, donor):
                return 1
            return 0
        if acc.residue.fixed:
            if self.try_donor(donor, acc):
                return 1
            return 0

        _LOGGER.debug("Working on %s %s (donor) to %s %s (acceptor)",
                      donor.residue, donor.name, acc.residue, acc.name)

        if self.is_hbond(donor, acc):
            if accobj.try_acceptor(acc, donor):
                self.fix(donor, acc)
                _LOGGER.debug("NET BOND SUCCESSFUL!")
                return 1
            return 0
        return 0

    def is_carboxylic_hbond(self, donor, acc):
        """Determine whether this donor acceptor pair is a hydrogen bond"""
        for donorhatom in donor.bonds:
            if not donorhatom.is_hydrogen:
                continue

            # Check the H(D)-A distance
            dist = util.distance(donorhatom.coords, acc.coords)
            if dist > DIST_CUTOFF:
                continue

            # Check the A-D-H(D) angle
            angle = self.get_hbond_angle(acc, donor, donorhatom)
            if angle <= ANGLE_CUTOFF:
                _LOGGER.debug("Found HBOND! %.4f %.4f", dist, angle)
                return 1

        # If we get here, no bond is formed
        return 0

    def try_acceptor(self, acc, donor):
        """The main driver for adding an LP to an optimizeable residue."""
        residue = acc.residue

        # Do some error checking
        if not donor.hdonor:
            return 0

        _LOGGER.debug("Working on %s %s (acceptor) to %s %s (donor)",
                      acc.residue, acc.name, donor.residue, donor.name)

        # We want to ignore the Hs on the acceptor
        if self.is_carboxylic_hbond(donor, acc):

            # Eliminate the closer hydrogen
            hyds = []
            dist = None
            donorhatom = None
            for hatom in self.hlist:
                if hatom.is_hydrogen:
                    hyds.append(hatom)

            if len(hyds) < 2:
                return 1

            dist = util.distance(hyds[0].coords, donor.coords)
            dist2 = util.distance(hyds[1].coords, donor.coords)
            # Eliminate hyds[0]
            if dist < dist2:
                self.hlist.remove(hyds[0])
                self.routines.cells.remove_cell(hyds[0])
                residue.remove_atom(hyds[0].name)
                donorhatom = residue.get_atom(hyds[1].name)
            elif hyds[1] in self.hlist:
                self.hlist.remove(hyds[1])
                self.routines.cells.remove_cell(hyds[1])
                residue.remove_atom(hyds[1].name)
                if residue.has_atom(hyds[0].name):
                    donorhatom = residue.get_atom(hyds[0].name)
                elif len(self.hlist) != 0 and residue.has_atom(self.hlist[0].name):
                    donorhatom = residue.get_atom(self.hlist[0].name)

            # If only one H is left, we're done
            if len(self.hlist) == 1:
                if donorhatom != None:
                    self.rename(donorhatom)
                residue.fixed = 1
            return 1

        else:
            return 0

    def try_donor(self, donor, acc):
        """The main driver for adding a hydrogen to an optimizeable residue."""
        # Do some error checking
        if not acc.hacceptor:
            return 0

        if self.is_hbond(donor, acc):
            self.fix(donor, acc)
            return 1
        return 0

    def fix(self, donor, acc):
        """Fix the carboxylic residue."""
        _LOGGER.debug("Fixing residue %s due to %s", donor.residue, donor.name)
        residue = donor.residue

        # Grab the H(D) that caused the bond
        for donorhatom in donor.bonds:
            if donorhatom.is_hydrogen:
                if self.get_hbond_angle(acc, donor, donorhatom) <= ANGLE_CUTOFF:
                    the_donorhatom = donorhatom
                    break

        # Remove all the other available bonded hydrogens
        hydrogens = self.hlist[:]
        for atom in hydrogens:
            if atom != the_donorhatom:
                self.routines.cells.remove_cell(atom)
                self.hlist.remove(atom)
                residue.remove_atom(atom.name)

        # Rename the atoms
        self.rename(the_donorhatom)
        residue.fixed = 1

    def finalize(self):
        """Finalize a protontated residue.  Try to minimize conflict with
        nearby atoms.
        """

        # Initialize some variables
        hydrogens = []
        bestatom = None
        residue = self.residue

        if residue.fixed:
            return

        # For each atom, get the closest atom
        bestenergy = 999.99
        for hydatom in self.hlist:
            energy = 0
            bondedatom = hydatom.bonds[0]
            closeatoms = self.routines.cells.get_near_cells(bondedatom)
            for catom in closeatoms:
                energy += self.get_pair_energy(bondedatom, catom) + \
                    self.get_pair_energy(catom, bondedatom)
            if energy < bestenergy:
                bestenergy = energy
                bestatom = hydatom

        # Keep the bestatom
        for hydatom in self.hlist:
            hydrogens.append(hydatom)
        for hydatom in hydrogens:
            if bestatom != hydatom:
                self.hlist.remove(hydatom)
                residue.remove_atom(hydatom.name)

        # Rename the atoms
        if bestatom != None and len(bestatom.name) == 4:
            self.rename(bestatom)
        else:
            pass
        residue.fixed = 1

    def rename(self, hydatom):
        """Rename the optimized atoms appropriately.

        This is done since the forcefields tend to require that the hydrogen
        is linked to a specific oxygen, and this atom may have different
        parameter values.

        Args:
            hydatom:  The hydrogen atom that was added. (atom)
        """
        residue = self.residue
        optinstance = self.optinstance

        # No need to rename if hydatom is not in residue.map
        if hydatom.name not in residue.map.keys():
            return

        # Take off the extension
        if len(hydatom.name) == 4:
            hname = hydatom.name[:-1]
            residue.rename_atom(hydatom.name, hname)

        # PATCHES.xml expects *2 - if it's *1 that left, flip names
        if len(self.atomlist) == 2:
            if hydatom.name.endswith("1"):
                residue.rename_atom(hydatom.name, "%s2" %  hydatom.name[:-1])
                bondname0 = self.atomlist[0].name
                bondname1 = self.atomlist[1].name
                tempname = "FLIP"
                residue.rename_atom(self.atomlist[0].name, tempname)
                residue.rename_atom(self.atomlist[1].name, bondname0)
                residue.rename_atom(tempname, bondname1)
            elif hydatom.name.endswith("OXT"):
                residue.rename_atom(hydatom.name, "HO")
                bondname0 = self.atomlist[0].name
                bondname1 = self.atomlist[1].name
                tempname = "FLIP"
                residue.rename_atom(self.atomlist[0].name, tempname)
                residue.rename_atom(self.atomlist[1].name, bondname0)
                residue.rename_atom(tempname, bondname1)

        elif len(self.atomlist) == 1:
            # Appending the other bondatom to self.atomlist
            hnames = [hname[:-1] + "1", hname[:-1] + "2"]
            for hname in hnames:
                bondatom = residue.get_atom(optinstance.map[hname].bond)
                if bondatom.name != self.atomlist[0].name:
                    self.atomlist.append(bondatom)
                else:
                    pass

            if hydatom.name.endswith("1"):
                if hydatom.name[:-1] + "2" in residue.map.keys():
                    residue.remove_atom("%s2" % hydatom.name[:-1])
                residue.rename_atom(hydatom.name, "%s2" % hydatom.name[:-1])
                bondname0 = self.atomlist[0].name
                bondname1 = self.atomlist[1].name
                tempname = "FLIP"
                residue.rename_atom(self.atomlist[0].name, tempname)
                residue.rename_atom(self.atomlist[1].name, bondname0)
                residue.rename_atom(tempname, bondname1)
            elif hydatom.name.endswith("OXT"):
                residue.rename_atom(hydatom.name, "HO")
                bondname0 = self.atomlist[0].name
                bondname1 = self.atomlist[1].name
                tempname = "FLIP"
                residue.rename_atom(self.atomlist[0].name, tempname)
                residue.rename_atom(self.atomlist[1].name, bondname0)
                residue.rename_atom(tempname, bondname1)

    def complete(self):
        """If not already fixed, finalize"""
        if (len(self.hlist) == 2) and (self.residue.fixed == 1):
            self.residue.fixed = 0
        if not self.residue.fixed:
            self.finalize()


class PotentialBond(object):
    """A small class containing the hbond structure"""

    def __init__(self, atom1, atom2, dist):
        """Initialize the class

        Args:
            atom1:  The first atom in the potential bond (Atom)
            atom2:  The second atom in the potential bond (Atom)
            dist:  The distance between the two atoms (float)
        """
        self.atom1 = atom1
        self.atom2 = atom2
        self.dist = dist

    def __str__(self):
        txt = "%s %s" % (self.atom1.name, self.atom1.residue)
        txt += " to "
        txt += "%s %s" % (self.atom2.name, self.atom2.residue)
        txt += " (%.2f A)" % self.dist
        return txt


class HydrogenDefinition(object):
    """HydrogenDefinition class

    The HydrogenDefinition class provides information on possible
    ambiguities in amino acid hydrogens.  It is essentially the hydrogen
    definition file in object form.
    """

    def __init__(self, name, opttype, optangle, map_):
        """Initialize the object with information from the definition file

        Args:
            name:          The name of the grouping (string)
            opttype:       The optimization type of the grouping (string)
            optangle:      The optimization angle of the grouping (string)
            map:           The map of Hydrogens.

            See HYDROGENS.XML for more information
        """
        self.name = name
        self.opttype = opttype
        self.optangle = optangle
        self.map = map_
        self.conformations = []

    def __str__(self):
        output = "Name:                  %s\n" % self.name
        output += "Opttype:               %s\n" % self.opttype
        output += "Optangle:              %s\n" % self.optangle
        output += "map:                   %s\n" % self.map
        output += "Conformations:\n"
        for conf in self.conformations:
            output += "\n%s" % conf
        output += "*****************************************\n"
        return output

    def add_conf(self, conf):
        """Add a HydrogenConformation to the list of conformations

        Args:
            conf:  The conformation to be added (HydrogenConformation)
        """
        self.conformations.append(conf)


class HydrogenConformation(object):
    """HydrogenConformation class

    The HydrogenConformation class contains data about possible
    hydrogen conformations as specified in the hydrogen data file.
    """

    def __init__(self, hname, boundatom, bondlength):
        """
        Args:
            hname      : The hydrogen name (string)
            boundatom  : The atom the hydrogen is bound to (string)
            bondlength : The bond length (float)
        """
        self.hname = hname
        self.boundatom = boundatom
        self.bondlength = bondlength
        self.atoms = []

    def __str__(self):
        output = "Hydrogen Name: %s\n" % self.hname
        output += "Bound Atom:    %s\n" % self.boundatom
        output += "Bond Length:   %.2f\n" % self.bondlength
        for atom in self.atoms:
            output += "\t%s\n" % atom
        return output

    def add_atom(self, atom):
        """Add an atom to the list of atoms

        Args:
            atom: The atom to be added (DefinitionAtom)
        """
        self.atoms.append(atom)


class HydrogenAmbiguity(object):
    """
        A class containing information about the ambiguity
    """
    def __init__(self, residue, hdef, routines):
        """If the residue has a rotateable hydrogen, remove it.

        If it can be flipped, pre-flip the residue by creating all additional
        atoms.

        Args:
            residue:  The residue in question (residue)
            hdef:     The hydrogen definition matching the residue
            routines: Pointer to the general routines object
        """
        self.residue = residue
        self.hdef = hdef
        self.routines = routines

    def __str__(self):
        text = "%s %i %s (%s)" % (self.residue.name, self.residue.res_seq, \
                                  self.residue.chain_id, self.hdef.opttype)
        return text


class HydrogenHandler(sax.ContentHandler):
    """Extends the SAX XML Parser to parse the Hydrogens.xml class"""
    def __init__(self):
        """
            Initalize the class.
        """

        self.curelement = ""
        self.curatom = None
        self.curobj = None
        self.curholder = None
        self.map = {}

    def startElement(self, name, _):
        """Create optimization holder objects or atoms"""
        if name == "class":
            obj = optimize.OptimizationHolder()
            self.curholder = obj
            self.curobj = obj
        elif name == "atom":
            obj = defns.DefinitionAtom()
            self.curatom = obj
            self.curobj = obj
        else:
            self.curelement = name
        return

    def endElement(self, name):
        """Complete whatever object is currently passed in by the name
        parameter
        """
        if name == "class": # Complete Residue object
            obj = self.curholder
            if not isinstance(obj, optimize.OptimizationHolder):
                raise ValueError("Internal error parsing XML!")

            self.map[obj.name] = obj
            self.curholder = None
            self.curobj = None

        elif name == "atom": # Complete atom object
            atom = self.curatom
            if not isinstance(atom, defns.DefinitionAtom):
                raise ValueError("Internal error parsing XML!")

            atomname = atom.name
            if atomname == "":
                raise ValueError("Atom name not set in XML!")
            else:
                self.curholder.map[atomname] = atom
                self.curatom = None
                self.curobj = self.curholder

        else: # Just free the current element namespace
            self.curelement = ""

        return self.map

    def characters(self, text):
        """Set a given attribute of the object to the text"""

        if text.isspace():
            return

        # If this is a float, make it so
        try:
            value = float(str(text))
        except ValueError:
            value = str(text)

        setattr(self.curobj, self.curelement, value)

