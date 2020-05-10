"""Hydrogen bond optimization routines."""
import logging
import math
from .. import quatfit as quat
from .. import utilities as util
from ..config import ANGLE_CUTOFF, DIST_CUTOFF


_LOGGER = logging.getLogger(__name__)


class Optimize(object):
    """The holder class for the hydrogen optimization routines.

    Individual optimization types inherit off of this class.  Any functions
    used by multiple types appear here.
    """
    def __init__(self):
        self.residue = None
        self.optinstance = None
        self.routines = None

    def __str__(self):
        txt = "%s (%s)" % (self.residue, self.optinstance.opttype)
        return txt

    @staticmethod
    def get_hbond_angle(atom1, atom2, atom3):
        """Get the angle between three atoms

        Args:
            atom1:  The first atom (atom)
            atom2:  The second (vertex) atom (atom)
            atom3:  The third atom (atom)

        Returns:
            angle:  The angle between the atoms (float)
        """
        atom_coords = atom2.coords
        coords1 = util.subtract(atom3.coords, atom_coords)
        coords2 = util.subtract(atom1.coords, atom_coords)
        norm1 = util.normalize(coords1)
        norm2 = util.normalize(coords2)
        dotted = util.dot(norm1, norm2)
        if dotted > 1.0: # If normalized, this is due to rounding error
            dotted = 1.0
        elif dotted < -1.0: # If normalized, this is due to rounding error
            dotted = -1.0
        rad = abs(math.acos(dotted))
        angle = rad*180.0/math.pi
        if angle > 180.0:
            angle = 360.0 - angle
        return angle

    def is_hbond(self, donor, acc):
        """Determine whether this donor acceptor pair is a hydrogen bond"""
        for donorhatom in donor.bonds:
            if not donorhatom.is_hydrogen:
                continue

            # Check the H(D)-A distance
            dist = util.distance(donorhatom.coords, acc.coords)
            if dist > DIST_CUTOFF:
                continue

            # Ensure no conflicts if H(A)s if present
            flag = 1
            for acchatom in acc.bonds:
                if not acchatom.is_hydrogen:
                    continue
                flag = 0

                # Check the H(D)-H(A) distance
                hdist = util.distance(donorhatom.coords, acchatom.coords)
                if hdist < 1.5:
                    continue

                # Check the H(D)-H(A)-A angle
                angle = self.get_hbond_angle(donorhatom, acchatom, acc)
                if angle < 110.0:
                    flag = 1

            if flag == 0:
                continue

            # Check the A-D-H(D) angle
            angle = self.get_hbond_angle(acc, donor, donorhatom)
            if angle <= ANGLE_CUTOFF:
                _LOGGER.debug("Found HBOND! %.4f %.4f", dist, angle)
                return 1

        # If we get here, no bond is formed
        return 0

    @staticmethod
    def get_pair_energy(donor, acceptor):
        """Get the energy between two atoms

        Args:
            donor:  The first atom in the pair (Atom)
            acceptor:  The second atom in the pair (Atom)
        Returns:
            energy:  The energy of the pair (float)
        """

        # TODO - lots of code in here that could be accelerated with numpy
        # Initialize some variables
        bump_energy = 10.0
        bump_distance = 1.5
        max_hbond_energy = -10.0
        max_ele_energy = -1.0
        adh_angle_cutoff = ANGLE_CUTOFF
        dhaha_angle_cutoff = 110.0
        max_dha_dist = DIST_CUTOFF
        max_ele_dist = 5.0
        energy = 0.0

        if not (donor.hdonor and acceptor.hacceptor):
            return energy

        # See if hydrogens are presently bonded to the acceptor and donor
        donorhs = (bond for bond in donor.bonds if bond.is_hydrogen)
        acceptorhs = [bond for bond in acceptor.bonds if bond.is_hydrogen]
        for donorhatom in donorhs:
            dist = util.distance(donorhatom.coords, acceptor.coords)
            if dist > max_dha_dist and dist < max_ele_dist:
                energy += max_ele_energy/(dist*dist)
                continue

            # Case 1: Both donor and acceptor hydrogens are present
            for acceptorhatom in acceptorhs:
                # Penalize if H(D) is too close to H(A)
                hdist = util.distance(donorhatom.coords, acceptorhatom.coords)
                if hdist < bump_distance:
                    energy += bump_energy
                    continue

                # Assign energies based on angles
                angle1 = Optimize.get_hbond_angle(acceptor, donor, donorhatom)
                if angle1 <= adh_angle_cutoff:
                    angle2 = Optimize.get_hbond_angle(donorhatom, acceptorhatom, acceptor)
                    if angle2 < dhaha_angle_cutoff:
                        angle2 = 1.0
                    else:
                        angle2 = (dhaha_angle_cutoff - angle2)/dhaha_angle_cutoff

                    angleterm = (adh_angle_cutoff - angle1)/adh_angle_cutoff
                    energy += max_hbond_energy/pow(dist, 3)*angleterm*angle2

            # Case 2: Only donor hydrogens are present
            if len(acceptorhs) == 0:
                # Assign energies based on A-D-H(D) angle alone
                angle1 = Optimize.get_hbond_angle(acceptor, donor, donorhatom)
                if angle1 <= adh_angle_cutoff:
                    angleterm = (adh_angle_cutoff - angle1)/adh_angle_cutoff
                    energy += max_hbond_energy/pow(dist, 2)*angleterm

        return energy

    def make_atom_with_no_bonds(self, atom, closeatom, addname):
        """Called for water oxygen atoms with no current bonds.

        Uses the closeatom to place the new atom directly colinear with the
        atom and the closeatom.

        Args:
            atom:      The oxygen atom of the water
            closeatom: The nearby atom (donor/acceptor)
            addname:   The name of the atom to add
        """
        newcoords = []
        residue = atom.residue

        # Place along line, 1 A away
        vec = util.subtract(closeatom.coords, atom.coords)
        dist = util.distance(atom.coords, closeatom.coords)

        for i in range(3):
            newcoords.append(vec[i]/dist + atom.coords[i])

        residue.create_atom(addname, newcoords)
        newatom = residue.get_atom(addname)
        self.routines.cells.add_cell(newatom)

        # Set the bonds (since not in reference structure)
        if newatom not in atom.bonds:
            atom.bonds.append(newatom)
        if atom not in newatom.bonds:
            newatom.bonds.append(atom)

    @classmethod
    def make_water_with_one_bond(cls, atom, addname):
        """Add an atom to a water residue that already has one bond.

        Uses the water reference structure to align the new atom.
        """
        residue = atom.residue
        nextatom = atom.bonds[0]
        coords = [atom.coords, nextatom.coords]
        refcoords = [residue.reference.map[atom.name].coords, \
                     residue.reference.map["H1"].coords]
        refatomcoords = residue.reference.map["H2"].coords

        # Make the atom
        newcoords = quat.find_coordinates(2, coords, refcoords, refatomcoords)
        residue.create_atom(addname, newcoords)

        # Set the bonds (since not in reference structure)
        newatom = residue.get_atom(addname)
        if newatom not in atom.bonds:
            atom.bonds.append(newatom)
        if atom not in newatom.bonds:
            newatom.bonds.append(atom)

    @classmethod
    def make_atom_with_one_bond_h(cls, atom, addname):
        """Add a hydrogen to an alcoholic donor with one existing bond."""

        residue = atom.residue
        nextatom = atom.bonds[0]
        coords = [atom.coords, nextatom.coords]
        refcoords = [residue.reference.map[atom.name].coords, \
                     residue.reference.map[nextatom.name].coords]
        refatomcoords = residue.reference.map[addname].coords

        # Make the atom
        newcoords = quat.find_coordinates(2, coords, refcoords, refatomcoords)
        residue.create_atom(addname, newcoords)

    @classmethod
    def make_atom_with_one_bond_lp(cls, atom, addname):
        """Add a lone pair to an alcoholic donor with one existing bond."""

        # Initialize some variables
        residue = atom.residue
        for refname in atom.reference.bonds:
            if refname.startswith("H"):
                the_refname = refname
                break

        nextatom = atom.bonds[0]
        coords = [atom.coords, nextatom.coords]
        refcoords = [residue.reference.map[atom.name].coords,
                     residue.reference.map[nextatom.name].coords]
        refatomcoords = residue.reference.map[the_refname].coords

        # Make the atom
        newcoords = quat.find_coordinates(2, coords, refcoords, refatomcoords)
        residue.create_atom(addname, newcoords)

        # Set the bonds (since not in reference structure)
        newatom = residue.get_atom(addname)
        if newatom not in atom.bonds:
            atom.bonds.append(newatom)
        if atom not in newatom.bonds:
            newatom.bonds.append(atom)

    def try_single_alcoholic_h(self, donor, acc, newatom):
        """After a new bond has been added using makeAtomWithOneBond*, try to
        find the best orientation by rotating to form a hydrogen bond.  If a
        bond cannot be formed, remove the newatom (thereby returning to a single
        bond). """

        # Initialize some variables
        besten = 999.99
        bestcoords = []
        residue = donor.residue
        pivot = donor.bonds[0]

        for _ in range(72):
            residue.rotate_tetrahedral(pivot, donor, 5.0)
            if self.is_hbond(donor, acc):
                energy = self.get_pair_energy(donor, acc)
                if energy < besten:
                    bestcoords = newatom.coords
                    besten = energy

        # If a hydrogen bond was made, set at best coordinates
        if bestcoords != []:
            newatom.x = bestcoords[0]
            newatom.y = bestcoords[1]
            newatom.z = bestcoords[2]
            self.routines.cells.add_cell(newatom)
            return 1
        residue.remove_atom(newatom.name)
        return 0

    def try_single_alcoholic_lp(self, acc, donor, newatom):
        """After a new bond has been added using makeAtomWithOneBond*, ensure
        that a hydrogen bond has been made.  If so, try to minimze the H(D)-A-LP
        angle.  If that cannot be minimized, ignore the bond and remove the
        atom.
        """

        # Initialize some variables
        residue = acc.residue
        pivot = acc.bonds[0]
        bestangle = 180.00
        bestcoords = []

        # If a hydrogen bond was made, set at best distance
        if not self.is_hbond(donor, acc):
            residue.remove_atom(newatom.name)
            return 0

        # Grab the H(D) that caused the bond

        for donorhatom in donor.bonds:
            if donorhatom.is_hydrogen:
                if self.get_hbond_angle(acc, donor, donorhatom) < ANGLE_CUTOFF:
                    the_donorhatom = donorhatom
                    break

        # TODO - where did 72 come from (= 360/5)
        for _ in range(72):
            residue.rotate_tetrahedral(pivot, acc, 5.0)
            angle = abs(self.get_hbond_angle(the_donorhatom, acc, newatom))
            if angle < bestangle:
                bestangle = angle
                bestcoords = newatom.coords

        # Remove if geometry does not work
        if bestangle > (ANGLE_CUTOFF * 2.0):
            _LOGGER.debug("Removing due to geometry %.2f > %.2f", bestangle, ANGLE_CUTOFF*2.0)
            residue.remove_atom(newatom.name)
            return 0

        # Otherwise set to best coordinates
        newatom.x = bestcoords[0]
        newatom.y = bestcoords[1]
        newatom.z = bestcoords[2]
        self.routines.cells.add_cell(newatom)

        return 1

    @classmethod
    def get_positions_with_two_bonds(cls, atom):
        """Given a tetrahedral geometry with two existing bonds, return the two
        potential sets of coordinates that are possible for a new bond."""

        # Initialize some variables
        residue = atom.residue
        fixed = atom.bonds[0]
        rotate = atom.bonds[1]

        # Rotate by 120 degrees twice
        residue.rotate_tetrahedral(fixed, atom, 120)
        loc1 = rotate.coords
        residue.rotate_tetrahedral(fixed, atom, 120)
        loc2 = rotate.coords

        # Set rotate back to original by one more rotation
        residue.rotate_tetrahedral(fixed, atom, 120)

        return loc1, loc2

    def try_positions_with_two_bonds_h(self, donor, acc, newname, loc1, loc2):
        """Try adding a new hydrogen two the two potential locations.
        If both form hydrogen bonds, place at whatever returns the best bond as
        determined by get_pair_energy.
        """

        # Initialize some variables
        besten = 999.99
        bestcoords = []
        residue = donor.residue

        # Try the first position
        residue.create_atom(newname, loc1)
        if self.is_hbond(donor, acc):
            besten = self.get_pair_energy(donor, acc)
            bestcoords = loc1

        # Try the second
        newatom = residue.get_atom(newname)
        newatom.x = loc2[0]
        newatom.y = loc2[1]
        newatom.z = loc2[2]
        if self.is_hbond(donor, acc):
            energy = self.get_pair_energy(donor, acc)
            if energy < besten:
                bestcoords = loc2

        # Set at best coords
        if bestcoords != []:
            newatom.x = bestcoords[0]
            newatom.y = bestcoords[1]
            newatom.z = bestcoords[2]
            self.routines.cells.add_cell(newatom)
            return 1
        residue.remove_atom(newname)
        return 0

    def try_positions_with_two_bonds_lp(self, acc, donor, newname, loc1, loc2):
        """Try placing an LP on a tetrahedral geometry with two existing bonds.
        If this isn't a hydrogen bond it can return - otherwise ensure that the
        H(D)-A-LP angle is minimized.
        """

        # Initialize some variables
        bestangle = 180.00
        bestcoords = []
        residue = acc.residue

        # If the donor/acceptor pair is not an hbond return
        if not self.is_hbond(donor, acc):
            return 0

        # Grab the H(D) that caused the bond
        for donorhatom in donor.bonds:
            if donorhatom.is_hydrogen:
                if self.get_hbond_angle(acc, donor, donorhatom) < ANGLE_CUTOFF:
                    the_donorhatom = donorhatom
                    break

        # Try the first position
        residue.create_atom(newname, loc1)
        newatom = residue.get_atom(newname)
        angle = abs(self.get_hbond_angle(the_donorhatom, acc, newatom))
        if angle < bestangle:
            bestangle = angle
            bestcoords = loc1

        # Try the second
        newatom.x = loc2[0]
        newatom.y = loc2[1]
        newatom.z = loc2[2]
        angle = self.get_hbond_angle(the_donorhatom, acc, newatom)
        if angle < bestangle:
            bestcoords = loc2

        # Remove if geometry does not work
        if bestangle > (ANGLE_CUTOFF * 2.0):
            residue.remove_atom(newname)
            return 0

        # Otherwise set at best coords
        newatom.x = bestcoords[0]
        newatom.y = bestcoords[1]
        newatom.z = bestcoords[2]
        self.routines.cells.add_cell(newatom)

        # Set the bonds (since not in reference structure)
        if newatom not in acc.bonds:
            acc.bonds.append(newatom)
        if acc not in newatom.bonds:
            newatom.bonds.append(acc)

        return 1

    @classmethod
    def get_position_with_three_bonds(cls, atom):
        """If there's three bonds in a tetrahedral geometry, there's only one
        available position.
        Find that position.
        """

        # Initialize some variables
        residue = atom.residue
        pivot = atom.bonds[0]
        rot1 = atom.bonds[1]
        rot2 = atom.bonds[2]

        # Find the two new positions
        residue.rotate_tetrahedral(pivot, atom, 120)
        newcoords1 = rot1.coords
        residue.rotate_tetrahedral(pivot, atom, 120)
        newcoords2 = rot1.coords
        residue.rotate_tetrahedral(pivot, atom, 120)

        # Determine which is unoccupied
        if util.distance(rot2.coords, newcoords1) > 0.1:
            return newcoords1
        return newcoords2

    def try_positions_three_bonds_h(self, donor, acc, newname, loc):
        """Try making a hydrogen bond with the lone available position."""
        residue = donor.residue
        residue.create_atom(newname, loc)
        if self.is_hbond(donor, acc):
            newatom = residue.get_atom(newname)
            self.routines.cells.add_cell(newatom)
            return 1
        residue.remove_atom(newname)
        return 0

    def try_positions_three_bonds_lp(self, acc, donor, newname, loc):
        """Try making a hydrogen bond using the lone available hydrogen
        position."""
        residue = acc.residue
        if not self.is_hbond(donor, acc):
            return 0

        # Grab the H(D) that caused the bond
        for donorhatom in donor.bonds:
            if donorhatom.is_hydrogen:
                if self.get_hbond_angle(acc, donor, donorhatom) < ANGLE_CUTOFF:
                    the_donorhatom = donorhatom
                    break

        residue.create_atom(newname, loc)
        newatom = residue.get_atom(newname)

        # Remove if geometry does not work
        angle = abs(self.get_hbond_angle(the_donorhatom, acc, newatom))
        if angle > (ANGLE_CUTOFF * 2.0):
            residue.remove_atom(newname)
            return 0

        # Otherwise keep it
        newatom = residue.get_atom(newname)
        self.routines.cells.add_cell(newatom)

        # Set the bonds (since not in reference structure)
        if newatom not in acc.bonds:
            acc.bonds.append(newatom)
        if acc not in newatom.bonds:
            newatom.bonds.append(acc)

        return 1


class OptimizationHolder(object):
    """A holder class for the XML parser."""

    def __init__(self):
        self.name = ""
        self.map = {}
        self.opttype = ""
        self.optangle = ""

    def __str__(self):
        text = "%s\n" % self.name
        text += "Type: %s\n" % self.opttype
        if self.optangle != "":
            text += "Optimization Angle: %s\n" % self.optangle
        text += "Atoms: \n"
        for atomname in self.map:
            text += "\t%s\n" % str(self.map[atomname])
        return text
