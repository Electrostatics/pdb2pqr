"""Hydrogen optimization for PDB2PQR

This is an module for hydrogen optimization routines.

Authors:  Todd Dolinsky, Jens Erik Nielsen, Yong Huang
"""
import logging
import math
from xml import sax
from . import topology
from .utilities import getDatFile, distance, subtract, normalize, dot, add
from .utilities import analyzeConnectivity, sortDictByValue, DuplicateFilter
from .quatfit import find_coordinates
from .definitions import DefinitionAtom
from .aa import Amino, WAT, HIS
from .routines import Cells, Routines

# TODO - This module is insane... so many lines!


_LOGGER = logging.getLogger(__name__)
_LOGGER.addFilter(DuplicateFilter())


HDEBUG = 0
HYDPATH = "dat/HYDROGENS.xml"
TOPOLOGYPATH = "dat/TOPOLOGY.xml"
ANGLE_CUTOFF = 20.0       # A - D - H(D) angle
DIST_CUTOFF = 3.3         # H(D) to A distance


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
            obj = OptimizationHolder()
            self.curholder = obj
            self.curobj = obj
        elif name == "atom":
            obj = DefinitionAtom()
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
            if not isinstance(obj, OptimizationHolder):
                raise ValueError("Internal error parsing XML!")

            self.map[obj.name] = obj
            self.curholder = None
            self.curobj = None

        elif name == "atom": # Complete atom object
            atom = self.curatom
            if not isinstance(atom, DefinitionAtom):
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
        coords1 = subtract(atom3.coords, atom_coords)
        coords2 = subtract(atom1.coords, atom_coords)
        norm1 = normalize(coords1)
        norm2 = normalize(coords2)
        dotted = dot(norm1, norm2)
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
            dist = distance(donorhatom.coords, acc.coords)
            if dist > DIST_CUTOFF:
                continue

            # Ensure no conflicts if H(A)s if present
            flag = 1
            for acchatom in acc.bonds:
                if not acchatom.is_hydrogen:
                    continue
                flag = 0

                # Check the H(D)-H(A) distance
                hdist = distance(donorhatom.coords, acchatom.coords)
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
            dist = distance(donorhatom.coords, acceptor.coords)
            if dist > max_dha_dist and dist < max_ele_dist:
                energy += max_ele_energy/(dist*dist)
                continue

            # Case 1: Both donor and acceptor hydrogens are present
            for acceptorhatom in acceptorhs:
                # Penalize if H(D) is too close to H(A)
                hdist = distance(donorhatom.coords, acceptorhatom.coords)
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
        vec = subtract(closeatom.coords, atom.coords)
        dist = distance(atom.coords, closeatom.coords)

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
        newcoords = find_coordinates(2, coords, refcoords, refatomcoords)
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
        newcoords = find_coordinates(2, coords, refcoords, refatomcoords)
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
        newcoords = find_coordinates(2, coords, refcoords, refatomcoords)
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
        if distance(rot2.coords, newcoords1) > 0.1:
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

class Flip(Optimize):
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
        moveablenames = self.routines.get_moveable_names(residue, pivot)
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
        if isinstance(residue, HIS):
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


class Alcoholic(Optimize):
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

# TODO - stopped de-linting here
class Water(Optimize):
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
                bestdist = distance(acc.coords, donor.coords)
                for donorh in donor.bonds:
                    dist = distance(acc.coords, donorh.coords)
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
                vec = subtract(atom.coords, closeatom.coords)
                dist = distance(atom.coords, closeatom.coords)

                for i in range(3):
                    newcoords.append(vec[i]/dist + atom.coords[i])

            else:
                newcoords = add(atom.coords, [1.0, 0.0, 0.0])

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

                dist = distance(nearatom.coords, newatom.coords)

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
                dist1 = distance(newatom.coords, nearatom.coords)

                # Place at other location
                self.routines.cells.remove_cell(atom)
                newatom.x = loc2[0]
                newatom.y = loc2[1]
                newatom.z = loc2[2]
                self.routines.cells.add_cell(atom)

                nearatom = self.routines.get_closest_atom(newatom)
                if nearatom != None:

                    # If this is worse, switch back
                    if distance(newatom.coords, nearatom.coords) < dist1:
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


class Carboxylic(Optimize):
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

        dist1 = distance(the_pivatom.coords, bondatom1.coords)
        dist2 = distance(the_pivatom.coords, bondatom2.coords)

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
            dist = distance(donorhatom.coords, acc.coords)
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

            dist = distance(hyds[0].coords, donor.coords)
            dist2 = distance(hyds[1].coords, donor.coords)
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


class Generic(Optimize):
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


class HydrogenRoutines(object):
    """The main routines for hydrogen optimization.
    This could potentially be extended from the routines object..."""

    def __init__(self, routines):
        self.routines = routines
        self.protein = routines.protein
        self.optlist = []
        self.atomlist = []
        self.resmap = {}
        self.hydrodefs = []

        handler = HydrogenHandler()
        sax.make_parser()

        # TODO - I don't think files should be loaded so deep in this module
        defpath = getDatFile(HYDPATH)
        if defpath == "":
            raise KeyError("Could not find %s!" % HYDPATH)

        hydrogen_file = open(defpath)
        sax.parseString(hydrogen_file.read(), handler)
        hydrogen_file.close()

        self.map = handler.map

    def switchstate(self, states, amb, state_id):
        """Switch a residue to a new state by first removing all hydrogens.

        Args:
            states: The list of states (list)
            amb   : The amibiguity to switch (tuple)
            state_id    : The state id to switch to (int)
        """
        if states == 'pKa':
            return self.pka_switchstate(amb, state_id)

        if state_id > len(states):
            raise IndexError("Invalid State ID!")

        # First Remove all Hs
        residue = getattr(amb, "residue")
        hdef = getattr(amb, "hdef")
        for conf in hdef.conformations:
            hname = conf.hname
            boundname = conf.boundatom
            if residue.get_atom(hname) != None:
                _LOGGER.debug('Removing %s %s %s', residue.name, residue.res_seq, hname)
                residue.remove_atom(hname)
            residue.get_atom(boundname).hacceptor = 1
            residue.get_atom(boundname).hdonor = 0

        # Update the IntraBonds
        name = residue.get("name")
        defresidue = self.routines.aadef.get_residue(name)
        residue.updateIntraBonds(defresidue)

        # Now build appropriate atoms
        state = states[state_id]
        for conf in state:
            _LOGGER.debug(conf)
            refcoords = []
            defcoords = []
            defatomcoords = []
            if conf == ():
                continue # Nothing to add
            hname = conf.hname
            for atom in conf.atoms:
                #print confatoms
                atomname = atom.get("name")
                resatom = residue.get_atom(atomname)
                if atomname == hname:
                    defatomcoords = atom.coords
                elif resatom != None:
                    refcoords.append(resatom.coords)
                    defcoords.append(atom.coords)
                else:
                    raise KeyError("Could not find necessary atom!")

            newcoords = find_coordinates(3, refcoords, defcoords, defatomcoords)
            boundname = conf.boundatom
            residue.create_atom(hname, newcoords, "ATOM")
            residue.addDebumpAtom(residue.get_atom(hname))
            residue.get_atom(boundname).addIntraBond(hname)
            residue.get_atom(boundname).hacceptor = 0
            residue.get_atom(boundname).hdonor = 1
            # Setting the SybylType for the newly built H
            residue.get_atom(hname).sybyl_type = 'H'
            # formal charge for PEOE_PB
            residue.get_atom(hname).formalcharge = 0.0
            # flag the added hydrogen
            residue.get_atom(hname).titratableH = True
            residue.get_atom(hname).addIntraBond(boundname)
        return None

    @classmethod
    def pka_switchstate(cls, amb, state_id_):
        """Switch a residue to a new state by first removing all hydrogens.
        This routine is used in pKa calculations only!

        Args:
            amb   : The amibiguity to switch (tuple)
            state_id    : The state id to switch to (list)
        """
        # TODO - This titration dictionary should not be buried in the code
        titrationdict = {'ASH1c': '1', 'ASH1t': '2', 'ASH2c': '3', 'ASH2t': '4', 'ASP': '0',
                         'GLH1c': '1', 'GLH1t': '2', 'GLH2c': '3', 'GLH2t': '4', 'GLU': '0',
                         'ARG0': '1+2+3+4', 'ARG': '1+2+3+4+5',
                         'LYS': '1', 'LYS0': '0',
                         'TYR': '1', 'TYR-': '0',
                         'HSD': '1', 'HSE': '2', 'HSP': '1+2',
                         'H3': '1', 'H2': '2', 'H3+H2': '1+2',
                         'CTR01c': '1', 'CTR01t': '2', 'CTR02c': '3', 'CTR02t': '4', 'CTR-': '0'}
        state_id = titrationdict[state_id_]
        state_id = state_id.split('+')
        new_state_id = []
        for i in state_id:
            new_state_id.append(int(i))
        residue = getattr(amb, "residue")
        hdef = getattr(amb, "hdef")
        for conf in hdef.conformations:
            hname = conf.hname
            boundname = conf.boundatom
            if residue.get_atom(hname) != None:
                residue.remove_atom(hname)
            residue.get_atom(boundname).hacceptor = 1
            residue.get_atom(boundname).hdonor = 0

        # Update the IntraBonds
        for state_id in new_state_id:
            if state_id == 0:
                continue
            conf = hdef.conformations[state_id-1]
            refcoords = []
            defcoords = []
            defatomcoords = []
            if conf == ():
                continue
            hname = conf.hname
            for atom in conf.atoms:
                if residue.is_n_term and residue.name == "PRO":
                    if atom.name == "H":
                        atom.name = "CD"
                        atom.x = 1.874
                        atom.y = 0.862
                        atom.z = 1.306

            if not Routines.rebuild_tetrahedral(residue, hname):
                for atom in conf.atoms:
                    atomname = atom.get("name")
                    resatom = residue.get_atom(atomname)
                    if atomname == hname:
                        defatomcoords = atom.coords
                    elif resatom != None:
                        refcoords.append(resatom.coords)
                        defcoords.append(atom.coords)
                    else:
                        raise KeyError("Could not find necessary atom!")

                newcoords = find_coordinates(3, refcoords, defcoords, defatomcoords)
                residue.create_atom(hname, newcoords)

            boundname = conf.boundatom
            residue.get_atom(boundname).hacceptor = 0
            residue.get_atom(boundname).hdonor = 1

            # Setting the SybylType for the newly built H
            residue.get_atom(hname).sybyl_type = 'H'

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
        return

    def cleanup(self):
        """If there are any extra carboxlyic *1 atoms, delete them.
        This may occur when no optimization is chosen
        """
        for residue in self.routines.protein.residues:
            if not isinstance(residue, Amino):
                continue
            if residue.name == "GLH" or "GLH" in residue.patches:
                if residue.has_atom("HE1") and residue.has_atom("HE2"):
                    residue.remove_atom("HE1")
            elif residue.name == "ASH" or "ASH" in residue.patches:
                if residue.has_atom("HD1") and residue.has_atom("HD2"):
                    residue.remove_atom("HD1")

    def is_optimizeable(self, residue):
        """Check to see if the given residue is optimizeable
        There are three ways to identify a residue:

        1.  By name (i.e. HIS)
        2.  By reference name - a PDB file HSP has
            a HIS reference name
        3.  By patch - applied by PropKa, terminal selection

        Args:
            residue:  The residue in question (Residue)
        Returns:
            optinstance: None if not optimizeable, otherwise
                            the OptimizationHolder instance that
                            corresponds to the residue.
        """
        optinstance = None
        if not isinstance(residue, (Amino, WAT)):
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
        if optinstance != None:
            if optinstance.opttype == "Alcoholic":
                atomname = list(optinstance.map.keys())[0]
                if not residue.reference.has_atom(atomname):
                    optinstance = None

        return optinstance

    def set_optimizeable_hydrogens(self):
        """Set any hydrogen listed in HYDROGENS.xml that is optimizeable.
        Used BEFORE hydrogen optimization to label atoms so that they won't be
        debumped - i.e. if SER HG is too close to another atom, don't debump
        but wait for optimization.  This function should not be used if full
        optimization is not taking place.
        """
        for residue in self.protein.residues:
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
        _LOGGER.info("Initializing full optimization...")

        # Do some setup
        self.routines.cells = Cells(5)
        self.routines.cells.assign_cells(self.protein)
        self.routines.calculate_dihedral_angles()
        self.routines.set_donors_acceptors()
        self.routines.update_internal_bonds()
        self.routines.set_reference_distance()
        self.optlist = []
        self.atomlist = []

        # First initialize the various types
        for residue in self.protein.residues:
            optinstance = self.is_optimizeable(residue)
            if isinstance(residue, Amino):
                if False in residue.stateboolean.values():
                    residue.fixed = 1
                else:
                    residue.fixed = 0
            if optinstance is None:
                continue

            type_ = optinstance.opttype
            if residue.fixed == 1:
                pass
            else:
                klass = globals()[type_]
                myobj = klass(residue, optinstance, self.routines)
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
        self.routines.cells = Cells(5)
        self.routines.cells.assign_cells(self.protein)
        self.routines.calculate_dihedral_angles()
        self.routines.set_donors_acceptors()
        self.routines.update_internal_bonds()
        self.routines.set_reference_distance()
        self.optlist = []

        # First initialize the various types
        for residue in self.protein.residues:
            optinstance = self.is_optimizeable(residue)
            if optinstance is None:
                continue

            type_ = optinstance.opttype
            if type_ == "Water":
                klass = globals()[type_]
                myobj = klass(residue, optinstance, self.routines)
                self.atomlist += myobj.atomlist
                self.optlist.append(myobj)
                self.resmap[residue] = myobj

        _LOGGER.debug("Done.")

    def optimize_hydrogens(self):
        """The main driver for the optimization.
        Should be called only after the optlist has been initialized.
        """
        _LOGGER.debug("Optimization progress:")

        optlist = self.optlist
        connectivity = {}

        # Initialize the detection progress
        if len(optlist) == 0:
            return

        _LOGGER.debug("  Detecting potential hydrogen bonds")
        progress = 0.0
        increment = 1.0/len(optlist)

        for obj in optlist:
            connectivity[obj] = []
            for atom in obj.atomlist:
                closeatoms = self.routines.cells.get_near_cells(atom)
                for closeatom in closeatoms:

                    # Conditions for continuing
                    if atom.residue == closeatom.residue:
                        continue
                    if not (closeatom.hacceptor or closeatom.hdonor):
                        continue
                    if atom.hdonor and not atom.hacceptor:
                        if not closeatom.hacceptor:
                            continue
                    if atom.hacceptor:
                        if not atom.hdonor and not closeatom.hdonor:
                            continue

                    dist = distance(atom.coords, closeatom.coords)
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
                _LOGGER.debug("%s has no nearby partners - fixing.", obj.residue)
                obj.finalize()

        # Determine the distinct networks
        networks = []
        seen = []
        for obj1 in optlist:
            if obj1.residue.fixed:
                continue
            if obj1 in seen:
                continue
            network = analyzeConnectivity(connectivity, obj1)
            for obj2 in network:
                if obj2 not in seen:
                    seen.append(obj2)
            networks.append(network)

        # Initialize the output progress
        if len(networks) > 0:
            _LOGGER.debug("Optimizing hydrogen bonds")
            progress = 0.0
            increment = 1.0/len(networks)

        # Work on the networks
        for network in networks:
            txt = ""
            for obj in network:
                txt += "%s, " % obj
            _LOGGER.debug("Starting network %s", txt[:-2])

            ###  FIRST:  Only optimizeable to backbone atoms
            _LOGGER.debug("* Optimizeable to backbone *")
            hbondmap = {}
            for obj in network:
                for hbond in obj.hbonds:
                    if hbond.atom2 not in self.atomlist:
                        hbondmap[hbond] = hbond.dist
            hbondlist = sortDictByValue(hbondmap)
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

            ### SECOND:  Non-dual water Optimizeable to Optimizeable
            _LOGGER.debug("* Optimizeable to optimizeable *")
            hbondmap = {}
            seenlist = []
            for obj in network:
                for hbond in obj.hbonds:
                    if hbond.atom2 in self.atomlist:
                        if not isinstance(hbond.atom1.residue, WAT):
                            if not isinstance(hbond.atom2.residue, WAT):
                                # Only get one hbond pair
                                if (hbond.atom2, hbond.atom1) not in seenlist:
                                    hbondmap[hbond] = hbond.dist
                                    seenlist.append((hbond.atom1, hbond.atom2))

            hbondlist = sortDictByValue(hbondmap)
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

            ### THIRD:  All water-water residues
            _LOGGER.debug("* Water to Water *")
            hbondmap = {}
            seenlist = []
            for obj in network:
                for hbond in obj.hbonds:
                    residue = hbond.atom1.residue
                    if isinstance(residue, WAT):
                        if isinstance(hbond.atom2.residue, WAT):
                            if (hbond.atom2, hbond.atom1) not in seenlist:
                                hbondmap[hbond] = hbond.dist
                                seenlist.append((hbond.atom1, hbond.atom2))

            hbondlist = sortDictByValue(hbondmap)
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

            ### FOURTH: Complete all residues
            for obj in network:
                obj.complete()

            # STEP 5:  Update progress meter
            progress += 100.0 * increment
            while progress >= 5.0:
                progress -= 5.0

    def parse_hydrogen(self, res):
        """Parse a list of lines in order to make a hydrogen definition

        Args:
            lines:  The lines to parse (list)
        Returns:
            mydef:  The hydrogen definition object (HydrogenDefinition)

        This is the current definition:  Name Ttyp  A R # Stdconf   HT Chi OPTm
        """
        # TODO - I don't think files should be loaded so deep in this module
        toppath = getDatFile(TOPOLOGYPATH)
        if toppath == "":
            raise KeyError("Could not find %s!" % TOPOLOGYPATH)

        with open(toppath) as topfile:
            top = topology.Topology(topfile)

        name = self.map[res].name
        opttype = self.map[res].opttype
        optangle = self.map[res].optangle
        map_ = self.map[res].map

        mydef = HydrogenDefinition(name, opttype, optangle, map_)
        patchmap = []
        refmap = {}
        titrationstatemap = {}
        tautomermap = {}
        conformermap = {}
        atommap = {}

        # reference map from TOPOLOGY.xml
        for res_ in top.residues:
            refmap[res_.name] = res_.reference
            for atom in refmap[res_.name].atoms:
                atommap[res_.name, atom.name] = atom
            for titrationstate in res_.titrationStates:
                titrationstatemap[titrationstate.name] = titrationstate
                for tautomer in titrationstate.tautomers:
                    tautomermap[tautomer.name] = tautomer
                    for conformer in tautomer.conformers:
                        conformermap[conformer.name] = conformer

        if name == 'CYS':
            _ = refmap['CYS']
            atoms = ['HG']
            refatoms = ['SG', 'CB']

        elif name == 'HIS':
            _ = refmap['HIS']
            atoms = ['HD1', 'HE2']
            for atom in atoms:
                refatoms = ['ND1', 'CG', 'CE1']

        elif name == 'LYS':
            _ = self.routines.protein.referencemap[name]
            patchmap = self.routines.protein.patchmap['LYN']
            atoms = patchmap.remove
            refatoms = ['HZ1', 'HZ2', 'NZ']

        elif name == 'TYR':
            _ = self.routines.protein.referencemap[name]
            patchmap = self.routines.protein.patchmap['TYM']
            atoms = patchmap.remove
            refatoms = ['OH', 'CZ', 'CE2']

        elif name == 'WAT':
            _ = self.routines.protein.referencemap[name]
            patchmap = self.routines.protein.patchmap['HOH']
            atoms = ['H1', 'H2']
            refatoms = None

        elif name == 'NTR':
            ntrmap = {}    # map for N-TERM
            for tautomer in titrationstatemap["NTER"].tautomers:
                for conformer in tautomermap[tautomer.name].conformers:
                    for conformeradds in conformermap[conformer.name].conformerAdds:
                        for atom in conformeradds.atoms:
                            ntrmap[atom.name] = atom
            atoms = ['H3', 'H2']
            refatoms = ['CA', 'H', 'N']

        elif name == 'CTR':
            hmap = {} # map for h atoms
            nonhmap = {} # map for refatoms
            conformernames = []
            for tautomer in titrationstatemap["CTER"].tautomers:
                for conformer in tautomermap[tautomer.name].conformers:
                    for conformeradds in conformermap[conformer.name].conformerAdds:
                        for atom in conformeradds.atoms:
                            nonhmap[atom.name] = atom
            for tautomer in titrationstatemap["CTER0"].tautomers:
                for conformer in tautomermap[tautomer.name].conformers:
                    conformernames.append(conformer.name)
                    for conformeradds in conformermap[conformer.name].conformerAdds:
                        for atom in conformeradds.atoms:
                            hmap[conformer.name, atom.name] = atom

            atoms = ['HO']
            refatoms = ['O', 'C', 'OXT']

        elif name in ['SER', 'GLN', 'THR', 'ARG', 'ASN']:
            _ = refmap[name]
            if name == 'SER':
                atoms = ['HG']
                refatoms = ['OG', 'CB']
            elif name == 'GLN':
                atoms = ['HE21']
                refatoms = ['NE2']
            elif name == 'THR':
                atoms = ['HG1']
                refatoms = ['OG1', 'CB']
            elif name == 'ARG':
                atoms = ['HH11', 'HH12', 'HH21', 'HH22', 'HE']
                for atom in atoms:
                    refatoms = ['NH1', 'NH2', 'CZ']
            elif name == 'ASN':
                atoms = ['HD21']
                refatoms = ['ND2']

        elif name == 'ASH':
            hmap = {}    # map for h atoms
            nonhmap = {}    # map for refatoms
            conformernames = []
            _ = refmap['ASP']
            for tautomer in titrationstatemap["ASH"].tautomers:
                for conformer in tautomermap[tautomer.name].conformers:
                    for conformeradds in conformermap[conformer.name].conformerAdds:
                        for atom in conformeradds.atoms:
                            hmap[conformer.name, atom.name] = atom
                            conformernames.append(conformer.name)
            atoms = ['HD1', 'HD2']
            refatoms = ['OD1', 'CG', 'OD2']

        elif name == 'GLH':
            hmap = {} # map for h atoms
            nonhmap = {} # map for refatoms
            conformernames = []
            _ = refmap['GLU']
            for tautomer in titrationstatemap["GLH"].tautomers:
                for conformer in tautomermap[tautomer.name].conformers:
                    for conformeradds in conformermap[conformer.name].conformerAdds:
                        for atom in conformeradds.atoms:
                            hmap[conformer.name, atom.name] = atom
                            conformernames.append(conformer.name)
            atoms = ['HE1', 'HE2']
            refatoms = ['OE1', 'CD', 'OE2']

        else:
            patchmap = self.routines.protein.patchmap[name]
            atoms = list(patchmap.map.keys())
            atoms.sort()

        if name in ['NTR']:
            for atom in atoms:
                hname = atom
                x = ntrmap[hname].x
                y = ntrmap[hname].y
                z = ntrmap[hname].z
                bondatom = ntrmap[hname].bonds[0]
                bondlength = 1.0
                myconf = HydrogenConformation(hname, bondatom, bondlength)
                atom = DefinitionAtom(hname, x, y, z)
                myconf.add_atom(atom)

                # TODO - lots of arbitrary undefined numbers in this section
                for atom_ in refatoms:
                    if atom_ == 'N':
                        natom = DefinitionAtom(atom_, 1.201, 0.847, 0.0)
                        myconf.add_atom(natom)
                    elif atom_ == 'CA':
                        caatom = DefinitionAtom(atom_, 0.0, 0.0, 0.0)
                        myconf.add_atom(caatom)
                    elif atom_ == 'H':
                        caatom = DefinitionAtom(atom_, 1.201, 1.847, 0.000)
                        myconf.add_atom(caatom)
                    else: pass
                mydef.add_conf(myconf)

        elif name in ['CTR']:
            for conformer in conformernames:
                for atom in atoms:
                    hname = atom
                    x = hmap[conformer, hname].x
                    y = hmap[conformer, hname].y
                    z = hmap[conformer, hname].z
                    bondatom = hmap[conformer, hname].bonds[0]
                    bondlength = 1.0
                    myconf = HydrogenConformation(hname, bondatom, bondlength)
                    atom = DefinitionAtom(hname, x, y, z)
                    myconf.add_atom(atom)

                    # TODO - the following code is almost nonsensical
                    for atom_ in refatoms:
                        if atom_ == 'C':
                            catom = DefinitionAtom(atom_, -1.250, 0.881, 0.000)
                            myconf.add_atom(catom)
                        else:
                            atomname = atom_
                            x = nonhmap[atom_].x
                            y = nonhmap[atom_].y
                            z = nonhmap[atom_].z
                            atom2 = DefinitionAtom(atomname, x, y, z)
                            myconf.add_atom(atom2)
                    mydef.add_conf(myconf)

        elif name in ['ASH', 'GLH']:
            for conformer in conformernames:
                for atom in atoms:
                    hname = atom
                    if ('1' in conformer and '1' in atom) or ('2' in conformer and '2' in atom):
                        x = hmap[conformer, hname].x
                        y = hmap[conformer, hname].y
                        z = hmap[conformer, hname].z
                        bondatom = hmap[conformer, hname].bonds[0]
                        bondlength = 1.0
                        myconf = HydrogenConformation(hname, bondatom, bondlength)
                        atom = DefinitionAtom(hname, x, y, z)
                        myconf.add_atom(atom)

                        for atom_ in refatoms:
                            atomname = atom_
                            if name == 'ASH':
                                refresname = 'ASP'
                            elif name == 'GLH':
                                refresname = 'GLU'
                            x = atommap[refresname, atom_].x
                            y = atommap[refresname, atom_].y
                            z = atommap[refresname, atom_].z
                            atom2 = DefinitionAtom(atomname, x, y, z)
                            myconf.add_atom(atom2)
                        mydef.add_conf(myconf)

        elif name in ['WAT']:
            pass

        else:
            for atom in atoms:
                hname = atom
                x = atommap[name, hname].x
                y = atommap[name, hname].y
                z = atommap[name, hname].z
                bondatom = atommap[name, hname].bonds[0]
                bondlength = 1.0
                myconf = HydrogenConformation(hname, bondatom, bondlength)
                atom = DefinitionAtom(hname, x, y, z)
                myconf.add_atom(atom)

                if refatoms != None:
                    if name == 'HIS' and atom.name == 'HE2':
                        refatoms = ['NE2', 'CE1', 'CD2']
                    if name == 'ARG' and atom.name == 'HE':
                        refatoms = ['NE', 'CZ', 'NH1']
                    for atom in refatoms:
                        atomname = atom
                        x = atommap[name, atomname].x
                        y = atommap[name, atomname].y
                        z = atommap[name, atomname].z
                        atom = DefinitionAtom(atomname, x, y, z)
                        myconf.add_atom(atom)
                    mydef.add_conf(myconf)
        return mydef

    def read_hydrogen_def(self):
        """Read the Hydrogen Definition file

        Returns
            hydrodef:  The hydrogen definition ()
        """
        self.hydrodefs = []
        for mapping in self.map:
            res = mapping
            mydef = self.parse_hydrogen(res)
            self.hydrodefs.append(mydef)
            res = ''


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
