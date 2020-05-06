"""This file contains classes associated with Amino Acid and Rotamer
definitions as used by PDB2PQR.

Authors:  Jens Erik Nielsen, Todd Dolinsky, Yong Huang
"""
import copy
import re
from xml import sax
from . import structures


class DefinitionHandler(sax.ContentHandler):
    """Handle definition XML file content."""

    def __init__(self):
        self.curelement = ""
        self.curatom = None
        self.curholder = None
        self.curobj = None
        self.map = {}
        self.patches = []
        return

    def startElement(self, name, _):
        if name == "residue":
            obj = DefinitionResidue()
            self.curholder = obj
            self.curobj = obj
        elif name == "patch":
            obj = Patch()
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
        if name == "residue": # Complete Residue object
            residue = self.curholder
            if not isinstance(residue, DefinitionResidue):
                raise RuntimeError("Internal error parsing XML!")
            resname = residue.name
            if resname == "":
                raise KeyError("Residue name not set in XML!")
            else:
                self.map[resname] = residue
                self.curholder = None
                self.curobj = None

        elif name == "patch": # Complete patch object
            patch = self.curholder
            if not isinstance(patch, Patch):
                raise RuntimeError("Internal error parsing XML!")
            patchname = patch.name
            if patchname == "":
                raise KeyError("Residue name not set in XML!")
            else:
                self.patches.append(patch)
                self.curholder = None
                self.curobj = None

        elif name == "atom": # Complete atom object
            atom = self.curatom
            if not isinstance(atom, DefinitionAtom):
                raise RuntimeError("Internal error parsing XML!")
            atomname = atom.name
            if atomname == "":
                raise KeyError("Atom name not set in XML!")
            else:
                self.curholder.map[atomname] = atom
                self.curatom = None
                self.curobj = self.curholder

        else: # Just free the current element namespace
            self.curelement = ""

        return self.map

    def characters(self, text):
        if text.isspace():
            return

        # If this is a float, make it so
        try:
            value = float(str(text))
        except ValueError:
            value = str(text)

        # Special cases - lists and dictionaries
        if self.curelement == "bond":
            self.curobj.bonds.append(value)
        elif self.curelement == "dihedral":
            self.curobj.dihedrals.append(value)
        elif self.curelement == "altname":
            self.curholder.altnames[value] = self.curatom.name
        elif self.curelement == "remove":
            self.curobj.remove.append(value)
        else:
            setattr(self.curobj, self.curelement, value)
        return


class Definition(object):
    """The Definition class contains the structured definitions found in the
    files and several mappings for easy access to the information.
    """

    def __init__(self, aa_file, na_file, patch_file):
        """Initialize object.

        Args:
            aa_file:  file object with amino acid definitions
            na_file:  file object with nucleic acid definitions
            patch_file:  file object with patch definitions
        """
        self.map = {}
        self.patches = {}

        handler = DefinitionHandler()
        sax.make_parser()

        for def_file in [aa_file, na_file]:
            sax.parseString(def_file.read(), handler)
            self.map.update(handler.map)

        handler.map = {}
        sax.parseString(patch_file.read(), handler)

        # Apply specific patches to the reference object, allowing users
        # to specify protonation states in the PDB file
        for patch in handler.patches:
            if patch.newname != "":

                # Find all residues matching applyto
                resnames = list(self.map.keys())
                for name in resnames:
                    regexp = re.compile(patch.applyto).match(name)
                    if not regexp:
                        continue
                    newname = patch.newname.replace("*", name)
                    self.add_patch(patch, name, newname)

            # Either way, make sure the main patch name is available
            self.add_patch(patch, patch.applyto, patch.name)

    def add_patch(self, patch, refname, newname):
        """Add a patch to a definition residue.

        Args:
            patch:  The patch object to add (Patch)
            refname:  The name of the object to add the patch to (string)
            newname:  The name of the new (patched) object (string)
        """
        try:
            aadef = self.map[refname] # The reference
            patch_residue = copy.deepcopy(aadef)

            # Add atoms from patch
            for atomname in patch.map:
                patch_residue.map[atomname] = patch.map[atomname]
                for bond in patch.map[atomname].bonds:
                    if bond not in patch_residue.map:
                        continue
                    if atomname not in patch_residue.map[bond].bonds:
                        patch_residue.map[bond].bonds.append(atomname)

            # Rename atoms as directed
            for key in patch.altnames:
                patch_residue.altnames[key] = patch.altnames[key]

            # Remove atoms as directed
            for remove in patch.remove:
                if not patch_residue.has_atom(remove):
                    continue
                removebonds = patch_residue.map[remove].bonds
                del patch_residue.map[remove]
                for bond in removebonds:
                    if remove in patch_residue.map[bond].bonds:
                        patch_residue.map[bond].bonds.remove(remove)

            # Add the new dihedrals
            for dihedral in patch.dihedrals:
                patch_residue.dihedrals.append(dihedral)

            # Point at the new reference
            self.map[newname] = patch_residue

            # Store the patch
            self.patches[newname] = patch

        except KeyError: # Just store the patch
            self.patches[newname] = patch


class Patch(object):
    """Patch the definitionResidue class"""
    def __init__(self):
        self.name = ""
        self.applyto = ""
        self.map = {}
        self.remove = []
        self.altnames = {}
        self.dihedrals = []
        self.newname = ""

    def __str__(self):
        """
            A basic string representation for debugging
        """
        text = "%s\n" % self.name
        text += "Apply to: %s\n" % self.applyto
        text += "Atoms to add: \n"
        for atom in self.map:
            text += "\t%s\n" % str(self.map[atom])
        text += "Atoms to remove: \n"
        for remove in self.remove:
            text += "\t%s\n" % remove
        text += "Alternate naming map: \n"
        text += "\t%s\n" % self.altnames
        return text


class DefinitionResidue(structures.Residue):
    """The DefinitionResidue class extends the Residue class to allow for a
    trimmed down initializing function."""
    def __init__(self):
        self.name = ""
        self.dihedrals = []
        self.map = {}
        self.altnames = {}

    def __str__(self):
        text = "%s\n" % self.name
        text += "Atoms: \n"
        for atom in self.map:
            text += "\t%s\n" % str(self.map[atom])
        text += "Dihedrals: \n"
        for dihedral in self.dihedrals:
            text += "\t%s\n" % dihedral
        text += "Alternate naming map: \n"
        text += "\t%s\n" % self.altnames
        return text

    def get_nearest_bonds(self, atomname):
        bonds = []
        lev2bonds = []
        atom = self.map[atomname]

        # Get directly bonded (length = 1) atoms
        for bondedatom in atom.bonds:
            if bondedatom not in bonds:
                bonds.append(bondedatom)

        # Get bonded atoms 2 bond lengths away
        for bondedatom in atom.bonds:
            for bond2 in self.map[bondedatom].bonds:
                if bond2 not in bonds and bond2 != atomname:
                    bonds.append(bond2)
                    lev2bonds.append(bond2)

        # Get bonded atoms 3 bond lengths away
        for lev2atom in lev2bonds:
            for bond3 in self.map[lev2atom].bonds:
                if bond3 not in bonds:
                    bonds.append(bond3)

        return bonds


class DefinitionAtom(structures.Atom):
    """A trimmed down version of the Atom class"""

    def __init__(self, name=None, x=None, y=None, z=None):
        self.name = name
        self.x = x
        self.y = y
        self.z = z
        if name is None:
            self.name = ""
        if x is None:
            self.x = 0.0
        if y is None:
            self.y = 0.0
        if z is None:
            self.z = 0.0
        self.bonds = []

    def __str__(self):
        text = "%s: %.3f %.3f %.3f" % (self.name, self.x, self.y, self.z)
        for bond in self.bonds:
            text += " %s" % bond
        return text

    @property
    def is_backbone(self):
        """Return true if atom name is in backbone, otherwise false"""
        if self.name in structures.BACKBONE:
            return True
        return False
