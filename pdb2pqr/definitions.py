"""XML handling for biomolecular residue topology definitions.

.. codeauthor::  Jens Erik Nielsen
.. codeauthor::  Todd Dolinsky
.. codeauthor::  Yong Huang
"""

import copy
import logging
import re
from xml import sax

from . import residue, structures

_LOGGER = logging.getLogger(__name__)


class DefinitionHandler(sax.ContentHandler):
    """Handle definition XML file content."""

    def __init__(self):
        self.curelement = ""
        self.curatom = None
        self.curholder = None
        self.curobj = None
        self.map = {}
        self.patches = []
        self.content = ""

    def startElement(self, name, _):
        """Start XML element parsing.

        :param name:  element name
        :type name:  str
        """
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

    def endElement(self, name):
        """End XML element parsing.

        :param name:  element name
        :type name:  str
        :raises KeyError:  for invalid atom or residue names
        :raises RuntimeError:  for unexpected XML features or format
        """
        if name == "residue":  # Complete Residue object
            residue_ = self.curholder
            if not isinstance(residue_, DefinitionResidue):
                raise RuntimeError("Internal error parsing XML!")
            resname = residue_.name
            if resname == "":
                raise KeyError("Residue name not set in XML!")
            else:
                self.map[resname] = residue_
                self.curholder = None
                self.curobj = None
        elif name == "patch":  # Complete patch object
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
        elif name == "atom":  # Complete atom object
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
        elif (self.curobj is not None) and (self.content.strip() != ""):
            self.content = self.content.strip()
            _LOGGER.debug(
                f"Got text for {self.curholder.name} <{name}>: {self.content}"
            )
            if name == "bond":
                self.curobj.bonds.append(self.content)
            elif name == "dihedral":
                self.curobj.dihedrals.append(self.content)
            elif name == "altname":
                self.curholder.altnames[self.content] = self.curatom.name
            elif name == "remove":
                self.curobj.remove.append(self.content)
            elif name == "name":
                self.curobj.name = self.content
            elif self.curelement != "":
                try:
                    self.content = float(self.content)
                except ValueError:
                    pass
                setattr(self.curobj, self.curelement, self.content)

        self.content = ""
        self.curelement = ""

        return self.map

    def characters(self, text):
        """Parse text data in XML.

        :param text:  text data to parse
        :type text:  str
        """
        self.content += text


class Definition:
    """Force field topology definitions.

    The Definition class contains the structured definitions found in the files
    and several mappings for easy access to the information.
    """

    def __init__(self, aa_file, na_file, patch_file):
        """Initialize object.

        :param aa_file:  file-like object with amino acid definitions
        :type aa_file:  file
        :param na_file:  file-like object with nucleic acid definitions
        :type na_file:  file
        :param patch_file:  file-like object with patch definitions
        :type patch_file:  file
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
        """Add a patch to a topology definition residue.

        :param patch:  the patch object to add
        :type patch:  Patch
        :param refname:  the name of the object to add the patch to
        :type refname:  str
        :param newname:  the name of the new (patched) object
        :type newname:  str
        """
        try:
            aadef = self.map[refname]  # The reference
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
        except KeyError:  # Just store the patch
            self.patches[newname] = patch


class Patch:
    """Residue patches for structure topologies."""

    def __init__(self):
        self.name = ""
        self.applyto = ""
        self.map = {}
        self.remove = []
        self.altnames = {}
        self.dihedrals = []
        self.newname = ""

    def __str__(self):
        text = f"{self.name}\n"
        text += f"Apply to: {self.applyto}\n"
        text += "Atoms to add: \n"
        for atom in self.map:
            text += f"\t{self.map[atom]!s}\n"
        text += "Atoms to remove: \n"
        for remove in self.remove:
            text += f"\t{remove}\n"
        text += "Alternate naming map: \n"
        text += f"\t{self.altnames}\n"
        return text


class DefinitionResidue(residue.Residue):
    """Force field toplogy representation for a residue."""

    def __init__(self):
        self.name = ""
        self.dihedrals = []
        self.map = {}
        self.altnames = {}

    def __str__(self):
        text = f"{self.name}\n"
        text += "Atoms: \n"
        for atom in self.map:
            text += f"\t{self.map[atom]!s}\n"
        text += "Dihedrals: \n"
        for dihedral in self.dihedrals:
            text += f"\t{dihedral}\n"
        text += "Alternate naming map: \n"
        text += f"\t{self.altnames}\n"
        return text

    def get_nearest_bonds(self, atomname):
        """Get bonded atoms near a given atom.

        :param atomname:  name of specific atom
        :type atomname:  str
        :return:  list of nearby bonded atom names
        :rtype:  [str]
        """
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
    """Store force field atom topology definitions."""

    def __init__(self, name=None, x=None, y=None, z=None):
        """Initialize class.

        :param name:  atom name
        :type name:  str
        :param x:  x-coordinate
        :type x:  float
        :param y:  y-coordinate
        :type y:  float
        :param z:  z-coordinate
        :type z:  float
        """
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
        text = f"{self.name}: {self.x:.3f} {self.y:.3f} {self.z:.3f}"
        for bond in self.bonds:
            text += f" {bond}"
        return text

    @property
    def is_backbone(self):
        """Identify whether atom is in backbone.

        :return:  true if atom name is in backbone, otherwise false
        :rtype:  bool
        """
        return self.name in structures.BACKBONE
