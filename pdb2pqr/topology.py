"""Parser for topology XML files.

.. todo:
   There are some very silly classes in this module that should be combined
   or removed.

.. codeauthor::  Nathan Baker
.. codeauthor::  Yong Huang
"""

import logging
from xml import sax

_LOGGER = logging.getLogger(__name__)
#: File name of topology XML data
TOPOLOGYPATH = "TOPOLOGY.xml"


class TopologyHandler(sax.ContentHandler):
    """Handler for XML-based topology files.

    Assumes the following hierarchy of tags:

    * topology

      * residue

        * reference
        * titrationstate

          * tautomer

            * conformer

    .. todo:  Why isn't ``super().__init__`` called here?
    """

    def __init__(self):
        self.curr_element = None
        self.curr_atom = None
        self.curr_dihedral = None
        self.curr_reference = None
        self.curr_residue = None
        self.curr_titration_state = None
        self.curr_tautomer = None
        self.curr_conformer = None
        self.curr_conformer_add = None
        self.curr_conformer_remove = None
        self.residues = []
        self.incomplete = 0

    def startElement(self, tag_name, _):
        """Start element parsing.

        .. note:  overrides super-class method invocation

        :param tag_name:  name of XML element tag to start parsing
        :type tag_name:  str
        """
        if self.incomplete:
            return
        if tag_name == "topology":
            pass
        elif tag_name == "residue":
            if self.curr_residue is not None:
                _LOGGER.info("** Overwriting current topology_residue object!")
            self.curr_residue = TopologyResidue(self)
        elif tag_name == "reference":
            if self.curr_reference is not None:
                _LOGGER.info(
                    "** Overwriting current TopologyReference object!"
                )
            self.curr_reference = TopologyReference(self.curr_residue)
        elif tag_name == "titrationstate":
            if self.curr_titration_state is not None:
                _LOGGER.info(
                    "** Overwriting current topology_titration_state "
                    "object!"
                )
            self.curr_titration_state = TopologyTitrationState(
                self.curr_residue
            )
        elif tag_name == "tautomer":
            if self.curr_tautomer is not None:
                _LOGGER.info("** Overwriting current Tautomer object!")
            self.curr_tautomer = TopologyTautomer(self.curr_titration_state)
        elif tag_name == "conformer":
            if self.curr_conformer is not None:
                _LOGGER.info("** Overwriting current Conformer object!")
            self.curr_conformer = TopologyConformer(self.curr_tautomer)
        elif tag_name == "name":
            self.curr_element = tag_name
        elif tag_name == "atom":
            if self.curr_conformer_add is not None:
                self.curr_atom = TopologyAtom(self.curr_conformer_add)
            elif self.curr_conformer_remove is not None:
                self.curr_atom = TopologyAtom(self.curr_conformer_remove)
            elif self.curr_reference is not None:
                self.curr_atom = TopologyAtom(self.curr_reference)
            else:
                _LOGGER.info("** Don't know what to do with this atom!")
        elif tag_name in ("x", "y", "z", "bond", "altname"):
            self.curr_element = tag_name
        elif tag_name == "dihedral":
            self.curr_element = tag_name
            if self.curr_conformer_add is not None:
                self.curr_dihedral = TopologyDihedral(self.curr_conformer_add)
            elif self.curr_conformer_remove is not None:
                self.curr_dihedral = TopologyDihedral(
                    self.curr_conformer_remove
                )
            elif self.curr_reference is not None:
                self.curr_dihedral = TopologyDihedral(self.curr_reference)
            else:
                _LOGGER.info("** Don't know what to do with this dihedral!")
        elif tag_name == "add":
            self.curr_conformer_add = TopologyConformerAdd(self.curr_conformer)
        elif tag_name == "remove":
            self.curr_conformer_remove = TopologyConformerRemove(
                self.curr_conformer
            )
        elif tag_name == "incomplete":
            self.incomplete = 1
        else:
            _LOGGER.info(f"** NOT handling {tag_name} start tag")

    def endElement(self, tag_name):
        """End parsing element.

        .. note: Overrides super-class method.

        :param tag_name:  name of XML element tag for which to end parsing
        :type tag_name:  str
        """
        if not self.incomplete:
            self.curr_element = None
            if tag_name in ("x", "y", "z", "name", "bond", "altname"):
                pass
            elif tag_name == "atom":
                self.curr_atom = None
            elif tag_name == "dihedral":
                self.curr_dihedral = None
            elif tag_name == "reference":
                self.curr_reference = None
            elif tag_name == "add":
                self.curr_conformer_add = None
            elif tag_name == "remove":
                self.curr_conformer_remove = None
            elif tag_name == "titrationstate":
                self.curr_titration_state = None
            elif tag_name == "conformer":
                self.curr_conformer = None
            elif tag_name == "tautomer":
                self.curr_tautomer = None
            elif tag_name == "residue":
                self.curr_residue = None
            elif tag_name != "topology":
                _LOGGER.info(f"** NOT handling {tag_name} end tag")
        elif tag_name == "incomplete":
            self.incomplete = 0

    def characters(self, text):
        """Parse characters in topology XML file.

        :param text:  XML character data
        :type text:  str
        """
        if text.isspace():
            return
        if not self.incomplete:
            if self.curr_element == "name":
                if self.curr_atom is not None:
                    self.curr_atom.name = text
                elif self.curr_conformer is not None:
                    self.curr_conformer.name = text
                elif self.curr_tautomer is not None:
                    self.curr_tautomer.name = text
                elif self.curr_titration_state is not None:
                    self.curr_titration_state.name = text
                elif self.curr_residue is not None:
                    self.curr_residue.name = text
                else:
                    _LOGGER.info(
                        f"    *** Don't know what to do with name {text}!"
                    )
            elif self.curr_element == "x":
                self.curr_atom.x = float(text)
            elif self.curr_element == "y":
                self.curr_atom.y = float(text)
            elif self.curr_element == "z":
                self.curr_atom.z = float(text)
            elif self.curr_element == "bond":
                self.curr_atom.bonds.append(text)
            elif self.curr_element == "altname":
                self.curr_atom.altname = text
            elif self.curr_element == "dihedral":
                self.curr_dihedral.atom_list = text
            else:
                _LOGGER.info(f"** NOT handling character text:  {text}")


class TopologyResidue:
    """A class for residue topology information.

    .. todo:
       This class (and other classes in this module) has (have) lots of
       repeated code that could be eliminated with better inheritance
       structures.
    """

    def __init__(self, topology_):
        """Initialize with a Topology object.

        :param topology_:  topology object
        """
        self.name = None
        self.reference = None
        self.titration_states = []
        self.topology = topology_
        self.topology.residues.append(self)

    def __str__(self):
        return self.name


class TopologyDihedral:
    """A class for dihedral angle topology information."""

    def __init__(self, parent):
        """Initialize with a parent that has a dihedral list.

        :param parent:  parent object with a dihedral list
        """
        self.parent = parent
        self.parent.dihedrals.append(self)
        self.atom_list = None

    def __str__(self):
        return self.atom_list


class TopologyAtom:
    """A class for atom topology information."""

    def __init__(self, parent):
        """Initialize with an upper-level class that contains an atom array.

        For example, initialize with :class:`TopologyReference` or
        :class:`TopologyConformerAddition`

        :param parent:  parent object with atom array
        """
        self.parent = parent
        self.parent.atoms.append(self)
        self.name = None
        self.x = None
        self.y = None
        self.z = None
        self.bonds = []
        self.altname = None

    def __str__(self):
        return self.name


class TopologyTitrationState:
    """A class for the titration state of a residue."""

    def __init__(self, topology_residue):
        """Initialize object.

        :param topology_residue:  residue topology
        :type topology_residue:  TopologyResidue
        """
        self.topology_residue = topology_residue
        self.topology_residue.titration_states.append(self)
        self.tautomers = []
        self.name = None

    def __str__(self):
        return self.name


class TopologyTautomer:
    """A class for topology tautomer information."""

    def __init__(self, topology_titration_state):
        """Initialize object.

        :param topology_titration_state:  titration state of this residue
        :type topology_titration_state:  TopologyTitrationState
        """
        self.topology_titration_state = topology_titration_state
        self.topology_titration_state.tautomers.append(self)
        self.conformers = []
        self.name = None

    def __str__(self):
        return self.name


class TopologyConformer:
    """A class for topology conformer information."""

    def __init__(self, topology_tautomer):
        """Initialize object.

        :param topology_tautomer:  tautomer for which to evaluate conformers
        :type topology_tautomer:  TopologyTautomer
        """
        self.topology_tautomer = topology_tautomer
        self.topology_tautomer.conformers.append(self)
        self.name = None
        self.conformer_adds = []
        self.conformer_removes = []

    def __str__(self):
        return self.name


class TopologyReference:
    """A class for the reference structure of a residue."""

    def __init__(self, topology_residue):
        """Initialize object.

        :param topology_residue:  reference residue
        :type topology_residue:  TopologyResidue
        """
        self.topology_residue = topology_residue
        self.topology_residue.reference = self
        self.atoms = []
        self.dihedrals = []


class TopologyConformerAdd:
    """A class for adding atoms to a conformer."""

    def __init__(self, topology_conformer):
        """Initialize object.

        :param topology_conformer:  conformer to which to add atoms
        :type topology_conformer:  TopologyConformer
        """
        self.topology_conformer = topology_conformer
        self.topology_conformer.conformer_adds.append(self)
        self.atoms = []
        self.name = None
        self.dihedrals = []


class TopologyConformerRemove:
    """A class for removing atoms from a conformer."""

    def __init__(self, topology_conformer):
        """Initialize object.

        :param topology_conformer:  conformer to which to add atoms
        :type topology_conformer:  TopologyConformer
        """
        self.topology_conformer = topology_conformer
        self.topology_conformer.conformer_removes.append(self)
        self.atoms = []
        self.name = None


class Topology:
    """Contains the structured definitions of residue reference coordinates
    as well as alternate titration, conformer, and tautomer states.
    """

    def __init__(self, topology_file):
        """Initialize object.

        :param topology_file:  topology file object ready for reading
        :type topology_file:  file-like object
        """
        handler = TopologyHandler()
        sax.make_parser()
        sax.parseString(topology_file.read(), handler)
        self.residues = handler.residues
