"""Parser for TOPOLOGY.xml

Authors:  Nathan Baker, Yong Huang
"""
from xml import sax
import logging


_LOGGER = logging.getLogger(__name__)
TOPOLOGYPATH = "TOPOLOGY.xml"


class TopologyHandler(sax.ContentHandler):
    """ Handler for XML-based topology files.  Assumes the following hierarchy of tags:
    topology
    -->residue
       |-->reference
       |-->titrationstate
           |-->tautomer
               |-->conformer
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

    def startElement(self, tagName, _):
        """Start element parsing."""
        if not self.incomplete:
            if tagName == "topology":
                pass
            elif tagName == "residue":
                if self.curr_residue != None:
                    _LOGGER.info("** Overwriting current topology_residue object!")
                self.curr_residue = TopologyResidue(self)
            elif tagName == "reference":
                if self.curr_reference != None:
                    _LOGGER.info("** Overwriting current TopologyReference object!")
                self.curr_reference = TopologyReference(self.curr_residue)
            elif tagName == "titrationstate":
                if self.curr_titration_state != None:
                    _LOGGER.info("** Overwriting current topology_titration_state object!")
                self.curr_titration_state = TopologyTitrationState(self.curr_residue)
            elif tagName == "tautomer":
                if self.curr_tautomer != None:
                    _LOGGER.info("** Overwriting current Tautomer object!")
                self.curr_tautomer = TopologyTautomer(self.curr_titration_state)
            elif tagName == "conformer":
                if self.curr_conformer != None:
                    _LOGGER.info("** Overwriting current Conformer object!")
                self.curr_conformer = TopologyConformer(self.curr_tautomer)
            elif tagName == "name":
                self.curr_element = tagName
            elif tagName == "atom":
                if self.curr_conformer_add != None:
                    self.curr_atom = TopologyAtom(self.curr_conformer_add)
                elif self.curr_conformer_remove != None:
                    self.curr_atom = TopologyAtom(self.curr_conformer_remove)
                elif self.curr_reference != None:
                    self.curr_atom = TopologyAtom(self.curr_reference)
                else:
                    _LOGGER.info("** Don't know what to do with this atom!")
            elif tagName == "x":
                self.curr_element = tagName
            elif tagName == "y":
                self.curr_element = tagName
            elif tagName == "z":
                self.curr_element = tagName
            elif tagName == "bond":
                self.curr_element = tagName
            elif tagName == "altname":
                self.curr_element = tagName
            elif tagName == "dihedral":
                self.curr_element = tagName
                if self.curr_conformer_add != None:
                    self.curr_dihedral = TopologyDihedral(self.curr_conformer_add)
                elif self.curr_conformer_remove != None:
                    self.curr_dihedral = TopologyDihedral(self.curr_conformer_remove)
                elif self.curr_reference != None:
                    self.curr_dihedral = TopologyDihedral(self.curr_reference)
                else:
                    _LOGGER.info("** Don't know what to do with this dihedral!")
            elif tagName == "add":
                self.curr_conformer_add = TopologyConformerAdd(self.curr_conformer)
            elif tagName == "remove":
                self.curr_conformer_remove = TopologyConformerRemove(self.curr_conformer)
            elif tagName == "incomplete":
                self.incomplete = 1
            else:
                _LOGGER.info("** NOT handling %s start tag", tagName)

    def endElement(self, tagName):
        """End parsing element"""
        if not self.incomplete:
            self.curr_element = None
            if tagName == "x":
                pass
            elif tagName == "y":
                pass
            elif tagName == "z":
                pass
            elif tagName == "name":
                pass
            elif tagName == "bond":
                pass
            elif tagName == "altname":
                pass
            elif tagName == "atom":
                self.curr_atom = None
            elif tagName == "dihedral":
                self.curr_dihedral = None
            elif tagName == "reference":
                self.curr_reference = None
            elif tagName == "add":
                self.curr_conformer_add = None
            elif tagName == "remove":
                self.curr_conformer_remove = None
            elif tagName == "titrationstate":
                self.curr_titration_state = None
            elif tagName == "conformer":
                self.curr_conformer = None
            elif tagName == "tautomer":
                self.curr_tautomer = None
            elif tagName == "residue":
                self.curr_residue = None
            elif tagName == "topology":
                pass
            else:
                _LOGGER.info("** NOT handling %s end tag", tagName)
        else:
            if tagName == "incomplete":
                self.incomplete = 0

    def characters(self, text):
        """Parse characters in topology file."""
        if text.isspace():
            return
        if not self.incomplete:
            if self.curr_element == "name":
                if self.curr_atom != None:
                    self.curr_atom.name = text
                elif self.curr_conformer != None:
                    self.curr_conformer.name = text
                elif self.curr_tautomer != None:
                    self.curr_tautomer.name = text
                elif self.curr_titration_state != None:
                    self.curr_titration_state.name = text
                elif self.curr_residue != None:
                    self.curr_residue.name = text
                else:
                    _LOGGER.info("    *** Don't know what to do with name %s!", text)
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
                _LOGGER.info("** NOT handling character text:  %s", text)


# TODO - lots of repeated code that could be eliminated with better inheritance
class TopologyResidue(object):
    """ A class for residue topology information """
    def __init__(self, topology_):
        """ Initialize with a Topology object """
        self.name = None
        self.reference = None
        self.titration_states = []
        self.topology = topology_
        self.topology.residues.append(self)

    def __str__(self):
        return self.name


class TopologyDihedral(object):
    """ A class for dihedral topology information.  """
    def __init__(self, parent):
        """ Needs a parent that has a dihedral list. """
        self.parent = parent
        self.parent.dihedrals.append(self)
        self.atom_list = None

    def __str__(self):
        return self.atom_list


class TopologyAtom(object):
    """ A class for atom topology information """
    def __init__(self, parent):
        """ Needs to be intialized with an upper-level class that contains an
        atoms array (e.g., TopologyReference or TopologyConformerAddition)"""
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


class TopologyTitrationState(object):
    """ A class for the titration state of a residue """
    def __init__(self, topology_residue):
        """ Initialize with a topology_residue object """
        self.topology_residue = topology_residue
        self.topology_residue.titration_states.append(self)
        self.tautomers = []
        self.name = None

    def __str__(self):
        return self.name


class TopologyTautomer(object):
    """ A class for topology tautomer information """
    def __init__(self, topology_titration_state):
        """ Initialize with a topology_titration_state object """
        self.topology_titration_state = topology_titration_state
        self.topology_titration_state.tautomers.append(self)
        self.conformers = []
        self.name = None

    def __str__(self):
        return self.name


class TopologyConformer(object):
    """ A class for topology conformer information """
    def __init__(self, topology_tautomer):
        """ Initialize with a TopologyTautomer object """
        self.topology_tautomer = topology_tautomer
        self.topology_tautomer.conformers.append(self)
        self.name = None
        self.conformer_adds = []
        self.conformer_removes = []

    def __str__(self):
        return self.name


class TopologyReference(object):
    """ A class for the reference structure of a residue """
    def __init__(self, topology_residue):
        """ Initialize with a topology_residue object """
        self.topology_residue = topology_residue
        self.topology_residue.reference = self
        self.atoms = []
        self.dihedrals = []


class TopologyConformerAdd(object):
    """ A class for adding atoms to a conformer """
    def __init__(self, topology_conformer):
        """ Initialize with a TopologyConformer object """
        self.topology_conformer = topology_conformer
        self.topology_conformer.conformer_adds.append(self)
        self.atoms = []
        self.name = None
        self.dihedrals = []


class TopologyConformerRemove(object):
    """ A class for removing atoms to a conformer """
    def __init__(self, topology_conformer):
        """ Initialize with a TopologyConformer object """
        self.topology_conformer = topology_conformer
        self.topology_conformer.conformer_removes.append(self)
        self.atoms = []
        self.name = None


class Topology(object):
    """ Contains the structured definitions of residue reference coordinates
    as well as alternate titration, conformer, and tautomer states.
    """
    def __init__(self, topology_file):
        """ Initialize with the topology file reference ready for reading """
        handler = TopologyHandler()
        sax.make_parser()
        sax.parseString(topology_file.read(), handler)
        self.residues = handler.residues


if __name__ == "__main__":
    with open(TOPOLOGYPATH, "r") as tfile_:
        Topology(tfile_)
