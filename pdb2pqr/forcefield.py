"""Force fields are fun.

The forcefield structure is modeled off of the structures.py file, where each
forcefield is considered a chain of residues of atoms.

.. codeauthor::  Todd Dolinsky
.. codeauthor::  Yong Huang
"""

import logging
import re
from xml import sax

from . import io

_LOGGER = logging.getLogger(__name__)


class ForcefieldHandler(sax.ContentHandler):
    """Process XML-format force field parameter files."""

    def __init__(self, map_, reference):
        """Initialize handler

        :param map_:  dictionary of paramaeter information
        :type map_:  dict
        :param reference:  reference map for force field
        :type reference:  dict
        """
        self.oldresname = None
        self.oldatomname = None
        self.curelement = None
        self.newatomname = None
        self.newresname = None
        self.atommap = {}
        self.map = map_
        self.reference = reference

    @classmethod
    def update_map(cls, toname, fromname, map_):
        """Update the given map by adding a pointer from a new name to an object.

        .. todo:: Should this be a staticmethod instead of classmethod?

        :param toname:  the new name for the object
        :type toname:  str
        :param fromname:  the old name for the object
        :type fromname:  str
        :param map_:  a dictionary of forcefield-related items
        :type map_: dict
        """
        fromobj = map_[fromname]
        if isinstance(fromobj, ForcefieldResidue):
            if toname not in map_:
                newres = ForcefieldResidue(fromname)
                map_[toname] = newres
            for atomname in fromobj.atoms:
                map_[toname].atoms[atomname] = fromobj.atoms[atomname]
        elif isinstance(fromobj, ForcefieldAtom):
            map_[toname] = fromobj

    @classmethod
    def find_matching_names(cls, regname, map_):
        """Find strings in the map that match the given regular expression.

        .. todo:: Should this be a staticmethod instead of classmethod?

        :param regname: the regular expression to search for in the map
        :type regname:  str
        :param map_:  the dictionary to search
        :type map_: dict
        :return:  a list of regular expression match objects for dictionary
            keys that match the regular expression.
        :rtype:  [re.Match]
        """
        name_list = []
        regname += "$"
        # Find the existing items that match this string
        for name in map_:
            regexp = re.compile(regname).match(name)
            if regexp:
                name_list.append(regexp)
        return name_list

    def startElement(self, name, _):
        """Start XML element parsing.

        :param name:  element name
        :type name:  str
        """
        if name != "name":
            self.curelement = name

    def endElement(self, name):
        """End XML element parsing.

        :param name:  the name of the element
        :type name:  str
        """
        if name == "residue":
            if self.oldresname is not None:  # Make a new residue hook
                newreslist = self.find_matching_names(
                    self.newresname, self.reference
                )
                # Multiple new residues
                if self.oldresname.find("$group") >= 0:
                    for resitem in newreslist:
                        resname = resitem.string
                        group = resitem.group(1)
                        fromname = self.oldresname.replace("$group", group)
                        if fromname in self.map:
                            self.update_map(resname, fromname, self.map)
                else:  # Work with a single new residue name
                    _ = self.find_matching_names(self.oldresname, self.map)
                    for resitem in newreslist:
                        resname = resitem.string
                        self.update_map(resname, self.oldresname, self.map)
            # If this was only a residue conversion, exit
            if self.atommap == {}:
                self.oldresname = None
                self.newresname = None
                return None
            # Apply atom conversions for all appropriate residues
            resmatchlist = self.find_matching_names(self.newresname, self.map)
            for resitem in resmatchlist:
                residue = self.map[resitem.string]
                for newname in self.atommap:
                    oldname = self.atommap[newname]
                    if oldname not in residue.atoms:
                        continue
                    self.update_map(newname, oldname, residue.atoms)
            # Clean up
            self.oldresname = None
            self.newresname = None
            self.atommap = {}
        elif name == "atom":
            self.atommap[self.newatomname] = self.oldatomname
            self.oldatomname = None
            self.newatomname = None
        else:  # Just free the current element namespace
            self.curelement = ""
        return self.map

    def characters(self, text):
        """Parse text information within XML tag.

        :param text:  the text value between the XML tags
        :type test:  str
        """
        if text.isspace():
            return
        text = str(text)
        if self.curelement == "residue":
            self.newresname = text
        elif self.curelement == "atom":
            self.newatomname = text
        elif self.curelement == "useatomname":
            self.oldatomname = text
        elif self.curelement == "useresname":
            self.oldresname = text


class Forcefield:
    """Parameter definitions for a given forcefield.

    .. note::
        Pass ff and ff names file like objects rather than sorting out whether
        to use user-created files here.

    The forcefield class contains definitions for a given forcefield.
    Each forcefield object contains a dictionary of residues, with each residue
    containing a dictionary of atoms.
    Dictionaries are used instead of lists as the ordering is not important.
    The forcefield definition files are unedited, directly from the
    forcefield - all transformations are done within.
    """

    def __init__(self, ff_name, definition, userff, usernames=None):
        """Initialize the class by parsing the definition file.

        .. todo:: Why are files being loaded so deep in this function?

        :param ff_name: the name of the forcefield (can be None)
        :type ff_name:  str
        :param definition:  the definition object for the Forcefield
        :type definition:  Definition
        :param userff:  path for user-defined forcefields file
        :type userff:  str
        :param usernames:  path to user-defined atom/residue names file
        :type usernames:  str
        :raises ValueError:  if invalid force field names specified
        """
        self.map = {}
        self.name = str(ff_name)
        defpath = ""
        defpath = io.test_dat_file(ff_name) if userff is None else userff
        with open(defpath, encoding="utf-8") as ff_file:
            lines = ff_file.readlines()
            for line in lines:
                if not line.startswith("#"):
                    fields = line.split()
                    if fields == []:
                        continue
                    try:
                        resname = fields[0]
                        atomname = fields[1]
                        charge = float(fields[2])
                        radius = float(fields[3])
                    except ValueError:
                        txt = (
                            "Unable to recognize user-defined forcefield file"
                        )
                        txt += f" {defpath}!" if defpath != "" else "!"
                        txt += " Please use a valid parameter file."
                        raise ValueError(txt)
                    try:
                        group = fields[4]
                        atom = ForcefieldAtom(
                            atomname, charge, radius, resname, group
                        )
                    except IndexError:
                        atom = ForcefieldAtom(
                            atomname, charge, radius, resname
                        )

                    my_residue = self.get_residue(resname)
                    if my_residue is None:
                        my_residue = ForcefieldResidue(resname)
                        self.map[resname] = my_residue
                    my_residue.add_atom(atom)
        # Now parse the XML file, associating with FF objects -
        # This is not necessary (if canonical names match ff names)
        try:
            defpath = io.test_names_file(ff_name)
        except FileNotFoundError:
            defpath = None
        if usernames:
            names_path = usernames
        elif defpath:
            names_path = defpath
        else:
            raise ValueError("Unable to identify .names file.")
        handler = ForcefieldHandler(self.map, definition.map)
        sax.make_parser()
        with open(names_path, encoding="utf-8") as namesfile:
            sax.parseString(namesfile.read(), handler)

    def has_residue(self, resname):
        """Check if the residue name is in the map or not.

        :param resname:  the residue name to search
        :type resname:  str
        :return:  indication of whether resname is in map
        :rtype:  bool
        """
        return resname in self.map

    def get_residue(self, resname):
        """Return the residue object with the given resname.

        :param resname:  the name of the residue
        :type resname:  str
        :return:  residue object
        :rtype:  ForcefieldResidue
        """
        if self.has_residue(resname):
            return self.map[resname]
        return None

    def get_names(self, resname, atomname):
        """Get forcefield names associated with residue and atom.

        The names passed in point to ForcefieldResidue and ForcefieldAtom
        objects which may have different names; grab these names and return.

        :param resname:  the residue name
        :type resname:  str
        :param atomname:  the atom name
        :type atomname:  str
        :return:  (forcefield's name for this residue, forcefield's name for
            this atom)
        :rtype:  (str, str)
        """
        rname = None
        aname = None
        if resname in self.map:
            res = self.map[resname]
            if res.has_atom(atomname):
                atom = res.atoms[atomname]
                aname = atom.name
                rname = atom.resname
        return rname, aname

    def get_group(self, resname, atomname):
        """Get the group/type associated with the input fields.

        :param resname:  the residue name
        :type resname:  str
        :param atomname:  the atom name
        :type atomname:  str
        :return:  group name or empty string
        :rtype:  str
        """
        group = ""
        if resname in self.map:
            resid = self.map[resname]
            if resid.has_atom(atomname):
                atom = resid.atoms[atomname]
                group = atom.group
        return group

    def get_params(self, resname, atomname):
        """Get the charge and radius parameters for an atom in a residue.

        The residue itself is needed instead of simply its name because the
        forcefield may use a different residue name than the standard amino
        acid name.

        .. todo::
           Why do both :func:`get_params` and :func:`get_params1` exist?

        :param resname:  the residue name
        :type resname:  str
        :param atomname:  the atom name
        :type atomname:  str
        :return:  (charge of the atom, radius of the atom)
        :rtype:  (float, float)
        """
        charge = None
        radius = None
        if resname in self.map:
            resid = self.map[resname]
            if resid.has_atom(atomname):
                atom = resid.atoms[atomname]
                charge = atom.charge
                radius = atom.radius
        return charge, radius

    def get_params1(self, residue, name):
        """Get the charge and radius parameters for an atom in a residue.

        The residue itself is needed instead of simply its name because the
        forcefield may use a different residue name than the standard amino
        acid name.

        .. todo::
           Why do both :func:`get_params` and :func:`get_params1` exist?

        :param resname:  the residue name
        :type resname:  str
        :param name:  the atom name
        :type name:  str
        :return:  (charge of the atom, radius of the atom)
        :rtype:  (float, float)
        """
        charge = None
        radius = None
        resname = ""
        atomname = ""
        if self.name == "amber":
            resname, atomname = self.get_amber_params(residue, name)
        elif self.name == "charmm":
            resname, atomname = self.get_charmm_params(residue, name)
        elif self.name == "parse":
            resname, atomname = self.get_parse_params(residue, name)
        else:
            resname = residue.name
            atomname = name
        defresidue = self.get_residue(resname)
        if defresidue is None:
            return charge, radius
        atom = defresidue.get_atom(atomname)
        if atom is not None:
            charge = atom.charge
            radius = atom.radius
        return charge, radius

    @classmethod
    def get_amber_params(cls, residue, name):
        """Get AMBER forcefield definitions.

        .. todo::  Should this be a staticmethod or a classmethod?

        .. todo::
           Figure out why :makevar:`residue.type` has :class:`int` values

        :param residue:  the residue
        :type residue:  Residue
        :param name:  the atom name
        :type name:  str
        :return:  (forcefield name of residue, forcefield name of atom)
        :rtype:  (str, str)
        """
        atomname = name
        res_type = residue.type
        resname = residue.naname if res_type == 4 else residue.name
        # Residue Substitutions
        if residue.name == "CYS" and "HG" not in residue.map:
            resname = "CYX"
        elif residue.name == "HIS":
            if "HD1" in residue.map and "HE2" in residue.map:
                resname = "HIP"
            elif "HD1" in residue.map:
                resname = "HID"
            elif "HE2" in residue.map:
                resname = "HIE"
            else:
                resname = "HID"  # Default for no hydrogens
        elif residue.name == "HSP":
            resname = "HIP"
        elif residue.name == "HSE":
            resname = "HIE"
        elif residue.name == "HSD":
            resname = "HID"
        elif residue.name in ["GLU", "GLH"]:
            if "HE1" in residue.map:
                resname = "GLH"
                if atomname == "HE1":
                    atomname = "HE2"
                elif atomname == "OE1":
                    atomname = "OE2"
                elif atomname == "OE2":
                    atomname = "OE1"
            elif "HE2" in residue.map:
                resname = "GLH"
        elif residue.name in ["ASP", "ASH"]:
            if "HD1" in residue.map:
                resname = "ASH"
                if atomname == "HD1":
                    atomname = "HD2"
                elif atomname == "OD1":
                    atomname = "OD2"
                elif atomname == "OD2":
                    atomname = "OD1"
            elif "HD2" in residue.map:
                resname = "ASH"
        if residue.is_c_term:
            resname = "C" + resname
        elif residue.is_n_term:
            resname = "N" + resname
        # Atom Substitutions
        if resname == "ILE":
            if atomname == "CD":
                atomname = "CD1"
        elif resname == "WAT":
            if atomname == "O":
                atomname = "OW"
            elif atomname in ["H1", "H2"]:
                atomname = "HW"
        # N-terminal
        if resname[0] == "N" and resname != "NME" and atomname == "H":
            atomname = "H1"
        if resname in ["CCYS", "NCYS"] and atomname == "HG":
            atomname = "HSG"
        if resname == "CYM" and atomname == "H":
            atomname = "HN"
        if residue.is_n_term and resname == "NPRO" and atomname == "HN2":
            atomname = "H2"
        if residue.is_n_term and resname == "NPRO" and atomname == "HN1":
            atomname = "H3"
        return resname, atomname

    @classmethod
    def get_parse_params(cls, residue, name):
        """Get PARSE forcefield definitions.

        .. todo::  Should this be a staticmethod or a classmethod?

        :param residue:  the residue
        :type residue:  Residue
        :param name:  the atom name
        :type name:  str
        :return:  (forcefield name of residue, forcefield name of atom)
        :rtype:  (str, str)
        """
        # TODO - this is a crazy list of conditionals!
        atomname = name
        resname = residue.name
        # Terminal/Water Substitutions
        nterm = residue.is_n_term
        cterm = residue.is_c_term
        if nterm and resname != "ACE":
            if resname == "PRO" and nterm == 2:
                resname = "PR+"
                if atomname == "H2":
                    atomname = "HN1"
                elif atomname == "H3":
                    atomname = "HN2"
            elif resname == "PRO" and nterm == 1:
                resname = "PRN"
                if atomname in ["H2", "H3"]:
                    atomname = "HN"
            elif nterm == 2:  # Neutral
                # TODO - there are a lot of hard-coded repeated lists like this
                # They should be replaced with module-level variables.
                if atomname in ["N", "H", "H2", "H3", "CA", "HA", "C", "O"]:
                    resname = "BKN"
                if atomname == "H":
                    atomname = "H1"
                if atomname == "H3":
                    atomname = "H2"
            elif nterm == 3:  # Positive
                if atomname in ["N", "H", "H2", "H3", "CA", "HA", "C", "O"]:
                    resname = "BK+"
                if atomname == "H":
                    atomname = "H1"
        elif cterm:
            if atomname == "O":
                atomname = "O1"
            elif atomname == "OXT":
                atomname = "O2"
            if cterm == 1 and atomname in [
                "N",
                "H",
                "HA",
                "CA",
                "C",
                "O1",
                "O2",
            ]:
                resname = "BK-"
            elif cterm == 2 and atomname in [
                "N",
                "H",
                "HA",
                "CA",
                "C",
                "O1",
                "O2",
                "HO",
            ]:
                if atomname == "HO":
                    atomname = "H2"
                resname = "BKC"
            # print 'Cterm resname is',resname
        elif residue.type == 3:
            resname = "H2O"
            if atomname == "O":
                atomname = "OH"
            elif atomname == "H1":
                atomname = "HH1"
            elif atomname == "H2":
                atomname = "HH2"
        # Residue Substitutions
        if resname == "HSD":
            resname = "HID"
        elif resname in ["HIE", "HSE"]:
            resname = "HIS"
        elif resname in ["HIP", "HSP"]:
            resname = "HI+"
        elif resname == "ILE":
            if atomname == "HG12":
                atomname = "HG11"
            elif atomname == "HG13":
                atomname = "HG12"
            elif atomname == "CD":
                atomname = "CD1"
        elif resname == "CYS" and "HG" not in residue.map:
            resname = "CSS"
        # Histidine
        elif resname == "HIS":
            if "HD1" in residue.map and "HE2" in residue.map:
                resname = "HI+"
            elif "HD1" in residue.map:
                resname = "HID"
            elif "HE2" in residue.map:
                resname = "HIS"
        elif resname in ["GLU", "GLH"]:
            if "HE1" in residue.map:
                resname = "GL0"
                if atomname == "HE1":
                    atomname = "HE2"
                elif atomname == "OE1":
                    atomname = "OE2"
                elif atomname == "OE2":
                    atomname = "OE1"
            elif "HE2" in residue.map:
                resname = "GL0"
        elif resname in ["ASP", "ASH"]:
            if "HD1" in residue.map:
                resname = "AS0"
                if atomname == "HD1":
                    atomname = "HD2"
                elif atomname == "OD1":
                    atomname = "OD2"
                elif atomname == "OD2":
                    atomname = "OD1"
            elif "HD2" in residue.map:
                resname = "AS0"
        elif resname == "ACE":
            if atomname == "HH31":
                atomname = "HA1"
            elif atomname == "HH32":
                atomname = "HA2"
            elif atomname == "HH33":
                atomname = "HA3"
            elif atomname == "CH3":
                atomname = "CA"
        elif resname == "TYR":
            if "HH" not in residue.map:
                resname = "TYM"
        elif resname == "TYM":
            resname = "TY-"
        elif resname == "CYM":
            resname = "CY-"
        elif resname == "LYN":
            resname = "LY0"
        # Neutral LYS and neutral ARG detection based on hydrogens - added by
        # Jens
        elif resname == "LYS":
            if "HZ3" not in residue.map:
                resname = "LY0"
        elif resname == "ARG":
            if "HE" not in residue.map:
                resname = "AR0"
        elif resname == "NME":
            resname = "N-M"
            if atomname == "CH3":
                atomname = "CA"
            elif atomname == "H":
                atomname = "H1"
            elif atomname.startswith("HH"):
                atomname = "HA" + atomname[-1]
        # Hydrogen Substitutions
        if atomname == "H":
            atomname = "HN"
        elif atomname == "HA2":
            atomname = "HA1"
        elif atomname == "HA3":
            atomname = "HA2"
        elif atomname == "HB2" and resname not in ["ALA"]:
            atomname = "HB1"
        elif atomname == "HB3" and resname not in ["ALA"]:
            atomname = "HB2"
        elif atomname == "HD2" and resname not in ["HIS", "HI+", "HID", "AS0"]:
            atomname = "HD1"
        elif atomname == "HD3" and resname not in ["HIS", "HI+", "HID"]:
            atomname = "HD2"
        elif atomname == "HE2" and resname not in [
            "TRP",
            "HIS",
            "HI+",
            "HID",
            "GL0",
        ]:
            atomname = "HE1"
        elif atomname == "HE3" and resname not in ["TRP", "HIS", "HI+", "HID"]:
            atomname = "HE2"
        elif atomname == "HG2":
            atomname = "HG1"
        elif atomname == "HG3":
            atomname = "HG2"
        elif atomname == "HZ2" and resname == "LY0":
            atomname = "HZ1"
        elif atomname == "HZ3" and resname == "LY0":
            atomname = "HZ2"
        return resname, atomname

    @classmethod
    def get_charmm_params(cls, residue, name):
        """Get CHARMM forcefield definitions.

        .. todo::  Should this be a staticmethod or a classmethod?

        .. todo::
           Figure out why :makevar:`residue.type` has :class:`int` values

        :param residue:  the residue
        :type residue:  Residue
        :param name:  the atom name
        :type name:  str
        :return:  (forcefield name of residue, forcefield name of atom)
        :rtype:  (str, str)
        """
        resname = residue.name
        atomname = name
        #  Nucleic Acid Substitutions
        if residue.type == 4:
            resname = resname[0]
            if resname == "A":
                resname = "ADE"
            elif resname == "C":
                resname = "CYT"
            elif resname == "G":
                resname = "GUA"
            elif resname == "T":
                resname = "THY"
                if atomname == "C7":
                    atomname = "C5M"
                elif atomname == "H71":
                    atomname = "H51"
                elif atomname == "H72":
                    atomname = "H52"
                elif atomname == "H73":
                    atomname = "H53"
            elif resname == "U":
                resname = "URA"
            if atomname == "H5'1":
                atomname = "H5'"
            elif atomname == "H5'2":
                atomname = "H5''"
            elif atomname == "H2'1":
                atomname = "H2'"
            elif atomname in ["H2'2", "HO'2"]:
                atomname = "H2''"
            if residue.get_atom("O2'") is None and atomname in [
                "C2'",
                "H2'",
                "H2''",
            ]:
                resname = "DEO1"
            if residue.get_atom("H5T") is not None and atomname in [
                "H5T",
                "O5'",
                "C5'",
            ]:
                resname = "5TER"
            if residue.get_atom("H3T") is not None and atomname in [
                "H3T",
                "O3'",
                "C3'",
            ]:
                resname = "3TER"
        # Terminal/Water Substitutions
        if residue.is_n_term:
            if resname == "GLY" and atomname in [
                "N",
                "H",
                "H2",
                "H3",
                "CA",
                "HA2",
                "HA3",
            ]:
                resname = "GLYP"
                if atomname == "H":
                    atomname = "HT1"
                elif atomname == "H2":
                    atomname = "HT2"
                elif atomname == "H3":
                    atomname = "HT3"
            elif resname == "PRO" and atomname in [
                "N",
                "HN1",
                "HN2",
                "CD",
                "CA",
                "HD1",
                "HD2",
                "HA",
                "H2",
                "H3",
            ]:
                resname = "PROP"
                if atomname == "H2":
                    atomname = "HN1"
                elif atomname == "H3":
                    atomname = "HN2"
            elif resname == "ACE":
                if atomname == "CH3":
                    atomname = "CAY"
                elif atomname == "HH31":
                    atomname = "HY1"
                elif atomname == "HH32":
                    atomname = "HY2"
                elif atomname == "HH33":
                    atomname = "HY3"
                elif atomname == "C":
                    atomname = "CY"
                elif atomname == "O":
                    atomname = "OY"
            else:
                if atomname in ["N", "H", "H2", "H3", "CA", "HA"]:
                    resname = "NTER"
                    if atomname == "H":
                        atomname = "HT1"
                    elif atomname == "H2":
                        atomname = "HT2"
                    elif atomname == "H3":
                        atomname = "HT3"
        elif residue.is_c_term:
            if atomname in ["O", "OXT", "C"]:
                resname = "CTER"
                if atomname == "O":
                    atomname = "OT1"
                elif atomname == "OXT":
                    atomname = "OT2"
        elif residue.type == 3:
            resname = "TP3M"
            if atomname == "O":
                atomname = "OH2"
        # Residue substitutions
        if resname == "ILE":
            if atomname == "CD1":
                atomname = "CD"
            elif atomname == "HD11":
                atomname = "HD1"
            elif atomname == "HD12":
                atomname = "HD2"
            elif atomname == "HD13":
                atomname = "HD3"
            elif atomname == "HG12":
                atomname = "HG11"
            elif atomname == "HG13":
                atomname = "HG12"
        elif resname == "CYS" and "HG" not in residue.map:
            resname = "CYS"
            if atomname == "CB":
                resname = "DISU"
                atomname = "1CB"
            elif atomname == "SG":
                resname = "DISU"
                atomname = "1SG"
        elif resname == "HIS":
            if "HD1" in residue.map and "HE2" in residue.map:
                resname = "HSP"
            elif "HD1" in residue.map:
                resname = "HSD"
            elif "HE2" in residue.map:
                resname = "HSE"
        elif resname in ["GLU", "GLH"]:
            if "HE1" in residue.map:
                if atomname == "HE1":
                    atomname = "HE2"
                elif atomname == "OE1":
                    atomname = "OE2"
                elif atomname == "OE2":
                    atomname = "OE1"
                if atomname in [
                    "CG",
                    "HG3",
                    "HG1",
                    "HG2",
                    "CD",
                    "OE1",
                    "OE2",
                    "HE2",
                ]:
                    resname = "GLUP"
                else:
                    resname = "GLU"
            elif "HE2" in residue.map:
                if atomname in [
                    "CG",
                    "HG3",
                    "HG1",
                    "HG2",
                    "CD",
                    "OE1",
                    "OE2",
                    "HE2",
                ]:
                    resname = "GLUP"
                else:
                    resname = "GLU"
        elif resname in ["ASP", "ASH"]:
            if "HD1" in residue.map:
                if atomname == "HD1":
                    atomname = "HD2"
                elif atomname == "OD1":
                    atomname = "OD2"
                elif atomname == "OD2":
                    atomname = "OD1"
                if atomname in [
                    "CB",
                    "HB3",
                    "HB1",
                    "HB2",
                    "CG",
                    "OD1",
                    "OD2",
                    "HD2",
                ]:
                    resname = "ASPP"
                else:
                    resname = "ASP"
            elif "HD2" in residue.map:
                if atomname in [
                    "CB",
                    "HB3",
                    "HB1",
                    "HB2",
                    "CG",
                    "OD1",
                    "OD2",
                    "HD2",
                ]:
                    resname = "ASPP"
                else:
                    resname = "ASP"
        # HETATM Substitutions
        if resname == "ACE":
            if atomname == "CH3":
                atomname = "CAY"
            elif atomname == "HH31":
                atomname = "HY1"
            elif atomname == "HH32":
                atomname = "HY2"
            elif atomname == "HH33":
                atomname = "HY3"
            elif atomname == "C":
                atomname = "CY"
            elif atomname == "O":
                atomname = "OY"
        elif resname == "ADP":
            atomname = atomname.replace("*", "'")
        elif resname == "NME":
            resname = "CT3"
            if atomname == "HH31":
                atomname = "HT1"
            elif atomname == "HH32":
                atomname = "HT2"
            elif atomname == "HH33":
                atomname = "HT3"
            elif atomname == "CH3":
                atomname = "CAT"
            elif atomname == "N":
                atomname = "NT"
            elif atomname == "H":
                atomname = "HNT"
        # Hydrogen Substitutions
        if atomname == "H":
            atomname = "HN"
        elif atomname == "HA2":
            atomname = "HA1"
        elif atomname == "HA3":
            atomname = "HA2"
        elif atomname == "HB2" and resname not in ["ALA"]:
            atomname = "HB1"
        elif atomname == "HB3" and resname not in ["ALA"]:
            atomname = "HB2"
        elif atomname == "HD2" and resname not in [
            "HSP",
            "HSE",
            "HSD",
            "ASPP",
        ]:
            atomname = "HD1"
        elif atomname == "HD3" and resname not in ["HIS", "HSE", "HSD"]:
            atomname = "HD2"
        elif atomname == "HE2" and resname not in [
            "TRP",
            "HSP",
            "HSE",
            "HSD",
            "GLUP",
        ]:
            atomname = "HE1"
        elif atomname == "HE3" and resname not in ["TRP", "HSP", "HSE", "HSD"]:
            atomname = "HE2"
        elif atomname == "HG2":
            atomname = "HG1"
        elif atomname == "HG3":
            atomname = "HG2"
        elif atomname == "HG" and resname in ["SER", "CYS"]:
            atomname = "HG1"
        return resname, atomname


class ForcefieldResidue:
    """ForcefieldResidue class

    The ForceFieldResidue class contains a mapping of all atoms within the
    residue for easy searching.
    """

    def __init__(self, name):
        """Initialize the ForceFieldResidue object.

        :param name:  the name of the residue
        :type name:  str
        """
        self.name = name
        self.atoms = {}

    def add_atom(self, atom):
        """Add an atom to the ForcefieldResidue object.

        :param atom:  the atom to be added
        :type atom:  Atom
        """
        atomname = atom.name
        self.atoms[atomname] = atom

    def has_atom(self, atomname):
        """Check to see if the named atom is in the current residue.

        :param atomname:  the name of the atom to search for
        :type atomname:  str
        :return:  indication of whether atmo is present
        :rtype:  bool
        """
        return atomname in self.atoms

    def get_atom(self, atomname):
        """Return the atom object with the given atomname.

        :param resname:  the name of the atom
        :type resname:  str
        :return: the atom object
        :rtype:  ForcefieldAtom
        """
        if self.has_atom(atomname):
            return self.atoms[atomname]
        return None


class ForcefieldAtom:
    """ForcefieldAtom class.

    Contains fields that are related to the forcefield at the atom level.
    """

    def __init__(self, name, charge, radius, resname, group=""):
        """Initialize the object.

        :param name:  atom name
        :type name:  str
        :param charge:  the charge on the atom
        :type charge:  float
        :param radius:  the radius of the atom
        :type radius:  float
        :param resname:  the residue name
        :type resname:  str
        :param group:  the group name
        :type group:  str
        """
        self.name = name
        self.charge = charge
        self.radius = radius
        self.resname = resname
        self.group = group

    def get(self, name):
        """Get a member of the ForcefieldAtom class.

        :param name:  the name of the member, including:

            * name:    The atom name (returns string)
            * charge:  The charge on the atom (returns float)
            * radius:  The radius of the atom (returns float)
            * epsilon: The epsilon assocaited with the atom (returns float)

        :type name:  str
        :return:  the value of the member
        :raises KeyError:  if member does not exist
        """
        try:
            return getattr(self, name)
        except AttributeError as e:
            message = (
                f'Unable to access object "{name}" in class ForcefieldAtom'
            )
            raise KeyError(message) from e

    def __str__(self):
        txt = f"{self.name}:\n"
        txt += f"  Charge: {self.charge:.4f}\n"
        txt += f"  Radius: {self.radius:.4f}"
        return txt
