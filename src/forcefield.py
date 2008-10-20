"""
    Forcefield.py

    This module takes a pdblist as input and replaces the occupancy and
    tempfactor fields with charge and radius fields, with values as defined
    by a particular forcefield.  The forcefield structure is modeled off of
    the structures.py file, where each forcefield is considered a chain of
    residues of atoms.

    ----------------------------
   
    PDB2PQR -- An automated pipeline for the setup, execution, and analysis of
    Poisson-Boltzmann electrostatics calculations

	Copyright (c) 2002-2007, Jens Erik Nielsen, University College Dublin; 
	Nathan A. Baker, Washington University in St. Louis; Paul Czodrowski & 
	Gerhard Klebe, University of Marburg

	All rights reserved.

	Redistribution and use in source and binary forms, with or without modification, 
	are permitted provided that the following conditions are met:

		* Redistributions of source code must retain the above copyright notice, 
		  this list of conditions and the following disclaimer.
		* Redistributions in binary form must reproduce the above copyright notice, 
		  this list of conditions and the following disclaimer in the documentation 
		  and/or other materials provided with the distribution.
		* Neither the names of University College Dublin, Washington University in 
		  St. Louis, or University of Marburg nor the names of its contributors may 
		  be used to endorse or promote products derived from this software without 
		  specific prior written permission.

	THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND 
	ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED 
	WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. 
	IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, 
	INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, 
	BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
	DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
	LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE 
	OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED 
	OF THE POSSIBILITY OF SUCH DAMAGE.

    ----------------------------
"""

__date__ = "6 November 2007"
__author__ = "Todd Dolinsky, Yong Huang"

import string
import sys
import getopt
import os
import re

from xml import sax
from utilities import *

class ForcefieldHandler(sax.ContentHandler):
   
    def __init__(self, map, reference):
        self.oldresname = None
        self.oldatomname = None
        self.curelement = None
        self.newatomname = None
        self.newresname = None
        self.atommap = {}
        self.map = map
        self.reference = reference           

    def updateMap(self, toname, fromname, map):
        """
            Update the given map by adding a pointer from a new
            name to an object.

            Parameters
                toname:  The new name for the object (string)
                fromname:  The old name for the object (string)
                map:  A dictionary of items (dict)
        """
        fromobj = map[fromname]
        if isinstance(fromobj, ForcefieldResidue):
            if toname not in map:
                newres = ForcefieldResidue(fromname)
                map[toname] = newres
            for atomname in fromobj.atoms:
                map[toname].atoms[atomname] = fromobj.atoms[atomname]
        elif isinstance(fromobj, ForcefieldAtom):
            map[toname] = fromobj
            
                        
    def findMatchingNames(self, regname, map):
        """
            Find a list of strings that match the given regular
            expression.

            Parameters
                regname: The regular expression (string)
                map:  The dictionary to search (dict)
                
            Returns
                list:  A list of regular expression objects that match
                       the regular expression.
        """
        list = [] 
        regname += "$"
     
        # Find the existing items that match this string

        for name in map:
            regexp = re.compile(regname).match(name)
            if regexp:
                list.append(regexp)

        return list
   
    def startElement(self, name, attributes):
        """
            Override the startElement function to keep track of the current
            element.
        """
        if name != "name": self.curelement = name
            
    def endElement(self, name):
        """
            At the end of the element, act on the stored information.

            Parameters
                name:  The name of the element (string)
        """
        if name == "residue":
            if self.oldresname != None:  # Make a new residue hook
              
                newreslist = self.findMatchingNames(self.newresname, self.reference)
                if self.oldresname.find("$group") >= 0:  # Multiple new residues
                    for resitem in newreslist:
                        resname = resitem.string
                        group = resitem.group(1)
                        fromname = self.oldresname.replace("$group", group)
                        if fromname in self.map:
                            self.updateMap(resname, fromname, self.map)
                              
                else: # Work with a single new residue name
                    oldreslist = self.findMatchingNames(self.oldresname, self.map)
                    for resitem in newreslist:
                        resname = resitem.string
                        self.updateMap(resname, self.oldresname, self.map)

            # If this was only a residue conversion, exit
                
            if self.atommap == {}:
                self.oldresname = None
                self.newresname = None
                return

            # Apply atom conversions for all appropriate residues
           
            resmatchlist = self.findMatchingNames(self.newresname, self.map)
            for resitem in resmatchlist:
                residue = self.map[resitem.string]        
                for newname in self.atommap:        
                    oldname = self.atommap[newname]
                    if oldname not in residue.atoms: continue
                    self.updateMap(newname, oldname, residue.atoms)
                    
            # Clean up

            self.oldresname = None
            self.newresname = None
            self.atommap = {}
            
        elif name == "atom": 

            self.atommap[self.newatomname] = self.oldatomname
            self.oldatomname = None
            self.newatomname = None

        else: # Just free the current element namespace
            self.curelement = ""

        return self.map

    def characters(self, text):
        """
            Store the information in the object for future use/

            Parameters
                text:  The text value between the XML tags
        """
        if text.isspace(): return
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
    """
        Forcefield class

        The forcefield class contains definitions for a given forcefield.
        Each forcefield object contains a dictionary of residues, with each
        residue containing a dictionary of atoms.  Dictionaries are used
        instead of lists as the ordering is not important. The forcefield
        definition files are unedited, directly from the forcefield - all
        transformations are done within.

    """

    def __init__(self, ff, definition, userff, usernames = None):
        """
            Initialize the class by parsing the definition file

            Parameters
                ff: The name of the forcefield (string)
                definition: The definition objects
                userff:  A link to the file for CGI based user-defined
                         forcefields
        """
        self.map = {}
        self.name = ff
        defpath = ""

        if userff == None:
            defpath = getFFfile(ff)
            if defpath == "":
                raise ValueError, "Unable to find forcefield parameter file %s!" % path
          
            file = open(defpath, 'rU')

        else: file = userff

        lines = file.readlines()
        for line in lines:
            if not line.startswith("#"):
                fields = string.split(line)
                if fields == []: continue  
                try:
                    resname = fields[0]
                    atomname = fields[1]
                    charge = float(fields[2])
                    radius = float(fields[3])
                except ValueError:
                    txt = "Unable to recognize user-defined forcefield file" 
                    if defpath != "": txt += " %s!" % defpath
                    else: txt += "!"
                    txt += " Please use a valid parameter file."
                    raise ValueError, txt
            
                try:
                    group = fields[4]
                    atom = ForcefieldAtom(atomname, charge, radius, resname, group)
                except:
                    atom = ForcefieldAtom(atomname, charge, radius, resname)

                myResidue = self.getResidue(resname)
                if myResidue == None:
                    myResidue = ForcefieldResidue(resname)
                    self.map[resname] = myResidue
                myResidue.addAtom(atom)

        file.close()

        # Now parse the XML file, associating with FF objects -
        # This is not necessary (if canonical names match ff names)
 

        defpath = getNamesFile(ff)
        if defpath != "":
        
            handler = ForcefieldHandler(self.map, definition.map)
            sax.make_parser()

            if usernames != None:
                namesfile = usernames
                sax.parseString(namesfile.read(), handler)
            else:
                namesfile = open(defpath)
                sax.parseString(namesfile.read(), handler)
            namesfile.close()

        # CGI based .names file handling
        else: 
            handler = ForcefieldHandler(self.map, definition.map)
            sax.make_parser()

            if usernames != None:
                namesfile = usernames            
                sax.parseString(namesfile.getvalue(), handler)
            else:
                raise ValueError, "Please provide a valid .names file!"
            namesfile.close()


    def hasResidue(self, resname):
        """
            Check if the residue name is in the map or not.

            Parameters
                resname:  The name to search for (string)

            Returns
                1 if the resname is in the map, 0 otherwise.
        """
        if resname in self.map: return 1
        else: return 0

    def getResidue(self, resname):
        """
            Return the residue object with the given resname

            Parameters
                resname: The name of the residue (string)
            Returns
                residue: The residue object (ForcefieldResidue)
        """
        if self.hasResidue(resname): return self.map[resname]
        else: return None

   
    def getNames(self, resname, atomname):
        """
            Get the actual names associated with the input fields.
            The names passed in point to ForcefieldResidue and
            ForcefieldAtom objects which may have different names;
            grab these names and return.

            Parameters
                resname:  The residue name (string)
                atomname: The atom name (string)
            Returns
                rname:    The forcefield's name for this residue (string)
                aname:    The forcefield's name for this atom (string)
        """
        rname = None
        aname = None
        if resname in self.map:
            res = self.map[resname]
            if res.hasAtom(atomname):
                atom = res.atoms[atomname]
                aname = atom.name
                rname = atom.resname
        return rname, aname

    def getGroup(self, resname, atomname):
        """
            Get the group/type associated with the input
            fields.  If not found, return a null string.
            
            Parameters:
                resname:  The residue name (string)
                atomname: The atom name (string)
        """
        group = ""
        if resname in self.map:
            resid = self.map[resname]
            if resid.hasAtom(atomname):
                atom = resid.atoms[atomname]
                group = atom.group
        return group

    def getParams(self, resname, atomname):
        """
            Get the parameters associated with the input fields.
            The residue itself is needed instead of simply its name
            because  the forcefield may use a different residue name
            than the standard amino acid name.

            Parameters
                resname:  The residue name (string)
                atomname: The atom name (string)
            Returns
                charge:   The charge on the atom (float)
                radius:   The radius of the atom (float)
        """
        charge = None
        radius = None

        #print self.map.keys()

        if resname in self.map:
            resid = self.map[resname]
            if resid.hasAtom(atomname):
                atom = resid.atoms[atomname]
                charge = atom.charge
                radius = atom.radius
                
        return charge, radius

class ForcefieldResidue:
    """
        ForcefieldResidue class

        The ForceFieldResidue class contains a mapping of all atoms within
        the residue for easy searching.
    """
    def __init__(self, name):
        """
            Initialize the ForceFieldResidue object

            Parameters
                name: The name of the residue (string)
        """
        self.name = name
        self.atoms = {}


    def addAtom(self, atom):
        """
            Add an atom to the ForcefieldResidue

            Parameters
                atom:  The atom to be added (atom)
        """
        atomname = atom.get("name")
        self.atoms[atomname] = atom

    def getAtoms(self):
        """
            Return the list of atoms in this residue.
        """
        return self.atoms

    def hasAtom(self, atomname):
        """
            Check to see if the atomname is in the current residue.

            Parameters
                atomname:  The name of the atom to search for
            Returns
                1 if the atom is present in the residue, 0 otherwise
        """
        if atomname in self.atoms: return 1
        else: return 0

    def getAtom(self, atomname):
        """
            Return the atom object with the given atomname

            Parameters
                resname: The name of the atom (string)
            Returns
                residue: The atom object (ForcefieldAtom)
        """
        if self.hasAtom(atomname): return self.atoms[atomname]
        else: return None

class ForcefieldAtom:
    """
        ForcefieldAtom class

        The ForcefieldAtom object contains fields that are related to the
        forcefield at the atom level
    """
    
    def __init__(self, name, charge, radius, resname, group=""):
        """
            Initialize the object

            Parameters
                name:    The atom name (string)
                charge:  The charge on the atom (float)
                radius:  The radius of the atom (float)
                resname: The residue name (string)
                group:   The group name (string)
        """
        self.name = name
        self.charge = charge
        self.radius = radius
        self.resname = resname
        self.group = group

    def get(self, name):
        """
            Get a member of the ForcefieldAtom class

            Parameters
                name:       The name of the member (string)
            Possible Values
                name:    The atom name (string)
                charge:  The charge on the atom (float)
                radius:  The radius of the atom (float)
                epsilon: The epsilon assocaited with the atom (float)
            Returns
                item:       The value of the member
        """
        try:
            item = getattr(self, name)
            return item
        except AttributeError:
            message = "Unable to access object \"%s\" in class ForcefieldAtom" % name
            raise ValueError, message

    def __str__(self):
        """
            String representation of the forcefield atom.
        """
        txt = "%s:\n"% self.name
        txt += "  Charge: %.4f\n" % self.charge
        txt += "  Radius: %.4f" % self.radius
        return txt
