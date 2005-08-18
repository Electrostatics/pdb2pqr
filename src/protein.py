"""
    Routines for PDB2PQR

    This module contains the protein object used in PDB2PQR and associated
    methods
    
    ----------------------------
   
    PDB2PQR -- An automated pipeline for the setup, execution, and analysis of
    Poisson-Boltzmann electrostatics calculations

    Nathan A. Baker (baker@biochem.wustl.edu)
    Todd Dolinsky (todd@ccb.wustl.edu)
    Dept. of Biochemistry and Molecular Biophysics
    Center for Computational Biology
    Washington University in St. Louis

    Jens Nielsen (Jens.Nielsen@ucd.ie)
    University College Dublin

    Additional contributing authors listed in documentation and supporting
    package licenses.

    Copyright (c) 2003-2005.  Washington University in St. Louis.  
    All Rights Reserved.

    This file is part of PDB2PQR.

    PDB2PQR is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 2 of the License, or
    (at your option) any later version.
 
    PDB2PQR is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
   
    You should have received a copy of the GNU General Public License
    along with PDB2PQR; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA

    ----------------------------

"""

__date__ = "14 November 2003"
__author__ = "Todd Dolinsky"

from pdb import *
from structures import *

class Protein:
    """
        Protein class

        The protein class represents the parsed PDB, and provides a
        hierarchy of information - each Protein contains a list of Chain
        objects as provided in the PDB file.  Each Chain then contains its
        associated list of Residue objects, and each Residue contains a list
        of Atom objects, completing the hierarchy.
    """

    def __init__(self, pdblist):
        """
            Initialize using parsed PDB file

            Parameters
                pdblist: List of Classes of PDB lines as created
                         by pdb.py->readPDB
        """

        self.chainmap, self.chains = self.createProtein(pdblist)

    def createProtein(self, pdblist):
        """
            Fill the Protein with chains, residues, and atoms

            Parameters
                pdblist: List of Classes of PDB lines as created
                         by pdb.py->readPDB (list)
            Returns
                dict:    Mapping of chain ID to chain object
                list:    List of chain objects sorted by chain ID (dict)
        """

        dict = {}
        list = []
        
        previousAtom = None
        residue = []
        numModels = 0
        numChains = 1
        count = 0

        for record in pdblist: # Find number of chains
            if isinstance(record, TER):
                numChains += 1

        for record in pdblist:
            if isinstance(record, ATOM) or isinstance(record, HETATM):

                if record.chainID == "" and numChains > 1 and record.resName not in ["WAT","HOH"]:
                    # Assign a chain ID
                    record.chainID = string.ascii_uppercase[count]
                
                chainID = record.chainID
                resSeq = record.resSeq
                resName = record.resName
                iCode = record.iCode

                if previousAtom is None:
                    previousAtom = record
                
                if chainID not in dict:
                    myChain = Chain(chainID)
                    dict[chainID] = myChain
                        
                if resSeq != previousAtom.resSeq or \
                       iCode != previousAtom.iCode:
                    myResidue = Residue(residue, previousAtom)
                    dict[previousAtom.chainID].addResidue(myResidue)
                    residue = []

                residue.append(record)
                previousAtom = record

            elif isinstance(record, END):
                myResidue = Residue(residue, previousAtom)
                dict[previousAtom.chainID].addResidue(myResidue)
                residue = []

            elif isinstance(record, MODEL):
                numModels += 1
                if residue == []: continue
                if numModels > 1:
                    myResidue = Residue(residue, previousAtom)
                    dict[previousAtom.chainID].addResidue(myResidue)
                    break

            elif isinstance(record, TER):
                count += 1

        if residue != [] and numModels <= 1:
            myResidue = Residue(residue, previousAtom)
            dict[previousAtom.chainID].addResidue(myResidue)

        chainmap = dict.copy()
        if dict.has_key(""):
            dict["ZZ"] = dict[""]
            del dict[""]

        keys = dict.keys()
        keys.sort()

        for key in keys:
            list.append(dict[key])
            
        return chainmap, list

    def printAtoms(self, atomlist):
        """
            Get the text for the entire protein
            Parameters
                atomlist: The list of atoms to include (list)
            Returns
                text:     The list of (stringed) atoms (list)
        """
        self.reSerialize()
        text = []
        for atom in atomlist:
            text.append("%s\n" % str(atom))
        return text

    def reSerialize(self):
        """
            Generate new serial numbers for atoms in the protein
        """
        count = 1
        for atom in self.getAtoms():
            atom.set("serial", count)
            count += 1
    
    def numResidues(self):
        """
            Get the number of residues for the entire protein (including
            multiple chains)

            Returns
                count:  Number of residues in the protein (int)
        """
        count = 0
        for chain in self.chains:
            count += chain.numResidues()
        return count

    def numAtoms(self):
        """
            Get the number of atoms for the entire protein(including
            multiple chains)

            Returns
                count:  Number of atoms in the protein (int)
        """
        count = len(self.getAtoms())
        return count

    def getAtoms(self):
        """
            Return all Atom objects in list format

            Returns
                atomlist:  List of Atom objects in the protein (list)
        """

        atomlist = []
        for chain in self.chains:
            for atom in chain.getAtoms():
                atomlist.append(atom)
        return atomlist

    def getCharge(self):
        """
            Get the total charge on the protein
            NOTE:  Since the misslist is used to identify incorrect
                   charge assignments, this routine does not list the
                   3 and 5 termini of nucleic acid chains as having
                   non-integer charge even though they are (correctly)
                   non-integer.
            Returns:
                misslist: List of residues with non-integer
                          charges (list)
                charge:   The total charge on the protein (float)
        """
        charge = 0.0
        misslist = []
        for chain in self.chains:
            for residue in chain.get("residues"):
                rescharge = residue.getCharge()
                charge = charge + rescharge
                if residue.get("is3term") or residue.get("is5term"):
                    continue
                if float("%i" % rescharge) != rescharge:
                    misslist.append(residue)
        return misslist, charge

    def getChains(self):
        """
            Get the chains object

            Returns
                chains: The list of chains in the protein (chain)
        """
        return self.chains
