"""Routines for PDB2PQR

This module contains the protein object used in PDB2PQR and associated methods

Authors:  Todd Dolinsky, Yong Huang
"""
# TODO - Remove import * from this module
# from .aa import Amino, WAT
from .aa import *
from .na import Nucleic
from .structures import Chain, Residue
from .pdb import TER, ATOM, HETATM, END, MODEL
from .forcefield import Forcefield


class Protein:
    """Protein class
    
    The protein class represents the parsed PDB, and provides a hierarchy of
    information - each Protein contains a list of Chain objects as provided in
    the PDB file.  Each Chain then contains its associated list of Residue
    objects, and each Residue contains a list of Atom objects, completing the
    hierarchy.
    """

    def __init__(self, pdblist, definition):
        """Initialize using parsed PDB file
        
        Args:
            pdblist: List of Classes of PDB lines as created
        """
        self.chainmap = {}
        self.chains = []
        self.residues = []
        self.referencemap = definition.map
        self.patchmap = definition.patches

        chainDict = {}
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

                if previousAtom == None:
                    previousAtom = record
                
                if chainID not in chainDict:
                    myChain = Chain(chainID)
                    chainDict[chainID] = myChain
                        
                if resSeq != previousAtom.resSeq or \
                      iCode != previousAtom.iCode or \
                      chainID != previousAtom.chainID:
                    myResidue = self.createResidue(residue, previousAtom.resName)
                    chainDict[previousAtom.chainID].addResidue(myResidue)
                    residue = []

                residue.append(record)
                previousAtom = record

            elif isinstance(record, END):
                myResidue = self.createResidue(residue, previousAtom.resName)
                chainDict[previousAtom.chainID].addResidue(myResidue)
                residue = []

            elif isinstance(record, MODEL):
                numModels += 1
                if residue == []: continue
                if numModels > 1:
                    myResidue = self.createResidue(residue, previousAtom.resName)    
                    chainDict[previousAtom.chainID].addResidue(myResidue)
                    break

            elif isinstance(record, TER):
                count += 1

        if residue != [] and numModels <= 1:
            myResidue = self.createResidue(residue, previousAtom.resName)
            chainDict[previousAtom.chainID].addResidue(myResidue)

        # Keep a map for accessing chains via chainID

        self.chainmap = chainDict.copy()

        # Make a list for sequential ordering of chains
        
        if "" in chainDict:
            chainDict["ZZ"] = chainDict[""]
            del chainDict[""]

        keys = list(chainDict.keys())
        keys.sort()

        for key in keys:
            self.chains.append(chainDict[key])

        for chain in self.chains:
            for residue in chain.getResidues():
                self.residues.append(residue)

    def createResidue(self, residue, resname):
        """Create a residue object.
        
        If the resname is a known residue type, try to make that specific
        object, otherwise just make a standard residue object.

        Args:
            residue:  A list of atoms (list)
            resname:  The name of the residue (string)

        Returns:
            residue:  The residue object (Residue)
        """
        try:
            refobj = self.referencemap[resname]
            if refobj.name != resname: #Patched!
                obj = "%s(residue, refobj)" % refobj.name
                residue = eval(obj)
                residue.reference = refobj
            else:
                obj = "%s(residue, refobj)" % resname
                residue = eval(obj)
        except (KeyError, NameError):
            residue = Residue(residue)
        return residue

    def printAtoms(self, atomlist, chainflag=False, pdbfile=False):
        """Get the text for the entire protein
        
        Args:
            atomlist:  The list of atoms to include (list)
            chainflag:  Flag whether to print chainid or not
        Returns:
            text:  list of (stringed) atoms (list)
        """
        self.reSerialize()
        text = []
        currentchainID = None
        for atom in atomlist:
            # Print the "TER" records between chains
            if currentchainID == None:
                currentchainID = atom.chainID
            elif atom.chainID != currentchainID:
                currentchainID = atom.chainID
                text.append("TER\n")
            
            if pdbfile == True:
                text.append("%s\n" % atom.getPDBString())
            else:
                text.append("%s\n" % atom.getPQRString(chainflag=chainflag))
        text.append("TER\nEND")
        return text

    def createHTMLTypeMap(self, definition, outfilename):
        """Create an HTML typemap file at the desired location.
        
        If a type cannot be found for an atom a blank is listed.
        
        Args:
            definition: The definition objects.
            outfilename:  The name of the file to write (string)
        """
        # Cache the initial atom numbers
        numcache = {}
        for atom in self.getAtoms():
            numcache[atom] = atom.serial
        self.reSerialize()

        amberff = Forcefield("amber", definition, None)
        charmmff = Forcefield("charmm", definition, None)

        file = open(outfilename, "w")
        file.write("<HTML>\n")
        file.write("<HEAD>\n")
        file.write("<TITLE>PQR Typemap (beta)</TITLE>\n")
        file.write("</HEAD>\n")
        file.write("<BODY>\n")
        file.write("<H3>This is a developmental page including the atom type for the atoms in the PQR file.</H3><P>\n")
        file.write("<TABLE CELLSPACING=2 CELLPADDING=2 BORDER=1>\n")
        file.write("<tr><th>Atom Number</th><th>Atom Name</th><th>Residue Name</th><th>Chain ID</th><th>AMBER Atom Type</th><th>CHARMM Atom Type</th></tr>\n")
       
        for atom in self.getAtoms():
            if isinstance(atom.residue, (Amino, WAT, Nucleic)):
                resname = atom.residue.ffname
            else:
                resname = atom.residue.name

            ambergroup = amberff.getGroup(resname, atom.name)
            charmmgroup  = charmmff.getGroup(resname, atom.name)
        
            
            file.write("<tr><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>\n" % (atom.serial, atom.name, resname, atom.chainID, ambergroup, charmmgroup))
        

        file.write("</table>\n")
        file.write("</BODY></HTML>\n")
        file.close()

        # Return the original numbers back
        for atom in self.getAtoms():
            atom.serial = numcache[atom]

    def reSerialize(self):
        """Generate new serial numbers for atoms in the protein"""
        count = 1
        for atom in self.getAtoms():
            atom.set("serial", count)
            count += 1

    def getResidues(self):
        """Return the list of residues in the entire protein"""
        return self.residues
    
    def numResidues(self):
        """Get the number of residues for the entire protein (including
        multiple chains)

        Returns:
            count:  Number of residues in the protein (int)
        """
        return len(self.getResidues())

    def numAtoms(self):
        """Get the number of atoms for the entire protein (including multiple
        chains)
        """
        return len(self.getAtoms())

    def getAtoms(self):
        """Return all Atom objects in list format

        Returns:
            atomlist:  List of Atom objects in the protein (list)
        """
        atomlist = []
        for chain in self.chains:
            for atom in chain.getAtoms():
                atomlist.append(atom)
        return atomlist

    def getCharge(self):
        """Get the total charge on the protein
        
        NOTE:  Since the misslist is used to identify incorrect charge
        assignments, this routine does not list the 3 and 5 termini of nucleic
        acid chains as having non-integer charge even though they are
        (correctly) non-integer.
 
        Returns:
            misslist:  List of residues with non-integer charges (list)
            charge:  The total charge on the protein (float)
        """
        charge = 0.0
        misslist = []
        for chain in self.chains:
            for residue in chain.get("residues"):
                rescharge = residue.getCharge()
                charge += rescharge
                if isinstance(residue, Nucleic):               
                    if residue.is3term or residue.is5term: continue
                if float("%i" % rescharge) != rescharge:
                    misslist.append(residue)
        return misslist, charge

    def getChains(self):
        """Get the chains object

        Returns
            chains:  The list of chains in the protein (chain)
        """
        return self.chains
    
    def getSummary(self):
        """Some sort of undefined text output."""
        output = []
        for chain in self.chains:
            output.append(chain.getSummary())
        return ' '.join(output)
