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

                if record.chain_id == "" and numChains > 1 and record.res_name not in ["WAT","HOH"]:
                    # Assign a chain ID
                    record.chain_id = string.ascii_uppercase[count]

                chain_id = record.chain_id
                res_seq = record.res_seq
                res_name = record.res_name
                ins_code = record.ins_code

                if previousAtom == None:
                    previousAtom = record
                
                if chain_id not in chainDict:
                    myChain = Chain(chain_id)
                    chainDict[chain_id] = myChain
                        
                if res_seq != previousAtom.res_seq or \
                      ins_code != previousAtom.ins_code or \
                      chain_id != previousAtom.chain_id:
                    my_residue = self.createResidue(residue, previousAtom.res_name)
                    chainDict[previousAtom.chain_id].addResidue(my_residue)
                    residue = []

                residue.append(record)
                previousAtom = record

            elif isinstance(record, END):
                my_residue = self.createResidue(residue, previousAtom.res_name)
                chainDict[previousAtom.chain_id].addResidue(my_residue)
                residue = []

            elif isinstance(record, MODEL):
                numModels += 1
                if residue == []: continue
                if numModels > 1:
                    my_residue = self.createResidue(residue, previousAtom.res_name)    
                    chainDict[previousAtom.chain_id].addResidue(my_residue)
                    break

            elif isinstance(record, TER):
                count += 1

        if residue != [] and numModels <= 1:
            my_residue = self.createResidue(residue, previousAtom.res_name)
            chainDict[previousAtom.chain_id].addResidue(my_residue)

        # Keep a map for accessing chains via chain_id

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
            for residue in chain.get_residues():
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
        currentchain_id = None
        for atom in atomlist:
            # Print the "TER" records between chains
            if currentchain_id == None:
                currentchain_id = atom.chain_id
            elif atom.chain_id != currentchain_id:
                currentchain_id = atom.chain_id
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
        for atom in self.get_atoms():
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
       
        for atom in self.get_atoms():
            if isinstance(atom.residue, (Amino, WAT, Nucleic)):
                resname = atom.residue.ffname
            else:
                resname = atom.residue.name

            ambergroup = amberff.get_group(resname, atom.name)
            charmmgroup  = charmmff.get_group(resname, atom.name)
        
            
            file.write("<tr><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>\n" % (atom.serial, atom.name, resname, atom.chain_id, ambergroup, charmmgroup))
        

        file.write("</table>\n")
        file.write("</BODY></HTML>\n")
        file.close()

        # Return the original numbers back
        for atom in self.get_atoms():
            atom.serial = numcache[atom]

    def reSerialize(self):
        """Generate new serial numbers for atoms in the protein"""
        count = 1
        for atom in self.get_atoms():
            atom.set("serial", count)
            count += 1

    def get_residues(self):
        """Return the list of residues in the entire protein"""
        return self.residues
    
    def num_residues(self):
        """Get the number of residues for the entire protein (including
        multiple chains)

        Returns:
            count:  Number of residues in the protein (int)
        """
        return len(self.get_residues())

    def numAtoms(self):
        """Get the number of atoms for the entire protein (including multiple
        chains)
        """
        return len(self.get_atoms())

    def get_atoms(self):
        """Return all Atom objects in list format

        Returns:
            atomlist:  List of Atom objects in the protein (list)
        """
        atomlist = []
        for chain in self.chains:
            for atom in chain.get_atoms():
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
