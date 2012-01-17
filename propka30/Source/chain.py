#
# * This library is free software; you can redistribute it and/or
# * modify it under the terms of the GNU Lesser General Public
# * License as published by the Free Software Foundation; either
# * version 2.1 of the License, or (at your option) any later version.
# *
# * This library is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# * Lesser General Public License for more details.
#

#propka3.0, revision 182                                                                      2011-08-09
#-------------------------------------------------------------------------------------------------------
#--                                                                                                   --
#--                                   PROPKA: A PROTEIN PKA PREDICTOR                                 --
#--                                                                                                   --
#--                              VERSION 3.0,  01/01/2011, COPENHAGEN                                 --
#--                              BY MATS H.M. OLSSON AND CHRESTEN R. SONDERGARD                       --
#--                                                                                                   --
#-------------------------------------------------------------------------------------------------------
#
#
#-------------------------------------------------------------------------------------------------------
# References:
#
#   Very Fast Empirical Prediction and Rationalization of Protein pKa Values
#   Hui Li, Andrew D. Robertson and Jan H. Jensen
#   PROTEINS: Structure, Function, and Bioinformatics 61:704-721 (2005)
#
#   Very Fast Prediction and Rationalization of pKa Values for Protein-Ligand Complexes
#   Delphine C. Bas, David M. Rogers and Jan H. Jensen
#   PROTEINS: Structure, Function, and Bioinformatics 73:765-783 (2008)
#
#   PROPKA3: Consistent Treatment of Internal and Surface Residues in Empirical pKa predictions
#   Mats H.M. Olsson, Chresten R. Sondergard, Michal Rostkowski, and Jan H. Jensen
#   Journal of Chemical Theory and Computation, 7, 525-537 (2011)
#-------------------------------------------------------------------------------------------------------
import math, sys
import lib
import mutate
from residue import Residue


class Chain:
    """
        Chain class - contains a chain and properties of chain
    """

    def __init__(self, atoms, resInfo=None, options=None):
        """
        Constructer of chain object
        """
        self.chainID        = atoms[ atoms["keys"][0] ][0].chainID
        self.alignment      = None
        self.residues       = []
        self.configurations = []
        self.last_residue   = None

        if options.verbose == True:
          print("constructing chain %c (atoms: configurations)" % (self.chainID))

        # creating the 'residues'
        for key in atoms["keys"]:
          myResidue = Residue(atoms[key], resInfo=resInfo)
          if myResidue.checked():
            self.last_residue = myResidue
            self.residues.append(myResidue)

        # setting up list of configuration keys in this chain
        for residue in self.residues:
          for key in residue.configurations:
            if key not in self.configurations:
              self.configurations.append(key)

        # checking and correcting C-terminus and N-terminus
        self.addNTerminus(resInfo=resInfo)
        self.addCTerminus(resInfo=resInfo)



    def addNTerminus(self, resInfo=None):
        """
        Creating a N-terminus residue for this chain and adding it to 'residues' list
        """
        atom = None
        for residue in self.residues:
          if residue.type == "amino-acid":
            atom = residue.getAtom(name="N")
            break
        if atom != None:
          Nterminus = Residue([atom], resName="N+ ", resInfo=resInfo)
          self.residues.append(Nterminus)


    def addCTerminus(self, resInfo=None):
        """
        Creating a C-terminus residue for this chain and adding it to 'residues' list
        """
        last_residue = None
        for residue in self.residues:
          if residue.type == "amino-acid":
            last_residue = residue
        if last_residue != None:
          atoms = last_residue.checkOXT()
          Cterminus = Residue(atoms, resName="C- ", resInfo=resInfo)
          self.residues.append(Cterminus)


    def fillUnknownConfigurations(self, keys=None, options=None):
        """
        Fills in  the configurations that have not been read
        """
        for residue in self.residues:
          residue.fillUnknownConfigurations(keys=keys, options=options)


    def checkResidues(self, options=None):
        """
        Checks that there are only known residues - I ignore ligands for now
        """
        resNumb = 0
        for residue in self.residues:
            resNumb += 1
            while resNumb < residue.resNumb:
               print("XXX%4d - missing residue" % (resNumb))
               resNumb += 1
            residue.checkResidue(options=options)


    def calculateDesolvation(self, atoms, version=None, options=None):
        """
        Calculates the desolvation for each residue in this chain. Note, 'atoms' contains ALL atoms, not only 
        atoms belonging to this chain. Thus, you will get the desolvation of this chain in the presence of all.
        """
        propka1_list = lib.residueList("propka1")
        for residue in self.residues:
          if   residue.location == "BONDED":
            do_it = False
          elif residue.type     == "ion":
            do_it = True
          elif residue.resName in propka1_list:
            do_it = True
          else:
            do_it = False

          if do_it == True:
              residue.calculateDesolvation(atoms, version=version, options=options)


    def appendToBackBoneLists(self, NHlist, COlist):
        """
        Creates and returns CO-list list for base-backbone interactions
        """
        for residue in self.residues:
            if residue == self.last_residue:
                """ do nothing """
            elif residue.type == "amino-acid":
                N, H, C, O = residue.extractBackBoneAtoms()
                if N == None or H == None:
                  """ do nothing, probably PRO; print("%s missing N-H" % (residue.label)) """
                else:
                  NHlist.append([N, H])
                if C == None or O == None:
                  """ do nothing, no idear what this is; print("%s missing N-H" % (residue.label)) """
                else:
                  COlist.append([C, O])
            else:
                """ do nothing """


    def appendToResidueDictionary(self, residue_dictionary):
        """
        Adds all propka interesting residues for this chain to list.
        """
        # list of residue names that should be included in the dictionary

        for residue in self.residues:
          if   residue.type in ["amino-acid", "N-terminus", "C-terminus"]:
            key = residue.resName
          elif residue.type == "ligand":
            key = "LIG"
          elif residue.type == "ion":
            key = "ION"
          else:
            print("don't know what I have here %s (%s)" % (residue.type, residue.resName))

          if residue.location != "BONDED":
            if key in residue_dictionary:
              residue_dictionary[key].append(residue)
            else:
              residue_dictionary[key] = [residue]


    def appendPropkaResidues(self, propka_residues):
        """
        Adds all propka interesting residues for this chain to list.
        """
        residue_interaction_list = lib.residueInteractionList("ALL")
        #ligand_interaction_list = version.ions.keys()
        ligand_interaction_list = []
        for residue in self.residues:
            if residue.resName in residue_interaction_list or residue.resName in ligand_interaction_list:
                if residue.location == "BONDED":
                    """ do nothing """
                    #print("%s %s" % (residue.label, residue.type))
                else:
                    propka_residues.append(residue)
                    

    def appendAcidicResidues(self, acidic_residues, pka_residues):
        """
        Creates and returns a list for acid-backbone interactions
        """
        acidic_residue_list = lib.residueList("acids")
        for residue in self.residues:
            if residue.resName in acidic_residue_list:
                if residue.location == "BONDED":
                    """ do nothing """
                else:
                    acidic_residues.append(residue)
                    pka_residues.append(residue)


    def appendBasicResidues(self, basic_residues, pka_residues):
        """
        Creates and returns a list for base-backbone interactions
        """
        basic_residue_list = lib.residueList("bases")
        for residue in self.residues:
            if residue.resName in basic_residue_list:
                basic_residues.append(residue)
                pka_residues.append(residue)


    def calculateTotalPKA(self):
        """
        Calculates the total pKa from pKa_mod and the determinants
        """
        residue_list = lib.residueList("propka1")
        for residue in self.residues:
            if residue.resName in residue_list:
                residue.calculateTotalPKA()


    def setConfiguration(self, key=None):
        """
        set the 'current possition' to a 'configuration'
        """
        #print( "switching to configuration %s (chain)" % (key) )
        for residue in self.residues:
            residue.setConfiguration(key=key)


    def cleanupResidues(self):
        """
        Initializing/cleaning residues from 'old' determinants etc.
        """
        for residue in self.residues:
            residue.cleanupPKA()


    def setChain(self, chainID):
        """
        Sets the chainID to a specific label, 'chainID'
        """
        self.chainID = chainID
        for residue in self.residues:
            residue.setChain(chainID)


    def makeSequence(self, mutation=None):
        """
        creates a sequence from the chain object to be used in Scwrl-mutations
        """
        sequence = ""
        for residue in self.residues:
          if residue.resName in lib.residueList("standard"):
            if   mutation == None:
              code, resName = lib.convertResidueCode(resName=residue.resName)
              sequence += code.lower()
            elif residue.label in mutation:
              new_label = mutation[residue.label]['label']
              code, resName = lib.convertResidueCode(resName=new_label[:3])
              sequence += code.upper()
            else:
              code, resName = lib.convertResidueCode(resName=residue.resName)
              sequence += code.lower()

        return sequence


    def setAlignment(self, alignment=None):
        """
        Sets the alignment information (from self)
        """
        if alignment == None:
          print("setting alignment according to protein")
          self.alignment = ""
          for residue in self.residues:
            if residue.resName != "N+ " and residue.resName != "C- ":
              code, resName = lib.convertResidueCode(resName=residue.resName)
              if code in "ARNDCQEGHILKMFPSTWYV":
                self.alignment += (code)
          print(self.alignment)
          print( len(self.alignment) )
        else:
          self.alignment = alignment


    def shiftResidueNumber(self, shift):
        """
        Shift the residue numbers with 'shift'
        """
        for residue in self.residues:
            residue.shiftResidueNumber(shift)


    def printDeterminants(self):
        """
        prints the resulting pKa values and determinants to stdout
        """
        # Get same order as propka2.0
        residue_list = lib.residueList("propka1")
        for residue_type in residue_list:
          for residue in self.residues:
            if residue.resName == residue_type:
                print("%s" % ( residue.getDeterminantString() ))


    def writeDeterminants(self, file, verbose=True):
        """
        prints the resulting pKa values and determinants to stdout
        """
        # Get same order as propka2.0
        residue_list = lib.residueList("propka1")
        for residue_type in residue_list:
          for residue in self.residues:
            if residue.resName == residue_type:
                str = "%s\n" % ( residue.getDeterminantString() )
                file.write(str)
                if verbose == True:
                  print(str)


    def writeSummary(self, file, verbose=True):
        """
        prints the resulting pKa values in summary form to stdout
        """

        residue_list = lib.residueList("propka1")
        for residue_type in residue_list:
          for residue in self.residues:
            if residue.resName == residue_type:
                str = "%s\n" % ( residue.getSummaryString() )
                file.write(str)
                if verbose == True:
                  print(str)


    def printSummary(self):
        """
        prints the resulting pKa values in summary form to stdout
        """

        residue_list = lib.residueList("propka1")
        for residue_type in residue_list:
          for residue in self.residues:
            if residue.resName == residue_type:
                str = residue.getSummaryString()
                print(str)


    def calculateCharge(self, pH):
        """
        Calculates the total charge of this chain at pH 'pH'
        """
        propka1_residue_labels = lib.residueList("propka1")
        Qpro = 0.00
        Qmod = 0.00
        for residue in self.residues:
          if residue.resName in propka1_residue_labels:
            Qpro += residue.getCharge(pH, "folded")
            Qmod += residue.getCharge(pH, "unfolded")

        return Qpro, Qmod


    def calculateFoldingEnergy(self, pH=None, reference=None, options=None):
        """
        Calculates the folding energy given the correct pKa values; given
        pKa values calculated for the entire protein 'all chains' will give 
        the total folding energy partitioned to chains, not chain folding 
        energies.
        """
        propka1_residue_labels = lib.residueList("propka1")
        dG = 0.00
        for residue in self.residues:
          if residue.resName in propka1_residue_labels:
            ddG = residue.calculateFoldingEnergy(pH=pH, reference=reference, options=options)
            dG += ddG
            #print "%s %6.2lf" % (residue.label, ddG)

        return dG


