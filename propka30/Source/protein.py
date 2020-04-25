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

# propka3.0, revision 182                                                                      2011-08-09
# -------------------------------------------------------------------------------------------------------
# --                                                                                                   --
# --                                   PROPKA: A PROTEIN PKA PREDICTOR                                 --
# --                                                                                                   --
# --                              VERSION 3.0,  01/01/2011, COPENHAGEN                                 --
# --                              BY MATS H.M. OLSSON AND CHRESTEN R. SONDERGARD                       --
# --                                                                                                   --
# -------------------------------------------------------------------------------------------------------
#
#
# -------------------------------------------------------------------------------------------------------
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
# -------------------------------------------------------------------------------------------------------
from . import pdb
from .chain import Chain
from . import calculator as calculate
from . import coupled_residues
from . import output
from . import determinants
import math
from sys import exit
import os
import time
import string
from . import lib
pka_print = lib.pka_print
#import debug


class Protein:
    """
        Protein class - contains chains and protein properties
    """

    def __init__(self, atoms=None, pdbfile=None, name=None, options=None):
        """
        constructer of the protein object.
        """
        self.name = None     # name of this protein
        self.chains = []       # array of chains
        self.atoms = atoms    # dictionary of atom objects, i.e. atoms[chainID][resLabel]
        self.property = {}       # a dictionary that will hold protein properties such as e.g. pI, folding-profile etc.
        self.configurations = None     # a list of configuration keys valid for this protein.
        self.propka_residues = []       # residues needed by propka, Have to be in resNumb order to break the pair-loop correctly
        self.NHlist = []
        self.COlist = []
        self.residue_dictionary = None     # dictionary with references to protein residues, for easy search and storage
        self.alignment = None
        self.pI = None
        self.coupled_residues = False
        self.status = {'protonated': False,
                       'done pka':   False}

        # seting protein name (in this case probably the pdbcode)
        self.setName(name=name, pdbfile=pdbfile)

        # reading pdbfile and generating 'atoms' dictionary
        self.setAtomsDictionary(atoms=atoms, pdbfile=pdbfile)

        # generating residue parameter information needed to make the protein
        resInfo = getResidueParameters()

        if options.verbose == True:
            pka_print("constructing protein \"%s\"" % (self.name))

        # creating chains from atom objects read from the pdbfile
        for chainID in sorted(self.atoms.keys()):

            # creating a new chain with the atoms dictionary
            if options.chains == None or chainID in options.chains:
                self.chains.append(Chain(self.atoms[chainID], resInfo=resInfo, options=options))
            else:
                # killing this chain, not strictly needed, but hey ...
                del self.atoms[chainID]

        # setting up sorted and available keys for this protein
        self.setupConfigurationKeys()

        # enforcing all atoms to have each residue configurations
        for chain in self.chains:
            chain.fillUnknownConfigurations(keys=self.configurations, options=options)
            chain.checkResidues(options=options)

        # printing available protein configurations
        if options.verbose == True:
            self.printAvailableConfigurations()

        # protonate protein
        self.protonate(scheme=options.protonation)

        # setting up residue reference lists for the protein
        self.setupReferenceLists()

        # checking and setting SS-bonds
        self.findSSbonds()

        # cleaning the 'bonded' residues from 'propka_list'
        self.cleanBondedResiduesFromList(self.propka_residues)

    def setName(self, name=None, pdbfile=None):
        """ 
        set protein name
        """
        if name != None:
            self.name = name
        elif pdbfile != None:
            self.name = lib.extractName(pdbfile)
        else:
            self.name = "1xxx"

        return

    def setAtomsDictionary(self, atoms=None, pdbfile=None):
        """ 
        set reference to the atoms dictionary
        atoms = dictionary with atom objects to make the protein from, i.e. atoms[chainID][resName]
        """
        if atoms != None:
            self.atoms = atoms
        elif pdbfile != None:
            self.atoms = pdb.readPDB(filename=pdbfile)
        else:
            pka_print("need either an atoms dictionary or pdbfile to create a protein")
            exit(9)

        return

    def setupConfigurationKeys(self, options=None):
        """ 
        sets protein configuration keys and makes sure there is an available key on each residue
        """
        available_keys = []
        # collecting all configuration keys in this protein
        for chain in self.chains:
            for key in chain.configurations:
                if key not in available_keys:
                    available_keys.append(key)

        # sort configuration keys
        self.configurations = lib.sortConfigurationKeys(available_keys)

        # enforcing all atoms to have each residue configuration, and setting default configuration
        for chain in self.chains:
            chain.fillUnknownConfigurations(keys=self.configurations, options=options)

        return

    def printAvailableConfigurations(self, options=None):
        """ 
        prints available keys for the atom configurations-dictionary
        """
        str = "configurations:"
        for key in self.configurations:
            str += "%6s" % (key)
        pka_print(str)

        return

    def setupReferenceLists(self):
        """ 
        setup lists with references to residues and atoms needed to calculate pKa values
        lists: residue_dictionary
               NHlist
               COlist
               propka_residues
        """
        # initializing empty dictionary
        self.residue_dictionary = {}
        for key in lib.residueList("propka1"):
            self.residue_dictionary[key] = []
        self.residue_dictionary['ION'] = []

        for chain in self.chains:
            # make list with N-H & C=O fragments
            chain.appendToBackBoneLists(self.NHlist, self.COlist)

            # setup the residue dictionary with references to all residues
            chain.appendToResidueDictionary(self.residue_dictionary)

            # residues with propka interactions,
            # IMPORTANT, this list is assumed to be ordered according to resNumb to loop ver pairs
            chain.appendPropkaResidues(self.propka_residues)

        return

    def findSSbonds(self):
        """
        Finds all SS-pairs and sets location to 'BONDED' if distance < 2.50
        """
        # create CYS-list if not on protein object
        if self.residue_dictionary == None:
            CYSlist = []
            for chain in self.chains:
                for residue in chain.residues:
                    if residue.resName == "CYS":
                        CYSlist.append(residue)
        else:
            CYSlist = self.residue_dictionary["CYS"]

        # Searching all CYS pairs
        for residue1 in CYSlist:
            S1 = residue1.getAtom("SG")
            for residue2 in CYSlist:
                if residue1 == residue2:
                    break
                else:
                    S2 = residue2.getAtom("SG")
                    distance = calculate.InterAtomDistance(S1, S2)
                    if distance < 2.50:
                        residue1.location = "BONDED"
                        residue2.location = "BONDED"
                        residue1.pKa_mod = 99.99
                        residue2.pKa_mod = 99.99

        # cleaning out 'bonded' CYS from dictionary
        self.cleanBondedResiduesFromList(self.residue_dictionary["CYS"])

        return

    def cleanBondedResiduesFromList(self, residue_list):
        """ 
        Cleans out "BONDED" residues, e.g. CYS from the CYS dictionary
        """
        killme = []
        index = 0
        for residue in residue_list:
            if residue.location == "BONDED":
                killme.append(index)
            index += 1
        killme.reverse()
        for kill in killme:
            del residue_list[kill]

        return

    def removeHydrogens(self, options=None):
        """ 
        removes all protons from the protein
        """
        for chain in self.chains:
            for residue in chain.residues:
                killme = []
                i = -1
                for atom in residue.atoms:
                    i += 1
                    if atom.element == 'H':
                        killme.insert(0, i)
                for i in killme:
                    del residue.atoms[i]

        self.status['protonated'] = False

        return

    def protonate(self, scheme=None):
        """ 
        protonates the protein according to given scheme
        """
        from .protonator import makeProtonator
        please = makeProtonator(scheme=scheme)
        self.removeHydrogens()

        please.protonate(protein=self)

        self.renumberAtoms()
        self.status['protonated'] = True

        return

    def getChain(self, chainID=None):
        """ 
        returns chain chainID=chainID
        """
        for chain in self.chains:
            if chain.chainID == chainID:
                return chain

    def calculatePKA(self, version=None, options=None):
        """ 
        Calculates the pKa values, average if there are more than one configuration
        """
        # create a default version if not provided
        if version == None:
            from . import version
            version = version.makeVersion(label=options.version_label)

        if len(self.configurations) == 1:
            # found only one configuration in pdbfile
            self.calculateConfigurationPKA(version=version, options=options)
        else:
            # found multiple configurations in pdbfile
            self.calculateAveragePKA(version=version, options=options)

        # Check for coupled residues
        coupled_residues.identify_coupled_residues(self, options=options)

        # signalling that the pKa values are available for this protein.
        self.status['done pka'] = True

        # printing our averaged total pKa values
        if options.verbose == True:
            output.printResult(self)

        return

    def calculateAveragePKA(self, version=None, options=None):
        """ 
        Calculates the pKa values of each configuration and averages them
        """
        # initializing dictionaries for making averages, key = 'residue.label'
        pkas = {}
        Nmass = {}
        Emass = {}
        Nlocl = {}
        Elocl = {}
        determinants = {}
        number_of_configurations = len(self.configurations)
        # setting up list for 'propka1 pka residues'
        pka_residues = []
        for resName in lib.residueList("propka1"):
            pka_residues.extend(self.residue_dictionary[resName])
        # initializing dictionaries for making averages, key = 'configuration label'
        for residue in pka_residues:
            key = residue.label
            pkas[key] = {}
            Nmass[key] = {}
            Emass[key] = {}
            Nlocl[key] = {}
            Elocl[key] = {}
            # initializing dictionary for making averages, key = 'determinant.label'
            determinants[key] = [{}, {}, {}]

        # perform pKa calulation on each configuration and add pKa properties to dictionaries
        for key in self.configurations:
            self.setConfiguration(key=key, options=options)
            self.calculateConfigurationPKA(version=version, options=options)
            # printing out the pKa section for this configuration
            if options.verbose == True:
                output.printPKASection(self)

            # transferring property to right dictionary
            for residue in pka_residues:
                residue_key = residue.label
                pkas[residue_key][key] = residue.pKa_pro
                Nmass[residue_key][key] = residue.Nmass
                Emass[residue_key][key] = residue.Emass
                Nlocl[residue_key][key] = residue.Nlocl
                Elocl[residue_key][key] = residue.Elocl
                for type in range(3):
                    for determinant in residue.determinants[type]:
                        if determinant.label in determinants[residue_key][type]:
                            determinants[residue_key][type][determinant.label] += determinant.value
                        else:
                            determinants[residue_key][type][determinant.label] = determinant.value

        # print each configuration-pKa in nice matrix
        residue_list = self.residue_dictionary["ASP"]
        residue_list.extend(self.residue_dictionary["GLU"])
        str = "%4s" % ("Res#")
        for residue in residue_list:
            str += "%6d" % (residue.resNumb)
        pka_print(str)
        index = 0
        for key in self.configurations:
            str = "%4s:" % (key)
            for residue in residue_list:
                reskey = residue.label
                str += "%6.2lf" % (pkas[reskey][key])
                #str += "%6.2lf" % (Emass[reskey][key])
                #str += "%6.2lf" % (Elocl[reskey][key])
            pka_print(str)

        # get the average pKa properties by dividing the sum with number of configurations, len(configuration keys)
        from determinants import Determinant
        for residue in pka_residues:
            residue_key = residue.label
            sum_pka = 0.00
            sum_Nmass = 0.00
            sum_Emass = 0.00
            sum_Nlocl = 0.00
            sum_Elocl = 0.00
            for key in pkas[residue_key].keys():
                sum_pka += pkas[residue_key][key]
                sum_Nmass += Nmass[residue_key][key]
                sum_Emass += Emass[residue_key][key]
                sum_Nlocl += Nlocl[residue_key][key]
                sum_Elocl += Elocl[residue_key][key]
            residue.pKa_pro = sum_pka/len(pkas[residue_key].keys())
            residue.Nmass = sum_Nmass/len(pkas[residue_key].keys())
            residue.Emass = sum_Emass/len(pkas[residue_key].keys())
            residue.Nlocl = sum_Nlocl/len(pkas[residue_key].keys())
            residue.Elocl = sum_Elocl/len(pkas[residue_key].keys())
            residue.determinants = [[], [], []]
            for type in range(3):
                for key in determinants[residue_key][type].keys():
                    value = determinants[residue_key][type][key] / len(pkas[residue_key].keys())
                    if abs(value) > 0.005:  # <-- removing determinant that appears as 0.00
                        newDeterminant = Determinant(key, value)
                        residue.determinants[type].append(newDeterminant)

        return

    def calculateConfigurationPKA(self, version=None, options=None):
        """ 
        Calculates the pKa values
        """
        # 1. Initialize/clean up residue properties (Nmass, Emass, Nlocl, Elocl, buried, determinants)
        for chain in self.chains:
            chain.cleanupResidues()

        # 2. calculating the desolvation contribution
        # using square distance for selecting close atoms
        version.desolv_cutoff_sqr = version.desolv_cutoff*version.desolv_cutoff
        version.buried_cutoff_sqr = version.buried_cutoff*version.buried_cutoff
        for chain in self.chains:
            chain.calculateDesolvation(self.atoms, version=version, options=options)

        # 3. calculating the pKa determinants/perturbations
        # setup interaction list for acid and base back-bone interactions
        backbone_interactions = self.makeBackBoneInteractionList()

        # setting non-iterative back-bone interaction determinants
        determinants.setBackBoneDeterminants(backbone_interactions, version=version)

        # setting ion determinants
        determinants.setIonDeterminants(self, version=version)

        # calculating the back-bone reorganization/desolvation term - testing new term
        version.calculateBackBoneReorganization(self)

        # setting remaining non-iterative and iterative side-chain & Coulomb interaction determinants
        determinants.setDeterminants(self.propka_residues, version=version, options=options)

        # 4. calculating the total pKa values
        for chain in self.chains:
            chain.calculateTotalPKA()

        # 6. printing out timings
        # pka_print("   desolvation   %lf s" % (t1-t0))
        # pka_print("   making lists  %lf s" % (t2-t1))
        # pka_print("   back-bone     %lf s" % (t3-t2))
        # pka_print("   determinants  %lf s" % (t4-t3))

        return

    def makeBackBoneInteractionList(self):
        """ 
        making a nice 'object' to make the acid/base interactions nice and easy to loop over in 'determinants'
        backbone_interactions = [[acids, NH], [bases, CO]]
        """
        backbone_interactions = []
        # --- ACIDS ---
        residues = []
        for resName in lib.residueList("acids"):
            residues.extend(self.residue_dictionary[resName])
        backbone_interactions.append([residues, self.NHlist])
        # --- BASES ---
        residues = []
        # famous propka2.0 exceptions
        for resName in ["HIS"]:
            residues.extend(self.residue_dictionary[resName])
        backbone_interactions.append([residues, self.COlist])

        return backbone_interactions

    def setAvailableConfigurations(self, configurations=None, options=None):
        """ 
        set the 'available configuration keys' for this protein
        """
        if configurations == None:
            pka_print("need to specify a list of available configuration keys")
            exit(8)
        self.configurations = configurations

        return

    def setConfiguration(self, key=None, options=None):
        """ 
        set the 'current possition' to a 'configuration'
        """
        if options.verbose == True:
            pka_print("switching to configuration %6s (protein)" % (key))
        for chain in self.chains:
            chain.setConfiguration(key=key)

        return

    def setChain(self, chainID):
        """ 
        force a chainID, e.g. 'A'
        """
        for chain in self.chains:
            chain.setChain(chainID)

    def renumberAtoms(self):
        """ 
        remunber all protein atoms
        """
        number = 0
        for chain in self.chains:
            for residue in chain.residues:
                for atom in residue.atoms:
                    number += 1
                    atom.setProperty(numb=number)

    def shiftResidueNumber(self, shift):
        """ 
        Shift the residue numbers with 'shift'
        """
        for chain in self.chains:
            chain.shiftResidueNumber(shift)

    def printPKA(self):
        """ 
        Prints the result for current configuration (determinants + summary)
        """
        output.printPKASection(self)

    def printResult(self):
        """ 
        Prints all the resulting output from determinants and down
        """
        output.printResult(self)

    def checkPKA(self, version=None, options=None):
        """ 
        checks calculated pKa values against tabulated
        """
        from test import testCalculatedPKA
        testCalculatedPKA(self, version=version, options=options)

    def writePKA(self, filename=None, reference="neutral", direction="folding", options=None):
        """ 
        Writes a pkafile based on what is in the protein object
        """
        output.writePKA(self, filename=filename, reference=reference, direction=direction, options=options)

    def writePQR(self, label=None, options=None):
        """ 
        Writes a new pqrfile based on what is in the protein
        """
        output.writePQR(self, label=label, options=options)

    def writePDB(self, file=None, filename=None, all_configuration=False, hydrogens=False, options=None):
        """ 
        Writes a new pdbfile based on what is in the protein
        """

        configurations = lib.get_sorted_configurations(self.configurations)

        if len(configurations) == 1 or all_configuration == False:
            output.writePDB(self, filename=filename, hydrogens=hydrogens)
        else:
            # file is opened here
            if file == None:
                if filename == None:
                    filename = "%s.pdb" % (protein.name)
                file = open(filename, 'w')
            # write configurations
            for configuration in configurations:
                self.setConfiguration(configuration)
                file.write("MODEL%9d\n" % (int(configuration[1])))
                output.writePDB(self, file=file, hydrogens=hydrogens)
                file.write('ENDMDL\n')
            # file is closed here
            if file == None:
                file.close()

    def getPHopt(self, reference="neutral", direction="folding", grid=[0., 14., 0.1], options=None):
        """ 
        returns the pH of optimum stability
        """
        if "pH-opt" in self.property:
            return self.property["pH-opt"]
        else:
            profile = self.getFoldingProfile(reference=reference, direction=direction, grid=grid, options=options)
            dG_opt = 999.
            for pH, dG in profile:
                if dG < dG_opt:
                    pH_opt = pH
                    dG_opt = dG
            self.property["pH-opt"] = [pH_opt, dG_opt]
            return self.property["pH-opt"]

    def getStabilityRange(self, reference="neutral", direction="folding", grid=[0., 14., 0.1], options=None):
        """ 
        returns the range where the protein is 'stable'
        """
        if "stable" in self.property:
            return self.property["stable"]
        else:
            profile = self.getFoldingProfile(reference=reference, direction=direction, grid=grid, options=options)
            pH_min = None
            pH_max = None
            for pH, dG in profile:
                if pH_min == None and dG < 0.00:
                    pH_min = pH
                if pH_min != None and dG < 0.00:
                    pH_max = pH
            self.property["stable"] = [pH_min, pH_max]
            return self.property["stable"]

    def getDG80(self, reference="neutral", direction="folding", grid=[0., 14., 0.1], options=None):
        """ 
        returns the dG 80 values
        """
        if "dG80" in self.property:
            return self.property["dG80"]
        else:
            profile = self.getFoldingProfile(reference=reference, direction=direction, grid=grid, options=options)
            pH_opt, dG_opt = self.getPHopt(reference=reference, direction=direction, grid=grid, options=options)
            dG80_min = None
            dG80_max = None
            pH_min = None
            pH_max = None
            for pH, dG in profile:
                if dG80_min == None and dG < dG_opt*0.80:
                    pH_min = pH
                    dG80_min = dG
                if dG80_min != None and dG < dG_opt*0.80:
                    pH_max = pH
                    dG80_max = dG
            self.property["dG80"] = [pH_min, pH_max]
            return self.property["dG80"]

    def getPI(self):
        """ 
        returns the protein isoelectric point (pI)
        """
        if "pI" in self.property:
            return self.property["pI"]
        else:
            pI = calculate.pI(self)
            self.property["pI"] = pI
            return self.property["pI"]

    def getFoldingProfile(self, reference="neutral", direction="folding", grid=[0., 14., 0.1], options=None):
        """ 
        returns the folding energy profile
        """
        if "folding-profile" not in self.property:
            profile = []
            pH, end, increment = grid
            while pH <= end:
                dG = self.calculateFoldingEnergy(pH=pH, options=options)
                profile.append([pH, dG])
                pH += increment

            self.property["folding-profile"] = profile

        return self.property["folding-profile"]

    def getTmProfile(self, reference="neutral", grid=[0., 14., 0.1], Tm=None, Tms=None, ref=None, options=None):
        """ 
        returns the Tm profile
        """
        if "Tm-profile" not in self.property:
            if False:
                Nres = 0
                for chain in self.chains:
                    Nres += len(chain.residues)
                dS = 0.0173*Nres
                pH, Tm_ref = Tm
                dG_ref = self.calculateFoldingEnergy(pH=pH, options=options)
                profile = []
                pH, end, increment = grid
                while pH <= end:
                    dG = self.calculateFoldingEnergy(pH=pH, options=options)
                    dTm = -4.187*(dG - dG_ref)/dS
                    profile.append([pH, Tm_ref+dTm])
                    pH += increment
            else:
                self.property["Tm-profile"] = calculate.TmProfile(self, reference=reference, grid=grid, Tm=Tm, Tms=Tms, ref=ref, options=options)

        return self.property["Tm-profile"]

    def getChargeProfile(self, grid=[0., 14., 1.], options=None):
        """ 
        returns the charge profile
        """
        if "charge-profile" not in self.property:
            profile = []
            pH, end, increment = grid
            while pH <= end:
                Q_pro, Q_mod = self.calculateCharge(pH)
                profile.append([pH, Q_pro, Q_mod])
                pH += increment

            self.property["charge-profile"] = profile

        return self.property["charge-profile"]

    def calculateCharge(self, pH):
        """ 
        Calculates the protein charge
        """
        Q_pro = 0.00
        Q_mod = 0.00
        for chain in self.chains:
            dQ_pro, dQ_mod = chain.calculateCharge(pH)
            Q_pro += dQ_pro
            Q_mod += dQ_mod

        return Q_pro, Q_mod

    def mutateToAla(self, mutation):
        """ 
        mutates a residue to 'ALA'
        """
        label = ""
        for i in range(len(mutation)):
            if mutation[i] in string.digits:
                label += mutation[i]
        residue_number = int(label)

        for chain in self.chains:
            for residue in chain.residues:
                if residue.resNumb == residue_number:
                    residue.mutateToAla()

    def makeMutant(self, mutation=None, atoms=None, version=None, options=None):
        """ 
        mutates the protein according to 'mutation' using 'method'
        """
        import mutate
        newProtein = mutate.makeMutatedProtein(self, mutation=mutation, atoms=atoms, options=options)

        return newProtein

    def optimizeMutationDeterminants(self, mutation=None, atoms=None, alignment=None, version=None, options=None):
        """ 
        permutes multiple mutations and determins the most stable combination; note, you need the version for stability calculations
        """
        import mutate
        best_mutation = mutate.optimizeMutationDeterminants(self, mutation=mutation, atoms=atoms,
                                                            alignment=alignment, version=version, options=options)

        return best_mutation

    def optimizeMultipleMutations(self, mutations=None, atoms=None, alignment=None, version=None, options=None):
        """ 
        permutes multiple mutations and determins the most stable combination; note, you need the version for stability calculations
        """
        import mutate
        best_mutation = mutate.optimizeMultipleMutations(self, mutations=mutations, atoms=atoms, alignment=None, version=version, options=options)

        return best_mutation

    def makeSequence(self, mutation=None, program=None, options=None):
        """ 
        making the protein sequence from chains - note, chains defined differetly depending on program convention
        """
        sequence = ""
        for chain in self.chains:
            sequence += chain.makeSequence(mutation=mutation)
            if program == "scwrl":
                """ do nothing """
            elif chain != self.chains[-1]:
                sequence += "/"
            else:
                """ do nothing ? """
                #pka_print("multiple-chain convention not verified in sequence")
                # exit(9)

        return sequence

    def setAlignment(self, filename=None, **argv):
        """ 
        sets the alignment information (from self for now)
        """
        if filename == None:
            alignments = None
        else:
            first, last, alignments = lib.readAlignments(filename, self.name)

        i = 0
        for chain in self.chains:
            if alignments == None:
                chain.setAlignment()
            else:
                chain.setAlignment(alignment=alignments[i])
                i += 1

    def getResidue(self, label=None):
        """ 
        returns the residue with label='label'
        """
        for chain in self.chains:
            for residue in chain.residues:
                if residue.label == label:
                    return residue
                    break

    def getPKA(self, label=None):
        """ 
        returns the pka value of 'label'
        """
        for chain in self.chains:
            for residue in chain.residues:
                if residue.label == label:
                    return residue.pKa_pro
                    break

    def calculateTitrationCurve(self, label=None, grid=[0., 14., 0.1], options=None):
        """ 
        Calculates the titration curve of residue 'label'
        """
        if label == None:
            pka_print("Must specify residue label, cannot calculate titration curve on whole protein")
            exit(8)
        else:
            residue = self.getResidue(label=label)

        return residue.calculateTitrationCurve(grid=grid)

    def printStabilityBars(self, options=None):
        """ 
        Calculates the pKa-dependant folding energy
        """
        pka_print("\n# stability bars for %s" % (self.name))
        targets = ["GLU", "ASP", "HIS", "LYS", "ARG"]
        targets = ["C- ", "ASP", "GLU", "HIS", "CYS", "TYR", "LYS", "ARG", "N+ "]
        for target in targets:
            for chain in self.chains:
                for residue in chain.residues:
                    if residue.resName == target:
                        dG = 1.36*residue.Q*(residue.pKa_pro - residue.pKa_mod)
                        str = "%6.2lf %6.2lf" % (residue.pKa_pro, dG)
                        str += "  %s" % (residue.label)
                        pka_print(str)
            pka_print("")

        return None

    def calculateFoldingEnergy(self, pH=None, reference=None, options=None):
        """ 
        Calculates the pKa-dependant folding energy
        """
        dG = 0.00
        for chain in self.chains:
            dG += chain.calculateFoldingEnergy(pH=pH, reference=reference, options=options)

        return dG

    def compareWithExperiment(self, file, list, set=None, labels=None, restype=None, experiment=None, options=None):
        """ 
        Compares the calculated pKa values with experiment if it is found in 'experiment'
        """

        str = "# %s\n" % (self.name)
        file.write(str)

        # set 'protein' based on pdbcode
        for protein in experiment.keys():
            str = "%s" % (protein)
            # pka_print(str)
            #str += " %s" % (experiment[protein]['pdb'])
            pka_print(str)
            if experiment[protein]['pdb'] == self.name:
                break

        # setting a list of labels to work with
        if labels == None or labels == "ALL":
            labels = []
            for label in experiment[protein].keys():
                labels.append(label)

        # iterating that list
        for label in lib.extractResidueType(labels, restype=restype, sort=False):
            residue = self.getResidue(label=label)
            pka_exp = experiment[protein][label]
            pka_print(label)
            dpka_exp = pka_exp - residue.pKa_mod
            dpka_clc = residue.pKa_pro - residue.pKa_mod
            diff = residue.pKa_pro - pka_exp
            #pka_print( "compare:%s%6.2lf%6.2lf (%6.2lf%6.2lf) %6.2lf" % (label, pka_exp, residue.pKa_pro, dpka_exp, dpka_clc, diff) )
            #pka_print( "xxx%6.2lf%6.2lf " % (dpka_exp, residue.Emass+residue.Elocl) )
            list.append(diff)
            str = "%8.2lf%8.2lf" % (dpka_exp, dpka_clc)
            #str += " %s" % (residue.label)
            str += "\n"
            #str = "%8.2lf%8.2lf\n" % (dpka_exp, residue.Emass+residue.Elocl)
            #str = "%8.2lf%8.2lf%6.2lf\n" % (dpka_exp, residue.Emass, residue.Elocl)
            #str = "%8.2lf%8.2lf\n" % (residue.buried, residue.Emass+residue.Elocl)
            if False:
                coulomb = 0.00
                for determinant in residue.determinants[2]:
                    coulomb += determinant.value
                str = "%8.2lf%8.2lf\n" % (dpka_exp, coulomb)
            file.write(str)


# ----- protein-related methods to make life easier -----


def getResidueParameters():
    """ 
    Reads necessary information about residues (includes ions)
    """
    from .parameters_new import resName2Type, getQs, pKa_mod
    resInfo = {}
    # reading residue information from parameters.py
    resInfo['resType'] = resName2Type()
    resInfo['Q'] = getQs()
    resInfo['pKa'] = pKa_mod()
    resInfo['type'] = {'C- ': "C-terminus",
                       'N+ ': "N-terminus"}
    # setting up 'type' = 'amino-acid' for residues
    for resName in resInfo['resType'].keys():
        if resName not in resInfo['type']:
            resInfo['type'][resName] = "amino-acid"

    # reading ion information
    filename = os.path.join(os.path.dirname(__file__), 'ions.list')
    file = open(filename, 'r')
    for line in file.readlines():
        if (line[0] == '#'):
            """ comment line """
        else:
            words = line.split()
            key = "%-3s" % (words[0])
            resInfo['resType'][key] = "ION"
            resInfo['type'][key] = "ion"
            resInfo['Q'][key] = float(words[1])
            resInfo['pKa'][key] = None
    file.close()

    return resInfo
