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
import sys
import string
import math
import copy
from . import lib
from .pdb import Atom
pka_print = lib.pka_print


class Residue:
    """
        Residue class - contains atoms and properties of the residue
    """

    def __init__(self, atoms, resName=None, resNumb=None, chainID=None, resInfo=None, options=None):
        """
        Constructer of the residue object.
        """
        self.label = None
        self.atoms = []
        if chainID == None:
            self.chainID = atoms[0].chainID
        else:
            self.chainID = chainID
        if resNumb == None:
            self.resNumb = atoms[0].resNumb
        else:
            self.resNumb = resNumb
        if resName == None:
            self.resName = atoms[0].resName
        else:
            self.resName = resName
        self.resType = None                         # determins the interaction parameters, e.g. 'COO'
        self.Q = None                         # residue charge
        self.type = None                         # 'amino-acid' / 'ion' / 'ligand' / 'N-terminus' / 'C-terminus'
        self.buried = None
        self.location = None
        self.determinants = [[], [], []]
        self.Nmass = 0
        self.Nlocl = 0
        self.Emass = 0.0
        self.Vmass = 0.0
        self.Elocl = 0.0
        self.Vlocl = 0.0
        self.x = 0.00
        self.y = 0.00
        self.z = 0.00
        self.pKa_mod = None
        self.pKa_pro = None
        self.pKas = []                           # list with several pKa-objects: not used yet
        self.center_atoms = []                    # list of atoms constituting the 'residue center', needed!
        self.default_key = None                  # the 'default' configuration if not all are there
        self.configurations = []                    # keys to the atom configurations belonging to this residue
        self.coupled_residues = []

        # setting residue information i.e. resType, type, Q, pKa_mod & pKa_pro from resInfo dictionary
        self.setResidueInformation(resInfo=resInfo)

        # setting residue label of the residue to create
        self.setResidueLabel()

        # setting up the atom labels used to calculate the 'residue center'
        residue_center_atom_labels = lib.residueCenterAtomList(self.resName)

        # adding/setting various atom-related properties such as atoms, configurations etc.
        for atom in atoms:

            # checking if this residue has 'HETATM' atoms; good chance it's a ligand then
            if atom.type == "hetatm" and self.type == None:
                self.type = "ligand"

            # setting the number of configurations for this residue
            for key in atom.configurations.keys():
                if key not in self.configurations:
                    self.configurations.append(key)

            # if atom.name[0] != 'H':
            if True:
                atom.set_residue(self)
                self.atoms.append(atom)

            # setting 'center atoms'; needs to be on self since it is reset when you switch configurations
            if len(residue_center_atom_labels) == 0:
                self.center_atoms.append(atom)
            elif atom.name in residue_center_atom_labels:
                self.center_atoms.append(atom)

        # set residue center, i.e. give 'x, y, z' values
        if len(self.center_atoms) > 0:
            self.setResidueCenter()

    def setResidueInformation(self, resInfo=None):
        """
        set residue information based on resName - it is set here since it is a convenience thing
        """
        # resType - determines interaction parameters
        if self.resName in resInfo['resType']:
            self.resType = resInfo['resType'][self.resName]

        # Q - group charge
        if self.resName in resInfo['Q']:
            self.Q = resInfo['Q'][self.resName]
        else:
            self.Q = 0.00

        # type - 'amino-acid' / 'ion' / 'ligand' / 'N-terminus' / 'C-terminus'
        if self.resName in lib.residueList("standard"):
            self.type = "amino-acid"
        elif self.resName in resInfo['type']:
            self.type = resInfo['type'][self.resName]

        # pKa_mod - model or water pKa value
        if self.resName in resInfo['pKa']:
            self.pKa_mod = resInfo['pKa'][self.resName]
            self.pKa_pro = resInfo['pKa'][self.resName]
        else:
            self.pKa_mod = resInfo['pKa']['default']
            self.pKa_pro = resInfo['pKa']['default']

    def setConfiguration(self, key=None):
        """
        set the 'current possition' to a specific 'configuration', or default if it doesn't exist
        """
        self.cleanupPKA()
        if key in self.configurations:
            configuration = key
        else:
            configuration = self.default_key

        for atom in self.atoms:
            atom.setConfiguration(key=configuration)
        self.setResidueCenter()

    def cleanupPKA(self):
        """
        Initializing/cleaning up residue!
        """
        self.Nmass = 0
        self.Nlocl = 0
        self.Emass = 0.00
        self.Elocl = 0.00
        self.buried = 0.00
        self.pKa_pro = self.pKa_mod
        self.determinants = [[], [], []]

    def setResidueCenter(self):
        """
        sets the center of the residue based on center_atoms
        """
        number_of_atoms = len(self.center_atoms)
        self.x = 0.00
        self.y = 0.00
        self.z = 0.00
        for atom in self.center_atoms:
            self.x += atom.x
            self.y += atom.y
            self.z += atom.z
        if number_of_atoms > 0:
            self.x = self.x/number_of_atoms
            self.y = self.y/number_of_atoms
            self.z = self.z/number_of_atoms

    def getThirdAtomInAngle(self, atom=None):
        """
        finds and returns the third atom in angular dependent interactions
        expecting one of ["HIS", "ARG", "AMD", "TRP"]
        """
        if self.resName == "HIS":
            if atom.name == "HD1":
                return self.getAtom(name="ND1")
            elif atom.name == "HE2":
                return self.getAtom(name="NE2")
        elif self.resName == "ARG":
            if atom.name in ["HE"]:
                return self.getAtom(name="NE")
            elif atom.name in ["1HH1", "2HH1"]:
                return self.getAtom(name="NH1")
            elif atom.name in ["1HH2", "2HH2"]:
                return self.getAtom(name="NH2")
        elif self.resName == "ASN":
            return self.getAtom(name="ND2")
        elif self.resName == "GLN":
            return self.getAtom(name="NE2")
        elif self.resName == "TRP":
            return self.getAtom(name="NE1")

    def getAtom(self, name=None):
        """
        finds and returns the specified atom in this residue.
        """
        for atom in self.atoms:
            if atom.name == name:
                return atom
                break

        return None

    def checkOXT(self):
        """
        Checks that OXT is present or creates it.
        """
        O = self.getAtom(name="O")
        OXT = self.getAtom(name="OXT")

        # NMR Xplor over-write
        if O == None:
            O = self.getAtom(name="OT1")
            if O != None:
                O.setProperty(name="O")
        if OXT == None:
            OXT = self.getAtom(name="OT2")
            if OXT != None:
                OXT.setProperty(name="OXT")

        # continuing after 'NMR exception'; creating OXT if not found
        if OXT == None:
            # did not find OXT, creating it 'on the fly'
            CA = self.getAtom(name="CA")
            C = self.getAtom(name="C")
            if O == None or CA == None or C == None:
                pka_print("ERROR: cannot create OXT atom - missing CA, C, or O atoms; please correct pdbfile")
                sys.exit(8)
            dX = -((CA.x-C.x) + (O.x-C.x))
            dY = -((CA.y-C.y) + (O.y-C.y))
            dZ = -((CA.z-C.z) + (O.z-C.z))
            distance = math.sqrt(dX*dX + dY*dY + dZ*dZ)
            x = C.x + 1.23*dX/distance
            y = C.y + 1.23*dY/distance
            z = C.z + 1.23*dZ/distance
            OXT = C.makeCopy(name="OXT", x=x, y=y, z=z)
            self.atoms.append(OXT)
            pka_print("creating %s atom" % (OXT.name))

        return [O, OXT]

    def fillUnknownConfigurations(self, keys=None, options=None):
        """
        Fills in  the configurations that have not been read
        """
        # getting default key as the first OK element in the sorted protein keys
        for key in keys:
            if key in self.configurations:
                self.default_key = key
                break

        # enforcing all residue keys on each atom
        for atom in self.atoms:
            for configuration in self.configurations:
                if key not in atom.configurations:
                    atom.configurations[configuration] = atom.configurations[self.default_key]

    def checked(self, options=None):
        """
        Checks that I understand all residues.
        """
        excluded_resNames = lib.residueList("excluded")
        if self.resName in excluded_resNames:
            return False
        else:
            return True

    def checkResidue(self, options=None):
        """
        Checks that I understand all residues. 
        """

        residue_list = lib.residueList("all")
        rename = {"LYP": "LYS",
                  "CYX": "CYS",
                  "HSD": "HIS",
                  }

        # renaming residue if in rename dictionary
        if self.resName in rename:
            newName = rename[self.resName]
            if options.verbose == True:
                pka_print("Warning: renaming %s to %s" % (self.resName, newName))
            self.resName = newName
            for atom in self.atoms:
                atom.resName = newName
        else:
            """OK - lets rock"""

        # after rename residue cases, check heavy atoms
        if self.resName in residue_list:
            # OK - lets rock
            self.checkAtoms(options=options)
        # Chresten's stuff
        elif self.type == "ion":
            outstr = "%s%4d -  OK %s" % (self.resName, self.resNumb, self.type)
            pka_print(outstr)
        # elif self.resName in version.ions.keys():
        #    str  = "%-3s%4d - %s with charge %+d" % (self.resName,
        #                                             self.resNumb,
        #                                             version.ions_long_names[self.resName],
        #                                             version.ions[self.resName])
        #    pka_print(str)
        else:
            outstr = "%s%4d - unidentified residue" % (self.resName, self.resNumb)
            pka_print(outstr)

    def checkAtoms(self, options=None):
        """
        Checks that all heavy atoms are there
        """
        outstr = "%s%4d - " % (self.resName, self.resNumb)
        atom_list = lib.atomList(self.resName)
        OK = True
        for name in atom_list:
            FOUND = False
            for atom in self.atoms:
                if atom.name == name:
                    FOUND = True
            if FOUND == False:
                outstr += " %s" % (name)
                OK = False

        if OK == True:
            self.checkConfigurations(verbose=False)
            outstr += " OK (%2d: %2d)" % (len(self.atoms), len(self.configurations))
            if options.verbose == True:
                pka_print(outstr)
        else:
            outstr += " missing"
            pka_print(outstr)

    def checkConfigurations(self, verbose=False):
        """
        checks that all atoms in this residue has the same number of configurations
        """
        for atom in self.atoms:
            for key in self.configurations:
                if key not in atom.configurations:
                    atom.configurations[key] = atom.configurations[self.default_key]

    def printLabel(self):
        """
        prints the residue ID
        """
        outstr = "%s%4d %s" % (self.resName, self.resNumb, self.chainID)
        pka_print(outstr)

    def __str__(self):
        return self.label

    def extractBackBoneAtoms(self):
        """
        Returns the back-bone atoms
        """
        N = None
        H = None
        C = None
        O = None
        for atom in self.atoms:
            if atom.name == "N":
                N = atom
            elif atom.name == "H":
                H = atom
            elif atom.name == "C":
                C = atom
            elif atom.name == "O":
                O = atom

        return N, H, C, O

    def makeDeterminantAtomList(self, resType=None, type=None, version=None):
        """
        Extracting reference to determinant atom - test stage still
        """
        if type == 'base' or type == 'acid':
            pair_type = type
        else:
            if resType in ["HIS", "LYS", "ARG", "N+ "]:
                pair_type = 'base'
            else:
                pair_type = 'acid'

        if self.resName not in version.atomInteractionList[pair_type]:
            pka_print("cannot find atomInteractionList for residue %s in residue.makeDeterminantAtomList()" % (self.resName))
            sys.exit(9)

        # Searching for determinant atom
        atoms = []
        for atom in self.atoms:
            if atom.name in version.atomInteractionList[pair_type][self.resName]:
                atoms.append(atom)

        return atoms

    def calculateTotalPKA(self):
        """
        Calculates the total pKa values from the desolvation and determinants
        """
        back_bone = 0.00
        for determinant in self.determinants[0]:
            value = determinant.value
            back_bone += value

        side_chain = 0.00
        for determinant in self.determinants[1]:
            value = determinant.value
            side_chain += value

        coulomb = 0.00
        for determinant in self.determinants[2]:
            value = determinant.value
            coulomb += value

        self.pKa_pro = self.pKa_mod + self.Emass + self.Elocl + back_bone + side_chain + coulomb

    def calculateIntrinsicPKA(self):
        """
        Calculates the intrinsic pKa values from the desolvation determinants, back-bone hydrogen bonds, 
        and side-chain hydrogen bond to non-titratable residues
        """
        back_bone = 0.00
        for determinant in self.determinants[1]:
            value = determinant.value
            back_bone += value

        side_chain = 0.00
        for determinant in self.determinants[0]:
            if determinant.label[0:3] not in ['ASP', 'GLU', 'LYS', 'ARG', 'HIS', 'CYS', 'TYR', 'C- ', 'N+ ']:
                value = determinant.value
                side_chain += value

        self.intrinsic_pKa = self.pKa_mod + self.Emass + self.Elocl + back_bone + side_chain

        return

    def calculateDesolvation(self, atoms, version=None, options=None):
        """
        Calculates the desolvation contribution
        """
        version.calculateDesolvation(self, atoms, options=options)

    def setChain(self, chainID):
        """
        Set a chainID
        """
        self.chainID = chainID
        for atom in self.atoms:
            atom.chainID = chainID

    def setResidueNumber(self, resNumb):
        """
        Set the residue numbers to 'resNumb'
        """
        self.resNumb = resNumb

    def setResidueLabel(self, label=None):
        """
        Set the residue label to e.g. 'GLU 145 A'
        """
        if label == None:
            self.label = lib.makeResidueLabel(self.resName, self.resNumb, self.chainID)
        else:
            self.label = label

    def shiftResidueNumber(self, shift):
        """
        Shift the residue numbers with 'shift'
        """
        self.resNumb = self.resNumb + shift

    def calculateFoldingEnergy(self, pH=None, reference=None, options=None):
        """
        returning the electrostatic energy of this residue at pH 'pH'
        """
        if pH == None:
            pH = options.pH
        if reference == None:
            reference = options.reference

        # calculating the ddG(neutral --> low-pH) contribution
        if self.resType not in ["COO", "HIS", "N+ ", "CYS", "TYR", "LYS", "ARG"]:
            ddG = 0.00
        else:
            if reference == "low-pH":
                ddG_neutral = 0.00
            else:
                if self.Q > 0.00:
                    pKa_prime = self.pKa_pro
                    coulomb_determinants = self.determinants[2]
                    for determinant in coulomb_determinants:
                        if determinant.value > 0.00:
                            pKa_prime -= determinant.value
                    ddG_neutral = -1.36*(pKa_prime - self.pKa_mod)
                else:
                    ddG_neutral = 0.00

            # calculating the ddG(low-pH --> pH) contribution
            # folded
            x = pH - self.pKa_pro
            y = 10**x
            Q_pro = math.log10(1+y)

            # unfolded
            x = pH - self.pKa_mod
            y = 10**x
            Q_mod = math.log10(1+y)

            ddG_low = -1.36*(Q_pro - Q_mod)
            ddG = ddG_neutral + ddG_low

        return ddG

    def getCharge(self, pH, state):
        """
        returning the charge of this residue at pH 'pH'
        """
        if state == "mod" or state == "unfolded":
            x = self.Q*(self.pKa_mod - pH)
        else:
            x = self.Q*(self.pKa_pro - pH)
            #x =  pH - self.pKa_pro
        y = 10**x
        charge = self.Q*(y/(1.0+y))
        #charge = math.log10(1+y)

        return charge

    def calculateTitrationCurve(self, grid=[0., 14., 0.10]):
        """
        calculates the titration curve of this residue
        """
        if grid == None:
            grid = [self.pKa_pro-2.5, self.pKa_pro+2.5, 0.10]
        state = "folded"
        titration_curve = []
        pH = grid[0]
        stop = grid[1] + grid[2]/2.0
        while pH < stop:
            Q = self.getCharge(pH, state)
            titration_curve.append([pH, Q])
            pH += grid[2]

        return titration_curve

    def getSummaryString(self):
        """
        Writing the summary string
        """
        outstr = "   %s%8.2lf%10.2lf" % (self.label, self.pKa_pro, self.pKa_mod)
        return outstr

    def mutateToAla(self):
        """
        mutating residue to 'ALA'
        Note, if you mutate a ionizable residue to Ala it will remain in 'ionizable_residues list'
        and propka will try to calcualte the pKa of it !!!
        """
        keep_atoms = lib.atomList("ALA")
        self.printLabel()
        self.resName = "ALA"
        self.setResidueLabel()
        self.printLabel()
        new_atoms = []
        for atom in self.atoms:
            if atom.name in keep_atoms:
                new_atoms.append(atom)
        pka_print(self.atoms)
        pka_print(new_atoms)
        self.cleanupResidue
        self.pKa_mod = pKa_mod(self.resName)
        self.pKa_pro = self.pKa_mod
        self.atoms = new_atoms

    def replaceWithResidue(self, new_residue):
        """
        replacing current residue with incoming residue
        """
        self.cleanupResidue()
        self.resName = new_residue.resName
        self.pKa_mod = new_residue.pKa_mod
        self.setResidueLabel()
        exclude_atoms = ["N", "CA", "C", "O"]
        tmp_atoms = []
        for atom in self.atoms:
            if atom.name in exclude_atoms:
                tmp_atoms.append(atom)
        self.atoms = tmp_atoms
        for new_atom in new_residue.atoms:
            if new_atom.name in exclude_atoms:
                """ do nothing """
            else:
                self.atoms.append(new_atom)

    def getDeterminantString(self):
        """
        Everything should be calculated, now, let's print the darn thing and be done!
        """
        # if self.location == "SURFACE" or self.location == "BURIED ":
        BURIED_RATIO = True

        empty_determinant = "%s%4d%2s" % ("XXX", 0, "X")
        number_of_sidechain = len(self.determinants[0])
        number_of_backbone = len(self.determinants[1])
        number_of_coulomb = len(self.determinants[2])
        number_of_determinants = number_of_sidechain + number_of_backbone + number_of_coulomb
        number_of_lines = max(1, number_of_sidechain, number_of_backbone, number_of_coulomb)
        outsting = ""
        #outsting += " number_of_sidechain = %d" % (number_of_sidechain)
        #outsting += " number_of_backbone  = %d" % (number_of_backbone )
        #outsting += " number_of_coulomb   = %d" % (number_of_coulomb  )
        #outsting += " number_of_lines     = %d" % (number_of_lines    )
        # print outsting

        for line_number in range(1, number_of_lines+1):
            outsting += "%s" % (self.label)
            if line_number == 1:
                outsting += " %6.2lf" % (self.pKa_pro)
                if len(self.coupled_residues) > 0:
                    outsting += '*'
                else:
                    outsting += ' '

                if BURIED_RATIO == True:
                    if self.type == "BONDED":
                        outsting += " BONDED "
                    else:
                        outsting += " %4d%2s " % (int(100.0*self.buried), "%")
                else:
                    outsting += "%8s" % (self.type)
                outsting += " %6.2lf %4d" % (self.Emass, self.Nmass)
                outsting += " %6.2lf %4d" % (self.Elocl, self.Nlocl)
            else:
                outsting += "%40s" % (" ")

            # Side-chain determinants
            if line_number > number_of_sidechain:
                outsting += "%8.2lf %s" % (0.0, empty_determinant)
            else:
                determinant = self.determinants[0][line_number-1]
                outsting += "%8.2lf %s" % (determinant.value, determinant.label)

            # Back-bone determinants
            if line_number > number_of_backbone:
                outsting += "%8.2lf %s" % (0.0, empty_determinant)
            else:
                determinant = self.determinants[1][line_number-1]
                outsting += "%8.2lf %s" % (determinant.value, determinant.label)

            # Coulomb determinants
            if line_number > number_of_coulomb:
                outsting += "%8.2lf %s" % (0.0, empty_determinant)
            else:
                determinant = self.determinants[2][line_number-1]
                outsting += "%8.2lf %s" % (determinant.value, determinant.label)

            # adding end-of-line
            outsting += "\n"

        return outsting

    def printResult(self):
        """
        Everything should be calculated, now, let's print the darn thing and be done!
        """
        # if self.location == "SURFACE" or self.location == "BURIED ":
        BURIED_RATIO = True

        empty_determinant = "%s%4d%2s" % ("XXX", 0, "X")
        number_of_sidechain = len(self.determinants[0])
        number_of_backbone = len(self.determinants[1])
        number_of_coulomb = len(self.determinants[2])
        number_of_determinants = number_of_sidechain + number_of_backbone + number_of_coulomb
        number_of_lines = max(1, number_of_sidechain, number_of_backbone, number_of_coulomb)
        outstr = ""
        outstr += " number_of_sidechain = %d" % (number_of_sidechain)
        outstr += " number_of_backbone  = %d" % (number_of_backbone)
        outstr += " number_of_coulomb   = %d" % (number_of_coulomb)
        outstr += " number_of_lines     = %d" % (number_of_lines)
        # print outstr

        if True:
            for line_number in range(1, number_of_lines+1):
                outstr = "%s%4d%2s" % (self.resName, self.resNumb, self.chainID)
                if line_number == 1:
                    outstr += " %6.2lf" % (self.pKa_pro)
                    if BURIED_RATIO == True:
                        outstr += "  %4d%2s " % (int(100.0*self.buried), "%")
                    else:
                        outstr += " %8s" % (self.type)
                    outstr += " %6.2lf %4d" % (self.Emass, self.Nmass)
                    outstr += " %6.2lf %4d" % (self.Elocl, self.Nlocl)
                else:
                    outstr += "%40s" % (" ")

                # Side-chain determinant
                if line_number > number_of_sidechain:
                    outstr += "%8.2lf %s" % (0.0, empty_determinant)
                else:
                    determinant = self.determinants[0][line_number-1]
                    outstr += "%8.2lf %s" % (determinant.value, determinant.label)

                # Back-bone determinant
                if line_number > number_of_backbone:
                    outstr += "%8.2lf %s" % (0.0, empty_determinant)
                else:
                    determinant = self.determinants[1][line_number-1]
                    outstr += "%8.2lf %s" % (determinant.value, determinant.label)

                # Coulomb determinant
                if line_number > number_of_coulomb:
                    outstr += "%8.2lf %s" % (0.0, empty_determinant)
                else:
                    determinant = self.determinants[2][line_number-1]
                    outstr += "%8.2lf %s" % (determinant.value, determinant.label)
                pka_print('%s' % (outstr))
        else:
            outstr = "%s%4d%2s%4d%2d" % (self.resName, self.resNumb, self.chainID, number_of_lines, number_of_determinants)
            pka_print('%s' % (outstr))

        pka_print('')

    def translate(self, translation):
        """
        translate residue according to 'translation'
        """
        for atom in self.atoms:
            atom.x += translation[0]
            atom.y += translation[1]
            atom.z += translation[2]
            for key in atom.configurations.keys():
                for i in range(3):
                    atom.configurations[key][i] += translation[i]

    def rotate(self, axis, theta, center=None):
        """
        rotate residue theta radians around axis with center=center
        """
        from rotate import generalRotationMatrix
        translate = [0.00, 0.00, 0.00]
        number_of_atoms = 0
        for atom in self.atoms:
            if atom.name in center or center == None:
                number_of_atoms += 1
                translate[0] += atom.x/len(self.atoms)
                translate[1] += atom.y/len(self.atoms)
                translate[2] += atom.z/len(self.atoms)
        for atom in self.atoms:
            for i in range(3):
                translate[i] = translate[i]/number_of_atoms

        # translate to rotation center
        for atom in self.atoms:
            atom.x -= translate[0]
            atom.y -= translate[1]
            atom.z -= translate[2]

        # get rotation matrix
        rotation_matrix = generalRotationMatrix(axis, theta)

        # rotating
        new_position = [None, None, None]
        for atom in self.atoms:

            # rotate actual position
            old_position = [atom.x, atom.y, atom.z]
            for xyz in range(3):
                new_position[xyz] = translate[xyz]
                for i in range(3):
                    new_position[xyz] += rotation_matrix[xyz][i]*old_position[i]
            # update position
            atom.x = new_position[0]
            atom.y = new_position[1]
            atom.z = new_position[2]

            # rotate configuration
            for key in atom.configurations.keys():
                for xyz in range(3):
                    new_position[xyz] = translate[xyz]
                    for i in range(3):
                        new_position[xyz] += rotation_matrix[xyz][i]*atom.configurations[key][i]
                for xyz in range(3):
                    atom.configurations[key][xyz] = new_position[xyz]

    def makeCopy(self,
                 chainID=None,
                 resNumb=None,
                 ):
        """
        making a copy of this residue
        """
        from protein import getResidueParameters
        if chainID == None:
            chainID = self.chainID
        if resNumb == None:
            resNumb = self.resNumb

        newAtoms = []
        for atom in self.atoms:
            newAtoms.append(atom.makeCopy(chainID=chainID, resNumb=resNumb))

        resInfo = getResidueParameters()
        newResidue = Residue(newAtoms, resInfo=resInfo)
        newResidue.resType = self.resType
        newResidue.Q = self.Q
        newResidue.type = self.type

        return newResidue
