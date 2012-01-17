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
import string, sys, copy
import lib

excluded_resNames = ["H2O", "HOH", "SO4", "PO4", "PEG", "EPE", "NAG", "TRS"]



def readPDB(filename, file=None, verbose=True, tags = ["ATOM"]):
    """
    Reads the pdb line and returns an atom-dictionary; atoms sorted according to chain and residue
    """
    atoms   = {}   # atoms[chainID][resLabel] = [A1, A2, A3 ...]

    # check the filename and open the file
    if file is None:
        file = openPdbFile(filename)

    # scan pdbfile for models and configurations
    number_of_configurations, models_configurations, models_lines = scanFileForConfigurations(file)

    file.close()

    # creating the atom objects. 
    # Note, only reading 'first' configuration, not adding anything to 'self.configurations'
    done_atoms = False
    for model in models_lines:
      if done_atoms == True:
        break
      for line in model:
        if line[16] == 'A' or line[16] == ' ' or line[16] == '1':
          newatom = Atom(line=line, verbose=verbose)
          key = newatom.makeResidueLabel()
          if newatom.chainID in atoms:
            if key in atoms[newatom.chainID]:
              # adding new atom to existing residue
              atoms[newatom.chainID][key].append(newatom)
            else:
              # adding new residue to existing chainID
              atoms[newatom.chainID]["keys"].append(key)
              atoms[newatom.chainID][key] = [newatom]
          else:
            # adding new chain to existing library
            atoms[newatom.chainID] = {"keys": [key]}
            atoms[newatom.chainID][key] = [newatom]
            done_atoms = True
        else:
          # printing alternative configurations
          """print(line)"""


    # creating and adding the existing configurations to atoms
    i_model = 0
    for current_model in models_lines:
      for line in current_model:
        chainID = getChainID(line)
        reskey = "%-3s%4d%2s" % (getResName(line), getResNumb(line), chainID)
        if chainID not in atoms:
          print("incorrect labeling in %s, please correct your pdbfile; could not find chain '%s'" % (filename, chainID)); sys.exit(8)
        if reskey not in atoms[chainID]:
          print("incorrect labeling in %s, please correct your pdbfile; could not find '%s'" % (filename, reskey)); sys.exit(8)
        for atom in atoms[chainID][reskey]:
          if atom.name == getAtomName(line):
            key = makeConfigurationKey(line, i_model)
            atom.configurations[key] = makeConfiguration(line)
      i_model += 1


    # some debugging printouts
    if False:
      for chainID in sorted( atoms.keys() ):
        for key in atoms[chainID]["keys"]:
          for atom in atoms[chainID][key]:
            str = "%s%4d  %4s%3d%7s" % (atom.resName, atom.resNumb, atom.name, len(atom.configurations.keys()), atom.type)
            for key in atom.configurations.keys():
              str += "%5s" % (key)
            print(str)

    #print(number_of_configurations, models_configurations)
    #sys.exit(9)

    return  atoms


def getResNumb(line):
    """
    reads the resNumb from the pdbline
    """
    if line == None:
      return 0
    else:
      return int( line[23:26].strip() )


def getResName(line):
    """
    reads the resName from the pdbline
    """
    if line == None:
      return ""
    else:
      return "%-3s" % (line[17:20].strip())


def getAtomName(line):
    """
    reads the name from the pdbline
    """
    if line == None:
      return ""
    else:
      return line[12:16].strip()


def getElement(line):
    """
    Chresten's stuff, have no idear why its like this
    """
    if line == None:
      return  ""
    elif len(line) > 75:
      element = line[76:78].strip()
    else:
      element = line[12:14].strip().strip(string.digits)
      # Xplor exception to HE, HD, HG, HH etc.
      if len(element) > 1 and element[0] == "H":
          element = element[0]

    return  element


def getAtomNumb(line):
    """
    reads the numb from the pdbline
    """
    if line == None:
      return 0
    else:
      return int( line[ 6:11].strip() )


def getChainID(line):
    """
    reads the chainID from the pdbline
    """
    if   line == None:
      return "A"
    elif line[21] == " ":
      return "A"
    else:
      return line[21]


def getOccupation(line):
    """
    reads the resName from the pdbline
    """
    if   line == None:
      return 1.0
    elif len(line) > 59:
      if line[56:60] != "    ":
        return float( line[56:60].strip() )


def getBeta(line):
    """
    reads the resName from the pdbline
    """
    if   line == None:
      return 0.0
    elif len(line) > 65:
      if line[60:66] != "      ":
        return float( line[60:66].strip() )


def getType(line):
    """
    reads the resName from the pdbline
    """
    if   line == None:
      return  ""
    else:
      return line[:6].strip().lower()


def makeConfiguration(line):
    """
    returns configuration based on the line
    """
    if line == None:
      configuration = [0.0, 0.0, 0.0]
    else:
      x = float( line[30:38].strip() )
      y = float( line[38:46].strip() )
      z = float( line[46:54].strip() )
      configuration = [x, y, z]

    return configuration


def makeConfigurationKey(line, i_model):
    """
    returns a configuration key
    """
    if line[16] == " " or line[16] == "1":
      return "M%dC%s" % (i_model, "A")
    else:
      return "M%dC%s" % (i_model, line[16])


def openPdbFile(filename):
    """
    check the name and open the pdbfile
    """
    root, extension = lib.splitFileName(filename)
    if   extension == "pdb":
        # all good, it's a pdbfile
        file = open(filename)
    elif extension == None:
        # probably trimmed name
        newname = "%s.pdb" % (filename)
        file = open(newname)
    elif extension == "ini":
        # this is a modeller ini pdbfile
        file = open(filename)
    elif extension[:-1] == "pdb":
        # this might be a biological unit from the ProteinDataBank
        file = open(filename)
    elif extension == "sa":
        # this might be a simulated annealing file from Xplor
        file = open(filename)
    else:
        print("trying to pass me a rotten pdbfile \"%s\" - no pdb extension" % (filename))
        print("check if there is a dot in full path ...")
        sys.exit(9)

    return  file


def scanFileForConfigurations(file, tags=["ATOM"], options=None):
    """
    Scanning through the pdbfile for models and configurations
    """
    acceptedAtomTypes = ["ATOM", "HETATM"]
    model_configurations = [[]]
    model_lines          = [[]]
    current_model_configurations = model_configurations[-1]
    current_model_lines  = model_lines[-1]

    while True:
        line = file.readline()
        if line == "":
          break
        line = line.strip()
        if line[:5] == "MODEL":
          model_configurations.append(["A"])
          current_model_configurations = model_configurations[-1]
          model_lines.append([])
          current_model_lines  = model_lines[-1]
        if line[:6] == "ENDMDL":
          # reset to 'first/0' model
          current_model_configurations = model_configurations[0]
          current_model_lines  = model_lines[0]
        record = line[0:6].strip()
        if record in acceptedAtomTypes:
          if getResName(line) not in excluded_resNames:
            current_model_lines.append(line)
            configuration = line[16]
            if configuration == " ":
              configuration = "A"
            if configuration not in current_model_configurations:
              current_model_configurations.append(configuration)

    number_of_configurations = 0
    for current_model_configurations in model_configurations:
      number_of_configurations += len(current_model_configurations)

    return number_of_configurations, model_configurations, model_lines



class Atom:
    """
      Atom class - contains all atom information found in the pdbfile
    """

    def __init__(self, line=None, verbose=False):

     
        self.name      =   getAtomName(line)
        self.numb      =   getAtomNumb(line)
        self.resName   =   getResName(line)
        self.resNumb   =   getResNumb(line)
        self.chainID   =   getChainID(line)
        self.occ       =   getOccupation(line)
        self.beta      =   getBeta(line)
        self.type      =   getType(line)
        self.element   =   getElement(line)
        self.x, self.y, self.z = makeConfiguration(line)
        self.configurations = {}
        
        self.bonded_atoms = []
        self.residue = None
        self.charge = 0
        self.steric_number = 0
        self.number_of_lone_pairs = 0
        self.number_of_protons_to_add = 0
        self.number_of_pi_electrons_in_double_and_triple_bonds = 0
        self.number_of_pi_electrons_in_conjugate_double_and_triple_bonds = 0


    def makeResidueLabel(self):
            """
            making a key = residue.label in readPDB()
            """
            return lib.makeResidueLabel(self.resName, self.resNumb, self.chainID)


    def translate(self, vector):
            """
            print Atom information
            """
            self.x += vector[0]
            self.y += vector[1]
            self.z += vector[2]
            for key in self.configurations.keys():
              for i in range(3):
                self.configurations[key][i] += vector[i]


    def printAtom(self):
            """
            print Atom information
            """
            str  = ""
            str += " %s" % (self.resName)
            str += " %6d" % (self.resNumb)
            str += " %s" % (self.chainID)
            str += "  %s" % (self.name)
            print(str)


    def trimConfigurations(self, configurations=None):
            """
            checks the configurations to make sure a 'align-mutated' residue contains the right configuration
            do not allow change in M value - not considered enough
            """

            return


    def setConfigurationPosition(self, key=None):
            """
            set the position of a 'configuration' to 'current position'
            """
            configuration = [self.x, self.y, self.z]
            self.configurations[key] = configuration


    def setConfiguration(self, key=None):
            """
            set the 'current possition' to a 'configuration'
            """
            if   key in self.configurations:
              self.x = self.configurations[key][0]
              self.y = self.configurations[key][1]
              self.z = self.configurations[key][2]
            elif len(self.configurations) == 1:
              # get single key if only one configuration: saving back-bone protonation when previous residue doesn't have 'key'
              for default_key in self.configurations.keys():
                break
              self.x = self.configurations[default_key][0]
              self.y = self.configurations[default_key][1]
              self.z = self.configurations[default_key][2]
            elif True:
              # get single key if only one configuration: saving back-bone protonation when previous residue doesn't have 'key'
              for default_key in self.configurations.keys():
                if key[:-2] == default_key[:-2]:
                  break
              self.x = self.configurations[default_key][0]
              self.y = self.configurations[default_key][1]
              self.z = self.configurations[default_key][2]
            else:
              resLabel = "%-3s%4d%2s" % (self.resName, self.resNumb, self.chainID)
              keys = ""
              for item in self.configurations.keys(): keys += "%5s" % (item)
              print("configuration '%s' not found in '%s' atom '%s' [%s]" % (key, resLabel, self.name, keys))
              sys.exit(8)



    def setProperty(self,
                    numb    = None, 
                    name    = None, 
                    resName = None, 
                    chainID = None,
                    resNumb = None,
                    x       = None,
                    y       = None,
                    z       = None,
                    occ     = None,
                    beta    = None,
                    element = None):
        """
        sets properties of the atom object
        """

        if numb    != None: self.numb    = numb
        if name    != None: self.name    = name
        if resName != None: self.resName = resName
        if chainID != None: self.chainID = chainID
        if resNumb != None: self.resNumb = resNumb
        if x       != None: self.x       = x
        if y       != None: self.y       = y
        if z       != None: self.z       = z
        if occ     != None: self.occ     = occ
        if beta    != None: self.beta    = beta
        if element != None: self.element = element



    def makeCopy(self, 
                    numb    = None,
                    name    = None,
                    resName = None,
                    chainID = None,
                    resNumb = None,
                    x       = None,
                    y       = None,
                    z       = None,
                    occ     = None,
                    beta    = None,
                    configs = None,
                    element = None):
        """
        making a copy of this atom
        """
        if numb    == None:  numb    = self.numb
        if name    == None:  name    = self.name
        if resName == None:  resName = self.resName
        if chainID == None:  chainID = self.chainID
        if resNumb == None:  resNumb = self.resNumb
        if x       == None:  x       = self.x
        if y       == None:  y       = self.y
        if z       == None:  z       = self.z
        if occ     == None:  occ     = self.occ
        if beta    == None:  beta    = self.beta
        if configs == None:  configs = sorted(self.configurations.keys())
        if element == None:  element = self.element

        line = self.makePDBLine(name=name, 
                                numb=numb, 
                                resName=resName, 
                                resNumb=resNumb, 
                                x=x, y=y, z=z, 
                                chainID=chainID, 
                                element=element)
        newAtom =  Atom(line)
        # copy configurations
        for key in configs:
          newAtom.configurations[key] = [x, y, z]

        return  newAtom


    def makePDBLine(self,
                    numb    = None,
                    name    = None,
                    resName = None,
                    chainID = None,
                    resNumb = None,
                    x       = None,
                    y       = None,
                    z       = None,
                    occ     = None,
                    beta    = None,
                    element = None):
        """
        returns a pdb ATOM-line for various purposes;
        specifying arguments over-writes.
        """
        if numb    == None: numb    = self.numb
        if name    == None: name    = self.name
        if resName == None: resName = self.resName
        if chainID == None: chainID = self.chainID
        if resNumb == None: resNumb = self.resNumb
        if x       == None: x       = self.x
        if y       == None: y       = self.y
        if z       == None: z       = self.z
        if occ     == None: occ     = self.occ
        if beta    == None: beta    = self.beta
        if element == None: element = self.element

        if len(name) > 4:
          name = name[:4]

        # making pdb-string
        str  = "ATOM "
        str += "%6d"      % (numb)
        if   len(element) == 2:
          str += " %-5s"  % (name)
        elif name[0] in "0123456789":
          str += " %-5s"  % (name)
        else:
          str += "  %-4s" % (name)
        str += "%s"       % (resName)
        str += "%2s"      % (chainID)
        str += "%4d"      % (resNumb)
        str += "%12.3lf"  % (x)
        str += "%8.3lf"   % (y)
        str += "%8.3lf"   % (z)
        str += "%6.2lf"   % (occ)
        str += "%6.2lf"   % (beta)
        str += "%12s"     % (element)

        return str


    def __str__(self):
        return '%5d-%4s %5d-%3s (%1s) [%8.3f %8.3f %8.3f]' %(self.numb, self.name, self.resNumb, self.resName, self.chainID, self.x, self.y, self.z)
            

    def get_element(self):
        """ try to extract element if not already done"""
        if self.element == '':
          if self.name[0] in "0123456789":
            self.element = self.name[1]
          else:
            self.element = self.name[0]
        return self.element
    

    def set_residue(self, residue):
        """ Makes a references to the parent residue"""
        if self.residue == None:
            self.residue = residue

 
