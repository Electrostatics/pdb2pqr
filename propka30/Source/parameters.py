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
import math
import lib
import sys, os
from calculator import calculate


def pKa_mod(resName):
    """ 
    definitions of min distance 
    """
    pKa_mod = {'C- ': 3.20, 
               'ASP': 3.80, 
               'GLU': 4.50,
               'HIS': 6.50,
               'CYS': 9.00,
               'TYR':10.00, 
               'LYS':10.50,
               'ARG':12.50,
               'N+ ': 8.00}

    if resName in pKa_mod:
        return pKa_mod[resName]
    else:
        return 20.00



def getAtomInteractionList():
    """ 
    storing label of interacting atoms
    """
    labels = {}
    labels['acid'] = {'SER': ['OG'],
                      'THR': ['OG1'],
                      'CYS': ['SG'],
                      'TYR': ['OH'],
                      'LYS': ['NZ'],
                      'GLN': ['HN1', 'HN2'],
                      'ASN': ['HN1', 'HN2'],
                      'ARG': ['HE', 'HN1', 'HN2', 'HN3', 'HN4'],
                      'TRP': ['HNE'],
                      'HIS': ['HND', 'HNE'],
                      'ASP': ['OD1', 'OD2'],
                      'GLU': ['OE1', 'OE2'],
                      'N+ ': ['N'],
                      'C- ': ['O', 'OXT'],
                     } 
    labels['base'] = {'SER': ['OG'],
                      'THR': ['OG1'],
                      'CYS': ['SG'],
                      'TYR': ['OH'],
                      'LYS': ['NZ'],
                      'GLN': ['OE1'],
                      'ASN': ['OD1'],
                      'ARG': ['HE', 'HN1', 'HN2', 'HN3', 'HN4'],
                      'TRP': ['HNE'],
                      'HIS': ['HND', 'HNE'],
                      'ASP': ['OD1', 'OD2'],
                      'GLU': ['OE1', 'OE2'],
                      'N+ ': ['N'],
                      'C- ': ['O', 'OXT'],
                     } 
    return labels


def getDesolvationDefault():
    """ 
    storage of desolvation default parameters
    """
    DesolvationDefault = {}
    DesolvationDefault['propka2']            = {'allowance': 400.00,
                                                'prefactor':  -0.01,
                                                'fudge':       0.00,
                                                'scaled':     False,
                                               }
    DesolvationDefault['ContactModel']       = {'allowance': 400.00,
                                                'prefactor':  -0.01,
                                                'fudge':       0.00,
                                                'scaled':     False,
                                               }
    DesolvationDefault['VolumeModel']        = {'allowance':   0.00,
                                                'prefactor': -13.50,
                                                'fudge':       0.00,
                                                'scaled':     True,
                                               }
    DesolvationDefault['ScaledVolumeModel']  = {'allowance':   0.00,
                                                'prefactor': -20.00,
                                                'scaled':     True,
                                               }

    return DesolvationDefault


def getVanDerWaalsVolumes():
    """ 
    storage of relative Van der Waals volumes for the volume desolvation model
    Note, uses a 'united-atom' description.
    """
    VanDerWaalsVolume = {'C':  1.40,  # 20.58 all 'C' and 'CA' atoms
                         'C4': 2.64,  # 38.79 hydrodphobic carbon atoms + unidentified atoms
                         'N':  1.06,  # 15.60 all nitrogen atoms
                         'O':  1.00,  # 14.71 all oxygen atoms
                         'S':  1.66,  # 24.43 all sulphur atoms
                        }

    return VanDerWaalsVolume


def getBackBoneDefault():
    """ 
    definitions of default back-bone interaction parameters
    """

    back_bone = {"COO": [-0.80, 2.00, 3.00],
                 "CYS": [-2.40, 3.50, 4.50],
                 "TYR": [-1.20, 3.50, 4.50],
                 "HIS": [ 1.20, 2.00, 3.50],
                 "N+ ": [ 1.20, 2.00, 3.50],
                 "LYS": [ 1.20, 2.00, 3.50],
                 "ARG": [ 1.20, 2.00, 3.50]}

    return back_bone


def getSideChainDefault():
    """ 
    definitions of default side-chain interaction parameters
    Note, parameters with 'None' are excluded since they are given by the reverse 
    (e.g. TYR-COO is given by COO-TYR). The sign is also worked out later from pKa_mod
    """
    side_chain = {}
    side_chain["COO"] = {"COO": [-0.80, [ 2.50, 3.50]],
                         "CYS": [-0.80, [ 3.00, 4.00]],
                         "TYR": [-0.80, [ 2.65, 3.65]],
                         "HIS": [-0.80, [ 2.00, 3.00]],
                         "N+ ": [-0.80, [ 2.85, 3.85]],
                         "LYS": [-0.80, [ 2.85, 3.85]],
                         "ARG": [-0.80, [ 1.85, 2.85]],
                         "ROH": [-0.80, [ 2.65, 3.65]],
                         "AMD": [-0.80, [ 2.00, 3.00]],
                         "TRP": [-0.80, [ 2.00, 3.00]]}
    side_chain["CYS"] = {"COO": None,
                         "CYS": [-1.60, [ 3.00, 5.00]],
                         "TYR": [-0.80, [ 3.50, 4.50]],
                         "HIS": [-1.60, [ 3.00, 4.00]],
                         "N+ ": [-2.40, [ 3.00, 4.50]],
                         "LYS": [-1.60, [ 3.00, 4.00]],
                         "ARG": [-1.60, [ 2.50, 4.00]],
                         "ROH": [-1.60, [ 3.50, 4.50]],
                         "AMD": [-1.60, [ 2.50, 3.50]],
                         "TRP": [-1.60, [ 2.50, 3.50]]}
    side_chain["TYR"] = {"COO": None,
                         "CYS": None,
                         "TYR": [ 0.80, [ 3.50, 4.50]],
                         "HIS": [-0.80, [ 2.00, 3.00]],
                         "N+ ": [-1.20, [ 3.00, 4.50]],
                         "LYS": [-0.80, [ 3.00, 4.00]],
                         "ARG": [-0.80, [ 2.50, 4.00]],
                         "ROH": [-0.80, [ 3.50, 4.50]],
                         "AMD": [-0.80, [ 2.50, 3.50]],
                         "TRP": [-0.80, [ 2.50, 3.50]]}
    side_chain["HIS"] = {"COO": None,
                         "CYS": None,
                         "TYR": None,
                         "HIS": [ 0.00, [ 0.00, 0.00]],
                         "N+ ": [ 0.00, [ 0.00, 0.00]],
                         "LYS": [ 0.00, [ 0.00, 0.00]],
                         "ARG": [ 0.00, [ 0.00, 0.00]],
                         "ROH": [ 0.00, [ 0.00, 0.00]],
                         "AMD": [ 0.80, [ 2.00, 3.00]],
                         "TRP": [ 0.00, [ 0.00, 0.00]]}
    side_chain["N+ "] = {"COO": None,
                         "CYS": None,
                         "TYR": None,
                         "HIS": None,
                         "N+ ": [ 0.00, [ 0.00, 0.00]],
                         "LYS": [ 0.00, [ 0.00, 0.00]],
                         "ARG": [ 0.00, [ 0.00, 0.00]],
                         "ROH": [ 0.00, [ 0.00, 0.00]],
                         "AMD": [ 0.00, [ 0.00, 0.00]],
                         "TRP": [ 0.00, [ 0.00, 0.00]]}
    side_chain["LYS"] = {"COO": None,
                         "CYS": None,
                         "TYR": None,
                         "HIS": None,
                         "N+ ": None,
                         "LYS": [ 0.00, [ 0.00, 0.00]],
                         "ARG": [ 0.00, [ 0.00, 0.00]],
                         "ROH": [ 0.00, [ 0.00, 0.00]],
                         "AMD": [ 0.00, [ 0.00, 0.00]],
                         "TRP": [ 0.00, [ 0.00, 0.00]]}
    side_chain["ARG"] = {"COO": None,
                         "CYS": None,
                         "TYR": None,
                         "HIS": None,
                         "N+ ": None,
                         "LYS": None,
                         "ARG": [ 0.00, [ 0.00, 0.00]],
                         "ROH": [ 0.00, [ 0.00, 0.00]],
                         "AMD": [ 0.00, [ 0.00, 0.00]],
                         "TRP": [ 0.00, [ 0.00, 0.00]]}

    return side_chain


def getCoulombDefault():
    """ 
    storage of Coulomb default parameters
    """
    CoulombDefault = {}
    CoulombDefault['Linear']                 =  {'max_dpka':     2.40,
                                                 'cutoff': [4.0, 7.0],
                                                 'diel':         None,
                                                 'scaled_diel':  False,
                                                 'scaled':       True,
                                                }
    CoulombDefault['Coulomb']                =  {'max_dpka':     None,
                                                 'cutoff': [4.0, 7.0],
                                                 'diel':        80.00,
                                                 'scaled_diel':  False,
                                                 'scaled':       True,
                                                }

    return CoulombDefault



def checkPair(resType1, resType2):
    """ 
    definitions of interaction parameters - side-chain hydrogen bonds & Coulomb
    """
    # If exactly one of the residues is ligand, do a simple interaction
    if [resType1, resType2].count("LIG") == 1:
        return "T"
    elif [resType1, resType2].count("LIG") == 2:
        print('Ligand-ligand interaction detected!!')
        return 'F'


    # matrix for propka interactions
    # T = True      - simple interaction
    # I = Iterative - iterative interaction
    # F = False     - no interaction
    #             COO  CYS  TYR  HIS   N+  LYS  ARG
    side_chain =[["I", "I", "T", "I", "T", "T", "T"], # COO
                 ["I", "I", "T", "I", "T", "T", "T"], # CYS
                 ["T", "T", "I", "I", "I", "I", "T"], # TYR
                 ["I", "I", "I", "I", "T", "T", "T"], # HIS
                 ["T", "T", "I", "T", "T", "T", "T"], # N+ 
                 ["T", "T", "I", "T", "T", "T", "T"], # LYS
                 ["T", "T", "T", "T", "T", "T", "I"], # ARG
                 ["T", "T", "T", "F", "F", "F", "F"], # ROH
                 ["T", "T", "T", "T", "F", "F", "F"], # AMD
                 ["T", "T", "T", "F", "F", "F", "F"]] # TRP

    i, j = getIndicies(resType1, resType2)

    if   i == None and j == None:
      do_this = "F"
    #elif indicies[0] < indicies[1]:
    #  do_this = side_chain[j][i]
    else:
      do_this = side_chain[i][j]

    return do_this


def getIndicies(name1, name2):
    """
    return indicies for parameter matricies - 'labels' are just to make it less abstract
    """
    indicies = []
    labels   = []
    names    = []
    names.append(name1)
    names.append(name2)

    
    for name in names:
        if   name == "COO":
            indicies.append(0)
        elif name == "CYS":
            indicies.append(1)
        elif name == "TYR":
            indicies.append(2)
        elif name == "HIS":
            indicies.append(3)
        elif name == "N+ ":
            indicies.append(4)
        elif name == "LYS":
            indicies.append(5)
        elif name == "ARG":
            indicies.append(6)
        elif name == "ROH":
            indicies.append(7)
        elif name == "AMD":
            indicies.append(8)
        elif name == "TRP":
            indicies.append(9)
        else:
            print("cannot find indicies for residue \"%s\" " % (name))
            sys.exit(9)

    if   indicies[0] > 6 and indicies[1] > 6:
      return None, None
    elif indicies[0] < indicies[1]:
      return indicies[1], indicies[0]
    else:
      return indicies[0], indicies[1]



# ----- methods for checking a number of exceptions -----


def checkCooArgException(residue_coo, residue_arg, dpka_max=None, cutoff=None, version=None):
    """
    checking Coo-Arg exception
    """
    # printing out all distances for debugging
    #import debug
    #debug.printCooArgAtomDistances(residue_coo, residue_arg)

    excluded_atoms = []
    exception = True
    value_tot = 0.00
    str = "xxx"
    # needs to be this way since you want to find shortest distance first
    for iter in ["shortest", "runner-up"]:
      distance = 999.
      closest_coo_atom = None
      closest_arg_atom = None
      for atom_coo in residue_coo.makeDeterminantAtomList(residue_arg.resName, version=version):
        for atom_arg in residue_arg.makeDeterminantAtomList(residue_coo.resName, version=version):
          fucked = False
          for atom in excluded_atoms:
            if atom == atom_coo or atom == atom_arg:
              fucked = True
          if fucked == False:
            current_distance = lib.calculateAtomDistance(atom_coo, atom_arg)
            if current_distance < distance:
              closest_arg_atom = atom_arg
              closest_coo_atom = atom_coo
              distance = current_distance
      atom3 = residue_arg.getThirdAtomInAngle(closest_arg_atom)
      distance, f_angle, nada = lib.calculateAngleFactor(closest_coo_atom, closest_arg_atom, atom3)
      value = calculate.SideChainEnergy(distance, dpka_max, cutoff, f_angle)
      value_tot += value
      #print(">>> %s %s: %6.2lf %6.2lf %s" % (closest_coo_atom.name, closest_arg_atom.name, distance, value, iter))
      str += "%6.2lf" % (value)
      #str += "%6.2lf" % (distance)
      excluded_atoms.append(closest_coo_atom)
      excluded_atoms.append(closest_arg_atom)
    #print(str)

    return exception, value


def checkCooArgException_old(residue_coo, residue_arg, version=None):
    """
    checking Coo-Arg exception
    """
    distances = []
    for atom_coo in residue_coo.makeDeterminantAtomList(residue_arg.resName):
      distance = 999.
      for atom_arg in residue_arg.makeDeterminantAtomList(residue_coo.resName):
        distance = min(lib.calculateAtomDistance(atom_coo, atom_arg), distance)
      distances.append(distance)
    exception = True
    for distance in distances:
      if distance > 2.2:
        exception = False

    # Value not set here anymore!!!
    if lib.checkBuried(residue_coo.Nmass, residue_arg.Nmass):
      return exception, None # 2.40 1.20
    else:
      return exception, None


def checkCooCooException(residue1, residue2):
    """
    checking Coo-Coo exception
    """
    exception = False
    if lib.checkBuried(residue1.Nmass, residue2.Nmass):
        exception = True

    # Value not set here anymore!!!
    return exception, 1.60 # 1.60 0.80


def checkCooHisException(residue1, residue2):
    """
    checking Coo-His exception
    """
    exception = False
    if lib.checkBuried(residue1.Nmass, residue2.Nmass):
        exception = True

    return exception, 1.60


def checkCysHisException(residue1, residue2):
    """
    checking Cys-His exception
    """
    exception = False
    if lib.checkBuried(residue1.Nmass, residue2.Nmass):
        exception = True

    return exception, 1.60


def checkCysCysException(residue1, residue2):
    """
    checking Cys-His exception
    """
    exception = False
    if lib.checkBuried(residue1.Nmass, residue2.Nmass):
        exception = True

    return exception, 3.60



# ----- methods for making a version -----


def makeVersion(label="Jan15", verbose=True):
    """
    return a version object for exceptions and weird stuff
    """
    if   label == "default":
      version = Version(verbose=verbose)
    elif label == "Jan01":
      version = Jan01(verbose=verbose)
    elif label == "Jan15":
      version = Jan15(verbose=verbose)
    elif label == "Aug24":
      version = Aug24(verbose=verbose)
    elif label == "Aug30":
      version = Aug30(verbose=verbose)
    elif label == "Aug31":
      version = Aug31(verbose=verbose)
    elif label == "Sep07":
      version = Sep07(verbose=verbose)
    elif label == "Sep08":
      version = Sep08(verbose=verbose)
    elif label == "Oct14":
      version = Oct14(verbose=verbose)
    elif label == "Dec18":
      version = Dec18(verbose=verbose)
    elif label == "Dec19":
      version = Dec19(verbose=verbose)
    elif label == "Oct13":
      version = Oct13(verbose=verbose)
    else:
      print("Could not find version %s" % (label))
      sys.exit(9)

    return version


def setVersion(label="default"):
    """
    return a version object for exceptions and weird stuff
    """
    if   label == "default":
      version = Version()
    elif label == "Jan01":
      version = Jan01()
    elif label == "Jan15":
      version = Jan15()
    elif label == "May13":
      version = May13()
    elif label == "Sep23":
      version = Sep23()
    elif label == "Dec18":
      version = Dec18()
    elif label == "Dec19":
      version = Dec19()
    else:
      print("Could not find version %s" % (label))
      sys.exit(9)

    return version




# === various classes for obtaining version for propka ===

class Version(object):
    """
        Version class - contains rules for calculating pKa values.
    """
    name = None
    buried_cutoff    =   15.50
    desolv_cutoff    =   20.00
    CoulombModel     = "Linear"
    CoulombDefault   = getCoulombDefault()
    coulomb_maxpka   = CoulombDefault[CoulombModel]['max_dpka']
    coulomb_cutoff   = CoulombDefault[CoulombModel]['cutoff']
    coulomb_diel     = CoulombDefault[CoulombModel]['diel']
    coulomb_scaled_diel = CoulombDefault[CoulombModel]['scaled_diel']
    coulomb_scaled   = CoulombDefault[CoulombModel]['scaled']
    DesolvationModel     = "propka2"
    DesolvationDefault   = getDesolvationDefault()
    desolvationPrefactor = DesolvationDefault[DesolvationModel]['prefactor']
    desolvationFudge     = DesolvationDefault[DesolvationModel]['fudge']
    desolvationAllowance = DesolvationDefault[DesolvationModel]['allowance']
    desolvationScaled    = DesolvationDefault[DesolvationModel]['scaled']
    Flocal               = -0.0700
    VanDerWaalsVolume    = getVanDerWaalsVolumes()
    LocalRadius          = {'ASP':4.5, 'GLU':4.5, 'HIS':4.5, 'CYS':3.5, 'TYR':3.5, 'LYS':4.5, 'ARG':5.0, 'C- ':4.5, 'N+ ':4.5}
    doingBackBoneReorganization = False
    scaleUpBuriedSideChain      = 0.00
    BackBoneParameters  = getBackBoneDefault()
    SideChainParameters = getSideChainDefault()
    sidechain_cutoff =  6.0
    Nmin             =  300
    Nmax             =  600
    valueCooArgException = 2.40
    valueCooCooException = 1.60
    exclude_list = ["H2O", "HOH", "SO4"]
    coulomb_list = ["COO", "CYS", "TYR", "HIS", "LYS", "ARG", "N+ "]
    exclude_sidechain_interactions = ["LYS", "ARG", "N+ "]
    # residue interactions that depends on the angle; note, back-bone interactions are explicitly angular dependent
    # also, not implemented for iterative interactions for now
    atomInteractionList = getAtomInteractionList()
    angularDependentSideChainInteractions = ["HIS", "ARG", "AMD", "TRP"] 


    def __init__(self, verbose=True):
        """
        constructer of the Default object.
        """
        str  = "WARNING: you are trying to create a 'default' version object. "
        str += " This object contains the cross section of all versions and is not a complete version itself."
        print(str)
        sys.exit(8)
        self.ions = {}
        self.name = "default"
        if verbose == True:
          print("creating propka version \"%s\"" % (self.name))


    def printVersion(self):
        """
        Print out the properties of this version - more a test and debug feature.
        """
        print("version = %s" % (self.name))
        print("  desolvation   = %s" % (self.DesolvationModel))
        print("    prefactor   = %6.2lf" % (self.desolvationPrefactor))
        print("    allowance   = %6.2lf" % (self.desolvationAllowance))
        print("    scaled      = %s" % (self.desolvationScaled))
        print("  Coulomb       = %s" % (self.CoulombModel))
        print("    prefactor   = %s" % (self.coulomb_maxpka))
        print("    cutoff      = %s" % (self.coulomb_cutoff))
        print("    diel        = %s" % (self.coulomb_diel))
        print("    scaled diel = %s" % (self.coulomb_scaled_diel))
        print("    scaled      = %s" % (self.coulomb_scaled))


    def calculateWeight(self, Nmass):
        """
        calculating the atom-based desolvation weight
        """
        weight = float(Nmass - self.Nmin)/float(self.Nmax - self.Nmin)
        weight = min(1.0, weight)
        weight = max(0.0, weight)

        return weight


    def calculatePairWeight(self, Nmass1, Nmass2):
        """
        calculating the atom-pair based desolvation weight
        """
        Nmass = Nmass1 + Nmass2
        Nmin  = 2*self.Nmin
        Nmax  = 2*self.Nmax
        weight = float(Nmass - Nmin)/float(Nmax - Nmin)
        weight = min(1.0, weight)
        weight = max(0.0, weight)

        return weight


    def getQ(self, resName):
        """
        Returns a residue charge
        """
        Q  =  {'C- ':-1.0,
               'ASP':-1.0,
               'GLU':-1.0,
               'CYS':-1.0,
               'TYR':-1.0,
               'HIS': 1.0,
               'LYS': 1.0,
               'ARG': 1.0,
               'N+ ': 1.0}
        
        if resName in Q:
            return Q[resName]
        else:
            return 0.00


    def resName2Type(self,resName):
        """ 
        definition of which parameter-group each residues belongs to
        """
        
        resType = {'C- ': "COO", 
                   'ASP': "COO", 
                   'GLU': "COO",
                   'HIS': "HIS",
                   'CYS': "CYS",
                   'TYR': "TYR", 
                   'LYS': "LYS",
                   'ARG': "ARG",
                   'N+ ': "N+ ",
                   'SER': "ROH",
                   'THR': "ROH",
                   'ASN': "AMD",
                   'GLN': "AMD",
                   'TRP': "TRP"}


        if resName in resType.keys():
            return resType[resName]
        else:
            return None

    def setCoulomb(self, label, max_dpka=None, cutoff=None, diel=None, scaled_diel=None, scaled=None):
        """
        setting the Coulomb model, and its parameters
        """
        if label not in self.CoulombDefault:
          print("do not accept Coulomb model \"%s\\n" % (label)); sys.exit(9)
        else:
          self.CoulombModel = label

        if max_dpka == None:
          self.coulomb_maxpka = self.CoulombDefault[label]['max_dpka']
        else:
          self.coulomb_maxpka = max_dpka
        if cutoff == None:
          self.coulomb_cutoff = self.CoulombDefault[label]['cutoff']
        else:
          self.coulomb_cutoff = cutoff
        if diel == None:
          self.coulomb_diel = self.CoulombDefault[label]['diel']
        else:
          self.coulomb_diel = diel
        if scaled_diel == None:
          self.coulomb_scaled_diel = self.CoulombDefault[label]['scaled_diel']
        else:
          self.coulomb_scaled_diel = scaled_diel
        if scaled == None:
          self.coulomb_scaled = self.CoulombDefault[label]['scaled']
        else:
          self.coulomb_scaled = scaled



    def setDesolvation(self, label, prefactor=None, allowance=None, fudge=None, scaled=None):
        """
        setting the desolvation model, and its parameters
        """
        desolvation_parameters = getDesolvationParameters()
        if label not in desolvation_parameters:
          print("do not accept solvation model \"%s\\n" % (label))
          sys.exit(9)
        else:
          self.DesolvationModel = label
          default_parameters = desolvation_parameters[label]

        if 'allawance' in default_parameters:
          if allowance == None:
            self.desolvationAllowance = default_parameters['allowance']
          else:
            self.desolvationAllowance = allowance
        if 'prefactor' in default_parameters:
          if prefactor == None:
            self.desolvationPrefactor = default_parameters['prefactor']
          else:
            self.desolvationPrefactor = prefactor
        if 'fudge' in default_parameters:
          if fudge == None:
            self.desolvationFudge     = default_parameters['fudge']
          else:
            self.desolvationFudge     = fudge
        if 'scaled' in default_parameters:
          if scaled    == None:
            self.desolvationScaled    = default_parameters['scaled']
          else:
            self.desolvationScaled    = scaled


    def calculateDesolvation(self, residue, atoms, verbose=True):
        """
        redirecting the desolvation calculation to the right model
        """
        if   self.DesolvationModel == "propka2":
          Nmass, Emass, Nlocl, Elocl = calculate.originalDesolvation(residue, atoms, self, verbose=verbose)
        elif self.DesolvationModel == "ContactModel":
          Nmass, Emass, Nlocl, Elocl = calculate.contactDesolvation(residue, atoms, self, verbose=verbose)
        elif self.DesolvationModel in ["VolumeModel", "ScaledVolumeModel"]:
          Nmass, Emass, Nlocl, Elocl = calculate.radialVolumeDesolvation(residue, atoms, self, verbose=verbose)
        else:
          print("Desolvation \"%s\" is not implemented" % (self.DesolvationModel))
          sys.exit(8)

        return Nmass, Emass, Nlocl, Elocl


    def calculateBackBoneReorganization(self, protein, verbose=True):
        """
        Testing new term, reorganization of back-bone CO groups
        """
        if self.doingBackBoneReorganization == True:
          calculate.BackBoneReorganization(protein, verbose=verbose)
        else:
          """ do nothing """


    def getSideChainParameters(self, residue1, residue2):
        """ 
        returns 'correct' side-chain interaction parameters, order key1 & key2
        """
        resTypes = ["COO", "CYS", "TYR", "HIS", "N+ ", "LYS", "ARG"]
        try:
          index1 = resTypes.index(residue1.resType)
        except ValueError:
          index1 = len(resTypes)
        try:
          index2 = resTypes.index(residue2.resType)
        except ValueError:
          index2 = len(resTypes)
        if index2 < index1:
          key1 = residue2.resType
          key2 = residue1.resType
        else:
          key1 = residue1.resType
          key2 = residue2.resType

        return self.SideChainParameters[key1][key2]


    def checkExceptions(self, residue1, residue2):
        """
        checks for exceptions for this version - using defaults
        """
        exception = False
        value = None
        resType1 = residue1.resType
        resType2 = residue2.resType
        if   (resType1 == "COO" and resType2 == "ARG"):
          exception, value = checkCooArgException_old(residue1, residue2, version=self)
          value = self.valueCooArgException
        elif (resType1 == "ARG" and resType2 == "COO"):
          exception, value = checkCooArgException_old(residue2, residue1, version=self)
          value = self.valueCooArgException
        elif (resType1 == "COO" and resType2 == "COO"):
          exception, value = checkCooCooException(residue1, residue2)
          value = self.valueCooCooException
        elif (resType1 == "CYS" and resType2 == "CYS"):
          exception, value = checkCysCysException(residue1, residue2)
        elif (resType1 == "COO" and resType2 == "HIS") or \
             (resType1 == "HIS" and resType2 == "COO"):
          exception, value = checkCooHisException(residue1, residue2)
        elif (resType1 == "CYS" and resType2 == "HIS") or \
             (resType1 == "HIS" and resType2 == "CYS"):
          exception, value = checkCysHisException(residue1, residue2)
        else:
          """ do nothing, no exception for this pair """

        return exception, value


    def checkCoulombPair(self, residue1, residue2, distance):
        """
        Checks if this Coulomb interaction should be done - a propka2.0 hack
        """
        Npair = residue1.Nmass + residue2.Nmass
        do_coulomb = True

        if residue1.resType in self.coulomb_list and residue2.resType in self.coulomb_list:
          # distance criteria
          if distance > self.coulomb_cutoff[1]:
            do_coulomb = False
          # famous COO-TYR exception
          if (residue1.resType == "COO" and residue2.resType == "TYR") or \
             (residue2.resType == "COO" and residue1.resType == "TYR"):
              """ do nothing """
          elif Npair < self.Nmin:
              do_coulomb = False
        else:
          do_coulomb = False

        #print "%s - %s Npair=%4d %s" % (residue1.label, residue2.label, Npair, do_coulomb)
        return do_coulomb



    def setCoulombCutOff(self, cutoff):
        """
        sets the cutoff for calculating Coulomb interactions
        """
        self.coulomb_cutoff = cutoff


    def calculateSideChainEnergy(self, distance, dpka_max, cutoff, weight, f_angle):
        """
        redirects to get the correct side-chain interaction
        """
        if True:
          prefactor = (1.0 + weight*self.scaleUpBuriedSideChain)
        else:
          prefactor = 1.00
        return calculate.SideChainEnergy(distance, prefactor*dpka_max, cutoff, f_angle)


    def calculateCoulombEnergy(self, distance, weight):
        """
        redirects to get the correct Coulomb interaction - linear for default
        """
        if   self.CoulombModel == "Linear":
          return calculate.linearCoulombEnergy(distance, weight, self, verbose=False)
        elif self.CoulombModel == "Coulomb":
          return calculate.CoulombEnergy(distance, weight, self, verbose=False)
        else:
          print("Coulomb \"%s\" is not implemented" % (self.CoulombModel))
          sys.exit(8)
          


    def calculateCoulombWeight(self, Nmass1, Nmass2):
        """
        calculates the weight for the Coulomb interaction - used for version "Dec18" & "Sep23"
        """
        N_pair = Nmass1 + Nmass2
        Nmin = 2*self.Nmin
        Nmax = 2*self.Nmax
        weight = float(N_pair - Nmin)/float(Nmax - Nmin)
        weight = min(1.0, weight)
        weight = max(0.0, weight)

        return weight


#   --- specific versions with different behaviour or initialization ---

class Jan01(Version):
    """ 
    This is a test to set up rules for different propka Jan15 version
    """

    def __init__(self, verbose=True):
        """
        Rules of action for version Jan01
        """
        Version.__init__(self, verbose=False)
    
        self.name = "Jan01"
        if verbose == True:
          print("creating propka version \"%s\"" % (self.name))

        self.Nmin             =  300
        self.Nmax             =  600
        self.buried_cutoff    =   15.50
        self.desolv_cutoff    =   15.50
        self.setDesolvation("propka2")
        self.setCoulomb("Linear", cutoff=[4.0, 7.0], diel=None, scaled_diel=None, scaled=False)
        self.doingBackBoneReorganization = False


    def checkCoulombPair(self, residue1, residue2, distance):
        """
        Checks if this Coulomb interaction should be done - a propka2.0 hack
        """
        do_coulomb = True

        # check all Coulomb criteria!
        if residue1.resType in self.coulomb_list and residue2.resType in self.coulomb_list:
          # distance criteria
          if distance > self.coulomb_cutoff[1]:
            do_coulomb = False
          # famous COO-TYR exception
          if (residue1.resType == "COO" and residue2.resType == "TYR") or \
             (residue2.resType == "COO" and residue1.resType == "TYR"):
              """ do nothing """
          elif lib.checkBuried(residue1.Nmass, residue2.Nmass) == False:
              do_coulomb = False
        else:
          do_coulomb = False

        return do_coulomb


    def calculateCoulombWeight(self, Nmass1, Nmass2):
        """
        calculates the weight for the Coulomb interaction - used for version "Dec18" & "Sep23"
        """
        if lib.checkBuried(Nmass1, Nmass2) == False:
          return 0.0
        else:
          return 1.0


class Oct13(Version):
    """ 
    To test ligand integration
    """

    def __init__(self, verbose=True):
        """
        Rules of action for version Oct13, based on Sep07
        """
        Version.__init__(self, verbose=False)

        self.name = "Oct13"
        if verbose == True:
          print("creating propka version \"%s\"" % (self.name))

        self.Nmin             =  280
        self.Nmax             =  560
        self.buried_cutoff    =   15.00
        self.desolv_cutoff    =   20.00
        self.setDesolvation("VolumeModel", prefactor=-13.12, fudge=0.40,  allowance=0.0, scaled=True)
        self.setCoulomb("Coulomb", cutoff=[4.0, 10.0], diel=80.0, scaled_diel=True, scaled=False)
        self.doingBackBoneReorganization = True
        
        coulomb_list = ["COO", "CYS", "TYR", "HIS", "LYS", "ARG", "N+ "]

        # ligand parameters
        self.read_ion_parameters()

        mono_atomic_ions = ['LIG']

        self.pl_coulomb_list = coulomb_list+mono_atomic_ions



    def read_ion_parameters(self):
        """ Reads in ions.list """

        self.ions = {}
        self.ions_long_names = {}
        file = os.path.join(os.path.dirname(__file__), 'ions.list')
        if not os.path.isfile(file):
            print('Error: Could not find ion parameter file:',file)
            exit(9)

        lines = open(file,'r').readlines()
        for line in lines:
            words = line[:line.find('#')].split()
            if len(words) == 2:
                self.ions[words[0]] = int(words[1])
                self.ions_long_names[words[0]] = line[line.find('#')+1:-1]
        return


    def getQ(self, resName):
        """
        Returns a residue charge
        """
        Q  =  {'C- ':-1.0,
               'ASP':-1.0,
               'GLU':-1.0,
               'CYS':-1.0,
               'TYR':-1.0,
               'HIS': 1.0,
               'LYS': 1.0,
               'ARG': 1.0,
               'N+ ': 1.0}
        
        if resName in Q:
            return Q[resName]
        elif resName in self.ions.keys():
            return self.ions[resName]
        else:
            return 0.00



    def resName2Type(self, resName):
        """ 
        Expands the standard resName2Type to make sure that ion names are included 
        """

        resType = {'C- ': "COO", 
                   'ASP': "COO", 
                   'GLU': "COO",
                   'HIS': "HIS",
                   'CYS': "CYS",
                   'TYR': "TYR", 
                   'LYS': "LYS",
                   'ARG': "ARG",
                   'N+ ': "N+ ",
                   'SER': "ROH",
                   'THR': "ROH",
                   'ASN': "AMD",
                   'GLN': "AMD",
                   'TRP': "TRP"}

        if resName in resType.keys():
            return resType[resName]
        elif resName in self.ions.keys():
            return 'LIG'
        else:
            return None


    def checkCoulombPair(self, residue1, residue2, distance):
        """
        Checks if this Coulomb interaction should be done - a propka2.0 hack
        """
        do_coulomb = True

        # check all Coulomb criteria!
        if residue1.resType in self.pl_coulomb_list and residue2.resType in self.pl_coulomb_list:
          # distance criteria
          if distance > self.coulomb_cutoff[1]:
            do_coulomb = False
          # famous COO-TYR exception
          if (residue1.resType == "COO" and residue2.resType == "TYR") or \
             (residue2.resType == "COO" and residue1.resType == "TYR"):
              """ do nothing """
          elif lib.checkBuried(residue1.Nmass, residue2.Nmass) == False:
              do_coulomb = False
        else:
          do_coulomb = False

        return do_coulomb





class Jan15(Version):
    """ 
    This is a test to set up rules for different propka Jan15 version
    """

    def __init__(self, verbose=True):
        """
        Rules of action for version Jan15
        """
        Version.__init__(self, verbose=False)

        self.name = "Jan15"
        if verbose == True:
          print("creating propka version \"%s\"" % (self.name))

        self.Nmin             =  300
        self.Nmax             =  600
        self.buried_cutoff    =   15.50
        self.desolv_cutoff    =   15.50
        self.setDesolvation("ContactModel", prefactor=-0.01, allowance=400.0, scaled=False)
        self.setCoulomb("Linear", cutoff=[4.0, 7.0], diel=None, scaled_diel=None, scaled=False)
        self.doingBackBoneReorganization = False


    def checkCoulombPair(self, residue1, residue2, distance):
        """
        Checks if this Coulomb interaction should be done - a propka2.0 hack
        """
        do_coulomb = True

        # check all Coulomb criteria!
        if residue1.resType in self.coulomb_list and residue2.resType in self.coulomb_list:
          # distance criteria
          if distance > self.coulomb_cutoff[1]:
            do_coulomb = False
          # famous COO-TYR exception
          if (residue1.resType == "COO" and residue2.resType == "TYR") or \
             (residue2.resType == "COO" and residue1.resType == "TYR"):
              """ do nothing """
          elif lib.checkBuried(residue1.Nmass, residue2.Nmass) == False:
              do_coulomb = False
        else:
          do_coulomb = False

        return do_coulomb


    def calculateCoulombWeight(self, Nmass1, Nmass2):
        """
        calculates the weight for the Coulomb interaction - used for version "Dec18" & "Sep23"
        """
        if lib.checkBuried(Nmass1, Nmass2) == False:
          return 0.0
        else:
          return 1.0



class May13(Version):
    """ 
    This is a test to set up rules for different propka May13 version
    """

    def __init__(self):
        """
        Rules of action for version May13
        """
        Version.__init__(self, verbose=False)

        self.name = "May13"
        self.coulomb_cutoff =   7.00
        print("creating propka version \"%s\"" % (self.name))


    def checkExceptions(self, residue1, residue2):
        """
        overwrites 'exceptions' from the default
        """
        exception = False
        value     = 0.00
        return exception, value


class Dec18(Version):
    """ 
    This is a test to set up rules for different propka versions
    """

    def __init__(self, verbose=True):
        """
        Rules of action for version Dec18
        """
        Version.__init__(self, verbose=False)

        self.name = "Dec18"
        if verbose == True:
          print("creating propka version \"%s\"" % (self.name))

        self.Nmin             =  300
        self.Nmax             =  600
        self.buried_cutoff    =   15.50
        self.desolv_cutoff    =   15.50
        self.setDesolvation("propka2", prefactor=-0.01, allowance=400.0, scaled=False)
        self.setCoulomb("Linear", cutoff=[4.0, 7.0], diel=None, scaled_diel=None, scaled=True)
        self.doingBackBoneReorganization = False


class Dec19(Version):
    """ 
    This is a test to set up rules for different propka versions
    """

    def __init__(self, verbose=True):
        """
        Rules of action for version Dec19
        """
        Version.__init__(self, verbose=False)
        
        self.name = "Dec19"
        if verbose == True:
          print("creating propka version \"%s\"" % (self.name))

        self.Nmin             =  300
        self.Nmax             =  600
        self.buried_cutoff    =   15.50
        self.desolv_cutoff    =   15.50
        self.setDesolvation("ContactModel", prefactor=-0.01, allowance=400.0, scaled=False)
        self.setCoulomb("Linear", cutoff=[4.0, 7.0], diel=None, scaled_diel=None, scaled=True)
        self.doingBackBoneReorganization = False


class Aug24(Version):
    """ 
    This is a test to set up rules for different propka versions
    """

    def __init__(self, verbose=True):
        """
        Rules of action for version Aug24
        """
        Version.__init__(self, verbose=False)

        self.name = "Aug24"
        if verbose == True:
          print("creating propka version \"%s\"" % (self.name))

        self.Nmin             =  280
        self.Nmax             =  560
        self.buried_cutoff    =   15.00
        self.desolv_cutoff    =   20.00
        self.setDesolvation("VolumeModel", prefactor=-6.75, fudge=0.00, allowance=0.0, scaled=False)
        self.setCoulomb("Linear", cutoff=[4.0, 7.0], diel=80.0, scaled_diel=None, scaled=True)
        self.doingBackBoneReorganization = True


class Aug30(Version):
    """ 
    This is a test to set up rules for different propka versions
    """

    def __init__(self, verbose=True):
        """
        Rules of action for version Aug30
        """
        Version.__init__(self, verbose=False)

        self.name = "Aug30"
        if verbose == True:
          print("creating propka version \"%s\"" % (self.name))

        self.Nmin             =  280
        self.Nmax             =  560
        self.buried_cutoff    =   15.00
        self.desolv_cutoff    =   20.00
        self.setDesolvation("VolumeModel", prefactor=-14.75, fudge=0.00, allowance=0.0, scaled=True)
        self.setCoulomb("Linear", cutoff=[4.0, 7.0], diel=80.0, scaled_diel=None, scaled=True)
        self.doingBackBoneReorganization = True


class Aug31(Version):
    """ 
    This is a test to set up rules for different propka versions
    """

    def __init__(self, verbose=True):
        """
        Rules of action for version Aug31
        """
        Version.__init__(self, verbose=False)

        self.name = "Aug31"
        if verbose == True:
          print("creating propka version \"%s\"" % (self.name))

        self.Nmin             =  280
        self.Nmax             =  560
        self.buried_cutoff    =   15.00
        self.desolv_cutoff    =   20.00
        self.setDesolvation("VolumeModel", prefactor=-13.50, fudge=0.25,  allowance=0.0, scaled=True)
        self.setCoulomb("Linear", cutoff=[4.0, 7.0], diel=80.0, scaled_diel=None, scaled=True)
        self.doingBackBoneReorganization = True


class Sep07(Version):
    """ 
    This is a test to set up rules for different propka versions
    """

    def __init__(self, verbose=True):
        """
        Rules of action for version Sep07
        """
        Version.__init__(self, verbose=False)

        self.name = "Sep07"
        if verbose == True:
          print("creating propka version \"%s\"" % (self.name))

        self.Nmin             =  280
        self.Nmax             =  560
        self.buried_cutoff    =   15.00
        self.desolv_cutoff    =   20.00
        self.setDesolvation("VolumeModel", prefactor=-13.12, fudge=0.40,  allowance=0.0, scaled=True)
        self.setCoulomb("Coulomb", cutoff=[4.0, 10.0], diel=80.0, scaled_diel=True, scaled=False)
        self.doingBackBoneReorganization = True


class Sep08(Version):
    """ 
    This is a test to set up rules for different propka versions
    """

    def __init__(self, verbose=True):
        """
        Rules of action for version Sep08
        """
        Version.__init__(self, verbose=False)

        self.name = "Sep08"
        if verbose == True:
          print("creating propka version \"%s\"" % (self.name))

        self.Nmin             =  280
        self.Nmax             =  560
        self.buried_cutoff    =   15.00
        self.desolv_cutoff    =   20.00
        self.setDesolvation("VolumeModel", prefactor=-13.75, fudge=0.20,  allowance=0.0, scaled=True)
        self.setCoulomb("Coulomb", cutoff=[4.0, 10.0], diel=20.0, scaled_diel=False, scaled=True)
        self.doingBackBoneReorganization = True


class Oct14(Version):
    """ 
    This is a test to set up rules for different propka versions
    """

    def __init__(self, verbose=True):
        """
        Rules of action for version Oct14
        """
        self.name = "Oct14"
        if verbose == True:
          print("creating propka version \"%s\"" % (self.name))

        self.Nmin             =  280
        self.Nmax             =  560
        self.buried_cutoff    =   15.00
        self.desolv_cutoff    =   20.00
        self.setDesolvation("VolumeModel", prefactor=-13.25, fudge=0.325,  allowance=0.0, scaled=True)
        self.setCoulomb("Coulomb", cutoff=[4.0, 10.0], diel=80.0, scaled_diel=True, scaled=False)
        self.doingBackBoneReorganization = True


    def checkExceptions(self, residue1, residue2):
        """
        checks for exceptions for this version - using defaults
        """
        exception = False
        value = None
        resType1 = residue1.resType
        resType2 = residue2.resType
        if   (resType1 == "COO" and resType2 == "ARG"):
          dpka_max, cutoff = self.getSideChainParameters(residue1, residue2)
          exception, value = checkCooArgException(residue1, residue2, dpka_max=dpka_max, cutoff=cutoff, version=self)
        elif (resType1 == "ARG" and resType2 == "COO"):
          dpka_max, cutoff = self.getSideChainParameters(residue1, residue2)
          exception, value = checkCooArgException(residue2, residue1, dpka_max=dpka_max, cutoff=cutoff, version=self)
        elif (resType1 == "COO" and resType2 == "COO"):
          exception, value = checkCooCooException(residue1, residue2)
          value = self.valueCooCooException
        elif (resType1 == "CYS" and resType2 == "CYS"):
          exception, value = checkCysCysException(residue1, residue2)
        elif (resType1 == "COO" and resType2 == "HIS") or \
             (resType1 == "HIS" and resType2 == "COO"):
          exception, value = checkCooHisException(residue1, residue2)
        elif (resType1 == "CYS" and resType2 == "HIS") or \
             (resType1 == "HIS" and resType2 == "CYS"):
          exception, value = checkCysHisException(residue1, residue2)
        else:
          """ do nothing, no exception for this pair """

        return exception, value


class Jan01_old:
    """ 
    This is a test to set up rules for different propka versions
    """

    def __init__(self):
        """
        Rules of action for version Jan01
        """
        self.name = "Jan01"
        self.coulomb_cutoff   = 7.00
        self.sidechain_cutoff = 7.00
        self.valueCooArgException = 2.40
        self.valueCooCooException = 1.60
        self.exclude_sidechain_interactions = ["LYS", "ARG", "N+ "]


    def printVersion(self):
        """
        Checks if this Coulomb interaction should be done - a propka2.0 hack
        """
        print("version = %s" % (self.name))


    def checkCoulombPair(self, residue1, residue2, distance):
        """
        Checks if this Coulomb interaction should be done - a propka2.0 hack
        """
        name1 = lib.groupName(residue1.resName)
        name2 = lib.groupName(residue2.resName)
        coulomb_list = ["COO", "ASP", "GLU", "CYS", "TYR", "HIS", "LYS", "ARG", "N+ ", "C- "]
        DIS2 = 7.0
        do_coulomb = True

        # check all Coulomb criteria - Needs to be changed for version Dec18
        if name1 in coulomb_list and name2 in coulomb_list:
          # distance criteria
          if distance > DIS2:
            do_coulomb = False
          # famous COO-TYR exception
          if (name1 == "COO" and name2 == "TYR") or \
             (name2 == "COO" and name1 == "TYR"):
              """ do nothing """
          elif lib.checkBuried(residue1.Nmass, residue2.Nmass) == False:
              do_coulomb = False
        else:
          do_coulomb = False

        return do_coulomb


    def setCooArgException(self, newValue):
        """
        Changing the value for the Coo - Arg exception
        """
        self.valueCooArgException = newValue


    def setCooCooException(self, newValue):
        """
        Changing the value for the Coo - Arg exception
        """
        self.valueCooCooException = newValue


    def checkExceptions(self, residue1, residue2):
        """
        checks for exceptions for this version - using defaults
        """
        exception = False
        value = None
        name1 = lib.groupName(residue1.resName)
        name2 = lib.groupName(residue2.resName)
        if   (name1 == "COO" and name2 == "ARG"):
          exception, value = checkCooArgException(residue1, residue2)
          value = self.valueCooArgException
        elif (name1 == "ARG" and name2 == "COO"):
          exception, value = checkCooArgException(residue2, residue1)
          value = self.valueCooArgException
        elif (name1 == "COO" and name2 == "COO"):
          exception, value = checkCooCooException(residue1, residue2)
          value = self.valueCooCooException
        elif (name1 == "CYS" and name2 == "CYS"):
          exception, value = checkCysCysException(residue1, residue2)
        elif (name1 == "COO" and name2 == "HIS") or \
             (name1 == "HIS" and name2 == "COO"):
          exception, value = checkCooHisException(residue1, residue2)
        elif (name1 == "CYS" and name2 == "HIS") or \
             (name1 == "HIS" and name2 == "CYS"):
          exception, value = checkCysHisException(residue1, residue2)
        else:
          """ do nothing, no exception for this pair """

        return exception, value


    def calculateCoulombWeight(self, Nmass1, Nmass2):
        """
        calculates the weight for the Coulomb interaction - propka version "Dec18"
        """
        return 1.0



class Dec18_old():
    """ 
    This is a test to set up rules for different propka versions
    """

    def __init__(self):
        """
        Rules for action for version Dec18
        """
        self.name = "Dec18"
        self.Nmax = 1200
        self.Nmin =  600
        self.coulomb_cutoff   = 7.00
        self.sidechain_cutoff = 7.00
        self.valueCooArgException = 2.40
        self.valueCooCooException = 1.60
        self.coulomb_list = ["COO", "CYS", "TYR", "HIS", "LYS", "ARG", "N+ "]
        self.exclude_sidechain_interactions = ["LYS", "ARG", "N+ "]


    def printVersion(self):
        """
        Print out the properties of this version - more a test and debug feature.
        """
        print("version = %s" % (self.name))


    def checkCoulombPair(self, residue1, residue2, distance):
        """
        Checks if this Coulomb interaction should be done - a propka2.0 hack
        """
        name1 = lib.groupName(residue1.resName)
        name2 = lib.groupName(residue2.resName)
        Npair = residue1.Nmass + residue2.Nmass
        do_coulomb = True

        if name1 in self.coulomb_list and name2 in self.coulomb_list:
          # distance criteria
          if distance > self.coulomb_cutoff:
            do_coulomb = False
          # famous COO-TYR exception
          if (name1 == "COO" and name2 == "TYR") or \
             (name2 == "COO" and name1 == "TYR"):
              """ do nothing """
          elif Npair < self.Nmin:
              do_coulomb = False
        else:
          do_coulomb = False

        #print "%s - %s Npair=%4d %s" % (residue1.label, residue2.label, Npair, do_coulomb)
        return do_coulomb


    def setCooArgException(self, newValue):
        """
        Changing the value for the Coo - Arg exception
        """
        self.valueCooArgException = newValue


    def setCooCooException(self, newValue):
        """
        Changing the value for the Coo - Coo exception
        """
        self.valueCooCooException = newValue


    def checkExceptions(self, residue1, residue2):
        """
        checks for exceptions for this version - using defaults
        """
        exception = False
        value = None
        name1 = lib.groupName(residue1.resName)
        name2 = lib.groupName(residue2.resName)
        if   (name1 == "COO" and name2 == "ARG"):
          exception, value = checkCooArgException(residue1, residue2)
          value = self.valueCooArgException
        elif (name1 == "ARG" and name2 == "COO"):
          exception, value = checkCooArgException(residue2, residue1)
          value = self.valueCooArgException
        elif (name1 == "COO" and name2 == "COO"):
          exception, value = checkCooCooException(residue1, residue2)
          value = self.valueCooCooException
        elif (name1 == "CYS" and name2 == "CYS"):
          exception, value = checkCysCysException(residue1, residue2)
        elif (name1 == "COO" and name2 == "HIS") or  \
             (name1 == "HIS" and name2 == "COO"):
          exception, value = checkCooHisException(residue1, residue2)
        elif (name1 == "CYS" and name2 == "HIS") or  \
             (name1 == "HIS" and name2 == "CYS"):
          exception, value = checkCysHisException(residue1, residue2)
        else:
          """ do nothing, no exception for this pair """

        return exception, value


    def getCoulombEnergy(self, distance):
        """
        redirects to get the correct Coulomb interaction - linear for default
        """
        return calculate.linearCoulombEnergy(distance)


    def calculateCoulombWeight(self, Nmass1, Nmass2):
        """
        calculates the weight for the Coulomb interaction - propka version "Dec18"
        """
        N_pair = Nmass1 + Nmass2
        #Npair = residue1.Nmass + residue2.Nmass
        weight = float(N_pair - self.Nmin)/float(self.Nmax - self.Nmin)
        weight = min(1.0, weight)
        weight = max(0.0, weight)

        return weight


