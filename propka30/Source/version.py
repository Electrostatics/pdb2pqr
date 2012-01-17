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
import calculator as calculate


def getAtomInteractionList():
    """ 
    storing label of interacting atoms
    """
    labels = {}
    # the atoms with which this residue interacts with acids
    labels['acid'] = {'SER': ['OG'],
                      'THR': ['OG1'],
                      'CYS': ['SG'],
                      'TYR': ['OH'],
                      'LYS': ['NZ'],
                      'GLN': ['1HE2', '2HE2'],
                      'ASN': ['1HD2', '2HD2'],
                      'ARG': ['HE', '1HH1', '2HH1', '1HH2', '2HH2'],
                      'TRP': ['HE1'],
                      'HIS': ['HD1', 'HE2'],
                      'ASP': ['OD1', 'OD2'],
                      'GLU': ['OE1', 'OE2'],
                      'N+ ': ['N'],
                      'C- ': ['O', 'OXT'],
                     } 
    # the atoms with which this residue interacts with bases
    labels['base'] = {'SER': ['OG'],
                      'THR': ['OG1'],
                      'CYS': ['SG'],
                      'TYR': ['OH'],
                      'LYS': ['NZ'],
                      'GLN': ['OE1'],
                      'ASN': ['OD1'],
                      'ARG': ['NE', 'NH1', 'NH2'],
                      'TRP': ['NE1'],
                      'HIS': ['ND1', 'NE2'],
                      'ASP': ['OD1', 'OD2'],
                      'GLU': ['OE1', 'OE2'],
                      'N+ ': ['N'],
                      'C- ': ['O', 'OXT'],
                     } 
    return labels


def getSideChainInteractionMatrix(side_chain=None):
    """ 
      setting up rules for iterative/non-iterative/non-interacting; 
      Note that only the LOWER part of the matrix is used!
     
      'N'   non-iterative interaction
      'I'   iterative interaction
      '-'   no interaction
    """

    # WARNING, this matrix has been moved to parameters_xxx since
    #          it is different for various versions!!!
    #                     COO  CYS  TYR  HIS   N+  LYS  ARG
    #side_chain = {'COO': ["I", "I", "N", "N", "N", "N", "N"],
    #              'CYS': ["I", "I", "N", "I", "N", "N", "N"],
    #              'TYR': ["N", "N", "I", "I", "I", "I", "N"],
    #              'HIS': ["I", "I", "I", "I", "N", "N", "N"],
    #              'N+ ': ["N", "N", "I", "N", "I", "N", "N"],
    #              'LYS': ["N", "N", "I", "N", "N", "I", "N"],
    #              'ARG': ["N", "N", "N", "N", "N", "N", "I"],
    #              'ROH': ["N", "N", "N", "-", "-", "-", "-"],
    #              'AMD': ["N", "N", "N", "N", "-", "-", "-"],
    #              'TRP': ["N", "N", "N", "-", "-", "-", "-"]}
    keys = ["COO", "CYS", "TYR", "HIS", "N+ ", "LYS", "ARG", "ROH", "AMD", "TRP"]

    # generating the true 'interaction' matrix
    interaction = {}
    for key1 in keys:
      i=0
      interaction[key1] = {}
      for key2 in keys:
        if   i >= len(side_chain[key1]):
          do_pair = False; iterative = None
        elif side_chain[key1][i] == "I":
          do_pair = True;  iterative = True
        elif side_chain[key1][i] == "N":
          do_pair = True;  iterative = False
        else:
          do_pair = False; iterative = None
          
        interaction[key1][key2] = [do_pair, iterative]
        interaction[key2][key1] = [do_pair, iterative]

        if key1 == key2:
          break
        else:
          i += 1

    return interaction




# ----- methods for checking a number of propka exceptions -----


def checkCooArgException(residue_coo, residue_arg, version=None):
    """
    checking Coo-Arg exception: uses the two shortes unique distances (involving 2+2 atoms)
    """
    # printing out all distances for debugging
    #import debug
    #debug.printCooArgAtomDistances(residue_coo, residue_arg)

    str = "xxx"
    excluded_atoms = []
    exception = True
    value_tot = 0.00
    dpka_max, cutoff = version.SideChainParameters[residue_coo.resType][residue_arg.resType]
    # needs to be this way since you want to find shortest distance first
    # print("--- exception for %s %s ---" % (residue_coo.label, residue_arg.label))
    for iter in ["shortest", "runner-up"]:
      distance = 999.
      closest_coo_atom = None
      closest_arg_atom = None
      for atom_coo in residue_coo.makeDeterminantAtomList(residue_arg.resName, version=version):
        for atom_arg in residue_arg.makeDeterminantAtomList(residue_coo.resName, version=version):
          fucked = False
          if atom_coo in excluded_atoms or atom_arg in excluded_atoms:
            fucked = True
          if fucked == False:
            current_distance = calculate.InterAtomDistance(atom_coo, atom_arg)
            #print("- %-4s %-4s: %6.2lf%6.2lf %s" % (atom_coo.name, atom_arg.name, current_distance, distance, iter))
            if current_distance < distance:
              closest_arg_atom = atom_arg
              closest_coo_atom = atom_coo
              distance = current_distance
      atom3 = residue_arg.getThirdAtomInAngle(closest_arg_atom)
      if   residue_arg.resType in version.angularDependentSideChainInteractions:
        distance, f_angle, nada = calculate.AngleFactorX(closest_coo_atom, closest_arg_atom, atom3)
      else:
        f_angle = 1.00
      value = calculate.HydrogenBondEnergy(distance, dpka_max, cutoff, f_angle)
      value_tot += value
      #print(">>> %s %s: %6.2lf %6.2lf %6.2lf %s" % (closest_coo_atom.name, closest_arg_atom.name, distance, f_angle, value, iter))
      str += "%6.2lf" % (value)
      #str += "%6.2lf" % (distance)
      excluded_atoms.append(closest_coo_atom)
      excluded_atoms.append(closest_arg_atom)
    #print(str)

    return exception, value_tot


def checkCooArgException_old(residue_coo, residue_arg, version=None):
    """
    checking Coo-Arg exception
    """
    distances = []
    for atom_coo in residue_coo.makeDeterminantAtomList(residue_arg.resName, version=version):
      distance = 999.
      for atom_arg in residue_arg.makeDeterminantAtomList(residue_coo.resName, version=version):
        distance = min(calculate.InterAtomDistance(atom_coo, atom_arg), distance)
      distances.append(distance)
    exception = True
    for distance in distances:
      if distance > 2.2:
        exception = False

    # Value not set here anymore!!!
    if lib.checkBuried(residue_coo.Nmass, residue_arg.Nmass):
      return exception, 2.40
    else:
      return exception, 2.40


def checkCooCooException(residue1, residue2, version=None):
    """
    checking Coo-Coo hydrogen-bond exception
    """
    exception = True
    distance = 999.
    for atom1 in residue1.makeDeterminantAtomList(residue2.resName, version=version):
      for atom2 in residue2.makeDeterminantAtomList(residue1.resName, version=version):
        distance = min(calculate.InterAtomDistance(atom1, atom2), distance)
    dpka_max, cutoff = version.SideChainParameters[residue1.resType][residue2.resType]
    f_angle = 1.00
    value = calculate.HydrogenBondEnergy(distance, dpka_max, cutoff, f_angle)
    weight = version.calculatePairWeight(residue1.Nmass, residue2.Nmass)
    value = value * (1.0 + weight)

    return exception, value


def checkCooCooException_old(residue1, residue2, version=None):
    """
    checking Coo-Coo exception
    """
    exception = False
    if lib.checkBuried(residue1.Nmass, residue2.Nmass):
        exception = True

    return exception, 1.60 # 1.60 0.80


def checkCooHisException(residue1, residue2, version=None):
    """
    checking Coo-His exception
    """
    exception = False
    if lib.checkBuried(residue1.Nmass, residue2.Nmass):
        exception = True

    return exception, 1.60


def checkCysHisException(residue1, residue2, version=None):
    """
    checking Cys-His exception
    """
    exception = False
    if lib.checkBuried(residue1.Nmass, residue2.Nmass):
        exception = True

    return exception, 1.60


def checkCysCysException(residue1, residue2, version=None):
    """
    checking Cys-Cys exception
    """
    exception = False
    if lib.checkBuried(residue1.Nmass, residue2.Nmass):
        exception = True

    return exception, 3.60



# ----- methods for making a version -----


def makeVersion(label="Nov30", options=None):
    """
    return a version object for exceptions and weird stuff
    """
    if   label == "Jan15":
      version = Jan15(options=options)
    elif label == "Jan01":
      version = Jan01(options=options)
    elif label == "Aug24":
      version = Aug24(options=options)
    elif label == "Aug30":
      version = Aug30(options=options)
    elif label == "Aug31":
      version = Aug31(options=options)
    elif label == "Sep05":
      version = Sep05(options=options)
    elif label == "Sep06":
      version = Sep06(options=options)
    elif label == "Sep07":
      version = Sep07(options=options)
    elif label == "Sep08":
      version = Sep08(options=options)
    elif label == "Dec18":
      version = Dec18(options=options)
    elif label == "Dec19":
      version = Dec19(options=options)
    elif label == "Oct13":
      version = Oct13(options=options)
    elif label == "Oct14":
      version = Oct14(options=options)
    elif label == "Nov28":
      version = Nov28(options=options)
    elif label == "Nov29":
      version = Nov29(options=options)
    elif label == "Nov30":
      version = Nov30(options=options)
    else:
      print("version \"%s\" not defined in makeVersion()" % (label))
      sys.exit(9)

    print("created version \"%s\"" % (version.name))

    return version




# --- definition of various classes for obtaining propka version ---


class Version(object):
    """
        Version class - contains rules for calculating pKa values.
    """
    name             =   None
    buried_cutoff    =   15.50
    desolv_cutoff    =   20.00
    Nmin             =  300
    Nmax             =  600
    sidechain_cutoff =  6.0
    CoulombModel     = "Linear"
    DesolvationModel = "ContactModel"
    BackBoneParameters  = None                    # set when class initialized
    SideChainParameters = None                    # set when class initialized
    angularDependentSideChainInteractions = None  # set when class initialized
    doingBackBoneReorganization = False
    scaleUpBuriedSideChain      = 0.00
    valueCooArgException = 2.40
    valueCooCooException = 1.60
    coulomb_list = ["COO", "CYS", "TYR", "HIS", "LYS", "ARG", "N+ "]
    exclude_sidechain_interactions = ["LYS", "ARG", "N+ "]
    atomInteractionList = getAtomInteractionList()


    def __init__(self, options=None):
        """
        constructer of the Default object.
        """
        str  = "WARNING: you are trying to create a 'default' version object. "
        str += " This object contains the cross section of all versions and is not a complete version itself."
        print(str)
        sys.exit(8)


    def printVersion(self):
        """
        Print out the properties of this version - more a test and debug feature.
        """
        print("version = %s" % (self.name))
        print("  desolvation   = %s" % (self.DesolvationModel))
        print("    prefactor   = %6.2lf" % (self.desolvationPrefactor))
        print("    allowance   = %6.2lf" % (self.desolvationAllowance))
        if   self.DesolvationModel in ["propka2", "ContactModel"]:
          print("    radii       = %s" % (self.desolvationRadii))
        if   self.DesolvationModel in ["VolumeModel", "ScaledVolumeModel"]:
          print("    surface     = %s" % (self.desolvationSurfaceScalingFactor))
        print("  Coulomb       = %s" % (self.CoulombModel))
        print("    cutoff      = %s" % (self.coulomb_cutoff))
        print("    scaled      = %s" % (self.coulomb_scaled))
        if   self.CoulombModel == "Linear":
          print("    prefactor   = %s" % (self.coulomb_maxpka))
        elif self.CoulombModel == "Coulomb":
          print("    diel        = %s" % (self.coulomb_diel))


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


    def setBackBoneMaxPKA(self, dpka=None, resType=None):
        """
        setting the max pKa shift due to all side-chain groups
        """
        if resType == None:
          # doing it for all residues
          #print("changing back-bone parameters")
          for key in self.BackBoneParameters.keys():
            self.BackBoneParameters[key][0] = dpka
            str  = " %s: %6.2lf " % (key, self.BackBoneParameters[key][0])
            str += "[%5.2lf,%5.2lf]" % (self.BackBoneParameters[key][1][0], self.BackBoneParameters[key][1][1])
            #print(str)
        else:
          # doing it just for residue 'key'
          key = resType
          #print("changing back-bone parameters")
          self.BackBoneParameters[key][0] = dpka
          str  = " %s: %6.2lf " % (key, self.BackBoneParameters[key][0])
          str += "[%5.2lf,%5.2lf]" % (self.BackBoneParameters[key][1][0], self.BackBoneParameters[key][1][1])
          #print(str)


    def setSideChainMaxPKA(self, dpka=None, resType=None):
        """
        setting the max pKa shift due to all side-chain groups
        """
        if resType == None:
          # doing it for all residues
          for key1 in self.SideChainParameters.keys():
            #print("changing side-chain parameters for resType \"%s\"" % (key1))
            for key2 in self.SideChainParameters[key1].keys():
              self.SideChainParameters[key1][key2][0] = dpka
              str  = " %s: %6.2lf " % (key2, self.SideChainParameters[key1][key2][0])
              str += "[%5.2lf,%5.2lf]" % (self.SideChainParameters[key1][key2][1][0], self.SideChainParameters[key1][key2][1][1])
              #print(str)
        else:
          # doing it just for residue 'key1'
          key1 = resType
          #print("changing side-chain parameters for resType \"%s\"" % (key1))
          for key2 in self.SideChainParameters[key1].keys():
            self.SideChainParameters[key1][key2][0] = dpka
            self.SideChainParameters[key2][key1][0] = dpka
            str  = " %s: %6.2lf " % (key2, self.SideChainParameters[key1][key2][0])
            str += "[%5.2lf,%5.2lf]" % (self.SideChainParameters[key1][key2][1][0], self.SideChainParameters[key1][key2][1][1])
            #print(str)


    def setCoulomb(self, label, max_dpka=None, cutoff=None, diel=None, scaled=None, mixed=False):
        """
        setting the Coulomb model, and its parameters. If parameter is not defined, it will take the default from 'parameters.py'
        """
        import parameters_new as parameters
        coulomb_parameters = parameters.getCoulombParameters()
        if label not in coulomb_parameters:
          print("do not accept Coulomb model \"%s\\n" % (label))
          sys.exit(9)
        else:
          self.CoulombModel = label
          default_parameters = coulomb_parameters[label]

        if 'max_dpka' in default_parameters:
          if max_dpka == None:
            self.coulomb_maxpka = default_parameters['max_dpka']
          else:
            self.coulomb_maxpka = max_dpka
        if 'cutoff' in default_parameters:
          if cutoff == None:
            self.coulomb_cutoff = default_parameters['cutoff']
          else:
            self.coulomb_cutoff = cutoff
        if 'diel' in default_parameters:
          if diel == None:
            self.coulomb_diel = default_parameters['diel']
          else:
            self.coulomb_diel = diel
        if 'scaled' in default_parameters:
          if scaled == None:
            self.coulomb_scaled = default_parameters['scaled']
          else:
            self.coulomb_scaled = scaled
        self.coulomb_mixed  = mixed



    def setDesolvation(self, label, prefactor=None, allowance=None, surface=None, volume=None, local=None, radii=None):
        """
        setting the desolvation model, and its parameters
        """
        import parameters_new as parameters
        desolvation_parameters = parameters.getDesolvationParameters()
        if label not in desolvation_parameters:
          print("do not accept solvation model \"%s\\n" % (label))
          sys.exit(9)
        else:
          self.DesolvationModel = label
          default_parameters = desolvation_parameters[label]

        if 'allowance' in default_parameters:
          if allowance == None:
            self.desolvationAllowance            = default_parameters['allowance']
          else:
            self.desolvationAllowance            = allowance
        if 'prefactor' in default_parameters:
          if prefactor == None:
            self.desolvationPrefactor            = default_parameters['prefactor']
          else:
            self.desolvationPrefactor            = prefactor
        if 'surface' in default_parameters:
          if surface == None:
            self.desolvationSurfaceScalingFactor = default_parameters['surface']
          else:
            self.desolvationSurfaceScalingFactor = surface
        if 'volume' in default_parameters:
          if volume    == None:
            self.desolvationVolume               = default_parameters['volume']
          else:
            self.desolvationVolume               = volume
        if 'local' in default_parameters:
          if local     == None:
            self.desolvationLocal                = default_parameters['local']
          else:
            self.desolvationLocal                = local
        if 'radii' in default_parameters:
          if radii     == None:
            self.desolvationRadii                = default_parameters['radii']
          else:
            self.desolvationRadii                = radii


    def calculateDesolvation(self, residue, atoms, options=None):
        """
        redirecting the desolvation calculation to the right model
        """
        if   self.DesolvationModel == "propka2":
          Nmass, Emass, Nlocl, Elocl = calculate.originalDesolvation(residue=residue, atoms=atoms, version=self, options=options)
        elif self.DesolvationModel == "ContactModel":
          Nmass, Emass, Nlocl, Elocl = calculate.contactDesolvation(residue, atoms, self, options=options)
        elif self.DesolvationModel in ["VolumeModel", "ScaledVolumeModel"]:
          Nmass, Emass, Nlocl, Elocl = calculate.radialVolumeDesolvation(residue, atoms, self, options=options)
        else:
          print("Desolvation \"%s\" is not implemented" % (self.DesolvationModel))
          sys.exit(8)

        return Nmass, Emass, Nlocl, Elocl


    def calculateBackBoneReorganization(self, protein, options=None):
        """
        Testing new term, reorganization of back-bone CO groups
        """
        if self.doingBackBoneReorganization == True:
          calculate.BackBoneReorganization(protein)
        else:
          """ do nothing """


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
        if False:
          prefactor = (1.0 + weight*self.scaleUpBuriedSideChain)
        else:
          prefactor = 1.00
        return calculate.HydrogenBondEnergy(distance, prefactor*dpka_max, cutoff, f_angle)


    def calculateCoulombEnergy(self, distance, weight, options=None):
        """
        redirects to get the correct Coulomb interaction - linear for default
        """
        if   self.CoulombModel == "Linear":
          return calculate.linearCoulombEnergy(distance, weight, self, options=options)
        elif self.CoulombModel == "Coulomb" and self.coulomb_mixed == True:
          return calculate.MixedCoulombEnergy(distance, weight, self, options=options)
        elif self.CoulombModel == "Coulomb":
          return calculate.CoulombEnergy(distance, weight, self, options=options)
        elif self.CoulombModel == "DistanceScaledCoulomb":
          return calculate.distanceScaledCoulombEnergy(distance, weight, self, options=options)
        else:
          print("Coulomb \"%s\" is not implemented" % (self.CoulombModel))
          sys.exit(8)
          

    def checkExceptions(self, residue1, residue2):
        """
        checks for exceptions for this version - using defaults
        """
        resType1 = residue1.resType
        resType2 = residue2.resType
        if   (resType1 == "COO" and resType2 == "ARG"):
          exception, value = checkCooArgException(residue1, residue2, version=self)
        elif (resType1 == "ARG" and resType2 == "COO"):
          exception, value = checkCooArgException(residue2, residue1, version=self)
        elif (resType1 == "COO" and resType2 == "COO"):
          exception, value = checkCooCooException(residue1, residue2, version=self)
        elif (resType1 == "CYS" and resType2 == "CYS"):
          exception, value = checkCysCysException(residue1, residue2, version=self)
        elif (resType1 == "COO" and resType2 == "HIS") or \
             (resType1 == "HIS" and resType2 == "COO"):
          exception, value = checkCooHisException(residue1, residue2, version=self)
        elif (resType1 == "CYS" and resType2 == "HIS") or \
             (resType1 == "HIS" and resType2 == "CYS"):
          exception, value = checkCysHisException(residue1, residue2, version=self)
        else:
          # do nothing, no exception for this pair
          exception = False; value = None

        return exception, value



#   --- specific versions with different behaviour or initialization ---

class Jan01(Version):
    """ 
    This is a test to set up rules for different propka Jan01 version
    """

    def __init__(self, options=None):
        """
        Rules of action for version Jan01
        """
        import parameters_std as parameters
    
        self.name             = "Jan01"
        self.Nmin             =  300
        self.Nmax             =  600
        self.buried_cutoff    =   15.50; self.buried_cutoff_sqr = self.buried_cutoff*self.buried_cutoff
        self.desolv_cutoff    =   15.50; self.desolv_cutoff_sqr = self.desolv_cutoff*self.desolv_cutoff
        self.setDesolvation("propka2")
        self.setCoulomb("Linear", cutoff=[4.0, 7.0], diel=None, scaled=False)
        self.doingBackBoneReorganization = False
        self.BackBoneParameters  = parameters.getHydrogenBondParameters(type='back-bone')
        self.SideChainParameters = parameters.getHydrogenBondParameters(type='side-chain')
        self.interaction         = getSideChainInteractionMatrix(side_chain=parameters.getInteraction()) # interaction rule matrix ['N'/'I'/'-']
        self.angularDependentSideChainInteractions = []
        self.exclude_sidechain_interactions = ["LYS", "ARG", "N+ "]


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


    def calculatePairWeight(self, Nmass1, Nmass2):
        """
        calculates the weight for the Coulomb interaction - used for version "Dec18" & "Sep23"
        """
        if lib.checkBuried(Nmass1, Nmass2) == False:
          return 0.0
        else:
          return 1.0


    def checkExceptions(self, residue1, residue2):
        """
        checks for exceptions for this version - using 'propka2' exceptions
        """
        resType1 = residue1.resType
        resType2 = residue2.resType
        if   (resType1 == "COO" and resType2 == "ARG"):
          exception, value = checkCooArgException_old(residue1, residue2, version=self)
        elif (resType1 == "ARG" and resType2 == "COO"):
          exception, value = checkCooArgException_old(residue2, residue1, version=self)
        elif (resType1 == "COO" and resType2 == "COO"):
          exception, value = checkCooCooException(residue1, residue2, version=self)
        elif (resType1 == "CYS" and resType2 == "CYS"):
          exception, value = checkCysCysException(residue1, residue2, version=self)
        elif (resType1 == "COO" and resType2 == "HIS") or \
             (resType1 == "HIS" and resType2 == "COO"):
          exception, value = checkCooHisException(residue1, residue2, version=self)
        elif (resType1 == "CYS" and resType2 == "HIS") or \
             (resType1 == "HIS" and resType2 == "CYS"):
          exception, value = checkCysHisException(residue1, residue2, version=self)
        else:
          # do nothing, no exception for this pair
          exception = False; value = None

        return exception, value





class Jan15(Version):
    """ 
    This is a test to set up rules for different propka Jan15 version
    """

    def __init__(self, options=None):
        """
        Rules of action for version Jan15
        """
        import parameters_std as parameters

        self.name             = "Jan15"
        self.Nmin             =  300
        self.Nmax             =  600
        self.buried_cutoff    =   15.50; self.buried_cutoff_sqr = self.buried_cutoff*self.buried_cutoff
        self.desolv_cutoff    =   15.50; self.desolv_cutoff_sqr = self.desolv_cutoff*self.desolv_cutoff
        self.setDesolvation("ContactModel", prefactor=-0.01, allowance=400.0)
        self.setCoulomb("Linear", cutoff=[4.0, 7.0], diel=None, scaled=False)
        self.doingBackBoneReorganization = False
        self.BackBoneParameters  = parameters.getHydrogenBondParameters(type='back-bone')
        self.SideChainParameters = parameters.getHydrogenBondParameters(type='side-chain')
        self.interaction         = getSideChainInteractionMatrix(side_chain=parameters.getInteraction()) # interaction rule matrix ['N'/'I'/'-']
        self.angularDependentSideChainInteractions = [] 
        self.exclude_sidechain_interactions = ["LYS", "ARG", "N+ "]


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


    def calculatePairWeight(self, Nmass1, Nmass2):
        """
        calculates the weight for the Coulomb interaction - used for version "Dec18" & "Sep23"
        """
        if lib.checkBuried(Nmass1, Nmass2) == False:
          return 0.0
        else:
          return 1.0


    def checkExceptions(self, residue1, residue2):
        """
        checks for exceptions for this version - using 'propka2' exceptions
        """
        resType1 = residue1.resType
        resType2 = residue2.resType
        if   (resType1 == "COO" and resType2 == "ARG"):
          exception, value = checkCooArgException_old(residue1, residue2, version=self)
        elif (resType1 == "ARG" and resType2 == "COO"):
          exception, value = checkCooArgException_old(residue2, residue1, version=self)
        elif (resType1 == "COO" and resType2 == "COO"):
          exception, value = checkCooCooException(residue1, residue2, version=self)
        elif (resType1 == "CYS" and resType2 == "CYS"):
          exception, value = checkCysCysException(residue1, residue2, version=self)
        elif (resType1 == "COO" and resType2 == "HIS") or \
             (resType1 == "HIS" and resType2 == "COO"):
          exception, value = checkCooHisException(residue1, residue2, version=self)
        elif (resType1 == "CYS" and resType2 == "HIS") or \
             (resType1 == "HIS" and resType2 == "CYS"):
          exception, value = checkCysHisException(residue1, residue2, version=self)
        else:
          # do nothing, no exception for this pair
          exception = False; value = None

        return exception, value





class May13(Version):
    """ 
    This is a test to set up rules for different propka May13 version
    """

    def __init__(self, options=None):
        """
        Rules of action for version May13
        """
        import parameters_std as parameters

        self.name             = "May13"
        self.coulomb_cutoff =   7.00
        print("creating propka version \"%s\"" % (self.name))
        self.BackBoneParameters  = parameters.getHydrogenBondParameters(type='back-bone')
        self.SideChainParameters = parameters.getHydrogenBondParameters(type='side-chain')
        self.interaction         = getSideChainInteractionMatrix(side_chain=parameters.getInteraction()) # interaction rule matrix ['N'/'I'/'-']
        self.angularDependentSideChainInteractions = [] 
        self.exclude_sidechain_interactions = ["LYS", "ARG", "N+ "]


    def checkExceptions(self, residue1, residue2):
        """
        overwrites 'exceptions' from the default - no exceptions here
        """
        exception = False
        value     = 0.00
        return exception, value


class Dec18(Version):
    """ 
    This is a test to set up rules for different propka versions
    """

    def __init__(self, options=None):
        """
        Rules of action for version Dec18
        """
        import parameters_std as parameters

        self.name             = "Dec18"
        self.Nmin             =  300
        self.Nmax             =  600
        self.buried_cutoff    =   15.50; self.buried_cutoff_sqr = self.buried_cutoff*self.buried_cutoff
        self.desolv_cutoff    =   15.50; self.desolv_cutoff_sqr = self.desolv_cutoff*self.desolv_cutoff
        self.setDesolvation("propka2")
        self.setCoulomb("Linear", cutoff=[4.0, 7.0], diel=None, scaled=True)
        self.doingBackBoneReorganization = False
        self.BackBoneParameters  = parameters.getHydrogenBondParameters(type='back-bone')
        self.SideChainParameters = parameters.getHydrogenBondParameters(type='side-chain')
        self.interaction         = getSideChainInteractionMatrix(side_chain=parameters.getInteraction()) # interaction rule matrix ['N'/'I'/'-']
        self.angularDependentSideChainInteractions = [] 
        self.exclude_sidechain_interactions = ["LYS", "ARG", "N+ "]


    def checkExceptions(self, residue1, residue2):
        """
        checks for exceptions for this version - using 'propka2' exceptions
        """
        resType1 = residue1.resType
        resType2 = residue2.resType
        if   (resType1 == "COO" and resType2 == "ARG"):
          exception, value = checkCooArgException_old(residue1, residue2, version=self)
        elif (resType1 == "ARG" and resType2 == "COO"):
          exception, value = checkCooArgException_old(residue2, residue1, version=self)
        elif (resType1 == "COO" and resType2 == "COO"):
          exception, value = checkCooCooException(residue1, residue2, version=self)
        elif (resType1 == "CYS" and resType2 == "CYS"):
          exception, value = checkCysCysException(residue1, residue2, version=self)
        elif (resType1 == "COO" and resType2 == "HIS") or \
             (resType1 == "HIS" and resType2 == "COO"):
          exception, value = checkCooHisException(residue1, residue2, version=self)
        elif (resType1 == "CYS" and resType2 == "HIS") or \
             (resType1 == "HIS" and resType2 == "CYS"):
          exception, value = checkCysHisException(residue1, residue2, version=self)
        else:
          # do nothing, no exception for this pair
          exception = False; value = None

        return exception, value





class Dec19(Version):
    """ 
    This is a test to set up rules for different propka versions
    """

    def __init__(self, options=None):
        """
        Rules of action for version Dec19
        """
        import parameters_std as parameters
        
        self.name             = "Dec19"
        self.Nmin             =  300
        self.Nmax             =  600
        self.buried_cutoff    =   15.50; self.buried_cutoff_sqr = self.buried_cutoff*self.buried_cutoff
        self.desolv_cutoff    =   15.50; self.desolv_cutoff_sqr = self.desolv_cutoff*self.desolv_cutoff
        self.setDesolvation("ContactModel", prefactor=-0.01, allowance=400.0)
        self.setCoulomb("Linear", cutoff=[4.0, 7.0], diel=None, scaled=True)
        self.doingBackBoneReorganization = False
        self.BackBoneParameters  = parameters.getHydrogenBondParameters(type='back-bone')
        self.SideChainParameters = parameters.getHydrogenBondParameters(type='side-chain')
        self.interaction         = getSideChainInteractionMatrix(side_chain=parameters.getInteraction()) # interaction rule matrix ['N'/'I'/'-']
        self.angularDependentSideChainInteractions = [] 
        self.exclude_sidechain_interactions = ["LYS", "ARG", "N+ "]


    def checkExceptions(self, residue1, residue2):
        """
        checks for exceptions for this version - using 'propka2' exceptions
        """
        resType1 = residue1.resType
        resType2 = residue2.resType
        if   (resType1 == "COO" and resType2 == "ARG"):
          exception, value = checkCooArgException(residue1, residue2, version=self)
        elif (resType1 == "ARG" and resType2 == "COO"):
          exception, value = checkCooArgException(residue2, residue1, version=self)
        elif (resType1 == "COO" and resType2 == "COO"):
          exception, value = checkCooCooException(residue1, residue2, version=self)
        elif (resType1 == "CYS" and resType2 == "CYS"):
          exception, value = checkCysCysException(residue1, residue2, version=self)
        elif (resType1 == "COO" and resType2 == "HIS") or \
             (resType1 == "HIS" and resType2 == "COO"):
          exception, value = checkCooHisException(residue1, residue2, version=self)
        elif (resType1 == "CYS" and resType2 == "HIS") or \
             (resType1 == "HIS" and resType2 == "CYS"):
          exception, value = checkCysHisException(residue1, residue2, version=self)
        else:
          # do nothing, no exception for this pair
          exception = False; value = None

        return exception, value





class Aug24(Version):
    """ 
    This is a test to set up rules for different propka versions
    """

    def __init__(self, options=None):
        """
        Rules of action for version Aug24
        """
        import parameters_new as parameters

        self.name             = "Aug24"
        self.Nmin             =  280
        self.Nmax             =  560
        self.buried_cutoff    =   15.00; self.buried_cutoff_sqr = self.buried_cutoff*self.buried_cutoff
        self.desolv_cutoff    =   20.00; self.desolv_cutoff_sqr = self.desolv_cutoff*self.desolv_cutoff
        self.setDesolvation("VolumeModel", prefactor=-7.00, surface=1.00)
        self.setCoulomb("Linear", cutoff=[4.0, 7.0], diel=80.0, scaled=True)
        self.doingBackBoneReorganization = True
        self.BackBoneParameters  = parameters.getHydrogenBondParameters(type='back-bone')
        self.SideChainParameters = parameters.getHydrogenBondParameters(type='side-chain')
        self.interaction         = getSideChainInteractionMatrix(side_chain=parameters.getInteraction()) # interaction rule matrix ['N'/'I'/'-']
        self.setBackBoneMaxPKA(resType=None, dpka=-1.15)
        self.setSideChainMaxPKA(resType=None, dpka=-1.15)
        self.angularDependentSideChainInteractions = ["HIS", "ARG", "AMD", "TRP"] 
        self.exclude_sidechain_interactions = []


class Aug30(Version):
    """ 
    This is a test to set up rules for different propka versions
    """

    def __init__(self, options=None):
        """
        Rules of action for version Aug30
        """
        import parameters_new as parameters

        self.name             = "Aug30"
        self.Nmin             =  280
        self.Nmax             =  560
        self.buried_cutoff    =   15.00; self.buried_cutoff_sqr = self.buried_cutoff*self.buried_cutoff
        self.desolv_cutoff    =   20.00; self.desolv_cutoff_sqr = self.desolv_cutoff*self.desolv_cutoff
        self.setDesolvation("VolumeModel", prefactor=-14.00, surface=0.00)
        self.setCoulomb("Linear", cutoff=[4.0, 7.0], diel=80.0, scaled=True)
        self.doingBackBoneReorganization = True
        self.BackBoneParameters  = parameters.getHydrogenBondParameters(type='back-bone')
        self.SideChainParameters = parameters.getHydrogenBondParameters(type='side-chain')
        self.interaction         = getSideChainInteractionMatrix(side_chain=parameters.getInteraction()) # interaction rule matrix ['N'/'I'/'-']
        self.setBackBoneMaxPKA(resType=None, dpka=-0.70)
        self.setSideChainMaxPKA(resType=None, dpka=-0.70)
        self.angularDependentSideChainInteractions = ["HIS", "ARG", "AMD", "TRP"] 
        self.exclude_sidechain_interactions = []


class Aug31(Version):
    """ 
    This is a test to set up rules for different propka versions
    """

    def __init__(self, options=None):
        """
        Rules of action for version Aug31
        """
        import parameters_new as parameters

        self.name             = "Aug31"
        self.Nmin             =  280
        self.Nmax             =  560
        self.buried_cutoff    =   15.00; self.buried_cutoff_sqr = self.buried_cutoff*self.buried_cutoff
        self.desolv_cutoff    =   20.00; self.desolv_cutoff_sqr = self.desolv_cutoff*self.desolv_cutoff0
        self.setDesolvation("VolumeModel", prefactor=-13.00, surface=0.25)
        self.setCoulomb("Linear", cutoff=[4.0, 7.0], diel=80.0, scaled=True)
        self.doingBackBoneReorganization = True
        self.BackBoneParameters  = parameters.getHydrogenBondParameters(type='back-bone')
        self.SideChainParameters = parameters.getHydrogenBondParameters(type='side-chain')
        self.interaction         = getSideChainInteractionMatrix(side_chain=parameters.getInteraction()) # interaction rule matrix ['N'/'I'/'-']
        self.setBackBoneMaxPKA(resType=None, dpka=-0.95)
        self.setSideChainMaxPKA(resType=None, dpka=-0.95)
        self.angularDependentSideChainInteractions = ["HIS", "ARG", "AMD", "TRP"] 
        self.exclude_sidechain_interactions = []


class Sep05(Version):
    """ 
    This is a test to set up rules for different propka versions
    """

    def __init__(self, options=None):
        """
        Rules of action for version Sep05
        """
        import parameters_new as parameters

        self.name             = "Sep05"
        self.Nmin             =  280
        self.Nmax             =  560
        self.buried_cutoff    =   15.00; self.buried_cutoff_sqr = self.buried_cutoff*self.buried_cutoff
        self.desolv_cutoff    =   20.00; self.desolv_cutoff_sqr = self.desolv_cutoff*self.desolv_cutoff
        self.setDesolvation("VolumeModel", prefactor=-7.50, surface=1.00)
        self.setCoulomb("Coulomb", cutoff=[4.0, 10.0], diel=[20.0, 80.0], scaled=False)
        self.doingBackBoneReorganization = True
        self.BackBoneParameters  = parameters.getHydrogenBondParameters(type='back-bone')
        self.SideChainParameters = parameters.getHydrogenBondParameters(type='side-chain')
        self.interaction         = getSideChainInteractionMatrix(side_chain=parameters.getInteraction()) # interaction rule matrix ['N'/'I'/'-']
        self.setBackBoneMaxPKA(resType=None, dpka=-1.00)
        self.setSideChainMaxPKA(resType=None, dpka=-1.00)
        self.angularDependentSideChainInteractions = ["HIS", "ARG", "AMD", "TRP"] 
        self.exclude_sidechain_interactions = []


class Sep06(Version):
    """ 
    This is a test to set up rules for different propka versions
    """

    def __init__(self, options=None):
        """
        Rules of action for version Sep06
        """
        import parameters_new as parameters

        self.name             = "Sep06"
        self.Nmin             =  280
        self.Nmax             =  560
        self.buried_cutoff    =   15.00; self.buried_cutoff_sqr = self.buried_cutoff*self.buried_cutoff
        self.desolv_cutoff    =   20.00; self.desolv_cutoff_sqr = self.desolv_cutoff*self.desolv_cutoff
        self.setDesolvation("VolumeModel", prefactor=-13.00, surface=0.00)
        self.setCoulomb("Coulomb", cutoff=[4.0, 10.0], diel=[20.0, 80.0], scaled=False)
        self.doingBackBoneReorganization = True
        self.BackBoneParameters  = parameters.getHydrogenBondParameters(type='back-bone')
        self.SideChainParameters = parameters.getHydrogenBondParameters(type='side-chain')
        self.interaction         = getSideChainInteractionMatrix(side_chain=parameters.getInteraction()) # interaction rule matrix ['N'/'I'/'-']
        self.setBackBoneMaxPKA(resType=None, dpka=-0.55)
        self.setSideChainMaxPKA(resType=None, dpka=-0.55)
        self.angularDependentSideChainInteractions = ["HIS", "ARG", "AMD", "TRP"] 
        self.exclude_sidechain_interactions = []


class Sep07(Version):
    """ 
    This is a test to set up rules for different propka versions
    """

    def __init__(self, options=None):
        """
        Rules of action for version Sep07
        """
        import parameters_new as parameters

        self.name             = "Sep07"
        self.Nmin             =  280
        self.Nmax             =  560
        self.buried_cutoff    =   15.00; self.buried_cutoff_sqr = self.buried_cutoff*self.buried_cutoff
        self.desolv_cutoff    =   20.00; self.desolv_cutoff_sqr = self.desolv_cutoff*self.desolv_cutoff
        self.setDesolvation("VolumeModel", prefactor=-12.50, surface=0.25)
        self.setCoulomb("Coulomb", cutoff=[4.0, 10.0], diel=[20.0, 80.0], scaled=False)
        self.doingBackBoneReorganization = True
        self.BackBoneParameters  = parameters.getHydrogenBondParameters(type='back-bone')
        self.SideChainParameters = parameters.getHydrogenBondParameters(type='side-chain')
        self.interaction         = getSideChainInteractionMatrix(side_chain=parameters.getInteraction()) # interaction rule matrix ['N'/'I'/'-']
        self.setBackBoneMaxPKA(resType=None, dpka=-0.75)
        self.setSideChainMaxPKA(resType=None, dpka=-0.75)
        self.angularDependentSideChainInteractions = ["HIS", "ARG", "AMD", "TRP"] 
        self.exclude_sidechain_interactions = []


class Sep08(Version):
    """ 
    This is a test to set up rules for different propka versions
    """

    def __init__(self, options=None):
        """
        Rules of action for version Sep08
        """
        import parameters_new as parameters

        self.name             = "Sep08"
        self.Nmin             =  280
        self.Nmax             =  560
        self.buried_cutoff    =   15.00; self.buried_cutoff_sqr = self.buried_cutoff*self.buried_cutoff
        self.desolv_cutoff    =   20.00; self.desolv_cutoff_sqr = self.desolv_cutoff*self.desolv_cutoff
        self.setDesolvation("VolumeModel", prefactor=-13.00, surface=0.25)
        self.setCoulomb("Coulomb", cutoff=[4.0, 10.0], diel=20.0, scaled=True)
        self.doingBackBoneReorganization = True
        self.BackBoneParameters  = parameters.getHydrogenBondParameters(type='back-bone')
        self.SideChainParameters = parameters.getHydrogenBondParameters(type='side-chain')
        self.interaction         = getSideChainInteractionMatrix(side_chain=parameters.getInteraction()) # interaction rule matrix ['N'/'I'/'-']
        self.setBackBoneMaxPKA(resType=None, dpka=-1.00)
        self.setSideChainMaxPKA(resType=None, dpka=-1.00)
        self.angularDependentSideChainInteractions = ["HIS", "ARG", "AMD", "TRP"] 
        self.exclude_sidechain_interactions = []


class Oct13(Version):
    """ 
    To test ligand integration
    """

    def __init__(self, options=None):
        """
        Rules of action for version Oct13, based on Sep07
        """
        import parameters_new as parameters

        self.name             = "Oct13"
        self.Nmin             =  280
        self.Nmax             =  560
        self.buried_cutoff    =   15.00; self.buried_cutoff_sqr = self.buried_cutoff*self.buried_cutoff
        self.desolv_cutoff    =   20.00; self.desolv_cutoff_sqr = self.desolv_cutoff*self.desolv_cutoff
        self.setDesolvation("VolumeModel", prefactor=-13.12, surface=0.40,  allowance=0.0)
        self.setCoulomb("Coulomb", cutoff=[4.0, 10.0], diel=80.0, scaled=False)
        self.doingBackBoneReorganization = True
        self.BackBoneParameters  = parameters.getHydrogenBondParameters(type='back-bone')
        self.SideChainParameters = parameters.getHydrogenBondParameters(type='side-chain')
        self.interaction         = getSideChainInteractionMatrix(side_chain=parameters.getInteraction()) # interaction rule matrix ['N'/'I'/'-']
        self.angularDependentSideChainInteractions = ["HIS", "ARG", "AMD", "TRP"] 
        self.exclude_sidechain_interactions = []
        
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





class Oct14(Version):
    """ 
    This is a test to set up rules for different propka versions
    """

    def __init__(self, options=None):
        """
        Rules of action for version Oct14
        """
        import parameters_new as parameters

        self.name             = "Oct14"
        self.Nmin             =  280
        self.Nmax             =  560
        self.buried_cutoff    =   15.00; self.buried_cutoff_sqr = self.buried_cutoff*self.buried_cutoff
        self.desolv_cutoff    =   20.00; self.desolv_cutoff_sqr = self.desolv_cutoff*self.desolv_cutoff
        self.setDesolvation("VolumeModel", prefactor=-13.00, surface=0.250)
        self.setCoulomb("Coulomb", cutoff=[3.5, 10.0], diel=80.0, scaled=False, mixed=True)
        self.doingBackBoneReorganization = True
        self.BackBoneParameters  = parameters.getHydrogenBondParameters(type='back-bone')
        self.SideChainParameters = parameters.getHydrogenBondParameters(type='side-chain')
        self.interaction         = getSideChainInteractionMatrix(side_chain=parameters.getInteraction()) # interaction rule matrix ['N'/'I'/'-']
        self.setBackBoneMaxPKA(resType=None, dpka=-0.75)
        self.setSideChainMaxPKA(resType=None, dpka=-0.75)
        self.angularDependentSideChainInteractions = ["HIS", "ARG", "AMD", "TRP"]
        self.exclude_sidechain_interactions = []


class Nov28(Version):
    """ 
    This is a test to set up rules for different propka versions
    """

    def __init__(self, options=None):
        """
        Rules of action for version Oct14
        """
        import parameters_new as parameters

        self.name             = "Nov28"
        self.Nmin             =  280
        self.Nmax             =  560
        self.buried_cutoff    =   15.00; self.buried_cutoff_sqr = self.buried_cutoff*self.buried_cutoff
        self.desolv_cutoff    =   20.00; self.desolv_cutoff_sqr = self.desolv_cutoff*self.desolv_cutoff
        self.setDesolvation("VolumeModel", prefactor=-8.00, surface=1.00)
        self.setCoulomb("DistanceScaledCoulomb", cutoff=[4.0, 10.0], diel=[20.0, 80.0])
        self.doingBackBoneReorganization = True
        self.BackBoneParameters  = parameters.getHydrogenBondParameters(type='back-bone')
        self.SideChainParameters = parameters.getHydrogenBondParameters(type='side-chain')
        self.interaction         = getSideChainInteractionMatrix(side_chain=parameters.getInteraction()) # interaction rule matrix ['N'/'I'/'-']
        self.setBackBoneMaxPKA(resType=None, dpka=-1.00)
        self.setSideChainMaxPKA(resType=None, dpka=-1.00)
        self.angularDependentSideChainInteractions = ["HIS", "ARG", "AMD", "TRP"]
        self.exclude_sidechain_interactions = []


class Nov29(Version):
    """ 
    This is a test to set up rules for different propka versions
    """

    def __init__(self, options=None):
        """
        Rules of action for version Nov29
        """
        import parameters_new as parameters

        self.name             = "Nov29"
        self.Nmin             =  280
        self.Nmax             =  560
        self.buried_cutoff    =   15.00; self.buried_cutoff_sqr = self.buried_cutoff*self.buried_cutoff
        self.desolv_cutoff    =   20.00; self.desolv_cutoff_sqr = self.desolv_cutoff*self.desolv_cutoff
        self.setDesolvation("VolumeModel", prefactor=-13.00, surface=0.00)
        self.setCoulomb("DistanceScaledCoulomb", cutoff=[4.0, 10.0], diel=[20.0, 80.0])
        self.doingBackBoneReorganization = True
        self.BackBoneParameters  = parameters.getHydrogenBondParameters(type='back-bone')
        self.SideChainParameters = parameters.getHydrogenBondParameters(type='side-chain')
        self.interaction         = getSideChainInteractionMatrix(side_chain=parameters.getInteraction()) # interaction rule matrix ['N'/'I'/'-']
        self.setBackBoneMaxPKA(resType=None, dpka=-0.40)
        self.setSideChainMaxPKA(resType=None, dpka=-0.40)
        self.angularDependentSideChainInteractions = ["HIS", "ARG", "AMD", "TRP"]
        self.exclude_sidechain_interactions = []


class Nov30(Version):
    """ 
    This is a test to set up rules for different propka versions
    """

    def __init__(self, options=None):
        """
        Rules of action for version Nov30
        """
        import parameters_new as parameters

        self.name             = "Nov30"
        self.Nmin             =  280
        self.Nmax             =  560
        self.buried_cutoff    =   15.00; self.buried_cutoff_sqr = self.buried_cutoff*self.buried_cutoff
        self.desolv_cutoff    =   20.00; self.desolv_cutoff_sqr = self.desolv_cutoff*self.desolv_cutoff
        self.setDesolvation("VolumeModel", prefactor=-13.00, surface=0.25)
        self.setCoulomb("DistanceScaledCoulomb", cutoff=[4.0, 10.0], diel=[30.0,160.0])
        self.doingBackBoneReorganization = True
        self.BackBoneParameters  = parameters.getHydrogenBondParameters(type='back-bone')
        self.SideChainParameters = parameters.getHydrogenBondParameters(type='side-chain')
        self.interaction         = getSideChainInteractionMatrix(side_chain=parameters.getInteraction()) # interaction rule matrix ['N'/'I'/'-']
        self.setBackBoneMaxPKA(resType=None, dpka=-0.85)
        self.setSideChainMaxPKA(resType=None, dpka=-0.85)
        self.angularDependentSideChainInteractions = ["HIS", "ARG", "AMD", "TRP"]
        self.exclude_sidechain_interactions = []


    def checkExceptions(self, residue1, residue2):
        """
        checks for exceptions for this version - using 'propka2' exceptions
        """
        resType1 = residue1.resType
        resType2 = residue2.resType
        if   (resType1 == "COO" and resType2 == "ARG"):
          exception, value = checkCooArgException(residue1, residue2, version=self)
        elif (resType1 == "ARG" and resType2 == "COO"):
          exception, value = checkCooArgException(residue2, residue1, version=self)
        elif (resType1 == "COO" and resType2 == "COO"):
          exception, value = checkCooCooException(residue1, residue2, version=self)
        elif (resType1 == "CYS" and resType2 == "CYS"):
          exception, value = checkCysCysException(residue1, residue2, version=self)
        #elif (resType1 == "COO" and resType2 == "HIS") or \
        #     (resType1 == "HIS" and resType2 == "COO"):
        #  exception, value = checkCooHisException(residue1, residue2, version=self)
        elif (resType1 == "CYS" and resType2 == "HIS") or \
             (resType1 == "HIS" and resType2 == "CYS"):
          exception, value = checkCysHisException(residue1, residue2, version=self)
        else:
          # do nothing, no exception for this pair
          exception = False; value = None

        return exception, value


