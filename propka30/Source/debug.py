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
import string
import lib
import calculator as calculate

from lib import pka_print


def interactionMatrix(interaction):
    """
    printing out all information in resInfo
    """
    keys = ["COO", "CYS", "TYR", "HIS", "N+ ", "LYS", "ARG", "ROH", "AMD", "TRP"]

    pka_print("interaction matrix:")
    for key1 in keys:
      str = "%6s:" % (key1)
      for key2 in keys:
        do_pair, iterative = interaction[key1][key2]
        if   do_pair == True and iterative == True:
          str += "%3s" % ("I")
        elif do_pair == True and iterative == False:
          str += "%3s" % ("N")
        else:
          str += "%3s" % ("-")

      pka_print(str)


def printResInfo(resInfo):
    """
    printing out all information in resInfo
    """
    pka_print("in resInfo:")
    for key1 in resInfo.keys():
      pka_print(" --- %s ---" % (key1))
      for key2 in resInfo[key1].keys():
        pka_print(key2, resInfo[key1][key2])


def printCooArgAtomDistances(residue_coo, residue_arg):
    """
    printing out all COO and ARG distances for debugging, picking closest + runner-up
    """
    for atom_coo in residue_coo.makeDeterminantAtomList(residue_arg.resName):
      for atom_arg in residue_arg.makeDeterminantAtomList(residue_coo.resName):
        distance = calculate.InterAtomDistance(atom_coo, atom_arg)
        pka_print("%3s %3s %6.2lf" % (atom_coo.name, atom_arg.name, distance))


def printBackBoneAtoms(list):
    """
    Prints out determinant information for debugging
    """
    pka_print("  --- debug back-bone atom list --- ")
    for atoms in list:
        label1 = "%s%4d%2s" % (atoms[0].resName, atoms[0].resNumb, atoms[0].chainID)
        label2 = "%s%4d%2s" % (atoms[1].resName, atoms[1].resNumb, atoms[1].chainID)
        pka_print("%s - %s" % (label1, label2))


def printResidues(residue_list):
    """
    Prints out determinant information for debugging
    """
    pka_print("  --- debug residue list --- ")
    for residue in residue_list:
        residue.printLabel()


def printIterativeDeterminants(all_determinants):
    """
    Prints out determinant information for debugging
    """

    if True:
      pka_print(" --- Iterative determinants ---")
      pka_print("%28s%8s" % ("H-bond", "Coulomb") )
      for determinant in all_determinants:
        pair   =  determinant[0]
        values =  determinant[1]
        annihilation = determinant[2]
        residue1 = pair[0]
        residue2 = pair[1]
        str  = ""
        str += " %s %s" % (residue1.label, residue2.label)
        str += " %6.2lf %6.2lf" % (values[0], values[1])
        pka_print(str)
    else:
      # priting out ALL types of determinants - mainly debugging
      pka_print("  --- debug iterative determinants --- ")
      pka_print("len(sidechain_determinants) = %d" % (len(all_determinants[0])))
      pka_print("len(backbone_determinants)  = %d" % (len(all_determinants[1])))
      pka_print("len(coulomb_determinants)   = %d" % (len(all_determinants[2])))
    
      sidechain_determinants = all_determinants[0]
      backbone_determinants  = all_determinants[1]
      coulomb_determinants   = all_determinants[2]
    
      if True:
        for type in range(0,3):
          determinants = all_determinants[type]
          if   type == 0:
            pka_print("Iterative side-chain interactions:")
          elif type == 1:
            pka_print("Iterative back-bone  interactions:")
          elif type == 2:
            pka_print("Iterative Coulomb    interactions:")
          for determinant in determinants:
            pair   =  determinant[0]
            values =  determinant[1]
    
            pair   = determinant[0]
            value  = determinant[1]
            residue1 = pair[0]
            residue2 = pair[1]
            str  = ""
            str += "'%s'  '%s'" % (residue1.label, residue2.label)
            str += "   "
            str += "%6.2lf  %6.2lf" % (value, value)
            pka_print(str)


def printAlignment(alignment):
    """
    Prints out alignment information for debugging
    """
    for key in alignment.keys():
      pka_print( " --- %s ---" % (key) )
      for key2 in alignment[key].keys():
        pka_print("%s %5d%2s" % (alignment[key][key2]["name"], alignment[key][key2]["resNumb"], alignment[key][key2]["chainID"]))
        pka_print("%s\n" % (alignment[key][key2]["sequence"]))


def printAllAtoms(protein):
    """
    Prints out determinant information for debugging
    """
    pka_print("Not implemented yet")
