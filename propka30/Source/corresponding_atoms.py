#!/usr/local/bin/python3.0

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

import sys
from lib import residueList, atomList


#    NOTE:
#    The sole purpose of this module it to set up a 'corresponding-atoms-dictionary',
#    which is used for the overlap procedure. The first loop sets up all back-bone 
#    atoms for all residue types and 'CB' for all residue-type pairs not involving GLY.
#    The second loop sets all atoms for self-residue overlap, i.e. ASP-ASP, VAL-VAL.
#    The final section sets up the remaining residue-pair specific corresponding atoms 
#    for the side-chain. The resulting dictionary might not be complete, or agree with 
#    what you expect for the moment; feel free to change it accordingly.


def makeCorrespondingAtomNames():
    """
    setting up a dictionary to define 'corresponding atoms' between two residues for 'overlap'
    """

    # getting list of all atoms
    resNames = residueList("standard")

    names = {}

    #   ----- back-bone & 'CB' section -----
    # simplifying the setup by including back-bone and 'CB' with this loop
    for resName1 in resNames:
      names[resName1] = {}
      for resName2 in resNames:
        names[resName1][resName2] = [['N',  'N'],
                                     ['CA', 'CA'],
                                     ['C',  'C'],
                                     ['O',  'O']]
        if resName1 != "GLY" and resName2 != "GLY":
          names[resName1][resName2].append(['CB', 'CB'])

    #   ----- self-overlap section -----
    # setting up all atoms for self-comparison
    for resName in resNames:
      atmNames = atomList(resName)
      for atmName in atmNames:
        if atmName not in ['N', 'CA', 'CB', 'C', 'O']:
          names[resName][resName].append([atmName, atmName])


    #   ----- side-chain section -----
    # side-chains left to consider (sorted alphabetically):
    # ['ARG', 'ASN', 'ASP', 'CYS', 'GLN', 'GLU', 'HIS', 'ILE', 'LEU', 'LYS', 'MET', 'PHE', 'PRO', 'SER', 'THR', 'TRP', 'TYR', 'VAL']


    # ARG
    # None 


    # ASN
    str1 = "ASN  CG"
    str2 = "ARG  CG"
    extendCorrespondingAtomsDictionary(names, str1, str2)


    # ASP
    str1 = "ASP  CG"
    str2 = "ARG  CG"
    extendCorrespondingAtomsDictionary(names, str1, str2)

    str1 = "ASP  OD1 OD2"
    str2 = "ASN  OD1 ND2"
    extendCorrespondingAtomsDictionary(names, str1, str2)


    # CYS
    # None 


    # GLN
    str1 = "GLN  CG CD"
    str2 = "ARG  CG CD"
    extendCorrespondingAtomsDictionary(names, str1, str2)

    str1 = "GLN  CG"
    str2 = "ASN  CG"
    extendCorrespondingAtomsDictionary(names, str1, str2)

    str1 = "GLN  CG"
    str2 = "ASP  CG"
    extendCorrespondingAtomsDictionary(names, str1, str2)


    # GLU


    # HIS


    # ILE


    # LEU


    # LYS


    # MET


    # PHE


    # PRO


    # SER


    # THR


    # TRP


    # TYR
    str1 = "TYR  CG"
    str2 = "LYS  CG"
    extendCorrespondingAtomsDictionary(names, str1, str2)


    # VAL


    return  names






def extendCorrespondingAtomsDictionary(names, str1, str2):
    """
    extends the pairs based on list1 & list2
    """
    list1 = str1.split()
    list2 = str2.split()
    for i in range(1, len(list1)):
      names[list1[0]][list2[0]].append([list1[i], list2[i]])
      names[list2[0]][list1[0]].append([list2[i], list1[i]])

    return  None



def main():
    """
    Simple check on the corresponding atoms-dictionary
    """

    corresponding_atoms = makeCorrespondingAtomNames()

    resNames = residueList("standard")

    for resName1 in resNames:
      for resName2 in resNames:
        str = "%s %s \n" % (resName1, resName2)
        for i in range(len(corresponding_atoms[resName1][resName2])):
          name1, name2 = corresponding_atoms[resName1][resName2][i]
          str += " %-3s %-3s" % (name1, name2)
          name1, name2 = corresponding_atoms[resName2][resName1][i]
          str += "%-5s %-3s %-3s\n" % (" ", name1, name2)
        print(str)
        if resName1 == resName2:
          break


if __name__ == '__main__': main()

