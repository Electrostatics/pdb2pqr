#!/usr/bin/env python
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


import string
import re
import sys
import os
import math

from .Source import lib
from .Source.protein import Protein


def main():
    """
    Simple main that creates the proteins, calculates pKa values, and prints pKa files
    """

    # loading options, flaggs and arguments
    options, pdbfiles = lib.loadOptions()

    for pdbfile in pdbfiles:
        # creating protein object
        myProtein = Protein(pdbfile=pdbfile, options=options)

        # calculating pKa values for ionizable residues
        myProtein.calculatePKA(options=options)
        # printing pka file
        myProtein.writePKA(options=options)


if __name__ == '__main__':
    main()
