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
from sys import exit
import copy
import math
from Source.lib import pka_print


def getCatchAllCompareWithExperiment(argv):
    """
    returns mutation label
    """
    import mats
    if "mutation" in argv:
        mutation = argv['mutation']
        if mutation.upper() == "NONE" or mutation.upper() == "ALL":
            mutation = None
    else:
        mutation = None

    if mutation == None:
        return "ALL"
    else:
        return mats.makeLabelFromMutationTag(mutation)


def getExperiment(fullname):
    """
    return experimental values
    """
    name = fullname[:4]
    if name == "1pga":
        # B1_ProteinG
        experiment = {
            "GLU  15 A":  4.40,
            "GLU  19 A":  3.70,
            "ASP  22 A":  2.90,
            "GLU  27 A":  4.50,
            "ASP  36 A":  3.80,
            "ASP  40 A":  4.00,
            "GLU  42 A":  4.40,
            "ASP  46 A":  3.60,
            "ASP  47 A":  3.40,
            "GLU  56 A":  4.00,
        }
    elif name == "1igd":
        # B2_ProteinG
        experiment = {
            "GLU  20 A":  4.30,
            "ASP  27 A":  2.90,
            "GLU  29 A":  4.20,
            "GLU  32 A":  4.60,
            "ASP  41 A":  3.90,
            "ASP  45 A":  4.40,
            "ASP  51 A":  3.60,
            "ASP  52 A":  3.40,
            "GLU  61 A":  4.20,
        }
    elif name == "1a2p":
        # Barnase
        experiment = {
            "ASP   8 A":  2.90,
            "ASP  12 A":  3.80,
            "ASP  22 A":  3.30,
            "GLU  29 A":  3.80,
            "ASP  44 A":  3.40,
            # "ASP  54 A":  2.20, # <= 2.2
            "GLU  60 A":  3.00,
            # "GLU  73 A":  2.10, # <= 2.1
            #   "ASP  75 A":  3.10, # difficult barnase pi-stacking
            "ASP  86 A":  4.20,
            # "ASP  93 A":  2.00, # <  2.0
            # "ASP 101 A":  2.00, # <= 2
        }
    elif name == "1rnz":
        # Bovine_RNase
        experiment = {
            "GLU   2 A":  2.70,
            "GLU   9 A":  4.00,
            # "ASP  14 A":  1.80, # < 2.0
            "ASP  38 A":  2.80,
            "GLU  49 A":  4.50,
            "ASP  53 A":  3.80,
            "ASP  83 A":  3.40,
            "GLU  86 A":  4.05,
            "GLU 111 A":  3.50,
            "ASP 121 A":  3.05,
        }
    elif name == "4icb":
        # Calbindin
        experiment = {
            "GLU   4 A":  3.80,
            "GLU   5 A":  3.40,
            "GLU  11 A":  4.70,
            "GLU  17 A":  3.60,
            "GLU  26 A":  4.10,
            "ASP  47 A":  3.00,
            "GLU  48 A":  4.60,
            "GLU  64 A":  3.80,
        }
    elif name == "1kxi":
        # Cardiotoxin
        experiment = {
            "GLU  17 A":  4.00,
            "ASP  42 A":  3.20,
            # "ASP  59 A":  2.30, # < 2.3
            "GLU  17 B":  4.00,
            "ASP  42 B":  3.20,
            # "ASP  59 B":  2.30, # < 2.3
        }
    elif name == "1kxi_A":
        # Cardiotoxin
        experiment = {"GLU  17 A":  4.00,
                      "ASP  42 A":  3.20,
                      "ASP  59 A":  2.30}
    elif name == "1kxi_B":
        # Cardiotoxin
        experiment = {"GLU  17 B":  4.00,
                      "ASP  42 B":  3.20,
                      "ASP  59 B":  2.30}
    elif name == "1cdc":
        # CD2d1
        experiment = {
            "ASP  25 A":  3.50,
            "ASP  26 A":  3.60,
            "ASP  28 A":  3.60,
            "GLU  29 A":  4.40,
            "GLU  33 A":  4.20,
            "GLU  41 A":  6.70,
            "GLU  56 A":  3.90,
            "ASP  62 A":  4.20,
            "ASP  71 A":  3.20,
            "ASP  72 A":  4.10,
            "ASP  94 A":  3.90,
            "GLU  99 A":  4.20,
            "ASP  25 B":  3.50,
            "ASP  26 B":  3.60,
            "ASP  28 B":  3.60,
            "GLU  29 B":  4.40,
            "GLU  33 B":  4.20,
            "GLU  41 B":  6.70,
            "GLU  56 B":  3.90,
            "ASP  62 B":  4.20,
            "ASP  71 B":  3.20,
            "ASP  72 B":  4.10,
            "ASP  94 B":  3.90,
            "GLU  99 B":  4.20,
        }
    elif name == "1beo":
        # Cryptogein
        experiment = {
            "ASP  21 A":  2.50,
            "ASP  30 A":  2.50,
            "ASP  72 A":  2.60,
        }
    elif name == "1hpx":
        # HIV_protease, note, homodimer, B is assumed to be A
        experiment = {
            # "ASP  25 B":  6.20, # > 6.2
            "ASP  29 B":  3.20,
            "ASP  30 B":  3.90,
            "ASP  60 B":  3.00,
            # "ASP  25 A":  2.50, # < 2.5
            "ASP  29 A":  3.70,
            "ASP  30 A":  3.80,
            "ASP  60 A":  3.00,
        }
    elif name == "1lys":
        # Hen lysozyme
        experiment = {
            "GLU   7 A":  2.90,
            "ASP  18 A":  2.70,
            "GLU  35 A":  6.20,
            # "ASP  48 A":  2.50, # < 2.5
            "ASP  52 A":  3.70,
            # "ASP  66 A":  2.00, # < 2.5
            "ASP  87 A":  2.10,
            "ASP 101 A":  4.10,
            "ASP 119 A":  3.20,
        }
    elif name == "135l":
        # Turkey lysozyme
        experiment = {
            "GLU   7 A":  2.70,
            "ASP  18 A":  2.70,
            "GLU  35 A":  6.10,
            # "ASP  48 A":  2.50, # < 2.5
            "ASP  52 A":  3.80,
            # "ASP  66 A":  2.00, # < 2.5
            "ASP  87 A":  2.10,
            "ASP 119 A":  3.40,
        }
    elif name == "2rn2":
        # Ribonuclease H1
        experiment = {
            "GLU   6 A":  4.50,
            "ASP  10 A":  6.10,
            "GLU  32 A":  3.60,
            "GLU  48 A":  4.40,
            "GLU  57 A":  3.20,
            "GLU  61 A":  3.90,
            "GLU  64 A":  4.40,
            "ASP  70 A":  2.60,
            "ASP  94 A":  3.20,
            # "ASP 102 A":  2.00, # < 2
            "ASP 108 A":  3.20,
            "GLU 119 A":  4.10,
            "GLU 129 A":  3.60,
            "GLU 131 A":  4.30,
            "ASP 134 A":  4.10,
            "GLU 135 A":  4.30,
            "GLU 147 A":  4.20,
            # "ASP 148 A":  2.00, # < 2
            "GLU 154 A":  4.40,
        }
    elif name == "2ovo":
        # Ovomucoid3Dom
        experiment = {
            "ASP   7 A":  2.50,
            "GLU  10 A":  4.10,
            "GLU  19 A":  3.20,
            "ASP  27 A":  2.20,
            "GLU  43 A":  4.80,
        }
    elif name == "1xnb":
        # Xylanase
        experiment = {
            "ASP   4 A":  3.00,
            "ASP  11 A":  2.50,
            "GLU  78 A":  4.60,
            # "ASP  83 A":  2.00, # < 2
            # "ASP 101 A":  2.00, # < 2
            "ASP 106 A":  2.70,
            "ASP 119 A":  3.20,
            "ASP 121 A":  3.60,
            "GLU 172 A":  6.70,
        }
    elif name == "1de3":
        # Sarcin
        experiment = {
            "ASP   9 A":  3.90,
            "GLU  19 A":  4.60,
            "GLU  31 A":  4.60,
            # "ASP  41 A":  3.00, # < 3
            "ASP  57 A":  4.30,
            "ASP  59 A":  4.10,
            "ASP  75 A":  3.90,
            # "ASP  77 A":  3.00, # < 3
            "ASP  85 A":  3.80,
            # "ASP  91 A":  3.00, # < 3
            "GLU  96 A":  5.10,
            # "ASP 102 A":  3.00, # < 3
            # "ASP 105 A":  3.00, # < 3
            "ASP 109 A":  3.70,
            "GLU 115 A":  4.90,
            "GLU 140 A":  4.30,
            "GLU 144 A":  4.30,
        }
    elif name == "2bus":
        # BullSeminalInhibitor
        experiment = {
            "ASP   6 A":  4.00,
            "GLU   9 A":  4.30,
            "ASP  12 A":  3.60,
            "GLU  20 A":  4.10,
        }
    elif name == "1egf":
        # Epidermal Growthfactor
        experiment = {
            "ASP  11 A":  3.90,
            "GLU  24 A":  4.10,
            "ASP  27 A":  4.00,
            "ASP  40 A":  3.60,
            "ASP  46 A":  3.80,
            "GLU  51 A":  4.00,
        }
    elif name == "1mhi":
        # Insulin
        experiment = {
            "GLU   4 A":  2.60,
            # "GLU  17 A":  3.70, # > 3.7
            "ASP   9 B":  2.60,
            "GLU  13 B":  2.20,
            "GLU  21 B":  3.70,
        }
    else:
        pka_print("could not find experimental data for \"%s\"" % (name))
        exit(9)

    return experiment


def makeErrorPlot(points):
    """
    do patented error graph data
    """
    error_list = []
    for i in range(0, 41):
        error_list.append(0.00)

    for point in points:
        abs_diff = abs(point)
        i = 0
        while abs_diff > float(i)/10.0 and i < 41:
            error_list[i] += 1.0
            i += 1

    number_of_points = float(len(points))

    for i in range(0, 41):
        error = float(i)/10.0
        fraction = error_list[i]/number_of_points
        pka_print("%6.2lf%6.2lf" % (error, fraction))
