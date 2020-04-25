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
from sys import exit
import os
import math
import Source.version as propka
import Source.lib as lib
from Source.protein import Protein
from Source.mutate import makeCompositeAtomsDictionary
pka_print = lib.pka_print


def main():
    """
    This is a tailor-made propka-executable for propka's GUI, but also useful for commandline execution
    """

    # I. preliminaries

    # loading options, flaggs and arguments; making a version object
    options, pdbfiles = lib.loadOptions()
    version = propka.makeVersion(label=options.version_label)

    # creating the protein object
    myProtein = Protein(pdbfile=pdbfiles[0], options=options)

    # creating a dictionary with atom objects needed for e.g. alignment mutations
    atoms = makeCompositeAtomsDictionary(protein=myProtein, pdbfiles=pdbfiles, options=options)

    # II. optimise 'single-site' mutation and selecting the good determinants
    mutations = []
    for mutation in options.mutations:
        # II.a. combinatorial search of 'determinants' for each site
        best_mutation = myProtein.optimizeMutationDeterminants(atoms=atoms, mutation=mutation, version=version, options=options)
        if best_mutation != None:
            mutations.append(best_mutation)

    # III. combinatorial search of single-site mutations with their resulting 'determinants'
    if len(mutations) > 0:
        best_mutation = myProtein.optimizeMultipleMutations(mutations=mutations, atoms=atoms, version=version, options=options)
    else:
        pka_print("Could not find any mutation combination more stable than WT\n")
        exit(8)


if __name__ == '__main__':
    main()
