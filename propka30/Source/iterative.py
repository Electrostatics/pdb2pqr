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
from .determinant import Determinant
import math
import time
from . import calculator as calculate
from . import lib
pka_print = lib.pka_print
#import debug


# Some library functions for the interative pKa determinants


def addtoDeterminantList(residue1, residue2, distance, iterative_interactions, version):
    """
    Adds 'iterative determinants' to list ..., [[R1, R2], [side-chain, coulomb], [A1, A2]], ...
    Note, the sign is determined when the interaction is added to the iterative object!
    Note, distance < coulomb_cutoff here
    """
    add_interaction = False
    hbond_value = 0.00
    coulomb_value = 0.00

    # Side-chain interactions
    if True:
        atoms1 = residue1.makeDeterminantAtomList(residue2.resName, version=version)
        atoms2 = residue2.makeDeterminantAtomList(residue1.resName, version=version)
        atom_distance = 999.
        for atom1 in atoms1:
            for atom2 in atoms2:
                # select the smallest inter-atom distance
                atom_distance = min(calculate.InterAtomDistance(atom1, atom2), atom_distance)
        dpka_max, cutoff = version.SideChainParameters[residue1.resType][residue2.resType]
        weight = version.calculatePairWeight(residue1.Nmass, residue2.Nmass)
        if atom_distance < cutoff[1]:
            add_interaction = True
            exception, hbond_value = version.checkExceptions(residue1, residue2)
            if residue1.resType == "COO" and residue2.resType == "COO":
                """ do nothing """
                #pka_print("xxx %6.2lf" % (atom_distance))
            # exception = False # circumventing exception
            if exception == True:
                """ do nothing, value should have been assigned """
                #pka_print(" exception for %s %s (I)" % (residue1.label, residue2.label))
            else:
                f_angle = 1.0
                hbond_value = version.calculateSideChainEnergy(atom_distance, dpka_max, cutoff, weight, f_angle)

    # Back-bone interactions
    """ Not done, never iterative """

    # Coulomb interactions
    do_coulomb = version.checkCoulombPair(residue1, residue2, distance)
    if do_coulomb == True:
        add_interaction = True
        weight = version.calculatePairWeight(residue1.Nmass, residue2.Nmass)
        coulomb_value = version.calculateCoulombEnergy(distance, weight)

    # adding the interaction to 'iterative_interactions'
    if add_interaction == True:
        interaction = []
        pair = [residue1, residue2]
        values = [hbond_value, coulomb_value]
        annihilation = [0., 0.]
        interaction = [pair, values, annihilation]
        iterative_interactions.append(interaction)


def addIterativeAcidPair(object1, object2, interaction):
    """ 
    Adding the Coulomb 'iterative' interaction (an acid pair):
    the higher pKa is raised  with QQ+HB
    the lower  pKa is lowered with HB
    """
    values = interaction[1]
    annihilation = interaction[2]
    hbond_value = values[0]
    coulomb_value = values[1]
    diff = coulomb_value + 2*hbond_value
    label1 = object1.label
    label2 = object2.label
    comp1 = object1.pKa_old + annihilation[0] + diff
    comp2 = object2.pKa_old + annihilation[1] + diff
    annihilation[0] = 0.
    annihilation[1] = 0.
    if comp1 > comp2:
        # side-chain
        determinant = [label2,  hbond_value]
        object1.determinants[0].append(determinant)
        determinant = [label1, -hbond_value]
        object2.determinants[0].append(determinant)
        # Coulomb
        determinant = [label2, coulomb_value]
        object1.determinants[2].append(determinant)
        annihilation[0] = -diff
    else:
        # side-chain
        determinant = [label1,  hbond_value]
        object2.determinants[0].append(determinant)
        determinant = [label2, -hbond_value]
        object1.determinants[0].append(determinant)
        # Coulomb
        determinant = [label1, coulomb_value]
        object2.determinants[2].append(determinant)
        annihilation[1] = -diff


def addIterativeBasePair(object1, object2, interaction):
    """ 
    Adding the Coulomb 'iterative' interaction (a base pair):
    the lower pKa is lowered
    """
    values = interaction[1]
    annihilation = interaction[2]
    hbond_value = values[0]
    coulomb_value = values[1]
    diff = coulomb_value + 2*hbond_value
    diff = -diff
    label1 = object1.label
    label2 = object2.label
    comp1 = object1.pKa_old + annihilation[0] + diff
    comp2 = object2.pKa_old + annihilation[1] + diff
    annihilation[0] = 0.
    annihilation[1] = 0.
    if comp1 < comp2:
        # side-chain
        determinant = [label2, -hbond_value]
        object1.determinants[0].append(determinant)
        determinant = [label1,  hbond_value]
        object2.determinants[0].append(determinant)
        # Coulomb
        determinant = [label2, -coulomb_value]
        object1.determinants[2].append(determinant)
        annihilation[0] = -diff
    else:
        # side-chain
        determinant = [label1, -hbond_value]
        object2.determinants[0].append(determinant)
        determinant = [label2,  hbond_value]
        object1.determinants[0].append(determinant)
        # Coulomb
        determinant = [label1, -coulomb_value]
        object2.determinants[2].append(determinant)
        annihilation[1] = -diff


def addIterativeIonPair(object1, object2, interaction, version):
    """ 
    Adding the Coulomb 'iterative' interaction (an acid-base pair):
    the pKa of the acid is lowered & the pKa of the base is raised
    """
    values = interaction[1]
    annihilation = interaction[2]
    hbond_value = values[0]
    coulomb_value = values[1]
    Q1 = object1.Q
    Q2 = object2.Q
    comp1 = object1.pKa_old + annihilation[0] + Q1*coulomb_value
    comp2 = object2.pKa_old + annihilation[1] + Q2*coulomb_value
    if object1.resName not in version.exclude_sidechain_interactions:
        comp1 += Q1*hbond_value
    if object2.resName not in version.exclude_sidechain_interactions:
        comp2 += Q2*hbond_value

    if Q1 == -1.0 and comp1 < comp2:
        add_term = True  # pKa(acid) < pKa(base)
    elif Q1 == 1.0 and comp1 > comp2:
        add_term = True  # pKa(base) > pKa(acid)
    else:
        add_term = False

    annihilation[0] = 0.00
    annihilation[1] = 0.00

    if add_term == True:

        # Coulomb
        if coulomb_value > 0.005:
            # residue1
            interaction = [object2.label, Q1*coulomb_value]
            annihilation[0] += -Q1*coulomb_value
            object1.determinants[2].append(interaction)
            # residue2
            interaction = [object1.label, Q2*coulomb_value]
            annihilation[1] += -Q2*coulomb_value
            object2.determinants[2].append(interaction)

        # Side-chain
        if hbond_value > 0.005:
            # residue1
            if object1.resName not in version.exclude_sidechain_interactions:
                interaction = [object2.label, Q1*hbond_value]
                annihilation[0] += -Q1*hbond_value
                object1.determinants[0].append(interaction)
            # residue2
            if object2.resName not in version.exclude_sidechain_interactions:
                interaction = [object1.label, Q2*hbond_value]
                annihilation[1] += -Q2*hbond_value
                object2.determinants[0].append(interaction)


def addDeterminants(iterative_interactions, version, options=None):
    """ 
    The iterative pKa scheme. Later it is all added in 'calculateTotalPKA'
    """
    # --- setup ---
    iteratives = []
    done_residue = []
    # debug.printIterativeDeterminants(iterative_interactions)
    # creating iterative objects with references to their real residue counterparts
    for interaction in iterative_interactions:
        pair = interaction[0]
        for residue in pair:
            if residue in done_residue:
                # print "done already"
                """ do nothing - already have an iterative object for this residue """
            else:
                newIterative = Iterative(residue)
                iteratives.append(newIterative)
                done_residue.append(residue)

    # Initialize iterative scheme
    if options.print_iterations == True:
        pka_print("\n   --- pKa iterations (%d residues, %d interactions) ---" % (len(iteratives), len(iterative_interactions)))
    converged = False
    iteration = 0
    for itres in iteratives:
        itres.pKa_iter.append(itres.pKa_NonIterative)

    # --- starting pKa iterations ---
    while converged == False:

        # initialize pKa_new
        iteration += 1
        for itres in iteratives:
            itres.determinants = [[], [], []]
            itres.pKa_new = itres.pKa_NonIterative

        # Adding interactions to temporary determinant container
        for interaction in iterative_interactions:
            pair = interaction[0]
            values = interaction[1]
            annihilation = interaction[2]
            # print "len(interaction) = %d" % (len(interaction))
            object1, object2 = findIterative(pair, iteratives)
            Q1 = object1.Q
            Q2 = object2.Q
            if Q1 < 0.0 and Q2 < 0.0:
                """ both are acids """
                addIterativeAcidPair(object1, object2, interaction)
            elif Q1 > 0.0 and Q2 > 0.0:
                """ both are bases """
                addIterativeBasePair(object1, object2, interaction)
            else:
                """ one of each """
                addIterativeIonPair(object1, object2, interaction, version)

        # Calculating pKa_new values
        for itres in iteratives:
            for type in range(0, 3):
                for determinant in itres.determinants[type]:
                    itres.pKa_new += determinant[1]

        # Check convergence
        converged = True
        for itres in iteratives:
            if itres.pKa_new == itres.pKa_old:
                itres.converged = True
            else:
                itres.converged = False
                converged = False

        # reset pKa_old & storing pKa_new in pKa_iter
        for itres in iteratives:
            itres.pKa_old = itres.pKa_new
            itres.pKa_iter.append(itres.pKa_new)

        if iteration == 10:
            pka_print("did not converge in %d iterations" % (iteration))
            break

    # --- Iterations finished ---

    # printing pKa iterations
    if options.print_iterations == True:
        str = "%12s" % (" ")
        for index in range(0, iteration+1):
            str += "%8d" % (index)
        pka_print(str)
        for itres in iteratives:
            str = "%s   " % (itres.label)
            for pKa in itres.pKa_iter:
                str += "%8.2lf" % (pKa)
            if itres.converged == False:
                str += " *"
            pka_print(str)

    # creating real determinants and adding them to residue object
    for itres in iteratives:
        for type in range(0, 3):
            for interaction in itres.determinants[type]:
                value = interaction[1]
                if value > 0.005 or value < -0.005:
                    label = interaction[0]
                    newDeterminant = Determinant(label, value)
                    itres.residue.determinants[type].append(newDeterminant)


def findIterative(pair, iteratives):
    """
    Function to find the two 'iteratives' that corresponds to the residues in 'pair'
    """
    for iterative in iteratives:
        if iterative.residue == pair[0]:
            iterative0 = iterative
        elif iterative.residue == pair[1]:
            iterative1 = iterative

    return iterative0, iterative1


class Iterative:
    """
        Iterative class - pKa values and references of iterative residues
        Note, this class has a fake determinant list, true determinants are
              made after the iterations are finished.
    """

    def __init__(self, residue):
        """
        Contructer of the iterative object
        """

        # print "creating 'iterative object' for %s" % (residue.label)

        self.label = residue.label
        self.resName = residue.resName
        self.Q = residue.Q
        self.pKa_old = None
        self.pKa_new = None
        self.pKa_iter = []
        self.pKa_NonIterative = 0.00
        self.determinants = [[], [], []]
        self.residue = residue
        self.converged = True

        # Calculate the Non-Iterative part of pKa from the residue object
        # Side chain
        side_chain = 0.00
        for determinant in residue.determinants[0]:
            value = determinant.value
            side_chain += value

        # Back bone
        back_bone = 0.00
        for determinant in residue.determinants[1]:
            value = determinant.value
            back_bone += value

        # Coulomb
        coulomb = 0.00
        for determinant in residue.determinants[2]:
            value = determinant.value
            coulomb += value

        self.pKa_NonIterative = residue.pKa_mod
        self.pKa_NonIterative += residue.Emass
        self.pKa_NonIterative += residue.Elocl
        self.pKa_NonIterative += side_chain
        self.pKa_NonIterative += back_bone
        self.pKa_NonIterative += coulomb

        self.pKa_old = self.pKa_NonIterative
