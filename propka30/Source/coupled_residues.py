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

import math
from .lib import pka_print

max_intrinsic_pKa_diff = 2.0
min_interaction_energy = 0.5

max_free_energy_diff = 1.0
min_swap_pka_shift = 1.0
#pH = 7.0
pH = 'variable'
reference = 'neutral'

min_pka = 0.0
max_pka = 10.0

do_intrinsic = False
do_pair_wise = False
do_prot_stat = True

# pkas 0-10 average?


def is_coupled_pairwise(residue1, residue2):
    interaction1 = get_interaction(residue1, residue2)
    interaction2 = get_interaction(residue2, residue1)

    interaction_free_pka1 = residue1.pKa_pro - interaction1
    interaction_free_pka2 = residue2.pKa_pro - interaction2
    max_interaction = max(interaction1, interaction2)

    factor_diff_intrinsic_pka = get_pka_diff_factor(interaction_free_pka1, interaction_free_pka2)
    factor_interaction = get_interaction_factor(max_interaction)
    res = factor_diff_intrinsic_pka*factor_interaction

    return (res, interaction_free_pka1, interaction_free_pka2, max_interaction)


def is_coupled_intrinsic_pka(residue1, residue2):
    """ Checks if residue1 and residue2 are coupled """

    # calculate intrinsic pKa's, if not already done
    for residue in [residue1, residue2]:
        if not hasattr(residue, 'intrinsic_pKa'):
            residue.calculateIntrinsicPKA()

    # check if intrinsic pKa's are similar
    factor_diff_intrinsic_pka = get_pka_diff_factor(residue1.intrinsic_pKa, residue2.intrinsic_pKa)

    # check if residues have electrostatic interaction
    interaction_energy = max(get_interaction(residue1, residue2), get_interaction(residue2, residue1))
    factor_interaction = get_interaction_factor(interaction_energy)

    res = factor_diff_intrinsic_pka * factor_interaction

    return (res, interaction_energy)


def is_coupled_protonation_state_probability(protein, residue1, residue2, options=None):

    interaction_energy = max(get_interaction(residue1, residue2), get_interaction(residue2, residue1))
    if interaction_energy <= min_interaction_energy:
        return {'coupling_factor': -1.0}

    # calculate intrinsic pKa's, if not already done
    for residue in [residue1, residue2]:
        if not hasattr(residue, 'intrinsic_pKa'):
            residue.calculateIntrinsicPKA()

    use_pH = pH
    if pH == 'variable':
        use_pH = min(residue1.pKa_pro, residue2.pKa_pro)

    #default_energy = protein.calculateFoldingEnergy(pH=use_pH, reference="neutral")
    default_energy = protein.calculateFoldingEnergy(pH=use_pH, options=options)
    default_pka1 = residue1.pKa_pro
    default_pka2 = residue2.pKa_pro

    # check that pka values are within relevant limits
    if max(default_pka1, default_pka2) < min_pka or min(default_pka1, default_pka2) > max_pka:
        return {'coupling_factor': -1.0}

    # Swap interactions and re-calculate pKa values
    swap_interactions(residue1, residue2, verbose=False)
    residue1.calculateTotalPKA()
    residue2.calculateTotalPKA()

    # store swapped energy and pka's
    swapped_energy = protein.calculateFoldingEnergy(pH=use_pH, options=options)
    swapped_pka1 = residue1.pKa_pro
    swapped_pka2 = residue2.pKa_pro

    pka_shift1 = swapped_pka1 - default_pka1
    pka_shift2 = swapped_pka2 - default_pka2

    # Swap back to original protonation state
    swap_interactions(residue1, residue2, verbose=False)
    residue1.calculateTotalPKA()
    residue2.calculateTotalPKA()

    if abs(default_energy - swapped_energy) <= max_free_energy_diff:
        if max(abs(pka_shift1), abs(pka_shift2)) >= min_swap_pka_shift:
            if abs(residue1.intrinsic_pKa - residue2.intrinsic_pKa) <= max_intrinsic_pKa_diff:
                factor = get_free_energy_diff_factor(default_energy, swapped_energy) *\
                    get_pka_diff_factor(residue1.intrinsic_pKa, residue2.intrinsic_pKa) *\
                    get_interaction_factor(interaction_energy)

                return {'coupling_factor': factor,
                        'default_energy': default_energy,
                        'swapped_energy': swapped_energy,
                        'interaction_energy': interaction_energy,
                        'swapped_pka1': swapped_pka1,
                        'swapped_pka2': swapped_pka2,
                        'pka_shift1': pka_shift1,
                        'pka_shift2': pka_shift2,
                        'pH': use_pH}

    return {'coupling_factor': -1.0}


def get_pka_diff_factor(pka1, pka2):
    intrinsic_pka_diff = abs(pka1-pka2)
    res = 0.0
    if intrinsic_pka_diff <= max_intrinsic_pKa_diff:
        res = 1-(intrinsic_pka_diff/max_intrinsic_pKa_diff)**2

    return res


def get_free_energy_diff_factor(energy1, energy2):
    free_energy_diff = abs(energy1-energy2)
    res = 0.0
    if free_energy_diff <= max_free_energy_diff:
        res = 1-(free_energy_diff/max_free_energy_diff)**2
    return res


def get_interaction_factor(interaction_energy):
    res = 0.0
    interaction_energy = abs(interaction_energy)
    if interaction_energy >= min_interaction_energy:
        res = (interaction_energy-min_interaction_energy)/(1.0+interaction_energy-min_interaction_energy)

    return res


def identify_coupled_residues(protein, options=None):
    """ Finds coupled residues in protein """

    verbose = options.display_coupled_residues

    if True:
        pka_print('')
        pka_print(' Detecting coupled residues')
        pka_print('   Maximum pKa difference:     %4.2f pKa units' % max_intrinsic_pKa_diff)
        pka_print('   Minimum interaction energy: %4.2f pKa units' % min_interaction_energy)
        pka_print('   Maximum free energy diff.:  %4.2f pKa units' % max_free_energy_diff)
        pka_print('   Minimum swap pKa shift:     %4.2f pKa units' % min_swap_pka_shift)
        pka_print('   pH:                         %6s ' % str(pH))
        pka_print('   Reference:                  %s' % reference)
        pka_print('   Min pKa:                    %4.2f' % min_pka)
        pka_print('   Max pKa:                    %4.2f' % max_pka)
        pka_print('')

    # make a single list of all residues in the protein
    all_residues = []
    for chain in protein.chains:
        for residue in chain.residues:
            if not residue in protein.residue_dictionary["ION"]:
                all_residues.append(residue)

    # find coupled residues
    for i in range(len(all_residues)):
        for j in range(len(all_residues)):
            if i == j:
                break
            swap = 0
            if do_intrinsic:
                (coupling_factor_intrinsic_pka, interaction_energy) = is_coupled_intrinsic_pka(all_residues[i], all_residues[j])
                if coupling_factor_intrinsic_pka > 0.0:
                    swap = 1
                    pka_print(' %s and %s coupled (intrinsic pKa): %4.2f | intrinsic pkas: %5.2f - %5.2f = %5.2f | int. energy: %4.2f' % (all_residues[i],
                                                                                                                                          all_residues[j],
                                                                                                                                          coupling_factor_intrinsic_pka,
                                                                                                                                          all_residues[i].intrinsic_pKa,
                                                                                                                                          all_residues[j].intrinsic_pKa,
                                                                                                                                          all_residues[i].intrinsic_pKa -
                                                                                                                                          all_residues[j].intrinsic_pKa,
                                                                                                                                          interaction_energy))

            if do_pair_wise:
                (coupling_factor_pairwise, int_free_pka1, int_free_pka2, interaction_energy) = is_coupled_pairwise(all_residues[i], all_residues[j])
                if coupling_factor_pairwise > 0.0:
                    swap = 1
                    pka_print(' %s and %s coupled (pair wise)    : %4.2f | int. free pkas: %5.2f - %5.2f = %5.2f | int. energy: %4.2f' % (all_residues[i],
                                                                                                                                          all_residues[j],
                                                                                                                                          coupling_factor_pairwise,
                                                                                                                                          int_free_pka1, int_free_pka2,
                                                                                                                                          int_free_pka1 - int_free_pka2,
                                                                                                                                          interaction_energy))

            if do_prot_stat:
                data = is_coupled_protonation_state_probability(protein, all_residues[i], all_residues[j], options=options)
                if data['coupling_factor'] > 0.0:
                    swap = 1
                    all_residues[i].coupled_residues.append(all_residues[j])
                    all_residues[j].coupled_residues.append(all_residues[i])
                    protein.coupled_residues = True

                    #pka_print('Coupled residue of', all_residues[i],':', all_residues[j])
                    #pka_print('Coupled residue of', all_residues[j],':', all_residues[i])

                    if verbose:
                        pka_print(make_data_to_string(data, all_residues[i], all_residues[j]))

            if swap and verbose:
                # swap...
                swap_interactions(all_residues[i], all_residues[j])
                # ...and swap back
                swap_interactions(all_residues[i], all_residues[j], verbose=False)

    return


def make_data_to_string(data, residue1, residue2):
    s = """ %s and %s coupled (prot.state): %5.2f
 Energy levels:       %6.2f, %6.2f  (difference: %6.2f) at pH %6.2f
 Interaction energy:  %6.2f 
 Intrinsic pka's:     %6.2f, %6.2f  (difference: %6.2f)
 Swapped pKa's:       %6.2f, %6.2f  (difference: %6.2f, %6.2f)""" % (residue1,
                                                                     residue2,
                                                                     data['coupling_factor'],
                                                                     data['default_energy'], data['swapped_energy'],
                                                                     data['default_energy'] - data['swapped_energy'],
                                                                     data['pH'],
                                                                     data['interaction_energy'],
                                                                     residue1.intrinsic_pKa,
                                                                     residue2.intrinsic_pKa,
                                                                     residue1.intrinsic_pKa-residue2.intrinsic_pKa,
                                                                     data['swapped_pka1'],
                                                                     data['swapped_pka2'],
                                                                     data['pka_shift1'],
                                                                     data['pka_shift2'])

    return s


def get_interaction(residue1, residue2, include_side_chain_hbs=True):
    determinants = residue1.determinants[2]
    if include_side_chain_hbs:
        determinants = residue1.determinants[0] + residue1.determinants[2]

    interaction_energy = 0.0
    for det in determinants:
        if residue2.label == det.label:
            interaction_energy += det.value

    pka_print(' '.join((str(residue1), str(residue2), str(interaction_energy))))

    return interaction_energy


def swap_interactions(residue1, residue2, include_side_chain_hbs=True, verbose=True):

    if verbose:
        pka_print(' '+'-'*113)
        tagged_pka_print(' Original|', residue1.getDeterminantString(), [residue1.label, residue2.label])
        tagged_pka_print(' Original|', residue2.getDeterminantString(), [residue1.label, residue2.label])

    # swap the interactions!
    transfer_determinant(residue1.determinants[2], residue2.determinants[2], residue1.label, residue2.label)
    if include_side_chain_hbs:
        transfer_determinant(residue1.determinants[0], residue2.determinants[0], residue1.label, residue2.label)

    # re-calculate pKa values
    residue1.calculateTotalPKA()
    residue2.calculateTotalPKA()

    if verbose:
        tagged_pka_print(' Swapped |', residue1.getDeterminantString(), [residue1.label, residue2.label])
        tagged_pka_print(' Swapped |', residue2.getDeterminantString(), [residue1.label, residue2.label])
        pka_print(' '+'='*113)
        pka_print('')
    return


def transfer_determinant(determinants1, determinants2, label1, label2):
    # find out what to transfer...
    from1to2 = []
    from2to1 = []
    for det in determinants1:
        if det.label == label2:
            from1to2.append(det)

    for det in determinants2:
        if det.label == label1:
            from2to1.append(det)

    # ...and transfer it!
    for det in from1to2:
        det.label = label1
        determinants2.append(det)
        determinants1.remove(det)

    for det in from2to1:
        det.label = label2
        determinants1.append(det)
        determinants2.remove(det)

    return


def tagged_pka_print(tag, s, labels):
    s = "%s %s" % (tag, s)
    s = s.replace('\n', '\n%s ' % tag)
    for label in labels:
        s = s.replace(label, '\033[31m%s\033[30m' % label)
    pka_print(s)
    return
