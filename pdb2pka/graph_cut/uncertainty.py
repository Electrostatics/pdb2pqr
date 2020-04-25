from __future__ import print_function
import sys
from itertools import product
import random
random.seed("Mmmmm, sandwiches...")
import logging


_LOGGER = logging.getLogger(__name__)


def resolve_uncertainty(protein_complex, labeling, uncertain, brute_force_limit=20, verbose=False):
    """Resolves all uncertain residue states using either brute force or a MC
       protein_complex - ProteinComplex instance for the protein
       labeling - ResidueVariable instance to ResidueInstance instance map for certain residue states
       uncertain - list of uncertain ResidueVariable instances
       brute_force_limit - limit of the number of uncertain residues before falling back to MC."""

    final_labeling = labeling.copy()

    if not uncertain:
        return final_labeling

    _LOGGER.debug("Uncertain count:", len(uncertain))

    if len(uncertain) > brute_force_limit:
        _LOGGER.debug("Using Monte Carlo")
        return monte_carlo(protein_complex, final_labeling, uncertain)
    else:
        _LOGGER.debug("Using brute force")
        return brute_force(protein_complex, final_labeling, uncertain)


def brute_force(pc, labeling, uncertain):
    """Iterate through all possible combinations of uncertain
       residue states picking the one that produces the best results."""
    result_labeling = None
    best_energy = sys.float_info.max

    state_pairs = ((x.instances["PROTONATED"], x.instances["DEPROTONATED"]) for x in uncertain)

    test_labeling = labeling.copy()
    for test_states in product(*state_pairs):
        for x, y in zip(uncertain, test_states):
            test_labeling.update([(x, y)])
        energy = pc.evaluate_energy(test_labeling, normal_form=True)
        if energy < best_energy:
            best_energy = energy
            result_labeling = test_labeling.copy()

    return result_labeling


def monte_carlo(pc, labeling, uncertain):
    """
    Pseudo code for the MC:
     Set starting best energy to the system max float value
     For each step in the MC
         Randomly assign states to the uncertain residues
         For each sub step in the MC
             Choose a random uncertain residue
             Choose the state of that residue that minimizes the total energy
             If both states result in the same energy
                 Choose a random state
         If the current total energy for the resulting is less than the best energy
             Save the current state as the new best state
             Save the current total energy at the new best energy
     Return the current best state.
    """
    result_labeling = None
    best_energy = sys.float_info.max

    state_pairs = [(x.instances["PROTONATED"], x.instances["DEPROTONATED"]) for x in uncertain]

    test_labeling = labeling.copy()
    #iterations = min(10000, 2 ** (len(uncertain) - 1))
    #iterations = max(iterations, 1)
    iterations = 1000
    sub_iterations = 500

    for _ in range(iterations):
        test_states = (random.choice(state_pair) for state_pair in state_pairs)
        test_labeling.update((x,y) for x,y in zip(uncertain, test_states))
        last_residue = None
        for _ in range(sub_iterations):
            random_residue = random.choice(uncertain)
            if random_residue is last_residue:
                continue
            last_residue = random_residue
            diff = pc.evaluate_energy_diff(random_residue, test_labeling, normal_form=True)
            if diff < 0.0:
                test_labeling[random_residue] = random_residue.instances["PROTONATED"]
            elif diff > 0.0:
                test_labeling[random_residue] = random_residue.instances["DEPROTONATED"]
            else:
                test_labeling[random_residue] = random.choice(list(random_residue.instances.values()))

        energy = pc.evaluate_energy(test_labeling, normal_form=True)
        if energy < best_energy:
            best_energy = energy
            result_labeling = test_labeling.copy()

    return result_labeling
