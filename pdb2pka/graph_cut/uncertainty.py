import sys
from itertools import product, izip
import random
random.seed("Mmmmm, sandwiches...")


def resolve_uncertainty(protein_complex, labeling, uncertain, brute_force_limit=12):

    final_labeling = labeling.copy()

    if not uncertain:
        return final_labeling

    print "Uncertain count:", len(uncertain)

    if len(uncertain) > brute_force_limit:
        print "Using Monte Carlo"
        return monte_carlo(protein_complex, final_labeling, uncertain)
    else:
        print "Using brute force"
        return brute_force(protein_complex, final_labeling, uncertain)


def brute_force(pc, labeling, uncertain):
    result_labeling = None
    best_energy = sys.float_info.max

    state_pairs = ((x.instances["PROTONATED"], x.instances["DEPROTONATED"]) for x in uncertain)

    test_labeling = labeling.copy()
    for test_states in product(*state_pairs):
        test_labeling.update((x,y) for x,y in izip(uncertain, test_states))
        energy = pc.evaluate_energy(test_labeling, normal_form=True)
        if energy < best_energy:
            best_energy = energy
            result_labeling = test_labeling.copy()

    return result_labeling


def monte_carlo(pc, labeling, uncertain):
    result_labeling = None
    best_energy = sys.float_info.max

    state_pairs = [(x.instances["PROTONATED"], x.instances["DEPROTONATED"]) for x in uncertain]

    test_labeling = labeling.copy()
    #iterations = min(10000, 2 ** (len(uncertain) - 1))
    #iterations = max(iterations, 1)
    iterations = 500
    sub_iterations = 200

    for _ in xrange(iterations):
        test_states = (random.choice(state_pair) for state_pair in state_pairs)
        test_labeling.update((x,y) for x,y in izip(uncertain, test_states))
        last_residue = None
        for _ in xrange(sub_iterations):
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
                test_labeling[random_residue] = random.choice(random_residue.instances.values())

        energy = pc.evaluate_energy(test_labeling, normal_form=True)
        if energy < best_energy:
            best_energy = energy
            result_labeling = test_labeling.copy()

    return result_labeling

