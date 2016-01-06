from __future__ import print_function
from graph import ProteinGraph
from uncertainty import resolve_uncertainty
from collections import defaultdict
import math
from pprint import pprint
import sys

# Use the number for R from https://en.wikipedia.org/wiki/Gas_constant
gas_constant = 8.3144621
Rln10 = math.log(10) * gas_constant
T = 300.0
Rln10_T = Rln10*T
RT = gas_constant * T

# TODO - figure out why modPkaHIP is hard-coded here!
modPkaHIP = 6.6
modPkaHIE = modPkaHIP
modPkaHID = modPkaHIP

def print_pc_state(pc, normal_form, out_file):
    """Dump protein_complex state to out_file
       normal_form - Dump the normal form part of the state."""
    rv = pc.residue_variables
    ie = pc.normalized_interaction_energies if normal_form else pc.interaction_energies_for_ph
    for v_residue in rv.values():
        for v_instance in v_residue.instances.values():
            for w_residue in rv.values():
                if v_residue == w_residue:
                    continue
                for w_instance in w_residue.instances.values():
                    out_file.write(str((v_instance, w_instance)) + " " + str(round(ie[v_instance, w_instance],4)) + '\n')
    keys = list(pc.residue_variables.keys())
    for key in keys:
        residue = pc.residue_variables[key]
        for instance in list(residue.instances.values()):
            if normal_form:
                out_file.write(str(instance) + " " + str(round(instance.energyNF,4)) + "\n")
            else:
                out_file.write(str(instance) + " " + str(round(instance.energy_with_ph,4)) + "\n")
    if normal_form:
        out_file.write("Normalized constant energy: " + str(round(pc.normalized_constant_energy,4)) + "\n")

def print_dg_state(dg, out_file):
    """Dump directed graph state to out_file"""
    out_file.write("Flow network:\nVertices:\n")

    nodes = list(dg.node.keys())
    nodes.sort()

    for node in nodes:
        out_file.write('_'.join(node)+"\n")

    out_file.write("\nEdges:\n")

    edges = []
    for edge in dg.edges_iter(data="capacity"):
        result = []
        if isinstance(edge[0], tuple):
            result.append('_'.join(edge[0]))
        else:
            result.append(edge[0])
        if isinstance(edge[1], tuple):
            result.append('_'.join(edge[1]))
        else:
            result.append(edge[1])

        result.append(edge[2])

        edges.append(result)

    edges.sort()

    for edge in edges:
        out_file.write("(")
        out_file.write(edge[0])
        out_file.write(", ")
        out_file.write(edge[1])
        out_file.write(")= ")
        out_file.write(str(round(edge[2],4))+"\n")


def get_titration_curves(protein_complex, state_file=None):
    """For each ph value:
           Get the normal form of the protein energies.
           Build a flow graph
           Get the min cut of the graph
           Find which state for each residue from the cut (labeling) and the unknown states (uncertain)
           Use brute force or MC to resolve the uncertain states.
           Calculate the curve value for each residue

        Returns results for all residues for each ph."""

    curves = defaultdict(list)
    pg = ProteinGraph(protein_complex)
    pH = 0.0
    step_size = 0.1
    end_ph = 20.0
    steps = int(end_ph / step_size) + 1

    for step in range(steps):
        pH = step * step_size
        print("pH", pH)

        if state_file is not None:
            state_file.write ("pH="+ str(pH)+"\n")
            state_file.write("REGULAR ENERGIES\n")
            protein_complex.energy_at_pH(pH)
            print_pc_state(protein_complex, False, state_file)
            state_file.write('\n')
            state_file.write("NORMAL FORM ENERGIES\n")

        protein_complex.normalize(pH)

        if state_file is not None:
            print_pc_state(protein_complex, True, state_file)
            state_file.write('\n')

        pg.update_graph()

        if state_file is not None:
            print_dg_state(pg.DG, state_file)
            state_file.write('\n')

        cv, s_nodes, t_nodes = pg.get_cut()
        labeling, uncertain = pg.get_labeling_from_cut(s_nodes, t_nodes)
        new_labeling = resolve_uncertainty(protein_complex, labeling, uncertain, verbose=True)

        curve_values = get_curve_values(protein_complex, new_labeling, pH)
        for key, value in curve_values.items():
            curves[key].append((pH, value))

    return curves

def get_curve_values(protein_complex, labeling, pH):
    """Using the given selected residue states (labeling) and pH get the
       current curve value for all titratable residues."""
    his_seen = set()
    results = {}

    aH = math.pow(10, -pH)

    for key, residue in protein_complex.residue_variables.items():
        name, chain, location = key

        if name in ("HId", "HIe"):
            #Do HIS stuff
            if (chain, location) in his_seen:
                continue
            his_seen.add((chain, location))

            if name == "HId":
                hid_residue = residue
                hie_residue = protein_complex.residue_variables["HIe", chain, location]
            else:
                hie_residue = residue
                hid_residue = protein_complex.residue_variables["HId", chain, location]

            # dge = HSP - HSE, dgd = HSP - HSD
            class Energies:
                pass
            energies = Energies()
            energies.aH = aH
            energies.dGdref = modPkaHID*math.log(10.0)
            energies.dGeref = modPkaHIE*math.log(10.0)
            energies.dGe, energies.dGd = protein_complex.evaluate_energy_diff_his(hie_residue, hid_residue, labeling,
                                                                normal_form=True)

            debug_craziness = False
            if debug_craziness:
                print("!!! DEBUG - SETTING dGd from %g to 0.0" % energies.dGd)
                energies.dGd = 0
                print("!!! DEBUG - SETTING dGe from %g to 0.0" % energies.dGe)
                energies.dGe = 0
            else:
                old = energies.dGd
                energies.dGd =  energies.dGd - math.log(aH) - energies.dGdref
                #print("Removed extra pH and pKa contributions from dGd: %g -> %g" % (old, energies.dGd))
                old = energies.dGe
                energies.dGe =  energies.dGe - math.log(aH) - energies.dGeref
                #print("Removing extra pH and pKa contributions from dGe: %g -> %g" % (old, energies.dGe))
            energies.ddG = (energies.dGe + energies.dGeref) - (energies.dGd + energies.dGdref)
            energies.dGp = energies.dGd - energies.dGdref
            #print(vars(energies))
            pHSD = 1.0
            pHSE = math.exp(-energies.ddG)
            pHSP = energies.aH*math.exp(-energies.dGp)
            Q = pHSD + pHSE + pHSP
            fracHSD = pHSD/Q
            fracHSE = pHSE/Q
            fracHSP = pHSP/Q

            # if not labeling[hie_residue].protonated and labeling[hid_residue].protonated:
            #     titration_value = fracHSD
            # elif labeling[hie_residue].protonated and not labeling[hid_residue].protonated:
            #     titration_value = fracHSE
            # elif labeling[hie_residue].protonated and labeling[hid_residue].protonated:
            #     titration_value = fracHSP
            # else:
            #     errstr = "How did we get here?"
            #     raise RuntimeError(errstr)
            results["HIS", chain, location] = fracHSP
            #results["HSE", chain, location] = fracHSE
            #results["HSD", chain, location] = fracHSD

        else:
            #Do not HIS stuff
            #energy_diff = protonated_energy - depotonated_energy
            energy_diff = protein_complex.evaluate_energy_diff(residue, labeling, normal_form=True)

            exp = -(energy_diff)
            #Handle case where there is an unresolved bump.
            try:
                e_exp = math.exp(exp)
                titration_value = e_exp/(1.0+e_exp)
            except OverflowError:
                titration_value = 1.0
            results[key] = titration_value

    return results
