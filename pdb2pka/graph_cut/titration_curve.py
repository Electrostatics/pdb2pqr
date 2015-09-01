from graph import ProteinGraph
from uncertainty import resolve_uncertainty
from collections import defaultdict
import math

gas_constant = 8.3144621
#ln(10) * gas constant
kln10 = math.log(10) * gas_constant
# Degrees Kelvin
T = 300
RT = 2479.0;
RT_gas_T = RT * gas_constant * T
kln10_T = kln10 * T

modPkaHIP = 6.6
modPkaHIE = modPkaHIP
modPkaHID = modPkaHIP


def get_titration_curves(protein_complex):
    curves = defaultdict(list)

    pg = ProteinGraph(protein_complex)

    pH = 0.0
    step_size = 0.1
    end_ph = 20.0
    steps = int(end_ph / step_size) + 1

    for step in xrange(steps):
        pH = step * 0.1
        print "Processing pH:", pH
        protein_complex.normalize(pH)
        pg.update_edges()

        cv, s_nodes, t_nodes = pg.get_cut()
        labeling, uncertain = pg.get_labeling_from_cut(s_nodes, t_nodes)
        new_labeling = resolve_uncertainty(protein_complex, labeling, uncertain)

        curve_values = get_curve_values(protein_complex, new_labeling, pH)
        for key, value in curve_values.iteritems():
            curves[key].append((pH, value))

    return curves

def get_curve_values(protein_complex, labeling, pH):
    his_seen = set()
    results = {}

    aH = math.pow(10, -pH)

    for key, residue in protein_complex.residue_variables.iteritems():
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

            dge, dgd = protein_complex.evaluate_energy_diff_his(hie_residue, hid_residue, labeling,
                                                                 normal_form=True)

            dge *= kln10_T
            dgd *= kln10_T

            dpkad = -math.log10(math.exp(dgd/RT_gas_T))
            dpkae = -math.log10(math.exp(dge/RT_gas_T))

            pkad = modPkaHIP + dpkad
            pkae = modPkaHIP + dpkae

            Gd = -math.log(math.pow(10, pkad))
            Ge = -math.log(math.pow(10, pkae))

            ThetaPEnerNumer = math.exp(-Gd)*aH
            ThetaPEnerDenom = 1.0+math.exp(-(Gd-Ge))+math.exp(-Gd)*aH

            titration_value = ThetaPEnerNumer/ThetaPEnerDenom
            results["HIS", chain, location] = titration_value

        else:
            #Do not HIS stuff
            #energy_diff = protonated_energy - depotonated_energy
            energy_diff = protein_complex.evaluate_energy_diff(residue, labeling, normal_form=True)

            exp = -(energy_diff*kln10_T)/RT
            e_exp = math.exp(exp)
            titration_value = e_exp/(1.0+e_exp)
            results[key] = titration_value

    return results





