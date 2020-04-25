from pdb2pka.graph_cut.protein_complex import ProteinComplex
from pdb2pka.graph_cut.titration_curve import get_titration_curves
from sys import version_info
import logging


_LOGGER = logging.getLogger(__name__)


def _add_state_pair(pc, inter_avg,
                    group1_type, group1_chain, group1_loc, group1_state,
                    group2_type, group2_chain, group2_loc, group2_state):

    pc.add_residue(group1_type, group1_chain, group1_loc)
    pc.add_residue(group2_type, group2_chain, group2_loc)

    if (group1_type, group1_chain, group1_loc) != (group2_type, group2_chain, group2_loc):
        instance1 = pc.get_instance(group1_type, group1_chain, group1_loc, group1_state)
        instance2 = pc.get_instance(group2_type, group2_chain, group2_loc, group2_state)

        pc.interaction_energies[instance1, instance2] = inter_avg

        flipped_inter_avg = pc.interaction_energies.get((instance2, instance1))
        if flipped_inter_avg is not None:
            diff = abs(inter_avg - flipped_inter_avg)
            if diff > 0.0:
                _LOGGER.info(group1_type, group1_chain, group1_loc, group1_state)
                _LOGGER.info(group2_type, group2_chain, group2_loc, group2_state)




def create_protein_complex_from_matrix(correct_matrix):
    """Builds a protein complex from a matrix returned from pKaRoutines.correct_matrix"""
    protein_complex = ProteinComplex()

    for pKa1 in correct_matrix:
        for titration1 in correct_matrix[pKa1]:
            for state1 in correct_matrix[pKa1][titration1]:
                for pKa2 in correct_matrix[pKa1][titration1][state1]:
                    for titration2 in correct_matrix[pKa1][titration1][state1][pKa2]:
                        for state2 in correct_matrix[pKa1][titration1][state1][pKa2][titration2]:
                            inter_avg = correct_matrix[pKa1][titration1][state1][pKa2][titration2][state2]
                            res1 = pKa1.residue
                            res2 = pKa2.residue
                            _add_state_pair(protein_complex, inter_avg,
                                            pKa1.pKaGroup.name, res1.chainID, str(res1.resSeq), state1,
                                            pKa2.pKaGroup.name, res2.chainID, str(res2.resSeq), state2)



    return protein_complex

def process_desolv_and_background(protein_complex, pKa):
    """Applies the background and desolvation energies from a pKa to a protein complex"""
    res_type = pKa.pKaGroup.name
    chain = pKa.residue.chainID
    location = str(pKa.residue.resSeq)

    if(version_info >= (3,0)):
        for state, energy in pKa.desolvation.items():
            _process_desolv_or_background_line(protein_complex, res_type, chain, location, state, energy)

        for state, energy in pKa.background.items():
            _process_desolv_or_background_line(protein_complex, res_type, chain, location, state, energy)

    else:
        for state, energy in pKa.desolvation.iteritems():
            _process_desolv_or_background_line(protein_complex, res_type, chain, location, state, energy)

        for state, energy in pKa.background.iteritems():
            _process_desolv_or_background_line(protein_complex, res_type, chain, location, state, energy)

def _process_desolv_or_background_line(protein_complex, res_type, chain, location, state_name, energy):
    instance = protein_complex.get_instance(res_type, chain, location, state_name)
    instance.energy += energy

def _create_protein_complex_from_pKa(pKa):
    protein_complex = ProteinComplex()
    res = pKa.residue
    group = pKa.pKaGroup

    #Grab all of the states
    states = []

    for titration1 in group.DefTitrations:
        states.extend(titration1.allstates)

    states.sort()

    for state1 in states:
        for state2 in states:
            _add_state_pair(protein_complex, 0.0,
                            group.name, res.chainID, str(res.resSeq), state1,
                            group.name, res.chainID, str(res.resSeq), state2)

    return protein_complex

def curve_for_one_group(pKa):
    """Roughly equivalent to pka_help.titrate_one_group
       Titrate a single group and return the titration curve for it.
       """
    protein_complex = _create_protein_complex_from_pKa(pKa)
    process_desolv_and_background(protein_complex, pKa)

    protein_complex.simplify()

    curve = get_titration_curves(protein_complex)

    return curve
