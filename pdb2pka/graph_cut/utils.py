from pdb2pka.graph_cut.protein_complex import ProteinComplex

def _add_state_pair(pc, inter_avg,
                    group1_type, group1_chain, group1_loc, group1_state,
                    group2_type, group2_chain, group2_loc, group2_state):

    pc.add_residue(group1_type, group1_chain, group1_loc)
    pc.add_residue(group2_type, group2_chain, group2_loc)

    if (group1_type, group1_chain, group1_loc) != (group2_type, group2_chain, group2_loc):
        instance1 = pc.get_instance(group1_type, group1_chain, group1_loc, group1_state)
        instance2 = pc.get_instance(group2_type, group2_chain, group2_loc, group2_state)

        pc.interaction_energies[instance1, instance2] = round(inter_avg,3)

        flipped_inter_avg = pc.interaction_energies.get((instance2, instance1))
        if flipped_inter_avg is not None:
            diff = abs(inter_avg - flipped_inter_avg)
            if diff > 0.0:
                print group1_type, group1_chain, group1_loc, group1_state
                print group2_type, group2_chain, group2_loc, group2_state




def create_protein_complex_from_matrix(correct_matrix):
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
    res_type = pKa.pKaGroup.name
    chain = pKa.residue.chainID
    location = str(pKa.residue.resSeq)
    for state, energy in pKa.desolvation.iteritems():
        _process_desolv_or_background_line(protein_complex, res_type, chain, location, state, energy)

    for state, energy in pKa.desolvation.iteritems():
        _process_desolv_or_background_line(protein_complex, res_type, chain, location, state, energy)

def _process_desolv_or_background_line(protein_complex, res_type, chain, location, state_name, energy):
    instance = protein_complex.get_instance(res_type, chain, location, state_name)
    instance.energy += round(energy,13)