import networkx as nx
from itertools import combinations

class ProteinGraph(object):
    def __init__(self, protein_complex):
        self.pc = protein_complex

    def _build_nodes(self):
        #Ditch existing graph and start over.
        self.DG = nx.DiGraph()
        self.DG.add_node("S")
        self.DG.add_node("T")

        #Create all state nodes.
        for key in self.pc.residue_variables:
            self.DG.add_node(key+("PROTONATED",))
            self.DG.add_node(key+("DEPROTONATED",))

    def update_graph(self):
        """Build a new graph based on the state of self.pc"""
        self._build_nodes()

        #Create edges going in and out of S and T.
        for key, v in self.pc.residue_variables.iteritems():
            prot_instance = v.instances["PROTONATED"]
            prot_capacity = prot_instance.energyNF / 2.0
            prot_node = key+("PROTONATED",)

            deprot_instance = v.instances["DEPROTONATED"]
            deprot_capacity = deprot_instance.energyNF / 2.0
            deprot_node = key+("DEPROTONATED",)

            if prot_capacity != 0.0:
                self.DG.add_edge("S", deprot_node, capacity=prot_capacity)
                self.DG.add_edge(prot_node, "T", capacity=prot_capacity)

            if deprot_capacity != 0.0:
                self.DG.add_edge("S", prot_node, capacity=deprot_capacity)
                self.DG.add_edge(deprot_node, "T", capacity=deprot_capacity)

        #Create all interaction energy edges.
        for p, q in combinations(self.pc.residue_variables.iteritems(),2):
            p_key, p_residue = p
            q_key, q_residue = q

            p_prot_instance = p_residue.instances["PROTONATED"]
            p_prot_node = p_key+("PROTONATED",)

            p_deprot_instance = p_residue.instances["DEPROTONATED"]
            p_deprot_node = p_key+("DEPROTONATED",)

            q_prot_instance = q_residue.instances["PROTONATED"]
            q_prot_node = q_key+("PROTONATED",)

            q_deprot_instance = q_residue.instances["DEPROTONATED"]
            q_deprot_node = q_key+("DEPROTONATED",)

            capacity = self.pc.normalized_interaction_energies[p_deprot_instance, q_deprot_instance] / 2.0

            if capacity != 0.0:
                self.DG.add_edge(p_deprot_node, q_prot_node, capacity=capacity)
                self.DG.add_edge(q_deprot_node, p_prot_node, capacity=capacity)

            capacity = self.pc.normalized_interaction_energies[p_prot_instance, q_deprot_instance] / 2.0

            if capacity != 0.0:
                self.DG.add_edge(p_prot_node, q_prot_node, capacity=capacity)
                self.DG.add_edge(q_deprot_node, p_deprot_node, capacity=capacity)

            capacity = self.pc.normalized_interaction_energies[p_deprot_instance, q_prot_instance] / 2.0

            if capacity != 0.0:
                self.DG.add_edge(p_deprot_node, q_deprot_node, capacity=capacity)
                self.DG.add_edge(q_prot_node, p_prot_node, capacity=capacity)

            capacity = self.pc.normalized_interaction_energies[p_prot_instance, q_prot_instance] / 2.0

            if capacity != 0.0:
                self.DG.add_edge(p_prot_node, q_deprot_node, capacity=capacity)
                self.DG.add_edge(q_prot_node, p_deprot_node, capacity=capacity)

    def get_cut(self):
        """Performs the min cut.
           Returns cut_value, s nodes, t nodes"""
        cut_value, partition = nx.minimum_cut(self.DG, "S", "T")
        s_nodes, t_nodes = partition

        return cut_value, set(s_nodes), set(t_nodes)

    def get_labeling_from_cut(self, s_nodes, t_nodes):
        """Creates a map of residues to instances based on the """
        labeling = {}
        uncertain = []
        for key, v in self.pc.residue_variables.iteritems():
            prot_node = key+("PROTONATED",)
            deprot_node = key+("DEPROTONATED",)

            if prot_node in s_nodes and deprot_node in t_nodes:
                labeling[v] = v.instances["PROTONATED"]
            elif deprot_node in s_nodes and prot_node in t_nodes:
                labeling[v] = v.instances["DEPROTONATED"]
            else:
                #Inconclusive
                uncertain.append(v)



        return labeling, uncertain





