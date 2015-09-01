import networkx as nx
from itertools import combinations

class ProteinGraph(object):
    def __init__(self, protein_complex):
        self.DG = nx.DiGraph()
        self.pc = protein_complex

        self._build_nodes()

    def _build_nodes(self):

        self.DG.add_node("S")
        self.DG.add_node("T")

        for key in self.pc.residue_variables:
            self.DG.add_node(key+("PROTONATED",))
            self.DG.add_node(key+("DEPROTONATED",))

    def update_edges(self):
        #Wipe existing edges.
        self.DG.succ.clear()
        self.DG.pred.clear()

        for key, v in self.pc.residue_variables.iteritems():
            prot_instance = v.instances["PROTONATED"]
            prot_capacity = prot_instance.energyNF / 2.0
            prot_node = key+("PROTONATED",)

            deprot_instance = v.instances["DEPROTONATED"]
            deprot_capacity = deprot_instance.energyNF / 2.0
            deprot_node = key+("DEPROTONATED",)

            if prot_capacity != 0.0:
                self.DG.add_edge("S", prot_node, capacity=prot_capacity)
                self.DG.add_edge(deprot_node, "T", capacity=prot_capacity)

            if deprot_capacity != 0.0:
                self.DG.add_edge("S", deprot_node, capacity=deprot_capacity)
                self.DG.add_edge(prot_node, "T", capacity=deprot_capacity)

        for v, w in combinations(self.pc.residue_variables.iteritems(),2):
            v_key, v_residue = v
            w_key, w_residue = w

            v_prot_instance = v_residue.instances["PROTONATED"]
            v_prot_node = v_key+("PROTONATED",)

            v_deprot_instance = v_residue.instances["DEPROTONATED"]
            v_deprot_node = w_key+("DEPROTONATED",)

            w_prot_instance = w_residue.instances["PROTONATED"]
            w_prot_node = w_key+("PROTONATED",)

            w_deprot_instance = w_residue.instances["DEPROTONATED"]
            w_deprot_node = w_key+("DEPROTONATED",)

            capacity = self.pc.normalized_interaction_energies[v_deprot_instance, w_deprot_instance] / 2.0
            self.DG.add_edge(v_prot_node, w_deprot_node, capacity=capacity)
            self.DG.add_edge(w_prot_node, v_deprot_node, capacity=capacity)

            capacity = self.pc.normalized_interaction_energies[v_prot_instance, w_deprot_instance] / 2.0
            self.DG.add_edge(v_prot_node, w_prot_node, capacity=capacity)
            self.DG.add_edge(w_deprot_node, v_deprot_node, capacity=capacity)

            capacity = self.pc.normalized_interaction_energies[v_deprot_instance, w_prot_instance] / 2.0
            self.DG.add_edge(v_deprot_node, w_deprot_node, capacity=capacity)
            self.DG.add_edge(w_prot_node, v_prot_node, capacity=capacity)

            capacity = self.pc.normalized_interaction_energies[v_prot_instance, w_prot_instance] / 2.0
            self.DG.add_edge(v_deprot_node, w_prot_node, capacity=capacity)
            self.DG.add_edge(w_deprot_node, v_prot_node, capacity=capacity)

    def get_cut(self):
        cut_value, partition = nx.minimum_cut(self.DG, "S", "T")
        s_nodes, t_nodes = partition

        return cut_value, set(s_nodes), set(t_nodes)

    def get_labeling_from_cut(self, s_nodes, t_nodes):
        labeling = {}
        uncertain = []
        for key, v in self.pc.residue_variables.iteritems():
            prot_node = key+("PROTONATED",)
            deprot_node = key+("DEPROTONATED",)

            if prot_node in s_nodes and deprot_node in t_nodes:
                labeling[v] = v.instances["DEPROTONATED"]
            elif deprot_node in s_nodes and prot_node in t_nodes:
                labeling[v] = v.instances["PROTONATED"]
            else:
                #Inconclusive
                uncertain.append(v)



        return labeling, uncertain





