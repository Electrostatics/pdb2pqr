
from itertools import combinations, product, permutations

from pprint import pprint, pformat
import sys

titratable_residue_data = {
                           "ARG": {"model_pka":13.0, "ionizable": 1, "tautomers":{"deprotonated":("1+2+3+4",), "protonated":("1+2+3+4+5",)}},
                           "ASP": {"model_pka": 3.9, "ionizable":-1, "tautomers":{"protonated":("1", "2", "3", "4"), "deprotonated":("0",)}},
                           "GLU": {"model_pka": 4.1, "ionizable":-1, "tautomers":{"protonated":("1", "2", "3", "4"), "deprotonated":("0",)}},
                           "LYS": {"model_pka":10.4, "ionizable": 1, "tautomers":{"protonated":("1",), "deprotonated":("0",)}},
                           "TYR": {"model_pka": 9.6, "ionizable": 1, "tautomers":{"protonated":("1",), "deprotonated":("0",)}},
                           "CTR": {"model_pka": 3.2, "ionizable":-1, "tautomers":{"protonated":("1", "2", "3", "4"), "deprotonated":("0",)}},
                           "NTR": {"model_pka": 8.3, "ionizable": 1, "tautomers":{"protonated":("1", "2"), "deprotonated":("1+2",)}},
                           "HIS": {"model_pka": 6.6, "ionizable": 1, "tautomers":{"protonated":("1", "2"), "deprotonated":("1+2",)}},
                          }


all_tautomers = set()
for residue_data in titratable_residue_data.itervalues():
    all_tautomers.update(residue_data["tautomers"]["deprotonated"])
    all_tautomers.update(residue_data["tautomers"]["protonated"])

state_to_tautomer = {}
for state in ("ASH1c", "GLH1c", "LYS", "TYR", "HSD", "H3", "CTR01c"):
    state_to_tautomer[state] = "1"
for state in ("ASH1t", "GLH1t", "HSE", "H2", "CTR01t"):
    state_to_tautomer[state] = "2"
for state in ("ASH2c", "GLH2c", "CTR02c"):
    state_to_tautomer[state] = "3"
for state in ("ASH2t", "GLH2t", "CTR02t"):
    state_to_tautomer[state] = "4"
for state in ("ASP", "GLU", "LYS0", "TYR-", "CTR-"):
    state_to_tautomer[state] = "0"

state_to_tautomer["ARG0"] = "1+2+3+4"
state_to_tautomer["ARG"] = "1+2+3+4+5"

for state in ("HSP", "H3+H2"):
    state_to_tautomer[state] = "1+2"

class ResidueInstance(object):
    def __init__(self, protonated, name, energy = 0.0):
        self.name = name
        self.protonated = protonated
        self.energy = energy
        self.energyNF = None

    def __repr__(self):
        return self.name
        #return "{}: {}, {}, {}".format(self.name, self.protonated, self.energy, self.energyNF)

    def __eq__(self, other):
        return self._get_self_comp_tuple() == other._get_self_comp_tuple()
    def __ne__(self, other):
        return self._get_self_comp_tuple() != other._get_self_comp_tuple()
    def __gt__(self, other):
        return self._get_self_comp_tuple() > other._get_self_comp_tuple()
    def __lt__(self, other):
        return self._get_self_comp_tuple() < other._get_self_comp_tuple()
    def __ge__(self, other):
        return self._get_self_comp_tuple() >= other._get_self_comp_tuple()
    def __le__(self, other):
        return self._get_self_comp_tuple() <= other._get_self_comp_tuple()

    def _get_self_comp_tuple(self):
        return (self.name, self.energy)


class ResidueVariable(object):
    def __init__(self, name):
        self.name = name
        self.instances = {}

    def get_instance(self, state):
        if state not in all_tautomers:
            state = state_to_tautomer.get(state)
        return self.instances.get(state)

    def get_prot_and_deprot_instances(self):
        #We don't want the placeholder instances
        prot = [value for key, value in self.instances.iteritems() if "PROTONATED" not in key and value.protonated]
        deprot = [value for key, value in self.instances.iteritems() if "PROTONATED" not in key and not value.protonated]

        return prot, deprot

    def __repr__(self):
        return self.name

    def __eq__(self, other):
        return self.name == other.name
    def __ne__(self, other):
        return self.name != other.name
    def __gt__(self, other):
        return self.name > other.name
    def __lt__(self, other):
        return self.name < other.name
    def __ge__(self, other):
        return self.name >= other.name
    def __le__(self, other):
        return self.name <= other.name


class ProteinComplex(object):
    def __init__(self, T=300):
        self.T = T
        self.normalized_constant_energy = 0.0
        self.residue_variables = {}
        self.interaction_energies = {}
        self.normalized_interaction_energies = None

    def add_residue(self, residue_type, chain, location):
        """Adds a residue to the protein if needed.
           Automatically creates all instances needed before simplification."""
        residue_tuple = (residue_type, chain, location)
        if residue_tuple not in self.residue_variables:
            residue_data = titratable_residue_data.get(residue_type)
            base_name = '_'.join(residue_tuple)+'_'
            res_var = ResidueVariable('_'.join(residue_tuple))

            tautomers = residue_data["tautomers"]
            model_pka = residue_data["model_pka"]


            for tautomer in tautomers["deprotonated"]:
                instance = ResidueInstance(False, base_name+tautomer)
                #instance.energy = 0
                res_var.instances[tautomer] = instance

            for tautomer in tautomers["protonated"]:
                instance = ResidueInstance(True, base_name+tautomer)
                instance.energy = -model_pka
                res_var.instances[tautomer] = instance

            if residue_type != "HIS":
                #Placeholder instances for later consolidation
                #HIS is split later so we don't add these as they are not needed.
                instance = ResidueInstance(False, base_name+"DEPROTONATED")
                #instance.energy = 0
                res_var.instances["DEPROTONATED"] = instance
                instance = ResidueInstance(True, base_name+"PROTONATED")
                #instance.energy = -model_pka
                res_var.instances["PROTONATED"] = instance

            self.residue_variables[residue_tuple] = res_var


    def get_residue(self, name, chain, location):
        return self.residue_variables.get((name, chain, location))

    def get_instance(self, residue_type, chain, location, state):
        res_var = self.residue_variables.get((residue_type, chain, location))
        if res_var is None:
            return None
        return res_var.get_instance(state)

    def drop_interaction_pairs(self, pair_list):
        for instance1, instance2 in pair_list:
            del self.interaction_energies[instance1, instance2]
            del self.interaction_energies[instance2, instance1]

    def get_interaction_combinations(self, pair1, pair2):
        product_list = [(x,y) for x,y in product(pair1, pair2)]
        return product_list

    def add_interaction_energy_pair(self, instance1, instance2, energy, normalized=False):
        ie = self.normalized_interaction_energies if normalized else self.interaction_energies
        ie[instance1, instance2] = energy
        ie[instance2, instance1] = energy


    def simplify(self):
        """Simplify the instances into protonated and deprotonated.
           Also divides HIS into HID and HIE."""
        self.update_special_instance_energies()
        self.consolidate()

        self.divide_his()

        #pprint(self.interaction_energies)

    def update_special_instance_energies(self):
        """After we have loaded the backgroind and desolvation files we need to
        update the consolidated instances to use the minimum energy of the
        instances they are meant to replace."""
        for key, residue in self.residue_variables.iteritems():
            name = key[0]
            if name == 'HIS':
                continue
            prot_consolidated = residue.instances["PROTONATED"]
            deprot_consolidated = residue.instances["DEPROTONATED"]
            protonated, deprotonated = residue.get_prot_and_deprot_instances()
            #Update consolidated instance energy
            prot_consolidated.energy = min(x.energy for x in protonated)
            deprot_consolidated.energy = min(x.energy for x in deprotonated)

    def consolidate(self):
        """Each residue has multiple protonated and deprotonated states. Here
           we consolidate those into two states for each residue, PROT and DEPROT.
           We take minimums of energies between states in each class. For example,
           assume we have two amino acids, A and B, where A has protonated states 1, 2, 3
           and deprotonated state 4, and B has protonated states 1, 2, and deprotonated
           state 3. Then
           E(A_PROT, B_PROT) = min{E(A1,B1), E(A1,B2), E(A2,B1), E(A2,B2), E(A3,B1), E(A3,B2)},
           E(A_PROT, B_DEPROT) = min{E(A1,B3), E(A2,B3), E(A3,B3)},
           E(A_DEPROT, B_PROT) = min{E(A4,B1), E(A4,B2)}, and
           E(A_DEPROT, B_DEPROT) = E(A4,B3).
           We do not deal with HIS here, it is kept in its 3 states for now.

           After this is finished all unused interaction energies will be removed from
           self.interaction_energies."""

        handled_interaction_pairs = set()

        for v, w in combinations(self.residue_variables.iteritems(),2):
            v_key, v_residue = v
            w_key, w_residue = w
            v_name = v_key[0]
            w_name = w_key[0]

            #Skip HIS for now
            if v_name == 'HIS' or w_name == 'HIS':
                continue

            v_protinated, v_unprotonated = v_residue.get_prot_and_deprot_instances()
            v_prot_consolidated = v_residue.instances["PROTONATED"]
            v_deprot_consolidated = v_residue.instances["DEPROTONATED"]

            w_protinated, w_unprotonated = w_residue.get_prot_and_deprot_instances()
            w_prot_consolidated = w_residue.instances["PROTONATED"]
            w_deprot_consolidated = w_residue.instances["DEPROTONATED"]

            #For every pairing of v and w (protonated and unprotonated) find the
            # minimum interaction energy and update the interaction map accordingly.
            v_stuff = ((v_protinated, v_prot_consolidated), (v_unprotonated, v_deprot_consolidated))
            w_stuff = ((w_protinated, w_prot_consolidated), (w_unprotonated, w_deprot_consolidated))

            for v_product, w_product in product(v_stuff, w_stuff):
                v_instances, v_consolidated = v_product
                w_instances, w_consolidated = w_product

                energies = []
                #Find all of the pairs
                for v_instance, w_instance in product(v_instances, w_instances):
                    energy = self.interaction_energies[v_instance, w_instance]
                    energies.append(energy)
                    #Mark for deletion.
                    handled_interaction_pairs.add((v_instance, w_instance))

                min_energy = min(energies)

                self.add_interaction_energy_pair(v_consolidated, w_consolidated, min_energy)
        #Now handle HIS.
        for v, w in permutations(self.residue_variables.iteritems(),2):
            his_key, his_residue = v
            other_key, other_residue = w
            his_name = his_key[0]
            other_name = other_key[0]

            #We only care about his this pass
            if his_name != 'HIS':
                continue
            #HIS - HIS is already correct for what this pass does.
            if  other_name == 'HIS':
                continue

            other_protinated, other_unprotonated = other_residue.get_prot_and_deprot_instances()
            other_prot_consolidated = other_residue.instances["PROTONATED"]
            other_deprot_consolidated = other_residue.instances["DEPROTONATED"]

            his_protinated, his_unprotonated = his_residue.get_prot_and_deprot_instances()
            his_stuff = his_protinated + his_unprotonated
            other_stuff = ((other_protinated, other_prot_consolidated), (other_unprotonated, other_deprot_consolidated))

            #For every pairing of a HIS instance and another non HIS residue find the
            # minimum interaction energy and update the interaction map accordingly.
            for his_instance, other_product in product(his_stuff, other_stuff):
                other_instances, other_consolidated = other_product

                energies = []

                for other_instance in other_instances:
                    energy = self.interaction_energies[his_instance, other_instance]
                    energies.append(energy)
                    handled_interaction_pairs.add((his_instance, other_instance))

                min_energy = min(energies)

                self.add_interaction_energy_pair(his_instance, other_consolidated, min_energy)

        #Clean out unused interaction energies
        self.drop_interaction_pairs(handled_interaction_pairs)

        for key, residue in self.residue_variables.iteritems():
            name = key[0]
            if name == 'HIS':
                continue

            residue.instances = dict((k,v) for k,v in residue.instances.iteritems() if "PROTONATED" in k)


    def divide_his(self):
        """Here we divide HIS into two residues - HID, HIE - each with half the pKa value. We
           have to set interaction energies between HIS and other residues and for HIS-HIS
           interactions. Do this based on the values given in the paper"""

        to_drop_his = set()
        new_to_old_his = {}
        handled_interaction_pairs = set()

        items = list(self.residue_variables.iteritems())

        for key, residue in items:
            name = key[0]
            if name != 'HIS':
                continue

            old_instance_1 = residue.get_instance('1')
            old_instance_2 = residue.get_instance('2')
            old_instance_12 = residue.get_instance('1+2')

            base_name = '_'+'_'.join(key[-2:])+'_'
            to_drop_his.add(key)

            hid = ResidueVariable("HId"+base_name)
            name = "HId"+base_name+"PROTONATED"
            hid_prot = ResidueInstance(True, name, energy=old_instance_1.energy)
            name = "HId"+base_name+"DEPROTONATED"
            hid_deprot = ResidueInstance(True, name)
            hid.instances["PROTONATED"] = hid_prot
            hid.instances["DEPROTONATED"] = hid_deprot

            res_tuple = ("HId",)+key[-2:]
            self.residue_variables[res_tuple] = hid

            new_to_old_his[hid] = residue

            hie = ResidueVariable("HIe"+base_name)
            name = "HIe"+base_name+"PROTONATED"
            hie_prot = ResidueInstance( True, name, energy=old_instance_2.energy)
            name = "HIe"+base_name+"DEPROTONATED"
            hie_deprot = ResidueInstance(True, name)
            hie.instances["PROTONATED"] = hie_prot
            hie.instances["DEPROTONATED"] = hie_deprot

            res_tuple = ("HIe",)+key[-2:]
            self.residue_variables[res_tuple] = hie

            new_to_old_his[hie] = residue

            energy = old_instance_12.energy - old_instance_1.energy - old_instance_2.energy
            self.add_interaction_energy_pair(hid_prot, hie_prot, energy)
            self.add_interaction_energy_pair(hid_prot, hie_deprot, 0.0)
            self.add_interaction_energy_pair(hid_deprot, hie_prot, 0.0)
            self.add_interaction_energy_pair(hid_deprot, hie_deprot, sys.float_info.max)

        for key in to_drop_his:
            del self.residue_variables[key]

        #This loop is order dependent for HIS/HIS interactions. That is why some of the code paths that
        #appear as though they should be duplicate or symmetrical are not
        #and why we compare chain locations.
        # (In case you were wondering why HIe <-> HId is not the same as
        #  HId <-> HIe.)
        for v, w in product(self.residue_variables.iteritems(), repeat=2):
            his_key, his_residue = v
            other_key, other_residue = w
            his_name, his_chain, his_location = his_key
            other_name, other_chain, other_location = other_key

            if his_name not in ('HId', 'HIe'):
                continue

            his_location = int(his_location)
            other_location = int(other_location)

            his_prot = his_residue.instances["PROTONATED"]
            his_deprot = his_residue.instances["DEPROTONATED"]

            old_his = new_to_old_his[his_residue]
            old_his_instance_1 = old_his.get_instance('1')
            old_his_instance_2 = old_his.get_instance('2')
            old_his_instance_12 = old_his.get_instance('1+2')

            other_prot = other_residue.instances["PROTONATED"]
            other_deprot = other_residue.instances["DEPROTONATED"]


            if other_name == 'HId':
                if his_chain > other_chain or his_location >= other_location:
                    continue

                old_other = new_to_old_his[other_residue]
                old_other_instance_1 = old_other.get_instance('1')
                old_other_instance_2 = old_other.get_instance('2')
                old_other_instance_12 = old_other.get_instance('1+2')

                if his_name == 'HIe':
                    energy = self.interaction_energies[old_his_instance_12, old_other_instance_12]

                    self.add_interaction_energy_pair(his_prot, other_prot, energy)

                    self.add_interaction_energy_pair(his_prot, other_deprot, 0.0)

                    self.add_interaction_energy_pair(his_deprot, other_prot, 0.0)

                    energy = (self.interaction_energies[old_his_instance_1, old_other_instance_2] -
                              self.interaction_energies[old_his_instance_1, old_other_instance_12] -
                              self.interaction_energies[old_his_instance_12, old_other_instance_2])

                    self.add_interaction_energy_pair(his_deprot, other_deprot, energy)

                elif his_name == 'HId':
                    self.add_interaction_energy_pair(his_prot, other_prot, 0.0)

                    energy = self.interaction_energies[old_his_instance_12, old_other_instance_2]
                    self.add_interaction_energy_pair(his_prot, other_deprot, energy)

                    energy = (self.interaction_energies[old_his_instance_2, old_other_instance_12] -
                              self.interaction_energies[old_his_instance_12, old_other_instance_12])
                    self.add_interaction_energy_pair(his_deprot, other_prot, energy)

                    energy = self.interaction_energies[old_his_instance_2, old_other_instance_2]
                    self.add_interaction_energy_pair(his_deprot, other_deprot, energy)

                combinations = self.get_interaction_combinations((old_his_instance_1, old_his_instance_2, old_his_instance_12),
                                                                 (old_other_instance_1, old_other_instance_2, old_other_instance_12))

                handled_interaction_pairs.update(combinations)

            elif other_name == 'HIe':
                if his_chain > other_chain or his_location >= other_location:
                    continue

                old_other = new_to_old_his[other_residue]
                old_other_instance_1 = old_other.get_instance('1')
                old_other_instance_2 = old_other.get_instance('2')
                old_other_instance_12 = old_other.get_instance('1+2')

                if his_name == 'HIe':
                    self.add_interaction_energy_pair(his_prot, other_prot, 0.0)

                    energy = (self.interaction_energies[old_his_instance_12, old_other_instance_1] -
                              self.interaction_energies[old_his_instance_12, old_other_instance_12])

                    self.add_interaction_energy_pair(his_prot, other_deprot, energy)

                    energy = self.interaction_energies[old_his_instance_1, old_other_instance_12]
                    self.add_interaction_energy_pair(his_deprot, other_prot, energy)

                    energy = self.interaction_energies[old_his_instance_1, old_other_instance_1]
                    self.add_interaction_energy_pair(his_deprot, other_deprot, energy)

                elif his_name == 'HId':
                    self.add_interaction_energy_pair(his_prot, other_prot, 0.0)
                    self.add_interaction_energy_pair(his_prot, other_deprot, 0.0)
                    self.add_interaction_energy_pair(his_deprot, other_prot, 0.0)

                    energy = (self.interaction_energies[old_his_instance_2, old_other_instance_1] +
                              self.interaction_energies[old_his_instance_12, old_other_instance_12] -
                              self.interaction_energies[old_his_instance_2, old_other_instance_12] -
                              self.interaction_energies[old_his_instance_12, old_other_instance_1])
                    self.add_interaction_energy_pair(his_deprot, other_deprot, energy)

                combinations = self.get_interaction_combinations((old_his_instance_1, old_his_instance_2, old_his_instance_12),
                                                                 (old_other_instance_1, old_other_instance_2, old_other_instance_12))

                handled_interaction_pairs.update(combinations)

            else:

                if his_name == 'HIe':
                    energy = (self.interaction_energies[old_his_instance_12, other_prot] -
                              self.interaction_energies[old_his_instance_1, other_prot])
                    self.add_interaction_energy_pair(his_prot, other_prot, energy)

                    energy = self.interaction_energies[old_his_instance_2, other_deprot]
                    self.add_interaction_energy_pair(his_prot, other_deprot, energy)

                    self.add_interaction_energy_pair(his_deprot, other_prot, 0.0)

                    energy = (self.interaction_energies[old_his_instance_1, other_deprot] +
                              self.interaction_energies[old_his_instance_2, other_deprot] -
                              self.interaction_energies[old_his_instance_12, other_deprot])
                    self.add_interaction_energy_pair(his_deprot, other_deprot, energy)

                elif his_name == 'HId':
                    energy = self.interaction_energies[old_his_instance_1, other_prot]
                    self.add_interaction_energy_pair(his_prot, other_prot, energy)

                    energy = (self.interaction_energies[old_his_instance_12, other_deprot] -
                              self.interaction_energies[old_his_instance_2, other_deprot])
                    self.add_interaction_energy_pair(his_prot, other_deprot, energy)

                    energy = (self.interaction_energies[old_his_instance_1, other_prot] +
                              self.interaction_energies[old_his_instance_2, other_prot] -
                              self.interaction_energies[old_his_instance_12, other_prot])
                    self.add_interaction_energy_pair(his_deprot, other_prot, energy)

                    self.add_interaction_energy_pair(his_deprot, other_deprot, 0.0)

                combinations = self.get_interaction_combinations((old_his_instance_1, old_his_instance_2, old_his_instance_12),
                                                                 (other_prot, other_deprot))
                handled_interaction_pairs.update(combinations)

        #Clean out unused interaction energies
        self.drop_interaction_pairs(handled_interaction_pairs)


    def evaluate_energy(self, labeling, normal_form = False):
        energy = 0.0
        for v in self.residue_variables.itervalues():
            v_instance = labeling[v]
            energy += v_instance.energyNF if normal_form else v_instance.energy

        ie = self.normalized_interaction_energies if normal_form else self.interaction_energies
        for v, w in combinations(self.residue_variables.itervalues(),2):
            v_instance = labeling[v]
            w_instance = labeling[w]
            energy += ie[v_instance, w_instance]

        if normal_form:
            energy += self.normalized_constant_energy
        return energy

    def instance_interaction_energy(self, instance, labeling, interaction_energy):
        return sum(interaction_energy.get((instance, labeling[x]),0) for x in self.residue_variables.itervalues())

    def evaluate_energy_diff(self, residue, labeling, normal_form = False):
        """Returns the total energy difference between the deprontated and protonated
           states for the supplied residue. (protonated total energy - deprotonated total energy)"""

        ie = self.normalized_interaction_energies if normal_form else self.interaction_energies

        prot_instance = residue.instances["PROTONATED"]
        prot_energy = sum(ie.get((prot_instance, labeling[x]),0.0) for x in self.residue_variables.itervalues())
        prot_energy += prot_instance.energyNF if normal_form else prot_instance.energy

        deprot_instance = residue.instances["DEPROTONATED"]
        deprot_energy = sum(ie.get((deprot_instance, labeling[x]),0.0) for x in self.residue_variables.itervalues())
        deprot_energy += deprot_instance.energyNF if normal_form else deprot_instance.energy

        return prot_energy - deprot_energy


    def evaluate_energy_diff_his(self, hie_residue, hid_residue, labeling, normal_form = False):
        """Returns the total energy differences between
             1. HSE total energy and HSP total energy (HSE total energy - HSP total energy)
             2. HSD total energy and HSP total energy (HSD total energy - HSP total energy)
            in a tuple."""

        labeling_copy = labeling.copy()
        ie = self.normalized_interaction_energies if normal_form else self.interaction_energies

        hie_prot_instance = hie_residue.instances["PROTONATED"]
        hie_deprot_instance = hie_residue.instances["DEPROTONATED"]
        hid_prot_instance = hid_residue.instances["PROTONATED"]
        hid_deprot_instance = hid_residue.instances["DEPROTONATED"]


        labeling_copy[hid_residue] = hid_deprot_instance
        e_energy = sum(ie.get((hie_prot_instance, labeling_copy[x]),0) for x in self.residue_variables.itervalues())
        labeling_copy[hid_residue] = hid_prot_instance
        p_energy = sum(ie.get((hie_prot_instance, labeling_copy[x]),0) for x in self.residue_variables.itervalues())

        hsp_hse_diff = p_energy - e_energy

        labeling_copy[hie_residue] = hie_deprot_instance
        d_energy = sum(ie.get((hid_prot_instance, labeling_copy[x]),0) for x in self.residue_variables.itervalues())
        labeling_copy[hie_residue] = hie_prot_instance
        p_energy = sum(ie.get((hid_prot_instance, labeling_copy[x]),0) for x in self.residue_variables.itervalues())

        hsp_hsd_diff = p_energy - d_energy

        return hsp_hse_diff, hsp_hsd_diff

    def normalize(self, pH):
        self.normalized_interaction_energies = self.interaction_energies.copy()
        self.normalized_constant_energy = 0.0

        for residue in self.residue_variables.itervalues():
            residue.instances["DEPROTONATED"].energyNF = residue.instances["DEPROTONATED"].energy
            residue.instances["PROTONATED"].energyNF = residue.instances["PROTONATED"].energy + pH


        for v, w in combinations(self.residue_variables.iteritems(), 2):
            v_key, v_residue = v
            w_key, w_residue = w
            is_same_his = (v_key[1:] == w_key[1:] and v_key[0] in ("HId", "HIe") and w_key[0] in ("HId", "HIe"))

            if is_same_his:
                v_prot = v_residue.instances["PROTONATED"]
                w_prot = w_residue.instances["PROTONATED"]
                self.normalized_interaction_energies[v_prot, w_prot] -= 2.0*pH
                self.normalized_interaction_energies[w_prot, v_prot] -= 2.0*pH

        for v, w in permutations(self.residue_variables.iteritems(), 2):
            v_key, v_residue = v
            w_key, w_residue = w

            for v_instance in v_residue.instances.itervalues():
                w_instances = w_residue.instances.values()
                min_energy = min(self.normalized_interaction_energies[v_instance, w_instance]
                                 for w_instance in w_instances)
                v_instance.energyNF += min_energy

                if min_energy != sys.float_info.max:
                    for w_instance in w_instances:
                        self.normalized_interaction_energies[v_instance, w_instance] -= min_energy
                        self.normalized_interaction_energies[w_instance, v_instance] -= min_energy

        for residue in self.residue_variables.itervalues():
            min_energy = min(instance.energyNF for instance in residue.instances.itervalues())
            self.normalized_constant_energy += min_energy
            if min_energy != sys.float_info.max:
                for instance in residue.instances.itervalues():
                    instance.energyNF -= min_energy


#         pprint(self.normalized_interaction_energies)
#         assert False
        #print "Normalized constant energy", self.normalized_constant_energy



