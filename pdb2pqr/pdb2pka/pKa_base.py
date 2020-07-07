"""
Helper classes for all the pKa stuff. I need to reorganise these into something
more logic at some point
"""

class pKa:
    """
        The main pKa object
    """
    def __init__(self, residue, group, amb):
        """
            Initialize the pKa object

            Parameters
                residue: The residue object (residue)
                group:   The pKaGroup object associated with the residue
                         (pKaGroup)
                amb:     The associated HydrogenAmbiguity object
                         (HydrogenAmbiguity)
        """
        self.residue = residue
        self.pka_group = group
        self.amb = amb
        self.desolvation = {}
        self.background = {}
        self.interaction_energies = {}
        self.intrinsic_pka = {}
        self.simulated_intrinsic_pka = None
        self.pka = None
        #
        # Unique identifier
        #
        self.unique_id = '%s_%s_%d_TITTYPE:%s' % (residue.name,
                                                  residue.chain_id,
                                                  residue.res_seq,
                                                  group.name)

    def __repr__(self):
        return self.unique_id


class pKaGroup:
    """
    pKaGroup holds the defintion on a single titratable entity. In most cases we will
    only have a single titration taking place in one group, but in some cases
    we might have to have several transitions in one group (e.g. His - -> His 0 -> His +)

    The name, resname and type should probably be in pKaTitrations, but for now I'll leave
    them here to get something working before we do complicated things...
    """


    def __init__(self, name, resname, type_, def_titrations):
        #
        #    Initialize the pKaGroup object
        #
        #    Parameters
        #        name:    The name of the group (string)
        #        resname: The residue name (string)
        #        type:    The type of group, acid or base (string)
        #        def_titrations: A list of DefTitration objects (list)
        #
        #
        self.name = name
        self.resname = resname
        self.type_ = type_
        self.def_titrations = def_titrations

    #
    # ------------------------
    #

    def __str__(self):
        """
            Print the pKa group object for debugging purposes

            Returns
                text:  The pKaGroup information (string)
        """
        text = "Group name:   %s\n" % self.name
        text += "Residue name: %s\n" % self.resname
        text += "Group type:   %s\n" % self.type_
        text += "Transitions:\n"
        for tran in self.def_titrations:
            text += str(tran)
        return text


class DefTitration:
    """
    pKa_Titration holds all the info on a specific titration
    We define a titration as a single group that has a number of
    startstates and a number of endstates which is modelled by a
    single model pKa value
    A single group can have several transitions depending on the
    number of startstates and endstates
    """

    def __init__(self, startstates, endstates, modelpka, name):
        #
        #    Initialize the pKaTransition object
        #
        #    Parameters
        #        startstates: A list of state numbers (list)
        #        endstates:   A list of state numbers (list)
        #        modelpka:    The model pKa associated with this titration
        #                     (float)
        #        transitions: A dictionary of the possible transitions for this group
        #                     (dictionary)
        #        interactions: A dictionary of the interaction energies with all other states
        #                      of all other titrations in the molecule. Only part of them
        #                      will be calculated
        #
        self.residue = None
        self.startstates = startstates
        self.endstates = endstates
        self.allstates = startstates + endstates
        self.modelpka = modelpka
        self.name = name
        #
        # Set transitions
        #
        self.transitions = {}
        count = 0
        for start_s in self.startstates:
            for end_s in self.endstates:
                count += 1
                self.transitions[count] = {'start':start_s, 'end':end_s}
        #
        # Interactions has to be set at a higher level
        #

    def __str__(self):
        """
            Print the pKa Transition object for debugging purposes

            Returns
                text:  The pKaTransition information (string)
        """
        text = "\tStartstates: %s\n" % self.startstates
        text += "\tEndstates:   %s\n" % self.endstates
        text += "\tmodelpKa:    %.1f\n" % self.modelpka
        text += "\tName:        %s\n" % self.name
        return text
