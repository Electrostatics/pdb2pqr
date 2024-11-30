"""Ligand support functions.

.. todo::
   Some of the definitions in this module belong in a configuration file other
   than here.

.. codeauthor:: Jens Erik Nielsen
"""

#: Radii for different atom types.
#: When using these tables, the most specific Sybyl atom type should be used
#: first and then the generic element should be used
RADII = {
    # NOTE - these are not the original PARSE radii but they are the ones
    # included in the previous version of PDB2PKA so I'm preserving them for
    # posterity. There's a claim they came from
    # http://amber.ch.ic.ac.uk/archive/ but that link no longer works.
    "not parse - do not use": {
        "C": 1.70,
        "N": 1.50,
        "O": 1.40,
        "S": 1.85,
        "H": 1.00,
        "Br": 2.50,
        "F": 1.20,
        "P": 1.90,
        "Cl": 1.75,
    },
    #: These are the PARSE radii from Table 4 of
    #: http://doi.org/10.1021/j100058a043
    "parse": {
        "C.1": 2.00,
        "C.2": 2.00,
        "C.3": 2.00,
        "C": 1.70,
        "H": 1.00,
        "O": 1.40,
        "N": 1.50,
        "S": 1.85,
    },
    #: These are the ZAP radii from Table 2 of
    #: http://doi.org/10.1021/jm070549%2B. Bondi radii should be used for
    #: atoms not found in this table.
    "zap9": {
        "C": 1.87,
        "H": 1.10,
        "O.co2": 1.76,
        "N": 1.40,
        "S": 2.15,
        "F": 2.40,
        "Cl": 1.82,
        "I": 2.65,
    },
    #: These are the Bondi radii from Table 2 of
    #: http://doi.org/10.1021/jm070549%2B
    "bondi-zap": {
        "C": 1.7,
        "H": 1.20,
        "O.co2": 1.52,
        "N": 1.55,
        "S": 1.80,
        "F": 1.47,
        "Cl": 1.75,
        "I": 1.98,
    },
    #: These are the Bondi radii from Table I of
    #: http://doi.org/10.1021/j100785a001. NOTE - there are some variations to
    #: the halogens in Table V that we might want to consider in the future.
    "bondi": {
        "H": 1.20,
        "He": 1.40,
        "C": 1.70,
        "N": 1.55,
        "O": 1.52,
        "F": 1.47,
        "Ne": 1.54,
        "Si": 2.10,
        "P": 1.80,
        "S": 1.80,
        "Cl": 1.75,
        "Ar": 1.88,
        "As": 1.85,
        "Se": 1.90,
        "Br": 1.85,
        "Kr": 2.02,
        "Te": 2.06,
        "I": 1.98,
        "Xe": 2.16,
    },
}


#: Numbers of valence electrons for the groups of the periodic table
VALENCE_BY_GROUP = {1: 1, 2: 2, 13: 3, 14: 4, 15: 5, 16: 6, 17: 7, 18: 8}


#: Groups of the periodic table
ELEMENT_BY_GROUP = {
    1: ["H", "Li", "Na", "K", "Rb", "Cs", "Fr"],
    2: ["Be", "Mg", "Ca", "Sr", "Ba", "Ra"],
    13: ["B", "Al", "Ga", "In", "Tl", "Nh"],
    14: ["C", "Si", "Ge", "Sn", "Pb", "Fl"],
    15: ["N", "P", "As", "Sb", "Bi", "Mc"],
    16: ["O", "S", "Se", "Te", "Po", "Lv"],
    17: ["F", "Cl", "Br", "I", "At", "Ts"],
    18: ["He", "Ne", "Ar", "Kr", "Xe", "Rn", "Og"],
}


#: Valence electrons by element
VALENCE_BY_ELEMENT = {}
for group, elem_list in ELEMENT_BY_GROUP.items():
    for elem in elem_list:
        VALENCE_BY_ELEMENT[elem] = VALENCE_BY_GROUP[group]


#: Numbers of non-bonded electrons for Sybyl-type atoms.  Adapted from
#: https://onlinelibrary.wiley.com/doi/abs/10.1002/jcc.540100804 (Table I).
NONBONDED_BY_TYPE = {
    "Al": 0,
    "Br": 6,
    "C.1": 0,
    "C.2": 0,
    "C.3": 0,
    "C.ar": 0,
    "Ca": 0,
    "Cl": 6,
    "F": 6,
    "H": 0,
    "I": 6,
    "K": 0,
    "Li": 0,
    "N.1": 2,
    "N.2": 2,
    "N.3": 2,
    "N.4": 0,
    "N.am": 0,
    "N.ar": 2,
    "N.pl3": 0,
    "Na": 0,
    "O.2": 4,
    "O.3": 4,
    "P.3": 0,
    "S.2": 4,
    "S.3": 4,
    "S.o": 2,
    "S.o2": 0,
    "Si": 0,
    "O.co2": 4.5,
}
