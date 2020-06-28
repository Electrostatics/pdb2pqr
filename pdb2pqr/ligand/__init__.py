"""Ligand support functions

Jens Erik Nielsen, University College Dublin 2004
"""
import sys
import pandas


assert sys.version_info >= (3, 5)


# TODO - this belongs in a configuration file somewhere other than here.
# PARSE radii data for C, N, O, S, H, Br, F, P are from Sitkoff et al's paper:
# 
#   Sitkoff D, Sharp KA, Honig B. Accurate Calculation of Hydration Free
#   Energies Using Macroscopic Solvent Models. J Phys Chem 98 (7) 1978-88,
#   1994. J. Phys. Chem. 1994, 98, 7, 1978–1988
#
# See also the AMBER mailing list: http://amber.ch.ic.ac.uk/archive/.
# 
# The van der Waals radius is used for chlorine.
PARSE_RADII = {
    "C": 1.70, "N": 1.50, "O": 1.40, "S": 1.85, "H": 1.00, "Br":2.50,
    "F": 1.20, "P": 1.90, "Cl": 1.75}


# Numbers of valence electrons for the groups of the periodic table
VALENCE_BY_GROUP = {1: 1, 2: 2, 13: 3, 14: 4, 15: 5, 16: 6, 17: 7, 18: 8}
# Groups of the periodic table
ELEMENT_BY_GROUP = {
    1: ["H", "Li", "Na", "K", "Rb", "Cs", "Fr"],
    2: ["Be", "Mg", "Ca", "Sr", "Ba", "Ra"],
    13: ["B", "Al", "Ga", "In", "Tl", "Nh"],
    14: ["C", "Si", "Ge", "Sn", "Pb", "Fl"],
    15: ["N", "P", "As", "Sb", "Bi", "Mc"],
    16: ["O", "S", "Se", "Te", "Po", "Lv"],
    17: ["F", "Cl", "Br", "I", "At", "Ts"],
    18: ["He", "Ne", "Ar", "Kr", "Xe", "Rn", "Og"]
}
# Valence electrons by element
VALENCE_BY_ELEMENT = {}
for group, elem_list in ELEMENT_BY_GROUP.items():
    for elem in elem_list:
        VALENCE_BY_ELEMENT[elem] = VALENCE_BY_GROUP[group]

# Numbers of non-bonded electrons for Sybyl-type atoms.  Adapted from
# https://onlinelibrary.wiley.com/doi/abs/10.1002/jcc.540100804 (Table I).
NONBONDED_BY_TYPE = {
    "Al": 0, "Br": 6, "C.1": 0, "C.2": 0, "C.3": 0, "C.ar": 0, "Ca": 0, 
    "Cl": 6, "F": 6, "H": 0, "I": 6, "K": 0, "Li": 0, "N.1": 2, "N.2": 2,
    "N.3": 2, "N.4": 0, "N.am": 0, "N.ar": 2, "N.pl3": 0, "Na": 0,
    "O.2": 4, "O.3": 4, "P.3": 0, "S.2": 4, "S.3": 4, "S.o": 2, "S.o2": 0,
    "Si": 0, "O.co2": 4.5
}
