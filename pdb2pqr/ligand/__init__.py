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
    "C": 1.70, "N": 1.50, "O": 1.40, "S": 1.85, "H": 1.00, "BR":2.50,
    "F": 1.20, "P": 1.90, "CL": 1.75}


# TODO - this belongs in a configuration file somewhere other than here.
# Bond lengths from
# http://www.chem.swin.edu.au/modules/mod2/bondlen.html
# We should get a better reference
_BOND_LENGTH_DICTS = [
    {"atom1": 'C', "atom2": 'C', "length": 1.54, "type": "single"},
    {"atom1": 'C', "atom2": 'C', "length": 1.34, "type": "double"},
    {"atom1": 'C', "atom2": 'C', "length": 1.20, "type": "triple"},
    {"atom1": 'C', "atom2": 'C', "length": 1.40, "type": "aromatic"},
    {"atom1": 'C', "atom2": 'O', "length": 1.43, "type": "single"},
    {"atom1": 'C', "atom2": 'O', "length": 1.21, "type": "double"},
    {"atom1": 'C', "atom2": 'N', "length": 1.47, "type": "single"},
    {"atom1": 'C', "atom2": 'N', "length": 1.25, "type": "double"},
    {"atom1": 'C', "atom2": 'N', "length": 1.16, "type": "triple"},
    {"atom1": 'C', "atom2": 'N', "length": 1.34, "type": "aromatic"},
    {"atom1": 'N', "atom2": 'N', "length": 1.35, "type": "aromatic"}
]
BOND_LENGTHS = pandas.DataFrame(_BOND_LENGTH_DICTS)
