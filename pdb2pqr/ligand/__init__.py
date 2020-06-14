"""Ligand support functions

Jens Erik Nielsen, University College Dublin 2004
"""
import sys
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
