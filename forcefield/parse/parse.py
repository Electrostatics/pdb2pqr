#!/usr/bin/python2 -O

"""
    Make the PDB2PQR parameter file for the PARSE force field

    This module generates the parameter file for the PARSE force field
    using two parameter files - one with the atom charges and another
    with the atom radii. The output file is in the PDB2PQR DAT file format.

    Original PARSE data from:

    Sitkoff, D., Sharp, K.A., and Honig, B.  (1994)
    Accurate Calculation of Hydration Free Energies Using Macroscopic 
    Solvent Models.
    J. Phys. Chem. 98:1978-1988.

    Todd Dolinsky (todd@ccb.wustl.edu)
    Washington University in St. Louis
"""

__date__ = "13 November 2003"
__author__ = "Todd Dolinsky"

import sys
import os
import string

CHARGEPATH = "parseres.crg"
RADPATH = "parseres.siz"
OUTPATH = "PARSE.DAT"

def makeparse():
    """
        Generate the new PARSE.DAT file
    """
    if not os.path.isfile(CHARGEPATH):
        print "Error: Unable to find PARSE charge file!"
        sys.exit(2)

    if not os.path.isfile(RADPATH):
        print "Error: Unable to find PARSE radius file!"
        sys.exit(2)

    outfile = open(OUTPATH,"w")

    cache = {}

    chargefile = open(CHARGEPATH)
    lines = chargefile.readlines()
    for line in lines:
        if line.startswith("!"):
            continue
        words = string.split(line)
        atom = words[0]
        res = words[1]
        charge = words[2]

        key = atom + res
        if key not in cache:
            cache[key] = charge

    chargefile.close()

    radfile = open(RADPATH)
    lines = radfile.readlines()
    for line in lines:
        if line.startswith("!"):
            continue
        words = string.split(line)
        atom = words[0]
        res = words[1]
        rad = words[2]

        key = atom + res
        try:
            charge = cache[key]
            outfile.write("%s\t%s\t%s\t%s\n" % (res, atom, charge, rad))
            
        except KeyError:
            pass

    radfile.close()
    outfile.close()
         

if __name__ == "__main__":
    makeparse()
