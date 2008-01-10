"""
    whitespace extension

    Put a white space between the atom name and residue name fields. Having
    no space between these two fields may cause problem when the pqr file is  
    read by APBS.

    Author:  Yong Huang
"""

__date__ = "10 January 2008"
__author__ = "Yong Huang"

from src.utilities import * 
from src.routines import * 

def usage():

    str =  "        --whitespace  :  Put an extra whitespace between \n"
    str += "                         atom name and residue name to \n"
    str += "                         {output-path}-whitespace.pqr \n"
    return str

def whitespace(routines, outroot):
    """
        Put an extra whitespace between atom name and residue name, even if
        this may break strict PDB formatting and cause problems for some
        visualization programs.

        Parameters
            routines:  A link to the routines object
            outroot:   The root of the output name
    """

    outname = outroot + "-whitespace.pqr" 
    file = open(outname, "w")

    routines.write("\nPutting a whitespace between atom name and residue name...\n")
    routines.write("----------------\n")

    protein = routines.protein

    for atom in protein.getAtoms():
        a = str(atom)
        if a[0:4] == 'ATOM':
            b = a[0:16] + ' ' + a[16:]
            file.write("%s\n" % (b))
        elif a[0:6] == 'HETATM':
            b = a[0:16] + ' ' + a[16:]
            file.write("%s\n" % (b))
        else:
            b = a
            file.write("%s\n" % (b))
            continue

    routines.write("\n")
    file.close()