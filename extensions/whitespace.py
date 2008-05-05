"""
    whitespace extension

    Put a white space between the atom name and residue name fields. Having
    no space between these two fields may cause problem when the pqr file is  
    read by APBS.

    Author:  Yong Huang
"""

__date__ = "5 May 2008"
__author__ = "Yong Huang"

from src.utilities import * 
from src.routines import * 

def usage():

    str =  "        --whitespace  :  Put whitespaces between atom name\n"
    str += "                         and residue name, x and y, y and\n"
    str += "                         z to {output-path}-whitespace.pqr \n"
    return str

def whitespace(routines, outroot):
    """
        Put extra whitespaces between atom name and residue name, x and 
        y, y and z; even if this may break strict PDB formatting and 
        cause problems for somevisualization programs.

        Parameters
            routines:  A link to the routines object
            outroot:   The root of the output name
    """

    outname = outroot + "-whitespace.pqr" 
    file = open(outname, "w")

    routines.write("\nPutting whitespaces between atom name and residue name, x and y, y and z...\n")
    routines.write("----------------\n")

    protein = routines.protein

    for atom in protein.getAtoms():
        a = str(atom)
        if a[0:4] == 'ATOM':
            b = a[0:16] + ' ' + a[16:38] + ' ' + a[38:46] + ' ' + a[46:]
            file.write("%s\n" % (b))
        elif a[0:6] == 'HETATM':
            b = a[0:16] + ' ' + a[16:38] + ' ' + a[38:46] + ' ' + a[46:]
            file.write("%s\n" % (b))
        else:
            b = a
            file.write("%s\n" % (b))
            continue

    routines.write("\n")
    file.close()