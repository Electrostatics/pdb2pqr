"""
    charmm.py - Driver for creating charmm data files for use with PDB2PQR

    This module takes parameter files from the charmm forcefield and creates
    a single CHARMM.DAT file containing the per atom charge and
    Van der Waals radii in the following order:

        Residue Name    Atom Name   Charge   Radius  Epsilon
"""

__date__ = "2 September 2004"
__author__ = "Todd Dolinsky"

import sys,string
CHARGEMAP = "top_all27_prot_na.rtf"
RADIUSMAP = "par_all27_prot_na.prm"
OUTFILE = "CHARMM.DAT"

header = """# Charmm Parameters from CHARMM22
#
# Radius and epsilon from par_all27_prot_na.prm
# Charges from par_all27_prot_na.prm
#
# V(Lennard-Jones) = Eps,i,j[(Rmin,i,j/ri,j)**12 - 2(Rmin,i,j/ri,j)**6]
#
# epsilon: kcal/mole, Eps,i,j = sqrt(eps,i * eps,j)
# Rmin/2: A, Rmin,i,j = Rmin/2,i + Rmin/2,j
#
"""


def makecharmm():
    """
        Create a CHARMM.DAT file for use with PDB2PQR
    """
    
    # Get Radii for each Group and store in Cache
    # Stored in form ID | Radius | Epsilon

    cache = {}
    count = 0
    radfile = open(RADIUSMAP, "r")

    line = radfile.readline()
    while string.find(line, "!V(Lennard-Jones)") == -1:
        line = radfile.readline()

    while string.find(line,"!DUM") != 0:
        if string.find(line,"!") != 0 and string.find(line," ") != 0 and string.find(line,"\t") != 0:
            entry = string.split(line)

            if len(entry) > 2:
                count += 1
                group = entry[0]
                rad = entry[3]
                eps = entry[2]
                cache[group] = [rad, eps]

        line = radfile.readline()

    # Special Case (in one file it is ON2b, while in the other it is ON2B):

    cache["ON2B"] = [cache["ON2b"][0], cache["ON2b"][1]]

    radfile.close()


    # Get residues, atoms, and charges

    chgfile = open(CHARGEMAP, "r")
    outfile = open(OUTFILE, "w")
    outfile.write(header)
    residue, atom, group, charge = "","","",""
    
    for line in chgfile.readlines():

        # Determine the Residue
        
        if string.find(line, "RESI") == 0 or string.find(line, "PRES") == 0:
            residue = string.split(line)[1]

        if string.find(line, "ATOM") == 0:
            entry = string.split(line)
            atom = entry[1]
            group = entry[2]
            charge = entry[3]
            
            if cache.has_key(group):
                rad = cache[group][0]
                eps = cache[group][1]
                outfile.write("%s\t%s\t%s\t%s\t%s\n" % \
                              (residue, atom, charge, rad, eps))
                
            else:
                print line
                stdout.write("Group %s not found\n" % group)

    outfile.close()
    chgfile.close()

if __name__ == "__main__": makecharmm()


