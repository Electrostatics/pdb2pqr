"""
    amber.py - Driver for creating amber data files for use with PDB2PQR

    This module takes parameter files from the amber forcefield and creates
    a single AMBER.DAT file containing the per atom charge and
    Van der Waals radii in the following order:

        Residue Name    Atom Name   Charge   Radius  Epsilon
"""

__date__ = "1 September 2004"
__author__ = "Todd Dolinsky"

import sys,string

CHARGEMAP1 = "all_amino94.lib"
CHARGEMAP2 = "all_aminoct94.lib"
CHARGEMAP3 = "all_aminont94.lib"
CHARGEMAP4 = "all_nucleic94.lib"

RADIUSMAP = "parm99.dat"
OUTFILE = "AMBER.DAT"

def makeamber():
    """
        Create an AMBER.DAT file for use with PDB2PQR
    """
    
    # Get Radii for each Group and store in Cache
    # Stored in form ID | Radius | Epsilon

    cache = {}
    radfile = open(RADIUSMAP, "r")

    line = radfile.readline()
    while string.find(line, "MOD4") == -1:
        line = radfile.readline()

    while string.find(line,"###") != 0:
        entry = string.split(line)
        if len(entry) > 2:
            group = entry[0]
            rad = entry[1]
            eps = entry[2]
            cache[group] = [rad, eps]
            
        line = radfile.readline()

    radfile.close()

    # Special Cases: C and N atoms
    # N   contains NA  N2  N*  NC  NB  N3  NP  NO
    # C   contains C*  CA  CB  CC  CN  CM  CK  CQ  CW  CV  CR  CA  CX  CY  CD

    clist = ["C*","CA","CB","CC","CD","CK","CM","CN","CQ","CR","CR","CV","CW","CX","CY","CZ"]
    nlist = ["NA","N2","N*","NC","NB","N3", "NT", "NP","NO", "NY"]

    # Get residues, atoms, and charges

    filelist = [CHARGEMAP1, CHARGEMAP2, CHARGEMAP3, CHARGEMAP4]

    outfile = open(OUTFILE, "w")    
    outfile.write("# Amber 99 parameters\n")
    
    residue, atom, group, charge = "","","",""
    for file in filelist:
        chgfile = open(file, "r")

        line = chgfile.readline()

        while line != "":

            # Determine the Residue
            
            if string.find(line, ".unit.atoms table") != -1:
                residue = string.split(line,".")[1]
                line = chgfile.readline()
                while string.find(line, "!") == -1:
                    entry = string.split(line)
                    atom = string.split(entry[0],"\"")[1]
                    group = string.split(entry[1],"\"")[1]
                    charge = entry[7]
            
                    if cache.has_key(group):
                        rad = cache[group][0]
                        eps = cache[group][1]
                        outfile.write("%s\t%s\t%s\t%s\t%s\n" % \
                                      (residue, atom, charge, rad, eps))
                    elif group in clist:
                        rad = cache["C"][0]
                        eps = cache["C"][1]
                        outfile.write("%s\t%s\t%s\t%s\t%s\n" % \
                                      (residue, atom, charge, rad, eps))

                    elif group in nlist:
                        rad = cache["N"][0]
                        eps = cache["N"][1]
                        outfile.write("%s\t%s\t%s\t%s\t%s\n" % \
                               (residue, atom, charge, rad, eps))
                    else:
                        stdout.write("Group %s not found\n" % group)

                    line = chgfile.readline()

            line = chgfile.readline()            
        chgfile.close()
        
    # Add charges for water molecules as per
    # J. Chem Phys Vol 79, p. 926 (1983)

    outfile.write("WAT\tHW\t0.417000\t0.0000\t0.0000\n")
    outfile.write("WAT\tOW\t-0.834000\t1.6612\t0.2100\n")

    outfile.close()

if __name__ == "__main__": makeamber()


