# CODE REUSED FROM https://github.com/isayevlab/pKa-ANI
# PLEASE CREDIT THE Isayev Lab if you use pKa-ANI to assign titration states!
#
# CITATION:
# Gokcan, H.; Isayev, O. Prediction of Protein pKa with Representation Learning.
# Chemical Science, 2022, 13, 2462–2474. https://doi.org/10.1039/d1sc05610g.

import sys
import os
import sys as sys
import subprocess as sp
import shutil

from parmed.tools import addPDB, parmout
from parmed.amber import AmberParm

# to download pdb file if not exist
from urllib.request import urlopen
from io import StringIO


def make_apo_pdb(pdbin, pdbout):
    Hatoms = ["H", "D"]
    proteinResidues = [
        "ALA",
        "ASN",
        "CYS",
        "GLU",
        "HIS",
        "LEU",
        "MET",
        "PRO",
        "THR",
        "TYR",
        "ARG",
        "ASP",
        "GLN",
        "GLY",
        "ILE",
        "LYS",
        "PHE",
        "SER",
        "TRP",
        "VAL",
        "HID",
        "HIE",
        "HIP",
        "LYN",
        "CYX",
        "CYM",
        "GLH",
        "ASH",
    ]

    altlocs = [" ", "A"]

    sschain = []
    ssres = []
    protein = []
    model_num = "1"

    with open(pdbin, "r+") as fd:
        contents = fd.readlines()

        # first find if there is an SSBOND
        for line in contents:
            if "SSBOND" in line:
                # SSBOND   1 CYS A    6    CYS A  127
                row = line.strip().split()
                if row[0] == "SSBOND":
                    sschain.append(row[3])
                    sschain.append(row[6])
                    ssres.append(row[4])
                    ssres.append(row[7])

        for line in contents:
            row = line.strip().split()

            if row[0] == "MODEL":
                model_num = row[-1].strip()
                if model_num == "1":
                    protein.append(line)

            if model_num == "1":
                # ATOM      1  N   MET A   1      14.274  15.403  10.777  1.00  0.00           N
                if (
                    row[0] == "ATOM"
                    and row[3] in proteinResidues
                    and row[-1] not in Hatoms
                    and line[16] in altlocs
                ):

                    newline = line[0:16] + " " + line[16 + 1 :]

                    ch = line[21]
                    rno = line[22:27].split()[0]
                    # make sure we have CYX if exists
                    for i, r in enumerate(ssres):
                        if rno == r and ch == sschain[i]:
                            newline = line[0:16] + " " + "CYX" + line[20:]

                    if row[3] == "HIS":
                        newline = line[0:16] + " " + "HIE" + line[20:]

                    protein.append(newline)

                    if row[0] == "TER":
                        protein.append("TER")

            if row[0] == "END":
                protein.append(line)

    with open(pdbout, "w") as f:
        for line in protein:
            f.write(line)


def tleap_vacuum(pdbfile):

    outname = pdbfile.rsplit(".", 1)[0]
    topfile = outname + "_vac.prmtop"
    crdfile = outname + "_vac.rst"

    model = "protein = loadpdb " + pdbfile + "\n"
    savetop = "saveamberparm protein " + topfile + " " + crdfile + "\n"

    if os.path.exists("tleap_vacuum.in"):
        os.remove("tleap_vacuum.in")

    # Append tleap_vacuum.in
    in_file = open("tleap_vacuum.in", "a")
    in_file.write("source leaprc.protein.ff14SB\n")
    in_file.write(model)
    in_file.write(savetop)
    in_file.write("quit\n")
    in_file.close()

    # Generate vacuum files
    p = sp.call("tleap -s -f tleap_vacuum.in > tleap_vacuum.log", shell=True)

    # Find out the charge of model for neutralizing
    # similar to AMBER APR example
    charge = 0
    f = open("tleap_vacuum.log", "r")
    for line in f:
        if "The unperturbed charge of the unit" in line:
            splitline = line.split()
            charge = splitline[6].strip("(").strip(")")
    f.close()

    # clean files
    os.remove("leap.log")
    os.remove("tleap_vacuum.log")
    os.remove("tleap_vacuum.in")

    outfiles = [topfile, crdfile]

    return outfiles, charge


def parm_top(pdbfile, topfile):

    base_name = topfile.rsplit(".", 1)[0]
    newtop = base_name + "_parmed.top"
    if os.path.exists(newtop):
        os.remove(newtop)

    parm = AmberParm(topfile)
    action = addPDB(parm, pdbfile)
    action.execute()
    action = parmout(parm, newtop)
    action.execute()
    return newtop


def add_missing_atoms(pdbin, pdbout):

    base_name = pdbin.rsplit(".", 1)[0]

    # add missing atoms
    outfiles, charge = tleap_vacuum(pdbin)

    # add pdb information to topology with parm
    newtop = parm_top(pdbin, outfiles[0])

    # generate pdb file with ambpdb
    pdbcall = (
        "ambpdb -ext -p " + newtop + " -c " + outfiles[1] + " > " + pdbout
    )
    p = sp.call(pdbcall, shell=True)

    # clean files
    os.remove(outfiles[0])

    return newtop, outfiles[1]


def run_sander_min(top, rst, pdbout):
    in_file = open("min.in", "a")
    in_file.write("MIN\n")
    in_file.write(" &cntrl\n")
    in_file.write("  imin   = 1,\n")
    in_file.write("  maxcyc = 500,\n")
    in_file.write("  ncyc   = 250,\n")
    in_file.write("  ntb    = 0,\n")
    in_file.write("  igb    = 0,\n")
    in_file.write("  cut    = 12\n")
    in_file.write(" /\n")
    in_file.close()

    base = top.rsplit(".", 1)[0]
    outf = base + "_min.out"
    rnow = base + "_min.rst"
    crd = base + "_min.crd"

    sander = "$AMBERHOME/bin/sander"
    mincall = (
        sander
        + " -O -i min.in -o %s -p %s -c %s -r %s -x %s -ref %s"
        % (outf, top, rst, rnow, crd, rst)
    )

    p = sp.call(mincall, shell=True)

    pdbcall = "ambpdb -ext -p " + top + " -c " + rnow + " > " + pdbout
    p = sp.call(pdbcall, shell=True)

    if os.path.exists(outf):
        os.remove(outf)
    if os.path.exists(rnow):
        os.remove(rnow)
    if os.path.exists(crd):
        os.remove(crd)
    if os.path.exists("mdinfo"):
        os.remove("mdinfo")
    if os.path.exists("min.in"):
        os.remove("min.in")


def prep_pdb(inputpdb):

    base_name = inputpdb.rsplit(".", 1)[0]
    print("Preparing %s for pKa calculations" % base_name)
    pdbin = base_name + ".pdb"
    pdb_cp = base_name + "_0.pdb"
    shutil.copyfile(pdbin, pdb_cp)

    pdbout0 = base_name + "_clean.pdb"
    make_apo_pdb(pdbin, pdbout0)
    pdbout1 = base_name + "_H.pdb"
    top, rst = add_missing_atoms(pdbout0, pdbout1)
    pdbout2 = base_name + ".pdb"
    run_sander_min(top, rst, pdbout2)
    if os.path.exists(top):
        os.remove(top)
    if os.path.exists(rst):
        os.remove(rst)
    os.remove(pdbout0)
    os.remove(pdbout1)
