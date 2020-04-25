#
# * This library is free software; you can redistribute it and/or
# * modify it under the terms of the GNU Lesser General Public
# * License as published by the Free Software Foundation; either
# * version 2.1 of the License, or (at your option) any later version.
# *
# * This library is distributed in the hope that it will be useful,
# * but WITHOUT ANY WARRANTY; without even the implied warranty of
# * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# * Lesser General Public License for more details.
#

# propka3.0, revision 182                                                                      2011-08-09
# -------------------------------------------------------------------------------------------------------
# --                                                                                                   --
# --                                   PROPKA: A PROTEIN PKA PREDICTOR                                 --
# --                                                                                                   --
# --                              VERSION 3.0,  01/01/2011, COPENHAGEN                                 --
# --                              BY MATS H.M. OLSSON AND CHRESTEN R. SONDERGARD                       --
# --                                                                                                   --
# -------------------------------------------------------------------------------------------------------
#
#
# -------------------------------------------------------------------------------------------------------
# References:
#
#   Very Fast Empirical Prediction and Rationalization of Protein pKa Values
#   Hui Li, Andrew D. Robertson and Jan H. Jensen
#   PROTEINS: Structure, Function, and Bioinformatics 61:704-721 (2005)
#
#   Very Fast Prediction and Rationalization of pKa Values for Protein-Ligand Complexes
#   Delphine C. Bas, David M. Rogers and Jan H. Jensen
#   PROTEINS: Structure, Function, and Bioinformatics 73:765-783 (2008)
#
#   PROPKA3: Consistent Treatment of Internal and Surface Residues in Empirical pKa predictions
#   Mats H.M. Olsson, Chresten R. Sondergard, Michal Rostkowski, and Jan H. Jensen
#   Journal of Chemical Theory and Computation, 7, 525-537 (2011)
# -------------------------------------------------------------------------------------------------------
import sys
from . import lib
revision = 182


def printHeader():
    """
    prints the header section
    """
    str = "%s\n" % (getPropkaHeader())
    str += "%s\n" % (getReferencesHeader())
    str += "%s\n" % (getWarningHeader())

    print(str)


def writePDB(protein, file=None, filename=None, hydrogens=False, options=None):
    """
    Write the residue to the new pdbfile
    """

    if file == None:
        # opening file if not given
        if filename == None:
            filename = "%s.pdb" % (protein.name)
        file = open(filename, 'w')
        if options.verbose:
            print("writing pdbfile %s" % (filename))
        close_file = True
    else:
        # don't close the file, it was opened in a different place
        close_file = False

    numb = 0
    for chain in protein.chains:
        for residue in chain.residues:
            if residue.resName not in ["N+ ", "C- "]:
                for atom in residue.atoms:
                    if hydrogens == False and atom.element == "H":
                        """ don't print """
                    else:
                        numb += 1
                        line = atom.makePDBLine(numb=numb)
                        line += "\n"
                        file.write(line)

    if close_file == True:
        file.close()


def writePQR(protein, label=None, hydrogens=False, options=None):
    """
    Write a quick pqr file for MEAD calculations - quick & dirty for my tests
    """
    import math

    filename = "%s.pqr" % (protein.name)
    solute = open(filename, 'w')
    filename = "%s_env.pqr" % (protein.name)
    environment = open(filename, 'w')
    radii = {'H': 1.00,
             'C': 1.70,
             'N': 1.50,
             'O': 1.40,
             'S': 1.85, }

    # setting center and making ogm-file
    residue = protein.getResidue(label=label)
    x = residue.x
    y = residue.y
    z = residue.z
    str = "%9.3lf%9.3lf%9.3lf     %3d%9.3lf\n" % (x, y, z, 39, 1.000)
    str += "%9.3lf%9.3lf%9.3lf     %3d%9.3lf\n" % (x, y, z, 75, 0.500)
    str += "%9.3lf%9.3lf%9.3lf     %3d%9.3lf\n" % (x, y, z, 147, 0.250)
    filename = "%s.ogm" % (protein.name)
    file = open(filename, 'w')
    file.write(str)
    file.close()

    numb = 0
    for chain in protein.chains:
        for residue in chain.residues:
            if residue.resName not in ["N+ ", "C- "]:
                for atom in residue.atoms:
                    dX = atom.x - x
                    dY = atom.y - y
                    dZ = atom.z - z
                    distance = math.sqrt(dX*dX + dY*dY + dZ*dZ)
                    if distance > 16.0:
                        """ don't print """
                    elif hydrogens == False and atom.element == "H":
                        """ don't print """
                    else:
                        numb += 1
                        pdbline = atom.makePDBLine(numb=numb, chainID=" ")
                        if residue.label == label:
                            if (residue.resName == "ASP" and atom.name in ["OD1", "OD2"]) or (residue.resName == "GLU" and atom.name in ["OE1", "OE2"]):
                                line = "%s %6.2lf%6.2lf\n" % (pdbline[:55], -0.70, radii[atom.element])
                                solute.write(line)
                            elif (residue.resName == "ASP" and atom.name in ["CG"]) or (residue.resName == "GLU" and atom.name in ["CD"]):
                                line = "%s %6.2lf%6.2lf\n" % (pdbline[:55],  0.40, radii[atom.element])
                                solute.write(line)
                            else:
                                line = "%s %6.2lf%6.2lf\n" % (pdbline[:55],  0.00, radii[atom.element])
                                environment.write(line)
                        else:
                            line = "%s %6.2lf%6.2lf\n" % (pdbline[:55], 0.00, radii[atom.element])
                            environment.write(line)

    solute.close()
    environment.close()


def writePKA(protein, filename=None, reference="neutral", direction="folding", options=None):
    """
    Write the pka-file based on the given protein
    """
    verbose = options.verbose if options is not None else False
    if filename == None:
        filename = "%s.pka" % (protein.name)
    file = open(filename, 'w')
    if verbose == True:
        print("writing pkafile %s" % (filename))

    # writing propka header
    str = "%s\n" % (getPropkaHeader())
    str += "%s\n" % (getReferencesHeader())
    str += "%s\n" % (getWarningHeader())

    # writing pKa determinant section
    str += getDeterminantSection(protein)

    # writing pKa summary section
    str += getSummarySection(protein)
    str += "%s\n" % (getTheLine())

    # printing Folding Profile
    str += getFoldingProfileSection(protein, reference=reference, direction=direction, window=[0., 14., 1.0], options=options)

    # printing Protein Charge Profile
    str += getChargeProfileSection(protein)
    str += "\n"

    # now, writing the pka text to file
    file.write(str)

    file.close()


def printTmProfile(protein, reference="neutral", window=[0., 14., 1.], Tm=[0., 0.], Tms=None, ref=None, verbose=False, options=None):
    """
    prints Tm profile
    """
    profile = protein.getTmProfile(reference=reference, grid=[0., 14., 0.1], Tms=Tms, ref=ref, options=options)
    if profile == None:
        str = "Could not determine Tm-profile\n"
    else:
        str = " suggested Tm-profile for %s\n" % (protein.name)
        for (pH, Tm) in profile:
            if pH >= window[0] and pH <= window[1] and (pH % window[2] < 0.01 or pH % window[2] > 0.99*window[2]):
                str += "%6.2lf%10.2lf\n" % (pH, Tm)
        if verbose:
            print(str)


def printResult(protein, verbose=False):
    """
    prints all resulting output from determinants and down
    """
    printPKASection(protein)


def printPKASection(protein, verbose=False):
    """
    prints out the pka-section of the result
    """
    # geting the determinants section
    str = getDeterminantSection(protein)
    if verbose:
        print(str[:-1])

    str = getSummarySection(protein)
    if verbose:
        print(str)


def getDeterminantSection(protein, verbose=False):
    """
    prints out the pka-section of the result
    """
    # getting the same order as in propka2.0
    residue_list = lib.residueList("propka1")
    str = "%s\n" % (getDeterminantsHeader())

    # printing determinants
    for chain in protein.chains:
        for residue_type in residue_list:
            for residue in chain.residues:
                if residue.resName == residue_type:
                    str += "%s\n" % (residue.getDeterminantString())

    # Add a warning in case of coupled residues
    if protein.coupled_residues:
        str += "%s\n" % (getTheLine())
        str += "  Residues that are found to be 'coupled', i.e. titrates together, has been marked by '*' in the above\n"
        str += "  section. Please rerun PropKa with the --display-coupled-residues option for detailed information.\n"

    return str


def getSummarySection(protein, verbose=False):
    """
    prints out the pka-section of the result
    """
    # getting the same order as in propka2.0
    residue_list = lib.residueList("propka1")
    str = "%s\n" % (getSummaryHeader())
    # printing pKa summary
    for chain in protein.chains:
        for residue_type in residue_list:
            for residue in chain.residues:
                if residue.resName == residue_type:
                    str += "%s\n" % (residue.getSummaryString())

    return str


def getFoldingProfileSection(protein, direction="folding", reference="neutral", window=[0., 14., 1.0], verbose=False, options=None):
    """
    returns the protein-folding-profile section
    """
    str = getTheLine()
    str += "\n"
    str += "Free energy of %9s (kcal/mol) as a function of pH (using %s reference)\n" % (direction, reference)

    profile = protein.getFoldingProfile(reference=reference, direction=direction, grid=[0., 14., 0.1], options=options)
    if profile == None:
        str += "Could not determine folding profile\n"
    else:
        for (pH, dG) in profile:
            if pH >= window[0] and pH <= window[1] and (pH % window[2] < 0.05 or pH % window[2] > 0.95):
                str += "%6.2lf%10.2lf\n" % (pH, dG)
        str += "\n"

    pH, dG = protein.getPHopt()
    if pH == None or dG == None:
        str += "Could not determine pH optimum\n"
    else:
        str += "The pH of optimum stability is %4.1lf for which the free energy is%6.1lf kcal/mol at 298K\n" % (pH, dG)

    dG_min, dG_max = protein.getDG80()
    if dG_min == None or dG_max == None:
        str += "Could not determine pH values where the free energy is within 80 %s of maximum\n" % ("%")
    else:
        str += "The free energy is within 80 %s of maximum at pH %4.1lf to %4.1lf\n" % ("%", dG_min, dG_max)

    pH_min, pH_max = protein.getStabilityRange()
    if pH_min == None or pH_max == None:
        str += "Could not determine where the free energy is positive\n\n"
    else:
        str += "The free energy is positive in the range %4.1lf - %4.1lf\n\n" % (dG_min, dG_max)

    return str


def getChargeProfileSection(protein, verbose=False, options=None):
    """
    returns the protein-folding-profile section
    """
    str = "Protein charge of folded and unfolded state as a function of pH\n"

    profile = protein.getChargeProfile(grid=[0., 14., 1.], options=options)
    if profile == None:
        str += "Could not determine charge profile\n"
    else:
        str += "%6s%10s%8s\n" % ("pH", "unfolded", "folded")
        for (pH, Q_pro, Q_mod) in profile:
            str += "%6.2lf%10.2lf%8.2lf\n" % (pH, Q_mod, Q_pro)

    pI_pro, pI_mod = protein.getPI()
    if pI_pro == None or pI_mod == None:
        str += "Could not determine the pI\n\n"
    else:
        str += "The pI is %5.2lf (folded) and %5.2lf (unfolded)" % (pI_pro, pI_mod)

    return str


def writeJackalScapFile(mutationData=None, filename="1xxx_scap.list", options=None):
    """
    writing a scap file for, i.e., generating a mutated protein
    """
    file = open(filename, 'w')

    for chainID, code1, resNumb, code2 in mutationData:
        str = "%s, %d, %s\n" % (chainID, resNumb, code2)
        file.write(str)
    file.close()


def writeScwrlSequenceFile(sequence, filename="x-ray.seq", options=None):
    """
    writing a scwrl sequence file for, e.g.,  generating a mutated protein
    """
    file = open(filename, 'w')

    start = 0
    while len(sequence[start:]) > 60:
        file.write("%s\n" % (sequence[start:start+60]))
        start += 60
    file.write("%s\n" % (sequence[start:]))

    file.close()


def writeWhatIfFile(mutationData=None, pdbfile=None, newfile=None, filename="whatif.sh", options=None):
    """
    writing a scap file for, i.e., generating a mutated protein
    """
    file = open(filename, 'w')

    debump_numbers = ""
    str = ""
    str += "#!/bin/sh\n"
    str += "if [ -f %s ]; then\n  rm %s\nfi\n" % (newfile, newfile)
    str += "whatif <<- EOF\n"
    str += "getmol %s\n" % (pdbfile)
    str += "\n"
    for chainID, code1, resNumb, code2 in mutationData:
        str += "%s %d %s\n" % ("mutate", resNumb, code2)
        debump_numbers += " %d" % (resNumb)
    str += "%s %s 0\n" % ("debump", debump_numbers)
    str += "\n"
    str += "soup\n"
    str += "makmol\n"
    str += "\n"
    str += "%s\n" % (newfile)
    str += "all\n"
    str += "0\n"
    str += "\n"
    str += "fullstp\n"
    str += "y\n"
    str += "EOF\n"
    str += "\n"
    str += "for file in ALTERR.LOG DSSPOUT fort.78 PDBFILE PDBFILE.PDB WHATIF.FIG\n"
    str += "do\n"
    str += "  if [ -f $file ]; then\n"
    str += "    rm $file\n"
    str += "  fi\n"
    str += "done\n"
    str += "\n"

    file.write(str)

    file.close()


# --- various header text --- #


def getPropkaHeader():
    """
    Creates the header
    """
    from datetime import date
    today = date.today()
    str = "propka3.0, revision %s %79s\n" % (revision, today)
    str += "-------------------------------------------------------------------------------------------------------\n"
    str += "--                                                                                                   --\n"
    str += "--                                   PROPKA: A PROTEIN PKA PREDICTOR                                 --\n"
    str += "--                                                                                                   --\n"
    str += "--                                VERSION 1.0,  04/25/2004, IOWA CITY                                --\n"
    str += "--                                             BY HUI LI                                             --\n"
    str += "--                                                                                                   --\n"
    str += "--                               VERSION 2.0,  11/05/2007, IOWA CITY/COPENHAGEN                      --\n"
    str += "--                                BY DELPHINE C. BAS AND DAVID M. ROGERS                             --\n"
    str += "--                                                                                                   --\n"
    str += "--                              VERSION 3.0,  xx/xx/2010, COPENHAGEN                                 --\n"
    str += "--                              BY MATS H.M. OLSSON AND CHRESTEN R. SONDERGARD                       --\n"
    str += "--                                                                                                   --\n"
    str += "-------------------------------------------------------------------------------------------------------\n"
    str += "\n"

    return str


def getReferencesHeader():
    """
    Returns the 'references' part in output file
    """

    str = ""
    str += "-------------------------------------------------------------------------------------------------------\n"
    str += " References:\n"
    str += "\n"
    str += "   Very Fast Empirical Prediction and Rationalization of Protein pKa Values\n"
    str += "   Hui Li, Andrew D. Robertson and Jan H. Jensen\n"
    str += "   PROTEINS: Structure, Function, and Bioinformatics 61:704-721 (2005)\n"
    str += "   \n"
    str += "   Very Fast Prediction and Rationalization of pKa Values for Protein-Ligand Complexes\n"
    str += "   Delphine C. Bas, David M. Rogers and Jan H. Jensen\n"
    str += "   PROTEINS: Structure, Function, and Bioinformatics 73:765-783 (2008)\n"
    str += "   \n"
    str += "   PROPKA3: Consistent Treatment of Internal and Surface Residues in Empirical pKa predictions\n"
    str += "   Mats H.M. Olsson, Chresten R. Sondergard, Michal Rostkowski, and Jan H. Jensen\n"
    str += "   Journal of Chemical Theory and Computation, to be submitted (2010)\n"
    str += "-------------------------------------------------------------------------------------------------------\n"

    return str


def getWarningHeader():
    """
    Returns the 'warning' part in output file
    """

    str = ""
    str += "-------------------------------------------------------------------------------------------------------\n"
    str += " WARNING !\n"
    str += "\n"
    str += "   Propka3.0 is not identical to propka2.0 and does not work with ligands\n"
    str += "-------------------------------------------------------------------------------------------------------\n"

    return str


def getDeterminantsHeader():
    """
    Creates the Determinant header
    """
    str = ""
    str += "---------  -----   ------   ---------------------    --------------    --------------    --------------\n"
    str += "                            DESOLVATION  EFFECTS       SIDECHAIN          BACKBONE        COULOMBIC\n"
    str += " RESIDUE    pKa    BURIED     REGULAR      RE        HYDROGEN BOND     HYDROGEN BOND      INTERACTION\n"
    str += "---------  -----   ------   ---------   ---------    --------------    --------------    --------------\n"

    return str


def getSummaryHeader():
    """
    returns the summary header
    """
    str = getTheLine()
    str += "\n"
    str += "SUMMARY OF THIS PREDICTION\n"
    str += "     RESIDUE    pKa   pKmodel   ligand atom-type"

    return str


def getTheLine():
    """
    draw the line - Johnny Cash would have been proud - or actually Aerosmith!
    """
    str = ""
    for i in range(0, 104):
        str += "-"

    return str
