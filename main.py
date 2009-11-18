"""
    Driver for PDB2PQR

    This module takes a PDB file as input and performs optimizations
    before yielding a new PDB-style file as output.

    Ported to Python by Todd Dolinsky (todd@ccb.wustl.edu)
    Washington University in St. Louis

    Parsing utilities provided by Nathan A. Baker (baker@biochem.wustl.edu)
    Washington University in St. Louis

	Copyright (c) 2002-2009, Jens Erik Nielsen, University College Dublin; 
	Nathan A. Baker, Washington University in St. Louis; Paul Czodrowski & 
	Gerhard Klebe, University of Marburg

	All rights reserved.

	Redistribution and use in source and binary forms, with or without modification, 
	are permitted provided that the following conditions are met:

		* Redistributions of source code must retain the above copyright notice, 
		  this list of conditions and the following disclaimer.
		* Redistributions in binary form must reproduce the above copyright notice, 
		  this list of conditions and the following disclaimer in the documentation 
		  and/or other materials provided with the distribution.
		* Neither the names of University College Dublin, Washington University in 
		  St. Louis, or University of Marburg nor the names of its contributors may 
		  be used to endorse or promote products derived from this software without 
		  specific prior written permission.

	THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND 
	ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED 
	WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. 
	IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, 
	INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, 
	BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, 
	DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF 
	LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE 
	OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED 
	OF THE POSSIBILITY OF SUCH DAMAGE.

"""

__date__  = "2 September 2009"
__author__ = "Todd Dolinsky, Nathan Baker, Jens Nielsen, Paul Czodrowski, Jan Jensen, Samir Unni, Yong Huang"
__version__ = "1.5"


import string
import sys
import getopt
import os
import time
from src import pdb
from src import utilities
from src import structures
from src import routines
from src import protein
from src.pdb import *
from src.utilities import *
from src.structures import *
from src.definitions import *
from src.forcefield import *
from src.routines import *
from src.protein import *
from src.server import *
from src.hydrogens import *
from src.aconf import *
from StringIO import *

def usage(rc):
    """
        Print usage for this script to stdout.

        Parameters
            rc:  Exit status (int)
    """

    str = "\n"
    str = str + "pdb2pqr  (Version %s)\n" % __version__
    str = str + "\n"
    str = str + "This module takes a PDB file as input and performs\n"
    str = str + "optimizations before yielding a new PDB-style file as\n"
    str = str + "output\n"
    str = str + "\n"
    str = str + "Usage: pdb2pqr.py [options] --ff=<forcefield> <path> <output-path>\n"
    str = str + "    Required Arguments:\n"
    str = str + "        <forcefield>  :  The forcefield to use - currently amber,\n"
    str = str + "                         charmm, parse, tyl06, peoepb and swanson\n"
    str = str + "                         are supported.\n"
    str = str + "        <path>        :  The path to the PDB file or an ID\n"
    str = str + "                         to obtain from the PDB archive\n"
    str = str + "        <output-path> :  The desired output name of the PQR file\n"
    str = str + "                         to be generated\n"
    str = str + "    Optional Arguments:\n"
    str = str + "        --nodebump    :  Do not perform the debumping operation\n"
    str = str + "        --noopt       :  Do not perform hydrogen optimization\n"
    str = str + "        --chain       :  Keep the chain ID in the output PQR file\n" 
    str = str + "        --assign-only :  Only assign charges and radii - do not add\n"
    str = str + "                         atoms, debump, or optimize.\n"
    str = str + "        --clean       :  Do no optimization, atom addition, or\n"
    str = str + "                         parameter assignment, just return the\n"
    str = str + "                         original PDB file in aligned format.\n"
    str = str + "        --ffout=<name>:  Instead of using the standard canonical\n"
    str = str + "                         naming scheme for residue and atom names,\n"
    str = str + "                         use the names from the given forcefield.\n"
    str = str + "        --with-ph=<ph>:  Use propka to calculate pKas and apply them\n"
    str = str + "                         to the molecule given the pH value. Actual\n"
    str = str + "                         PropKa results will be output to \n"
    str = str + "                         <output-path>.propka.\n"
    str = str + "        --apbs-input  :  Create a template APBS input file based on\n"
    str = str + "                         the generated PQR file.  Also creates a Python\n"
    str = str + "                         pickle for using these parameters in other programs.\n"
    str = str + "        --ligand=<path>: Calculate the parameters for the ligand in\n"
    str = str + "                         mol2 format at the given path. Pdb2pka must\n"
    str = str + "                         be compiled\n"
    str = str + "        --whitespace  :  Insert whitespaces between atom name and residue\n"
    str = str + "                         name, between x and y, and between y and z\n"
    str = str + "        --typemap :      Create Typemap output\n"
    str = str + "        --neutraln  :    Make the N-terminus of this protein neutral \n"
    str = str + "                         (default is charged)\n"
    str = str + "        --neutralc  :    Make the C-terminus of this protein neutral \n"
    str = str + "                         (default is charged)\n"
    str = str + "        --verbose (-v):  Print information to stdout\n"
    str = str + "        --help    (-h):  Display the usage information\n"

    # Check to see if there are usage statements from the
    # extensions directory

    extensions = getAvailableExtensions()
    if len(extensions) > 0:
        str = str + "\n    Optional Arguments from Extensions Directory:\n"
        for ext in extensions:
            str += extensions[ext].usage()
    
    str = str + "\n"
    sys.stderr.write(str)
    sys.exit(rc)

def printPQRHeader(atomlist, reslist, charge, ff, warnings, options):
    """
        Print the header for the PQR file

        Parameters:
            atomlist: A list of atoms that were unable to have
                      charges assigned (list)
            reslist:  A list of residues with non-integral charges
                      (list)
            charge:   The total charge on the protein (float)
            ff:       The forcefield name (string)
            warnings: A list of warnings generated from routines (list)
            options:  A dictionary of command lnie options (float)
        Returns
            header:   The header for the PQR file (string)
    """
    header = "REMARK   1 PQR file generated by PDB2PQR (Version %s)\n" % __version__
    header = header + "REMARK   1\n"
    header = header + "REMARK   1 Forcefield Used: %s\n" % ff
    if "ffout" in options:
        header = header + "REMARK   1 Naming Scheme Used: %s\n" % options["ffout"]
    header = header + "REMARK   1\n"
    
    if "ph" in options:
        header = header + "REMARK   1 pKas calculated by propka and assigned using pH %.2f\n" % options["ph"]
        header = header + "REMARK   1\n"

    for warning in warnings:
        header = header + "REMARK   5 " + warning 
    header = header + "REMARK   5\n"
    
    if len(atomlist) != 0:
        header += "REMARK   5 WARNING: PDB2PQR was unable to assign charges\n"
        header += "REMARK   5          to the following atoms (omitted below):\n"
        for atom in atomlist:
            header += "REMARK   5              %i %s in %s %i\n" % \
                      (atom.get("serial"), atom.get("name"), \
                       atom.get("residue").get("name"), \
                       atom.get("residue").get("resSeq"))
        header += "REMARK   5 This is usually due to the fact that this residue is not\n"
        header += "REMARK   5 an amino acid or nucleic acid; or, there are no parameters\n" 
        header += "REMARK   5 available for the specific protonation state of this\n" 
        header += "REMARK   5 residue in the selected forcefield.\n"
        header += "REMARK   5\n"
    if len(reslist) != 0:
        header += "REMARK   5 WARNING: Non-integral net charges were found in\n"
        header += "REMARK   5          the following residues:\n"
        for residue in reslist:
            header += "REMARK   5              %s - Residue Charge: %.4f\n" % \
                      (residue, residue.getCharge())
        header += "REMARK   5\n"
    header += "REMARK   6 Total charge on this protein: %.4f e\n" % charge
    header += "REMARK   6\n"

    return header

def runPDB2PQR(pdblist, ff, options):
    """
        Run the PDB2PQR Suite

        Parameters
            pdblist: The list of objects that was read from the PDB file
                     given as input (list)
            ff:      The name of the forcefield (string)
            options: A dictionary of PDB2PQR options, including:
                     verbose: When 1, script will print information to stdout
                              When 0, no detailed information will be printed (int)
                     debump:  When 1, debump heavy atoms (int)
                     opt:     When 1, run hydrogen optimization (int)
                     ph:      The desired ph of the system (float)
                     outname: The name of the desired output file
        Returns
            header:  The PQR file header (string)
            lines:   The PQR file atoms (list)
            missedligandresidues:  A list of ligand residue names whose charges could
                     not be assigned (ligand)
    """
    ph = None
    pkaname = ""
    outname = ""
    outroot = ""
    typemapname = ""
    neutraln = None
    neutralc = None
    lines = []
    Lig = None
    atomcount = 0   # Count the number of ATOM records in pdb

    # userff is CGI-based User Forcefield file object

    if "userff" in options: userff = options["userff"]
    else: userff = None

    if "usernames" in options: usernames = options["usernames"]
    else: usernames = None

    if "verbose" in options: verbose = 1
    else: verbose = 0

    if "opt" in options: optflag = 1
    else: optflag = 0

    if "typemap" in options: typemapflag = 1
    else: typemapflag = 0

    if "chain" in options: chainflag = 1
    else: chainflag = 0

    if "outname" not in options or options["outname"] == None:
        text = "Error: Output name not set!"
        raise ValueError, text
    else:
        outname = options["outname"]
        period = string.rfind(outname,".")
        if period > 0: outroot = outname[0:period]
        else: outroot = outname

    if "ph" in options:
        pka = 1
        ph = options["ph"]
        pkaname = outroot + ".propka"
        if os.path.isfile(pkaname): os.remove(pkaname)
    else: pka = 0

    typemapname = "%s-typemap.html" % outroot

    extmap = options["extensions"]
    
    start = time.time()

    if verbose:
        print "Beginning PDB2PQR...\n"

    myDefinition = Definition()
    if verbose:
        print "Parsed Amino Acid definition file."   

    # Check for the presence of a ligand!  This code is taken from pdb2pka/pka.py

    if "ligand" in options:
        from pdb2pka.ligandclean import ligff
        myProtein, myDefinition, Lig = ligff.initialize(myDefinition, options["ligand"], pdblist, verbose)        
        for atom in myProtein.getAtoms():
            if atom.type == "ATOM": 
                atomcount += 1
    else:
        myProtein = Protein(pdblist, myDefinition)

    if verbose:
        print "Created protein object -"
        print "\tNumber of residues in protein: %s" % myProtein.numResidues()
        print "\tNumber of atoms in protein   : %s" % myProtein.numAtoms()
        
    myRoutines = Routines(myProtein, verbose)

    for residue in myProtein.getResidues():
        multoccupancy = 0
        for atom in residue.getAtoms():
            if atom.altLoc != "":
                multoccupancy = 1
                txt = "Warning: multiple occupancies found: %s in %s\n" % (atom.name, residue)
                sys.stderr.write(txt)
        if multoccupancy == 1:
            myRoutines.warnings.append("WARNING: multiple occupancies found in %s,\n" % (residue))
            myRoutines.warnings.append("         at least one of the instances is being ignored.\n")

    if "neutraln" in options: neutraln = 1
    if "neutralc" in options: neutralc = 1

    myRoutines.setTermini(neutraln, neutralc)
    myRoutines.updateBonds()

    if "clean" in options:
        header = ""
        lines = myProtein.printAtoms(myProtein.getAtoms(), chainflag)
      
        # Process the extensions
        for ext in extmap:
            module = extmap[ext]
            call = "module.%s(myRoutines, outroot)" % ext
            eval(call)  
    
        if verbose:
            print "Total time taken: %.2f seconds\n" % (time.time() - start)
        return header, lines

    if not "assign-only" in options:
        # It is OK to process ligands with no ATOM records in the pdb
        if atomcount == 0 and Lig != None:
            pass
        else:
            myRoutines.findMissingHeavy()
        myRoutines.updateSSbridges()

        if "debump" in options:
            myRoutines.debumpProtein()  

        if pka:
            myRoutines.runPROPKA(ph, ff, pkaname)

        myRoutines.addHydrogens()

        myhydRoutines = hydrogenRoutines(myRoutines)

        if "debump" in options:
            myRoutines.debumpProtein()  

        if optflag:
            myhydRoutines.setOptimizeableHydrogens()
            myhydRoutines.initializeFullOptimization()
            myhydRoutines.optimizeHydrogens()
        else:
            myhydRoutines = hydrogenRoutines(myRoutines)
            myhydRoutines.initializeWaterOptimization()
            myhydRoutines.optimizeHydrogens()

        # Special for GLH/ASH, since both conformations were added
        myhydRoutines.cleanup()


    else:  # Special case for HIS if using assign-only
        for residue in myProtein.getResidues():
            if isinstance(residue, HIS):
                myRoutines.applyPatch("HIP", residue)

    myRoutines.setStates()

    myForcefield = Forcefield(ff, myDefinition, userff, usernames)
    hitlist, misslist = myRoutines.applyForcefield(myForcefield)
  
    ligsuccess = 0
    if "ligand" in options:

        # If this is independent, we can assign charges and radii here
 
        for residue in myProtein.getResidues():
            if isinstance(residue, LIG):
                templist = []
                Lig.make_up2date(residue)
                for atom in residue.getAtoms():
                    atom.ffcharge = Lig.ligand_props[atom.name]["charge"]
                    atom.radius = Lig.ligand_props[atom.name]["radius"]
                    if atom in misslist:
                        misslist.pop(misslist.index(atom))
                        templist.append(atom)

                charge = residue.getCharge()
                if abs(charge - int(charge)) > 0.001:
                    # Ligand parameterization failed
                    myRoutines.warnings.append("WARNING: PDB2PQR could not successfully parameterize\n")
                    myRoutines.warnings.append("         the desired ligand; it has been left out of\n")
                    myRoutines.warnings.append("         the PQR file.\n")
                    myRoutines.warnings.append("\n")
                    
                    # remove the ligand
                    myProtein.residues.remove(residue) 
                    for chain in myProtein.chains:
                        if residue in chain.residues: chain.residues.remove(residue)
                else:
                    ligsuccess = 1
                    # Mark these atoms as hits
                    hitlist = hitlist + templist
    
    # Temporary fix; if ligand was successful, pull all ligands from misslist
    if ligsuccess:
        templist = misslist[:]
        for atom in templist:
            if isinstance(atom.residue, Amino) or isinstance(atom.residue, Nucleic): continue
            misslist.remove(atom)

    # Creat the Typemap
    if typemapflag:
        myProtein.createHTMLTypeMap(myDefinition, typemapname)

    # Grab the protein charge

    reslist, charge = myProtein.getCharge()

    # If we want a different naming scheme, use that

    if "ffout" in options:
        scheme = options["ffout"]
        userff = None # Currently not supported
        if scheme != ff: myNameScheme = Forcefield(scheme, myDefinition, userff)
        else: myNameScheme = myForcefield
        myRoutines.applyNameScheme(myNameScheme)

    header = printPQRHeader(misslist, reslist, charge, ff, myRoutines.getWarnings(), options)
    lines = myProtein.printAtoms(hitlist, chainflag)

    # Determine if any of the atoms in misslist were ligands
    missedligandresidues = []
    for atom in misslist:
        if isinstance(atom.residue, Amino) or isinstance(atom.residue, Nucleic): continue
        if atom.resName not in missedligandresidues:
            missedligandresidues.append(atom.resName)

    # Process the extensions
 
    for ext in extmap:
        module = extmap[ext]
        call = "module.%s(myRoutines, outroot)" % ext
        eval(call)

    if verbose:
        print "Total time taken: %.2f seconds\n" % (time.time() - start)

    return header, lines, missedligandresidues

def getAvailableExtensions(displayflag=0):
    """
        Grab available extensions from the extensions directory

        Parameters
            displayflag: Display the error message if 1
        Returns
            extensions: A map containing the extensions name and
                        the module instance.
    """
    extensions = {}
    dir = "%s" % os.path.dirname(sys.argv[0])
    if dir == "": extdir = "extensions"
    else: extdir = "%s/extensions" % dir
    for filename in os.listdir(extdir):
        if filename.endswith(".py"):
            if filename == "__init__.py": continue
            
            # Test to see if we can find the function

            name = filename[:-3]
            try:
                e = __import__("extensions.%s" % name, globals(), locals(), name)
                if callable(eval("e.%s" % name)) and \
                   callable(eval("e.usage")):
                    extensions[name] = e
            except (AttributeError, ImportError):
                txt = "\nWarning: Missing either \"%s\" or \"usage\" functions in %s!" %\
                      (name, filename)
                txt += "\nThis extension will not be included.\n\n"
                if displayflag:
                    sys.stderr.write(txt)

    return extensions

def mainCommand(argv):
    """
        Main driver for running program from the command line.
    """

    # Append Numeric/Numpy path to sys.path if the user specified a non-standard location during configuration
    sys.argv=argv
    package_path = PACKAGE_PATH
    if package_path != "":
        sys.path.extend(package_path.split(":"))


    shortOptlist = "h,v"
    longOptlist = ["help","verbose","ff=","ffout=","nodebump","noopt","with-ph=","apbs-input","chain","clean","assign-only", "ligand=", "whitespace", "typemap", "neutraln", "neutralc"]

    extensions = getAvailableExtensions(1)
    longOptlist += extensions.keys()

    try: opts, args = getopt.getopt(sys.argv[1:], shortOptlist, longOptlist)
    except getopt.GetoptError, details:
        sys.stderr.write("GetoptError:  %s\n" % details)
        usage(2)

    if len(args) != 2:
        sys.stderr.write("Incorrect number (%d) of arguments!\n" % len(args))
        usage(2)

    options = {"debump":1,"opt":1,"extensions":{}}
 
    outpath = None
    ff = None
    for o,a in opts:
        undashed = o[2:]
        if o in ("-v","--verbose"):
            options["verbose"] = 1
        elif o in ("-h","--help"):
            usage(2)
            sys.exit()
        elif o == "--nodebump":  del options["debump"]
        elif o == "--noopt":    del options["opt"]
        elif o == "--apbs-input": options["input"] = 1
        elif o == "--whitespace": options["whitespace"]  = 1
        elif o == "--typemap": options["typemap"] = 1
        elif o == "--with-ph":
            try:
                ph = float(a)
                options["ph"] = ph
                if ph < 0.0 or ph > 14.0: raise ValueError
            except ValueError:
                text = "%s is not a valid pH!  " % a
                text += "Please choose a pH between 0.0 and 14.0."
                raise ValueError, text
        elif o == "--assign-only":
            del options["debump"]
            del options["opt"]
            options["assign-only"] = 1
        elif o == "--clean":
            del options["debump"]
            del options["opt"]
            options["clean"] = 1
        elif o == "--ff":      
            ff = a
            
            # Check to make sure forcefield file is available

            defpath = getFFfile(ff)
            if defpath == "":
                raise ValueError, "Unable to find parameter files for forcefield %s!" % ff
        elif o == "--userff" or o == "--usernames":
            pass
        elif o == "--neutraln": 
            if ff not in ["parse", "PARSE"]:
                raise ValueError, "neutraln option only works with PARSE forcefield!"
            options["neutraln"]  = 1

        elif o == "--neutralc": 
            if ff not in ["parse", "PARSE"]:
                raise ValueError, "neutralc option only works with PARSE forcefield!"
            options["neutralc"]  = 1

        elif o == "--chain": options["chain"] = 1
        elif o == "--ffout":
            if a.lower() in ["amber","charmm","parse","tyl06","peoepb","swanson"]:
                options["ffout"] = a
            else:
                raise ValueError, "Invalid forcefield naming scheme %s!" % a
        elif o == "--ligand":
            if os.path.isfile(a):
                options["ligand"] = open(a, 'rU')
            else:
                raise ValueError, "Unable to find ligand file %s!\n" % a
        elif undashed in extensions.keys():
            options["extensions"][undashed] = extensions[undashed]
            
    if ff == None and "clean" not in options:
        raise ValueError, "Forcefield not specified!"

    text =  "\n--------------------------\n"
    text += "PDB2PQR - a Python-based structural conversion utility\n"
    text += "--------------------------\n"
    text += "Please cite your use of PDB2PQR as:\n"
    text += "  Dolinsky TJ, Nielsen JE, McCammon JA, Baker NA.\n"
    text += "  PDB2PQR: an automated pipeline for the setup, execution,\n"
    text += "  and analysis of Poisson-Boltzmann electrostatics calculations.\n"
    text += "  Nucleic Acids Research 32 W665-W667 (2004).\n\n"
    sys.stdout.write(text)
            
    path = args[0]
    file = getPDBFile(path)
    pdblist, errlist = readPDB(file)
    
    if len(pdblist) == 0 and len(errlist) == 0:
        try: os.remove(path)
        except OSError: pass
        raise ValueError, "Unable to find file %s!\n" % path

    if len(errlist) != 0 and "verbose" in options:
        print "Warning: %s is a non-standard PDB file.\n" % path
        print errlist

    outpath = args[1]
    options["outname"] = outpath

    if "clean" not in options:
        header, lines, missedligands = runPDB2PQR(pdblist, ff, options)
    else:
        header, lines = runPDB2PQR(pdblist, ff, options)
        missedligands = None

    # Print the PQR file
    outfile = open(outpath,"w")
    outfile.write(header)
    # Adding whitespaces if --whitespace is in the options
    for line in lines:
        if "whitespace" in options: 
            if line[0:4] == 'ATOM':
                newline = line[0:16] + ' ' + line[16:38] + ' ' + line[38:46] + ' ' + line[46:]
                outfile.write(newline)
            elif line[0:6] == 'HETATM':
                newline = line[0:16] + ' ' + line[16:38] + ' ' + line[38:46] + ' ' + line[46:]
                outfile.write(newline)
        else: 
            outfile.write(line)
    outfile.close()

    if "input" in options:
        from src import inputgen
        from src import psize
        method = "mg-auto"
        size = psize.Psize()
        size.parseInput(outpath)
        size.runPsize(outpath)
        async = 0 # No async files here!
        input = inputgen.Input(outpath, size, method, async)
        input.printInputFiles()
        input.dumpPickle()


if __name__ == "__main__":
    mainCommand(sys.argv)
