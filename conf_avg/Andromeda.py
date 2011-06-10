#!/usr/bin/env python
#
# This is a module in development for making conformationally averaged PBE maps
#
debug=False
import sys, os

print __file__
import os
try:
    file_name=__file__
    if file_name[:2]=='./':
        scriptpath=os.getcwd()
    else:
        scriptpath=os.path.join(os.getcwd(),os.path.split(file_name)[0])
        if scriptpath[-1] == "/":
            scriptpath=scriptpath[:-1]
except:
    scriptpath=os.path.split(sys.argv[0])[0]
    if scriptpath=='.':
        scriptpath=os.getcwd()
#
# Add to import path
#
pdb2pqr_path=os.path.split(scriptpath)[0]
sys.path.append(pdb2pqr_path)

import string
import math
import string
import getopt
import time
from src.pdb import *
from src.utilities import *
from src.structures import *
from src.definitions import *
from src.forcefield import *  
from src.routines import *
from src.protein import *
from src.server import *
from StringIO import *
from src.hydrogens import *

#
# Read the pdb file
#
verbose=True
ff='parse'

#
# Read the PDB file
#
pdbfile_name='2lzt.pka.pdb'
pdbfile = getPDBFile(pdbfile_name)
pdblist, errlist = readPDB(pdbfile)
               
#
# Instantiate pdb2pqr
#
myDefinition = Definition()
myProtein = Protein(pdblist, myDefinition)

#
# Setup everything
#
myRoutines = Routines(myProtein, verbose)
myRoutines.updateResidueTypes()
myRoutines.updateSSbridges()
myRoutines.updateBonds()
myRoutines.setTermini()
myRoutines.updateInternalBonds()

myforcefield=Forcefield(ff, myDefinition, None)
myRoutines.applyNameScheme(myforcefield)

myRoutines.findMissingHeavy()
myRoutines.addHydrogens()
myRoutines.debumpProtein()
myProtein.reSerialize()

#
# Addn and optimze hydrogens:
# 
from src.hydrogens import hydrogenRoutines
myRoutines.updateInternalBonds()
myRoutines.calculateDihedralAngles()
myhydRoutines = hydrogenRoutines(myRoutines)
#
# Now optimize hydrogens
#
myhydRoutines.setOptimizeableHydrogens()
myhydRoutines.initializeFullOptimization()
myhydRoutines.optimizeHydrogens()
myhydRoutines.cleanup()
myRoutines.setStates()

print "Created protein object (after processing myRoutines) -"
print "\tNumber of residues in protein: %s" % myProtein.numResidues()
print "\tNumber of atoms in protein   : %s" % myProtein.numAtoms()

#
# Assign charges
#
for chain in myProtein.getChains():
    for residue in chain.get("residues"):
        for atom in residue.get("atoms"):
            atomname = atom.get("name")
            charge, radius = myforcefield.getParams1(residue, atomname)
            atom.set("radius", radius)
            atom.set("ffcharge", charge)

#import src.psize
#size=src.psize.Psize()

method=""
async=0
split=0
import pdb2pka.inputgen_pKa as IP
igen = IP.inputGen(pdbfile_name)
igen.maps=None
igen.set_type('intene')
igen.pdie=8.0
igen.sdie=80.0
all_center,extent=igen.getCenter()
igen.setfineCenter(all_center)
print 'Center: %5.1fA %5.1fA %5.1fA' %(all_center[0],all_center[1],all_center[2])
print 'Extent: %5.1fA %5.1fA %5.1fA'  %(extent[0],extent[1],extent[2])


apbs_inputfile=igen.printInput()
from pdb2pka.apbs import *
APBS=runAPBS()
potentials = APBS.runAPBS(myProtein, apbs_inputfile)
APBS.cleanup()

print potentials
