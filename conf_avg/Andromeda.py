#!/usr/bin/env python
print 'Hello'
import sys
sys.path.append('/home/mag_local/development/pdb2pqr')
import string
import math
import string
import getopt
import time
from src import pdb
from src import utilities
from src import structures
from src import routines
from src import protein
from src import server
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
myDefinition = Definition()
#
pdbfile_name='2lzt.pka.pdb'
pdbfile = getPDBFile(pdbfile_name)
pdblist, errlist = readPDB(pdbfile)
               

# Here, self.protein is the PDB2PQR instance. you will get this with something like:
myProtein = Protein(pdblist, myDefinition)

# then you do:
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

#myRoutines.randomizeWaters()
myProtein.reSerialize()

# to clean up and add hydrogens.

# optimze hydrogens:
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

# To print the charges do:
if False:
    for chain in myProtein.getChains():
        print 'Chain',chain
        for residue in chain.get("residues"):
            print residue
            for atom in residue.get("atoms"):
                atomname = atom.get("name")
                charge, radius = myforcefield.getParams1(residue, atomname)
                print atomname,charge,radius

#import src.psize
#size=src.psize.Psize()

method=""
async=0
split=0
import pdb2pka.inputgen_pKa as IP
igen = IP.inputGen(pdbfile_name)
igen.maps=0
igen.set_type('desolv')
igen.pdie=8

apbs_inputfile=igen.printInput()
import pdb2pka.apbs as RAPBS
APBS=RAPBS.runAPBS()
potentials = APBS.runAPBS(myProtein, apbs_inputfile)
APBS.cleanup()

print potentials
