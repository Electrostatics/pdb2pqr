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

class conf_avg:

	def __init__(self,options):
		
		potentials=[]
		
		listOfFiles=os.listdir(options.directoryPath)

#options.pdbfile_name

		for currentPDB in listOfFiles:
			print ".........................................................................................WORKING ON " + currentPDB
			pdbfile = getPDBFile(currentPDB)
			pdblist, errlist = readPDB(pdbfile)

			pots=self.get_potentials(currentPDB,pdblist)
			potentials.append(pots[0])
		#
		# Average potentials
		#
		avg_pots=self.average_potentials(potentials)
		return

	def get_potentials(self,currentPDB,pdblist):
		"""Get the potentials by first running pdb2pqr and then apbs"""
		myProtein,apbs_inputfile=self.run_pdb2pqr(currentPDB,pdblist)
		potentials=self.run_apbs(myProtein,apbs_inputfile)
		return potentials
	
	def run_pdb2pqr(self,currentPDB,pdblist):
		"""Run pdb2pqr, prepare input for apbs"""
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
		# Add and optimze hydrogens:
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
		igen = IP.inputGen(currentPDB)
		igen.maps=None
		igen.set_type('intene')
		igen.pdie=8.0
		igen.sdie=80.0
		all_center,extent=igen.getCenter()
		igen.setfineCenter(all_center)
		print 'Center: %5.1fA %5.1fA %5.1fA' %(all_center[0],all_center[1],all_center[2])
		print 'Extent: %5.1fA %5.1fA %5.1fA'  %(extent[0],extent[1],extent[2])

		apbs_inputfile=igen.printInput()

		return myProtein, apbs_inputfile

	def run_apbs(self,myProtein,apbs_inputfile):
		"""runs apbs"""
		from pdb2pka.apbs import *
		APBS=runAPBS()
		potentials = APBS.runAPBS(myProtein, apbs_inputfile)
		APBS.cleanup()

		return potentials

	def average_potentials(self,potentials):
		"""This function averages many potential maps"""
		avg_pots=[]
		for i in range(0,len(potentials[0])):
			currSum=0
			for j in range(0,len(potentials)):
				currSum+=potentials[j][i]
			currAvg=currSum/len(potentials)
			avg_pots.append(currAvg)

		print avg_pots
		return avg_pots


if __name__=='__main__':
	from optparse import OptionParser
	parser = OptionParser(usage='%prog [options]',version='%prog 1.0')
#	parser.add_option('-p','--pdb',dest='pdbfile_name',action='store',type='string',nargs=2,default='2lzt.pka.pdb',
#                  help='The PDB file. Default: %default')
	parser.add_option('-d','--dir',dest='directoryPath',action='store',type='string',default='',
				  help='Direcotry of the PDB files. Default: %default')

# I have to fix the path thing

	(options, args) = parser.parse_args()

	verbose=True
	ff='parse'

	I=conf_avg(options)

