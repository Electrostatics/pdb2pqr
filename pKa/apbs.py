"""
    APBS interface for PDB2PQR

    Todd Dolinsky (todd@ccb.wustl.edu)
    Washington University in St. Louis

    Jens Erik Nielsen

"""

__date__  = "23 September 2004"
__author__ = "Todd Dolinsky, Jens Erik Nielsen"

import sys
import time
from apbslib import *

Python_kb = 1.3806581e-23
Python_Na = 6.0221367e+23
NOSH_MAXMOL = 20
NOSH_MAXCALC = 20

class APBSError(Exception):
    """ APBSError class

        The APBSError class inherits off the Exception module and returns
        a string defining the nature of the error. 
    """
    
    def __init__(self, value):
        """
            Initialize with error message

            Parameters
                value:  Error Message (string)
        """
        self.value = value
        
    def __str__(self):
        """
            Return the error message
        """
        return `self.value`

def getUnitConversion():
    """
        Get the unit conversion from kT to kJ/mol

        Returns
            factor: The conversion factor (float)
    """
    temp = 298.15
    factor = Python_kb/1000.0 * temp * Python_Na
    return factor

def runAPBS(protein, inputpath):
    """
        Run APBS, using the protein instead of a pqr file

        Parameters
            protein:    The protein object (protein)
            inputpath:  The path to the APBS input file (string)
        Returns
            potentials: A list of lists of potentials at atom
                        locations - one list for each APBS
                        calculation
    """

    # Initialize the MALOC library
    startVio()

    # Initialize variables, arrays
    com = Vcom_ctor(1)
    rank = Vcom_rank(com)
    size = Vcom_size(com)
    mgparm = MGparm()
    pbeparm = PBEparm()
    mem = Vmem_ctor("Main")
    pbe = new_pbelist(NOSH_MAXMOL)
    pmg = new_pmglist(NOSH_MAXMOL)
    pmgp = new_pmgplist(NOSH_MAXMOL)
    realCenter = double_array(3)
    totEnergy = double_array(NOSH_MAXCALC)
    qfEnergy = double_array(NOSH_MAXCALC)
    qmEnergy = double_array(NOSH_MAXCALC)
    dielEnergy = double_array(NOSH_MAXCALC)
    npEnergy = double_array(NOSH_MAXCALC)
    nenergy = int_array(NOSH_MAXCALC)
    nforce = int_array(NOSH_MAXCALC)
    atomforce = new_atomforcelist(NOSH_MAXCALC)
    
    # Start the main timer
    main_timer_start = time.clock()

    # Check invocation
    #stdout.write(getHeader())

    # Parse the input file
    nosh = NOsh()
    NOsh_ctor2(nosh, rank, size)
    sys.stdout.write("Parsing input file %s...\n" % inputpath)
    if NOsh_parseFile(nosh, inputpath) != 1:
        sys.stderr.write("main:  Error while parsing input file.\n")
        raise APBSError, "Error while parsing input file!"

    # Load the molecules using Valist_load routine

    alist = new_valist(NOSH_MAXMOL)
    atoms = protein.getAtoms()
    protsize = len(atoms)
    x = double_array(protsize)
    y = double_array(protsize)
    z = double_array(protsize)
    chg = double_array(protsize)
    rad = double_array(protsize)

    for i in range(protsize):
        atom = atoms[i]
        set_entry(x,i,atom.get("x"))
        set_entry(y,i,atom.get("y"))
        set_entry(z,i,atom.get("z"))
        set_entry(chg,i,atom.get("ffcharge"))
        set_entry(rad,i,atom.get("radius"))
  
    myAlist = make_Valist(alist,0)
    Valist_load(myAlist, protsize, x,y,z,chg,rad) 

    potList = []
    potentials = double_array(protsize)
  
    # Load the dieletric maps

    dielXMap = new_gridlist(NOSH_MAXMOL)
    dielYMap = new_gridlist(NOSH_MAXMOL)
    dielZMap = new_gridlist(NOSH_MAXMOL)
  
    if loadDielMaps(nosh, dielXMap, dielYMap, dielZMap) != 1:
        sys.stderr.write("Error reading dielectric maps!\n")
        raise APBSError, "Error reading dielectric maps!"
 
    # Load the kappa maps
    kappaMap = new_gridlist(NOSH_MAXMOL)
    if loadKappaMaps(nosh, kappaMap) != 1:
        sys.stderr.write("Error reading kappa maps!\n")
        raise APBSError, "Error reading kappa maps!"

    # Load the charge maps
    chargeMap = new_gridlist(NOSH_MAXMOL)
    if loadChargeMaps(nosh, chargeMap) != 1:
        sys.stderr.write("Error reading charge maps!\n")
        raise APBSError, "Error reading charge maps!"
    
    # Do the calculations
 
    sys.stdout.write("Preparing to run %d PBE calculations. \n" % nosh.ncalc)
   
    for icalc in xrange(nosh.ncalc):
        sys.stdout.write("---------------------------------------------\n")
        calc = NOsh_getCalc(nosh, icalc)
        mgparm = calc.mgparm
        pbeparm = calc.pbeparm
        if calc.calctype != 0:
            sys.stderr.write("main:  Only multigrid calculations supported!\n")
            raise APBSError, "Only multigrid calculations supported!"

        for k in range(0, nosh.nelec):
            if NOsh_elec2calc(nosh,k) >= icalc:
                break

        name = NOsh_elecname(nosh, k+1)
        if name == "":
            sys.stdout.write("CALCULATION #%d:  MULTIGRID\n" % (icalc+1))
        else:
            sys.stdout.write("CALCULATION #%d (%s): MULTIGRID\n" % ((icalc+1),name))
        sys.stdout.write("Setting up problem...\n")
	
        # Routine initMG
	
        if initMG(icalc, nosh, mgparm, pbeparm, realCenter, pbe, 
              alist, dielXMap, dielYMap, dielZMap, kappaMap, chargeMap, 
              pmgp, pmg) != 1:
            sys.stderr.write("Error setting up MG calculation!\n")
            raise APBSError, "Error setting up MG calculation!"
	
        # Print problem parameters 
	
        printMGPARM(mgparm, realCenter)
        printPBEPARM(pbeparm)
        Python_kbT = Python_kb*Python_Na*(pbeparm.temp)/1000.0
	
        # Solve the problem : Routine solveMG
	
        thispmg = get_Vpmg(pmg,icalc)

        if solveMG(nosh, thispmg, mgparm.type) != 1:
            stderr.write("Error solving PDE! \n")
            raise APBSError, "Error Solving PDE!"

        # Set partition information : Routine setPartMG

        if setPartMG(nosh, mgparm, thispmg) != 1:
            sys.stderr.write("Error setting partition info!\n")
            raise APBSError, "Error setting partition info!"
	
        # Write out energies : Routine energyMG, npenergyMG
        # Create pointers to variables that store energies,
        # place in appropriate array, and delete pointers
        # See SWIG documentation for pointer function info
	
        totEng = ptrcreate("double",0.0)
        qfEng = ptrcreate("double",0.0)
        qmEng = ptrcreate("double",0.0)
        dielEng = ptrcreate("double",0.0)
        npEng = ptrcreate("double",0.0)
        neng = ptrcreate("int",0.0)

        energyMG(nosh, icalc, thispmg, nenergy,  totEng, qfEng, qmEng, dielEng)
        #npenergyMG(nosh, icalc, thispmg, nenergy, npEng)
	
        ptrset(totEnergy,ptrvalue(totEng),icalc)
        ptrset(qfEnergy,ptrvalue(qfEng),icalc)
        ptrset(qmEnergy,ptrvalue(qmEng),icalc)
        ptrset(dielEnergy,ptrvalue(dielEng),icalc)
        ptrset(npEnergy,ptrvalue(npEng),icalc)
        ptrset(nenergy,ptrvalue(neng),icalc)
        ptrfree(totEng)
        ptrfree(qfEng)
        ptrfree(qmEng)
        ptrfree(dielEng)
        ptrfree(npEng)
        ptrfree(neng)
	
        # Set partition information

        nfor = ptrcreate("int",0.0)
        forceMG(mem, nosh, pbeparm, mgparm, thispmg, nfor, atomforce, alist)
        ptrset(nforce,ptrvalue(nfor), icalc)
        ptrfree(nfor)
	
        # Write out data from MG calculations : Routine writedataMG	
        writedataMG(rank, nosh, pbeparm, thispmg)
	
        # Write out matrix from MG calculations	
        writematMG(rank, nosh, pbeparm, thispmg)

        # GET THE POTENTIALS
              
        potentials = getPotentials(nosh, pbeparm, thispmg, myAlist)
        potList.append(potentials)
        
    # Handle print statements

    if nosh.nprint > 0:
        sys.stdout.write("---------------------------------------------\n")
        sys.stdout.write("PRINT STATEMENTS\n")
    for iprint in xrange(nosh.nprint):
        if NOsh_printWhat(nosh, iprint) == NPT_ENERGY:
            printEnergy(com, nosh, totEnergy, iprint)
        elif NOsh_printWhat(nosh, iprint) == NPT_FORCE:
            printForce(com, nosh, nforce, atomforce, iprint)
        else:
            sys.stdout.write("Undefined PRINT keyword!\n")
            break

    # Put the potentials into Python readable arrays

    mypotlist = []
    for i in range(len(potList)):
        mypotlist.append([])
        for j in range(protsize):
            mypotlist[i].append(get_entry(potList[i], j))
	
    sys.stdout.write("----------------------------------------\n")
    sys.stdout.write("CLEANING UP AND SHUTTING DOWN...\n")

    # Clean up APBS structures
    killForce(mem, nosh, nforce, atomforce)
    killEnergy()
    killMG(nosh, pbe, pmgp, pmg)
    killChargeMaps(nosh, chargeMap)
    killKappaMaps(nosh, kappaMap)
    killDielMaps(nosh, dielXMap, dielYMap, dielZMap)
    killMolecules(nosh, alist)
    del nosh
    
    # Clean up MALOC structures
    del com
    del mem
    sys.stdout.write("\n")
    sys.stdout.write("Thanks for using APBS!\n\n")

    # Stop the main timer
    main_timer_stop = time.clock()
    sys.stdout.write("Total execution time:  %1.6e sec\n" % (main_timer_stop - main_timer_start))

    #Return potentials
    return mypotlist
