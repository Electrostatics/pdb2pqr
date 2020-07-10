"""APBS interface for PDB2PQR

Authors:  Todd Dolinsky and Jens Erik Nielsen
"""
import logging
from sys import version_info
import os
import time
import copy


_LOGGER = logging.getLogger(__name__)


try:
    from inspect import currentframe, getframeinfo
    # TODO - emliminate import *
    from .apbslib import *
except:
    ERRSTR = "Missing Python libraries for APBS interface. "
    ERRSTR += "This often means you need to build APBS with Python support "
    ERRSTR += "using the Cmake argument -DENABLE_PYTHON=ON"
    _LOGGER.error(ERRSTR)
    raise RuntimeError(ERRSTR)


PYTHON_KB = 1.3806581e-23
PYTHON_NA = 6.0221367e+23
NOSH_MAXMOL = 20
NOSH_MAXCALC = 20
ACD_ERROR = 2    # < Error setting up calculation>


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
        #return `self.value`
        return repr(self.value)

class runAPBS:

    def __init__(self):
        return


    def getUnitConversion(self):
        """
            Get the unit conversion from kT to kJ/mol

            Returns
                factor: The conversion factor (float)
        """
        temp = 298.15
        factor = PYTHON_KB / 1000.0 * temp * PYTHON_NA
        return factor

    def runAPBS(self, protein, inputpath, routines, cm_val=None):
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

        old_cwd = os.getcwd()
        new_cwd = os.path.dirname(os.path.abspath(inputpath))

        os.chdir(new_cwd)

        try:
            #
            # Initialize the MALOC library
            startVio()

            # Initialize variables, arrays
            self.com = Vcom_ctor(1)
            self.rank = Vcom_rank(self.com)
            self.size = Vcom_size(self.com)
            self.mgparm = MGparm()
            self.pbeparm = PBEparm()
            self.mem = Vmem_ctor("Main")
            self.pbe = new_pbelist(NOSH_MAXMOL)
            self.pmg = new_pmglist(NOSH_MAXMOL)
            self.pmgp = new_pmgplist(NOSH_MAXMOL)
            self.real_center = double_array(3)
            self.tot_energy = []
            self.x = []
            self.y = []
            self.z = []
            self.chg = []
            self.rad = []
            #nforce = int_array(NOSH_MAXCALC)
            #atomforce = new_atomforcelist(NOSH_MAXCALC)
            #nfor = ptrcreate("int",0)

            # Start the main timer
            self.main_timer_start = time.clock()

            # Check invocation
            #stdout.write(getHeader())

            # Parse the input file
            self.nosh = NOsh_ctor(self.rank, self.size)
            #nosh = NOsh()
            #NOsh_ctor2(nosh, rank, size)
            _LOGGER.info("Parsing input file %s...\n", inputpath)
            if NOsh_parseInputFile(self.nosh, inputpath) != 1:
                _LOGGER.error("main:  Error while parsing input file.\n")
                raise APBSError("Error while parsing input file!")

            # Load the molecules using Valist_load routine

            self.alist = new_valist(NOSH_MAXMOL)
            self.atoms = protein.atoms
            self.protsize = len(self.atoms)

            # SETUP CALCULATIONS

            if NOsh_setupElecCalc(self.nosh, self.alist) != 1:
                _LOGGER.error("Error setting up ELEC calculations\n")
                raise APBSError("Error while setting up calculations!")

            if NOsh_setupApolCalc(self.nosh, self.alist) == ACD_ERROR:
                _LOGGER.error("Error setting up APOL calculations\n")
                raise APBSError("Error while setting up calculations!")

            #
            # DEBUGGING
            #
            self.cm_list = []

            #print 'These are the charges in the PQR file'
            #print 'atom name\tresnum\tresname\tcharge\tradius'
            #res_charge = 0.0
            #old_res = -1
            #old_res_name = ''
            for i in range(len(self.atoms)):
                atom = self.atoms[i]
                self.x.append(atom.get("x"))
                self.y.append(atom.get("y"))
                self.z.append(atom.get("z"))
                self.chg.append(atom.get("ffcharge"))
                self.rad.append(atom.get("radius"))
                self.cm_list.append([atom.res_seq, atom.name, atom.get("ffcharge")])
            #
            # DEBUG
            #
            if cm_val:
                cm_val.display_charges(self.cm_list)

            self.my_alist = make_Valist(self.alist, 0)
            Valist_load(self.my_alist, self.protsize, self.x, self.y, self.z, self.chg, self.rad)

            # Initialize the energy holders
            self.tot_energy = [0.0] * int(self.nosh.ncalc)

            pot_list = []

            # Initialize the force holders
            force_list = []

            # Load the dieletric maps

            self.diel_xmap = new_gridlist(NOSH_MAXMOL)
            self.diel_ymap = new_gridlist(NOSH_MAXMOL)
            self.diel_zmap = new_gridlist(NOSH_MAXMOL)

            if loadDielMaps(self.nosh, self.diel_xmap, self.diel_ymap, self.diel_zmap) != 1:
                _LOGGER.error("Error reading dielectric maps!\n")
                raise APBSError("Error reading dielectric maps!")

            # Load the kappa maps
            self.kappa_map = new_gridlist(NOSH_MAXMOL)
            if loadKappaMaps(self.nosh, self.kappa_map) != 1:
                _LOGGER.error("Error reading kappa maps!\n")
                raise APBSError("Error reading kappa maps!")

            # Load the potential maps
            self.pot_map = new_gridlist(NOSH_MAXMOL)
            if loadPotMaps(self.nosh, self.pot_map) != 1:
                _LOGGER.error("Error reading potential maps!\n")
                raise APBSError("Error reading potential maps!")

            # Load the charge maps
            self.charge_map = new_gridlist(NOSH_MAXMOL)
            if loadChargeMaps(self.nosh, self.charge_map) != 1:
                _LOGGER.error("Error reading charge maps!\n")
                raise APBSError("Error reading charge maps!")

            # Do the calculations

            _LOGGER.info("Preparing to run %d PBE calculations. \n" % self.nosh.ncalc)

            for icalc in range(self.nosh.ncalc):
                _LOGGER.infowrite("---------------------------------------------\n")
                self.calc = NOsh_getCalc(self.nosh, icalc)
                self.mgparm = self.calc.mgparm
                self.pbeparm = self.calc.pbeparm
                if self.calc.calctype != 0:
                    _LOGGER.error("main:  Only multigrid calculations supported!\n")
                    raise APBSError("Only multigrid calculations supported!")

                for k in range(0, self.nosh.nelec):
                    if NOsh_elec2calc(self.nosh, k) >= icalc:
                        break

                name = NOsh_elecname(self.nosh, k + 1)

                # Routine initMG

                if initMG(icalc, self.nosh, self.mgparm, self.pbeparm,
                          self.real_center, self.pbe, self.alist, self.diel_xmap,
                          self.diel_ymap, self.diel_zmap, self.kappa_map,
                          self.charge_map, self.pmgp, self.pmg, self.pot_map) != 1:
                    _LOGGER.error("Error setting up MG calculation!\n")
                    raise APBSError("Error setting up MG calculation!")

                # Print problem parameters

                #printMGPARM(self.mgparm, self.real_center)
                #printPBEPARM(self.pbeparm)

                # Solve the problem : Routine solveMG

                self.thispmg = get_Vpmg(self.pmg, icalc)

                if solveMG(self.nosh, self.thispmg, self.mgparm.type) != 1:
                    _LOGGER.error("Error solving PDE! \n")
                    raise APBSError("Error Solving PDE!")

                # Set partition information : Routine setPartMG

                if setPartMG(self.nosh, self.mgparm, self.thispmg) != 1:
                    _LOGGER.error("Error setting partition info!\n")
                    raise APBSError("Error setting partition info!")

                ret, self.tot_energy[icalc] = energyMG(self.nosh, icalc, self.thispmg, 0,
                                                       self.tot_energy[icalc], 0.0, 0.0, 0.0)

                # Set partition information

                #aforce = get_AtomForce(atomforce, icalc)
                #forceMG(mem, nosh, pbeparm, mgparm, thispmg, nfor, aforce, alist)
                #ptrset(nforce,ptrvalue(nfor), icalc)

                # Write out data from MG calculations : Routine writedataMG
                writedataMG(self.rank, self.nosh, self.pbeparm, self.thispmg)

                # Write out matrix from MG calculations
                writematMG(self.rank, self.nosh, self.pbeparm, self.thispmg)

                # GET THE POTENTIALS
                potentials = getPotentials(self.nosh, self.pbeparm, self.thispmg, self.my_alist)
                pot_list.append(potentials)
        finally:
            os.chdir(old_cwd)

        #
        # Cleanup
        #
        return pot_list

    #
    # ------
    #

    def get_potentials(self, protein):

        #import copy

        #sys.setrecursionlimit(10000)

        delete_valist(self.alist)
        self.alist = new_valist(NOSH_MAXMOL)
        self.atoms = protein.atoms
        self.protsize = len(self.atoms)
        proteincopy = copy.copy(protein)

        for i in range(len(self.atoms)):
            atom = self.atoms[i]
            self.x.append(atom.get("x"))
            self.y.append(atom.get("y"))
            self.z.append(atom.get("z"))
            self.chg.append(atom.get("ffcharge"))
            self.rad.append(atom.get("radius"))
            #self.cm_list.append([atom.res_seq,atom.name,atom.get("ffcharge")])
        #
        # DEBUG
        #
        #if cm_val:
        #    cm_val.display_charges(self.cm_list)
        #
        self.my_alist = make_Valist(self.alist, 0)

        xlist = self.x[-1 * (self.protsize):]
        ylist = self.y[-1 * (self.protsize):]
        zlist = self.z[-1 * (self.protsize):]
        chglist = self.chg[-1 * (self.protsize):]
        radlist = self.rad[-1 * (self.protsize):]

        try:
            Valist_load(self.my_alist, self.protsize, xlist, ylist, zlist, chglist, radlist)
        except:
            frameinfo = getframeinfo(currentframe())
            _LOGGER.warning("%s[%d]: Valist_load Warning." % (frameinfo.filename, frameinfo.lineno))
        potentials = getPotentials(self.nosh, self.pbeparm, self.thispmg, self.my_alist)

        protein = copy.copy(proteincopy)

        # Free up the memory allocated for self.my_alist
        remove_Valist(self.my_alist)

        return potentials


    #
    # ------
    #

    def cleanup(self):

        # Handle print statements

        # Clean up APBS structures

        #killForce(mem, nosh, nforce, atomforce)
        killEnergy()
        killMG(self.nosh, self.pbe, self.pmgp, self.pmg)
        killChargeMaps(self.nosh, self.charge_map)
        killKappaMaps(self.nosh, self.kappa_map)
        killPotMaps(self.nosh, self.pot_map)
        killDielMaps(self.nosh, self.diel_xmap, self.diel_ymap, self.diel_zmap)

        if self.my_alist.number == 0:
            self.my_alist = make_Valist(self.alist, 0)
            Valist_load(self.my_alist, self.protsize, self.x, self.y, self.z, self.chg, self.rad)
        killMolecules(self.nosh, self.alist)


        del self.nosh

        # Clean up Python structures

        #ptrfree(nfor)
        delete_double_array(self.real_center)
        #delete_int_array(nforce)
        #delete_atomforcelist(atomforce)
        delete_valist(self.alist)
        delete_gridlist(self.diel_xmap)
        delete_gridlist(self.diel_ymap)
        delete_gridlist(self.diel_zmap)
        delete_gridlist(self.kappa_map)
        delete_gridlist(self.pot_map)
        delete_gridlist(self.charge_map)
        delete_pmglist(self.pmg)
        delete_pmgplist(self.pmgp)
        delete_pbelist(self.pbe)

        # Clean up MALOC structures
        del self.com
        del self.mem

        return
