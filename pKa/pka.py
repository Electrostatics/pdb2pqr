#!/bin/env python
#
# pKa calculations with APBS
#
# Jens Erik Nielsen & Todd Dolinsky
#
# University College Dublin & Washington University St. Louis
#
__date__="16 August 2005"
__author__="Jens Erik Nielsen, Todd Dolinsky"

import getopt
import sys, os
#
# find the path to the script, and to pdb2pqr
#
scriptpath=os.path.split(sys.argv[0])[0]
pdb2pqr_path=os.path.split(scriptpath)[0]
sys.path.append(pdb2pqr_path)

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

#from utilities import *
#from pdb import *
#from forcefield import *
#from protein import *
#from definitions import *
#from hydrogens import *
from apbs import *

TITRATIONFILE = os.path.join(scriptpath,"TITRATION.DAT")

class pKaRoutines:
    """
        Class for running all pKa related functions
    """
    def __init__(self, protein, routines, forcefield,apbs_inputfile):
        """
            Initialize the class using needed objects

            Parameters
                protein:    The PDB2PQR protein object
                routines:   The PDB2PQR routines object
                forcefield: The PDB2PQR forcefield object
        """
        self.protein = protein
        self.routines = routines
        self.forcefield = forcefield
        self.apbs_inputfile=apbs_inputfile
        self.pKagroups = self.readTitrationDefinition() 
        self.pKas = []

        myHydrogenRoutines = hydrogenRoutines(routines)
        myHydrogenRoutines.readHydrogenDefinition()
        self.hydrogenRoutines = myHydrogenRoutines
        #
        # Not sure this is the best place for the interaction energies...
        #
        self.matrix={}

        return

    #
    # -----------------------------------------
    #
        
    def runpKa(self):
        """
            Main driver for running pKa calculations
        """
        self.pKas = self.findTitratableGroups()

        self.calculateIntrinsicpKa()

        """ Calculate Pairwise Interactions """
        #self.calculatePairwiseInteractions()

        """ Calculate Full pKa Value """
        self.calculatepKaValue()

    #
    # -----------------------------------
    #

    def calculatePairwiseInteractions(self):
        #
        # Calculate the pairwise interaction energies
        #
        for pKa in self.pKas:
            residue = pKa.residue
            pKaGroup = pKa.pKaGroup
            ambiguity = pKa.amb
            atomnames = self.getAtomsForPotential(pKa)
            #
            # TODD: What does this one do?
            #
            fixedstates = self.hydrogenRoutines.getstates(ambiguity)
            #
            # Loop over each titration
            #
            if not self.matrix.has_key(pKa):
                self.matrix[pKa]={}
            #
            for titration in pKaGroup.pKaTitrations:
                if not self.matrix[pKa].has_key(titration):
                    self.matrix[pKa][titration]={}
                #
                # Get all states
                #
                possiblestates = titration.allstates
                for state in possiblestates:
                    self.hydrogenRoutines.switchstate(fixedstates, ambiguity, state) 
                    self.zeroAllRadiiCharges()
                    self.setAllRadii()
                    self.setCharges(residue, atomnames)
                    #
                    self.matrix[pKa][titration][state]=self.get_interaction_energies()
        return

    def get_interaction_energies(self):
        #
        # Run APBS and get the interaction with all other states
        #
        potentials=self.getAPBSPotentials()
        #
        # construct this side
        #
        energies={}
        #
        # Loop over all groups
        #
        for pKa in self.pKas:
            residue = pKa.residue
            pKaGroup = pKa.pKaGroup
            ambiguity = pKa.amb
            #
            # TODD: What does this one do?
            #
            fixedstates = self.hydrogenRoutines.getstates(ambiguity)
            #
            # Loop over each titration
            #
            if not energies.has_key(pKa):
                energies[pKa]={}
            #
            for titration in pKaGroup.pKaTitrations:
                if not energies[pKa].has_key(titration):
                    energies[pKa][titration]={}
                #
                # Get all states
                #
                possiblestates = titration.allstates
                for state in possiblestates:
                    #
                    # Switch to the particular state
                    #
                    atomnames = self.getAtomsForPotential(pKa)
                    self.hydrogenRoutines.switchstate(fixedstates, ambiguity, state) 
                    self.zeroAllRadiiCharges()
                    self.setAllRadii()
                    self.setCharges(residue, atomnames)
                    #
                    # Get atoms for potential
                    #
                    #print 'Atoms for measuring potential',atomnames
                    atomlist=[]
                    for atomname in atomnames:
                        atomlist.append(residue.getAtom(atomname))
                    energy=0.0
                    count=0
                    for atom in protein.getAtoms():
                        if atom in atomlist:
                            energy=energy+(potentials[3][count] - potentials[1][count])*atom.get("ffcharge")
                            #print 'Getting potential',residue.get('name'),atom.name,atom.get('ffcharge')
                        count=count+1
                    energies[pKa][titration][state]=energy
                
        return energies

    #
    # ----------------------------------
    #

    def calculatepKaValue(self):
        #
        #  Calculate the pKa Value
        #
        # We use a c++ class for the MC steps since it's a lot faster..
        #
        import pMC
        #
        # Matrix needs to be linear for transport to c++
        # TODD: Do you know how to pass a full dictionary or list of lists to c++?
        #
        linear=[]
        intpkas=[]
        acidbase=[]
        print
        print
        print 'Residue\tIntrinsic pKa'
        for pKa in self.pKas:
            pKaGroup = pKa.pKaGroup
            for titration in pKaGroup.pKaTitrations:
               
                # Acid/Base
                if pKaGroup.type=='acid':
                    acidbase.append(-1)
                else:
                    acidbase.append(1)
                # Intrinsic pKa value
                intpkas.append(titration.intrinsic_pKa)
                print pKa.residue.get('name'),titration.intrinsic_pKa
                #intpkas.append(4.0)
                possiblestates = titration.allstates
                for state in [possiblestates[1]]: #TODO: Modify so we use all states this means
                    # modifying the c++ module too....Correct for extra entropy etc
                    #
                    #
                    # Now for the states that this state interacts with
                    #
                    for pKa2 in self.pKas:
                        pKaGroup2 = pKa2.pKaGroup
                        for titration2 in pKaGroup2.pKaTitrations:
                            states2=titration2.allstates
                            for state2 in [states2[1]]: #TODO: Modify to use all states as above
                                linear.append(self.matrix[pKa][titration][state][pKa2][titration2][state2])
        #print intpkas
        print 'Linearized matrix',linear
        mcsteps=50000
        phstart=2.0
        phend=12.0
        phstep=0.1
        FAST=pMC.MC(intpkas,linear,acidbase)
        FAST.set_MCsteps(int(mcsteps))
        pKavals=FAST.calc_pKa(phstart,phend,phstep)
        count=0
        print
        print
        print 'Final pKa values'
        print
        for pKa in self.pKas:
            pKaGroup = pKa.pKaGroup
            for titration in pKaGroup.pKaTitrations:
                titration.pKa=pKavals[count]
                count=count+1
                print pKa.residue.get('name'),titration.pKa
        return

    #
    # ----------------------------------
    #
        
    def calculateIntrinsicpKa(self):
        #
        #  Calculate the intrinsic pKa
        #
        self.calculateDesolvation()
        self.calculateBackground()
        #
        # Calculate the intrinsic pKas
        #
        # Print what we got
        #
        for pKa in self.pKas:
            print 'State\tModel pKa\tDesolvation\tBackground'
            for titration in pKa.pKaGroup.pKaTitrations:
                for state in titration.allstates:
                    print '%3d\t%5.3f\t\t%5.3f\t\t%5.3f' %(state,titration.modelpKa,pKa.desolvation[state],pKa.background[state])
        print
        print
        #
        # We calculate an intrinsic pKa for every possible <startstate> -> <endstate> transition
        #
        import math
        ln10=math.log(10)
        for pKa in self.pKas:
            pKaGroup=pKa.pKaGroup
            for titration in pKaGroup.pKaTitrations:
                for transition_num in titration.transitions.keys():
                    transition=titration.transitions[transition_num]
                    #
                    # Get the start and end states for this transition
                    #
                    start_s=transition['start']
                    end_s=transition['end']
                    print 'Transition: %3d, start: %3d, end: %3d' %(transition_num,start_s,end_s)
                    intpKa=titration.modelpKa+(pKa.desolvation[end_s]-pKa.desolvation[start_s])+(pKa.background[end_s]-pKa.background[start_s])/ln10
                    print 'Intrinsic pKa: %5.3f' %intpKa
                    titration.intrinsic_pKa=intpKa
        return

    #
    # --------------------------------
    #

    def calculateBackground(self):
        """
            Calculate background interaction energies
        """
        for pKa in self.pKas:
            residue = pKa.residue
            pKaGroup = pKa.pKaGroup
            ambiguity = pKa.amb
            
            print "Finding Background Interaction Energy for %s %s" %(residue.name, residue.resSeq)
            #
            # Loop over all titrations in this group
            #
            for titration in pKaGroup.pKaTitrations:
                #
                # Get all the states for this titration
                #
                possiblestates = titration.startstates + titration.endstates        
                fixedstates = self.hydrogenRoutines.getstates(ambiguity)
                #
                # Loop over all states and calculate the Backgound Interaction energy for each
                #
                for state in possiblestates:
                    print "Calculating Backgound for state %i" % state
                    self.hydrogenRoutines.switchstate(fixedstates, ambiguity, state) 
                    self.zeroAllRadiiCharges()
                    self.setAllRadii()
                    #
                    """ JENS: SET CHARGES ON ALL BUT THE CURRENT RESIDUE? """
                    """Todd: Yes"""
                    #
                    # Set charges on all other residues
                    #
                    # - TODO: Make sure all other groups are in their neutral state
                    # - TODO: Here we also need to explore the number of neutral states
                    # - TODO: that all other groups can be in - if we want to do things 100%
                    # - TODO: correctly.
                    # 
                    for chain in self.protein.getChains():
                        for otherresidue in chain.get("residues"):
                            if residue == otherresidue: continue
                            otherlist = []
                            for atom in otherresidue.atoms:
                                otherlist.append(atom.get("name"))
                            self.setCharges(otherresidue, otherlist)
                    #
                    # Run APBS
                    #
                    potentials=self.getAPBSPotentials()
                    #
                    # Get the atoms where we will measure the potential
                    #
                    firststate = possiblestates[0]
                    atomnames = self.getAtomsForPotential(pKa)
                    print 'Atoms for measuring potential',atomnames
                    atomlist=[]
                    for atomname in atomnames:
                        atomlist.append(residue.getAtom(atomname))
                    #
                    # Assign charges to our residue
                    #
                    self.setCharges(residue, atomnames)
                    #
                    # Return the potentials - same order as in atomnames
                    #
                    energy=self.get_elec_energy(potentials,atomlist)
                    #
                    # Done with Background calc for this state
                    #
                    pKa.background[state] = energy
        return

    #
    # --------------------------------
    #
                    
    def calculateDesolvation(self):
        #
        #   Calculate the Desolvation Energies
        #
        for pKa in self.pKas:
            residue = pKa.residue
            pKaGroup = pKa.pKaGroup
            ambiguity = pKa.amb
            print "Calculating Desolvation Energy for %s %s" %(residue.name, residue.resSeq)
            for titration in pKaGroup.pKaTitrations:
                possiblestates = titration.allstates
                #
                # Get atoms for potential
                #
                atomnames = self.getAtomsForPotential(pKa)
                print 'Atoms for measuring potential',atomnames
                atomlist=[]
                for atomname in atomnames:
                    atomlist.append(residue.getAtom(atomname))
                #
                # Calculate the self energy for each state
                #        
                fixedstates = self.hydrogenRoutines.getstates(ambiguity)
                for state in possiblestates:
                    print "Calculating self energy for state %i" % state
                    #
                    # Switch to the state
                    # Assign, radii, charges
                    self.hydrogenRoutines.switchstate(fixedstates, ambiguity, state) 
                    self.zeroAllRadiiCharges()
                    self.setCharges(residue, atomnames)
                    self.setRadii(residue, atomnames)
                    #
                    # Run APBS first time for the state in solvent
                    #
                    solutionEnergy=self.get_elec_energy(self.getAPBSPotentials(),atomlist)
                    print "SOLUTION ENERGY FROM APBS: ", solutionEnergy
                    #
                    # Now we set all radii (= in protein)
                    #
                    self.setAllRadii()
                    #
                    # Run APBS again, - this time for the state in the protein
                    #
                    print
                    print 'Calculating self energy for %s in the protein' %(residue.name)
                    proteinEnergy = self.get_elec_energy(self.getAPBSPotentials(),atomlist)
                    print "PROTEIN ENERGY FROM APBS: ", proteinEnergy
                    #
                    # Calculate the difference in self energy for this state
                    #
                    desolvation = proteinEnergy - solutionEnergy
                    print 'Desolvation for %s %d in state %d is %5.3f'  \
                          %(residue.name,residue.resSeq,state,desolvation)
                    pKa.desolvation[state] = desolvation
        return

    #
    # -------------------------------
    #

    def get_elec_energy(self,potentials,atomlist):
        #
        # Given the electrostatic potential from getAPBSPotentials and a list
        # of atoms, this routine returns the energy in kT
        #
        # Todd: What happens if we try to get the charge for a non-existing atom? - no problem or is
        # there a "ghost" in protein.getAtoms()?
        #
        energy=0.0
        count=0
        for atom in protein.getAtoms():
            if atom in atomlist:
                energy=energy+(potentials[3][count] - potentials[1][count])*atom.get("ffcharge")
                #print atom.name,atom.get('ffcharge')
            count=count+1
        return energy

    #
    # ----------------------------------
    #

    def getAPBSPotentials(self):
        #
        #    Run APBS and get the potentials from atom
        #
        #    Parameters
        #        residue:   The residue to examine (residue)
        #        atomnames: A list of atomnames (list)
        #    Returns
        #        list of potentials (list of floats)
        #
        # Todd: I changed this since we in some cases need to assign
        # charges to a residue that was not charged in the PBE run
        #
        protein = self.protein
        #
        # Do it
        #
        potentials = runAPBS(protein, self.apbs_inputfile)
        #
        #for atomname in atomnames:
        #    atomlist.append(residue.getAtom(atomname))
        #
        # Return the potentials - same order as in atomnames
        #
        #potentials=[]
        #for atom in protein.getAtoms():
        #    if atom in atomlist:
        #        potentials.append(potentials[3][count] - potentials[1][count])
        #energy = 0.5 * energy * getUnitConversion()
        #
        # Todd: We keep energies in kT since it's convenient. 
        #
        return potentials



    #
    # ----------------------
    #
                
    def setRadii(self, residue, atomlist):
        """
            Set the radii for specific atoms in a residue

            Parameters
                residue:  The residue to set (residue)
                atomlist: A list of atomnames (list)
        """
        for atom in residue.get("atoms"):
            atomname = atom.get("name")
            if atomname not in atomlist: continue
            charge, radius = self.forcefield.getParams(residue, atomname)
            if radius != None:
                atom.set("radius", radius)
            else:
                text = "Could not find radius for atom %s" % atomname
                text += " in residue %s %i" % (residue.name, residue.resSeq)
                text += " while attempting to set radius!"
                raise ValueError, text

    #
    # ------------------------------------
    #
            
    def setCharges(self, residue, atomlist):
        """
            Set the charges for specific atoms in a residue
            
            Parameters
                residue:  The residue to set (residue)
                atomlist: A list of atomnames (list)
        """
        #print 'Setting charges'
        for atom in residue.get("atoms"):
            atomname = atom.get("name")
            if atomname not in atomlist: continue
            charge, radius = self.forcefield.getParams(residue, atomname)
            if charge != None:
                atom.set("ffcharge", charge)
                #print atom.name,charge
            else:
                text = "Could not find charge for atom %s" % atomname
                text += " in residue %s %i" % (residue.name, residue.resSeq)
                text += " while attempting to set charge!"
                raise ValueError, text
        #print 'done'
        return
    #
    # ----------------------------
    #

    def setAllRadii(self):
        """
            Set all radii for the entire protein
        """
        for chain in self.protein.getChains():
            for residue in chain.get("residues"):
                for atom in residue.get("atoms"):
                    atomname = atom.get("name")
                    charge, radius = self.forcefield.getParams(residue, atomname)
                    if radius != None:
                        atom.set("radius", radius)
                    else:
                        if residue.type != 2:
                            text = "Could not find radius for atom %s " % atomname
                            text +="in residue %s %i" % (residue.name, residue.resSeq)
                            text += " while attempting to set all radii!"
                            raise ValueError, text
    #
    # -------------------------------
    #
                        
    def zeroAllRadiiCharges(self):
        """
            Set all charges and radii for the protein to zero
        """
        for chain in self.protein.getChains():
            for residue in chain.get("residues"):
                for atom in residue.get("atoms"):
                    atom.set("ffcharge",0.0)
                    atom.set("radius",0.0)

    #
    # --------------------------------
    #
                
    def getAtomsForPotential(self, pKa):
        """
            Find the atoms that are needed for measuring the potential,
            only selecting atoms where the charge changes.

            JENS: For AMBER, this includes backbone atoms?
            TODD: Yep - if they've defined it that way then lets go with it. 

            Parameters
                pKa:  The pKa object (pKa)
            Returns:
                atomnames:  A list of atomnames to measure (list)
        """
        #
        # TODO: This routine needs to become independent of the state
        # TODO: and we need to pass it a titration instance, not a pKa instance
        # TODO: I've hacked a temporary solution below...
        #
        atomnames = []
        initialmap = {}
        residue = pKa.residue
        pKaGroup = pKa.pKaGroup
        ambiguity = pKa.amb
        states = self.hydrogenRoutines.getstates(ambiguity)
        #
        # Change to the start state
        #
        start_state=0
        self.hydrogenRoutines.switchstate(states, ambiguity, start_state)
        for atom in residue.atoms:
            atomname = atom.get("name")
            charge, radius = self.forcefield.getParams(residue, atomname)
            initialmap[atomname] = charge
        #
        # Change to another state
        #
        end_state=1
        self.hydrogenRoutines.switchstate(states, ambiguity, end_state)
        for atom in residue.atoms:
            atomname = atom.get("name")
            charge, radius = self.forcefield.getParams(residue, atomname)
            try:
                initcharge = initialmap[atomname]
                if charge != initcharge:
                    atomnames.append(atomname)
            except KeyError:     
                atomnames.append(atomname)
        return atomnames

    #
    # -------------------------------
    #


    def findTitratableGroups(self):
        """
            Find all titratable groups in the protein based on the definition

            Returns
                pKalist:  A list of pKa objects (list)
        """
        pKalist = []
        
        print "Finding Titratable residues:"
        for chain in self.protein.getChains():
            for residue in chain.get("residues"):
                resname = residue.get("name")
                for group in self.pKagroups:
                    if resname == group:
                        amb = None
                        for hydrodef in self.hydrogenRoutines.hydrodefs:
                            hydname = hydrodef.name
                            if hydname == group:
                                amb = hydrogenAmbiguity(residue, hydrodef,self.routines)
                        if amb == None:
                            text = "Could not find hydrogen ambiguity "
                            text += "for titratable group %s!" % group
                            raise ValueError, text
                             
                        thispKa = pKa(residue, self.pKagroups[group], amb)
                        pKalist.append(thispKa)
                        print "\t%s %s" % (resname, residue.resSeq)
        #
        # Print the residues that we have selected
        #
        print
        print
        print 'Titratable residues'
        for pKa_v in pKalist:
            print pKa_v.residue.name,pKa_v.residue.resSeq
        print
        print


        return pKalist

    #
    # ----------------------------------
    #
        
    def readTitrationDefinition(self):
        """
            Read the Titration Definition

            Returns:
               mygroups: A dictionary of pKaGroups
        """
        mygroups = {}
        filename = TITRATIONFILE
        if not os.path.isfile(TITRATIONFILE):
            raise ValueError, "Could not find TITRATION.DAT!"
        file = open(filename)
        
        while 1:
            line=file.readline()
            if line.startswith("//"): pass
            elif line == '': break
            elif line[0]=='*':
                name = ""
                resname = ""
                type = ""
                titrations = []

                name = string.strip(line[1:])
                line = file.readline()
                if line[:8] != 'Residue:':
                    text = "Wrong line found when looking for 'Residue'"
                    raise ValueError, "%s: %s" % (text, line)
                
                resname = string.strip(string.split(line)[1])
            
                line = file.readline()
                if line[:10] != 'Grouptype:':
                    text = "Wrong line found when looking for 'Grouptype'"
                    raise ValueError, "%s: %s" % (text, line)
                
                type = string.lower(string.strip(string.split(line)[1]))
                if type != 'acid' and type != 'base':
                    raise ValueError, 'Group type must be acid or base!'

                line = file.readline()
                while 1:  
                    """ Find next transition """
                    
                    startstates = []
                    endstates = []
                    modelpKa = None

                    if line[:11] != 'Transition:':
                        text = "Wrong line found when looking for 'Transition'"
                        raise ValueError, "%s: %s" % (text, line)
                    
                    split=string.split(line[11:],'->')
                    for number in string.split(split[0], ','):
                        startstates.append(int(number))                        
                    for number in string.split(split[1], ','):
                        endstates.append(int(number))

                    line = file.readline()
                    if line[:10]!='Model_pKa:':
                        text = "Wrong line found when looking for 'Model_pKa'"
                        raise ValueError, "%s: %s" % (text, line)
                       
                    modelpKa = float(string.split(line)[1])

                    thisTitration = pKaTitration(startstates, endstates,modelpKa)
                    titrations.append(thisTitration)
                    
                    line = file.readline()
                    if string.strip(line) == 'END': break

                thisGroup = pKaGroup(name, resname, type, titrations)
                mygroups[name] = thisGroup

                line = file.readline()
                if string.strip(line) == 'END OF FILE': break
        
        return mygroups

#
# -------------------------------------------
#

class pKa:
    """
        The main pKa object
    """
    def __init__(self, residue, group, amb):
        """
            Initialize the pKa object

            Parameters
                residue: The residue object (residue)
                group:   The pKaGroup object associated with the residue
                         (pKaGroup)
                amb:     The associated hydrogenAmbiguity object
                         (hydrogenAmbiguity)
        """
        self.residue = residue
        self.pKaGroup = group
        self.amb = amb
        self.desolvation = {}
        self.background = {}
        self.interactionEnergies = {}
        self.intrinsicpKa = {}
        self.calculatedpKa = {}

#
# -------------------------------------------
#
        
class pKaGroup:
    #
    # pKaGroup holds info on a single titratable entity. In most cases we will
    # only have a single titration taking place in one group, but in some cases
    # we might have to have several transitions in one group (e.g. His - -> His 0 -> His +)
    #
    # The name, resname and type should probably be in pKaTitrations, but for now I'll leave
    # them here to get something working before we do complicated things...
    #

    
    def __init__(self, name, resname, type, pKaTitrations):
        """
            Initialize the pKaGroup object

            Parameters
                name:    The name of the group (string)
                resname: The residue name (string)
                type:    The type of group, acid or base (string)
                pKaTitrations: A list of pKaTitration objects (list)

        """
        self.name = name
        self.resname = resname
        self.type = type
        self.pKaTitrations = pKaTitrations
        return

    #
    # ------------------------
    #
       
    def __str__(self):
        """
            Print the pKa group object for debugging purposes

            Returns
                text:  The pKaGroup information (string)
        """
        text  = "Group name:   %s\n" % self.name
        text += "Residue name: %s\n" % self.resname
        text += "Group type:   %s\n" % self.type
        text += "Transitions:\n"
        for tran in self.pKaTitrations:
            text += str(tran)
        return text


#
# -----------------------------------------------
#

class pKaTitration:
    #
    # pKa_Titration holds all the info on a specific titration
    # We define a titration as a single group that has a number of
    # startstates and a number of endstates which is modelled by a
    # single model pKa value
    # A single group can have several transitions depending on the
    # number of startstates and endstates
    #
    
    def __init__(self, startstates, endstates, modelpKa):
        #
        #    Initialize the pKaTransition object
        #
        #    Parameters
        #        startstates: A list of state numbers (list)
        #        endstates:   A list of state numbers (list)
        #        modelpKa:    The model pKa associated with this titration
        #                     (float)
        #        transitions: A dictionary of the possible transitions for this group
        #                     (dictionary)
        #        interactions: A dictionary of the interaction energies with all other states
        #                      of all other titrations in the molecule. Only part of them
        #                      will be calculated
        #
        self.startstates = startstates
        self.endstates = endstates
        self.allstates=startstates+endstates
        self.modelpKa = modelpKa
        self.intrinsic_pKa=0.0
        #
        # Set transitions
        #
        self.transitions={}
        count=0
        for start_s in self.startstates:
            for end_s in self.endstates:
                count=count+1
                self.transitions[count]={'start':start_s,'end':end_s}
        #
        # Interactions has to be set at a higher level
        #
        return

    #
    # ---------------------
    #

    def __str__(self):
        """
            Print the pKa Transition object for debugging purposes

            Returns
                text:  The pKaTransition information (string)
        """
        text  = "\tStartstates: %s\n" % self.startstates
        text += "\tEndstates:   %s\n" % self.endstates
        text += "\tmodelpKa:    %.1f\n" % self.modelpKa
        return text


#
# -------------
#

def usage(x):
    #
    # Print usage guidelines
    #
    print 'Usage: pka.py --ff <forcefield> <pdbfile>'
    print 'Force field can be amber, charmm and parse'
    print
    return

#
# --------------------------------------------------
#

def startpKa():
    """
        Function for starting pKa script from the command line.

        Returns
            protein:    The protein object as generated by PDB2PQR
            routines:   The routines object as generated by PDB2PQR
            forcefield: The forcefield object as generated by PDB2PQR
    """
    print
    print 'PDB2PQR pKa calculations'
    print
    shortOptlist = "h,v"
    longOptlist = ["help","verbose","ff="]

    try:
        opts, args = getopt.getopt(sys.argv[1:], shortOptlist, longOptlist)
    except getopt.GetoptError, details:
        sys.stderr.write("GetoptError:  %s\n" % details)
        usage(2)
        sys.exit(0)
    #
    #
    #
    if len(args) < 1 or len(args) > 2:
        sys.stderr.write("Incorrect number (%d) of arguments!\n" % len(args))
        usage(2)
        sys.exit(0)

    verbose = 0
    ff = None
    for o,a in opts:
        if o in ("-v","--verbose"):
            verbose = 1
        elif o in ("-h","--help"):
            usage(2)
            sys.exit()
        elif o == "--ff":
            if a in ["amber","AMBER","charmm","CHARMM","parse","PARSE"]:
                ff = string.lower(a)
            else:
                raise ValueError, "Invalid forcefield %s!" % a
    #
    # No forcefield?
    #
    if ff == None:
        raise ValueError, "Forcefield not specified!"
    #
    # Get the PDBfile
    #
    path = args[0]
    file = getFile(path)
    pdblist, errlist = readPDB(file)
    
    if len(pdblist) == 0 and len(errlist) == 0:
        print "Unable to find file %s!\n" % path
        os.remove(path)
        sys.exit(2)

    if len(errlist) != 0 and verbose:
        print "Warning: %s is a non-standard PDB file.\n" % path
        print errlist

    if verbose:
        print "Beginning PDB2PQR...\n"

    myProtein = Protein(pdblist)
    if verbose:
        print "Created protein object -"
        print "\tNumber of residues in protein: %s" % myProtein.numResidues()
        print "\tNumber of atoms in protein   : %s" % myProtein.numAtoms()

    myDefinition = Definition()
    if verbose:
        print "Parsed Amino Acid definition file."

    myRoutines = Routines(myProtein, verbose, myDefinition)                  
    myRoutines.updateResidueTypes()
    myRoutines.updateSSbridges()
    myRoutines.updateExtraBonds()
    myRoutines.correctNames() 
    myRoutines.findMissingHeavy()
    myRoutines.addHydrogens()
    #myRoutines.randomizeWaters()
    myProtein.reSerialize()
    myForcefield = Forcefield(ff)
    #
    # Create the APBS input file
    #
    import src.psize
    size=src.psize.Psize()

    method=""
    async=0
    split=0
    import src.inputgen
    igen = src.inputgen.inputGen(path, size, method, async)
    apbs_inputfile=igen.printInput()
    return myProtein, myRoutines, myForcefield,apbs_inputfile

#
# -----------------------------------------------
#
 
            
if __name__ == "__main__":
    protein, routines, forcefield,apbs_inputfile = startpKa()
    mypkaRoutines = pKaRoutines(protein, routines, forcefield,apbs_inputfile)
    mypkaRoutines.runpKa()
