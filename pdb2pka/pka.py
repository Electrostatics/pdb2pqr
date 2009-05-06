#!/usr/bin/python
#
# pKa calculations with APBS
#
# Jens Erik Nielsen & Todd Dolinsky
#
# Copyright University College Dublin & Washington University St. Louis 2004-2007
# All rights reserved
#
__date__="22 April, 2009"
__author__="Jens Erik Nielsen, Todd Dolinsky, Yong Huang, Tommy Carstensen"

debug=None
import getopt
import sys, os
from pKa_base import *

if debug:
    from Tkinter import *
    from charge_mon import *

    CM=charge_mon()
else:
    CM=None
#
# find the path to the script, and to pdb2pqr
#
print __file__
try:
    file_name=__file__
    if file_name[:2]=='./':
        scriptpath=os.getcwd()
    else:
        scriptpath=os.path.join(os.getcwd(),os.path.split(file_name)[0])
except:
    scriptpath=os.path.split(sys.argv[0])[0]
    if scriptpath=='.':
        scriptpath=os.getcwd()
print 'SC',scriptpath
                 
pdb2pqr_path=os.path.split(scriptpath)[0]
print 'PDB2PQRpath',pdb2pqr_path
sys.path.append(pdb2pqr_path)
#
# Set the phidir - where results of apbscalcs are stored
#
phidir=os.path.join(os.getcwd(),'phidir')
if not os.path.isdir(phidir):
    os.mkdir(phidir)

#
# A dirty global variable
#
pdbfile_name='Not set'
#
# --
#
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

import ligandclean.ligff
from apbs import *

##
# ---
# Function for comparing two atoms
#

def is_sameatom(atom1,atom2):
    #
    # Compare atom1 and atom2
    #
    properties=['name','resSeq','chainID']
    for attr in properties:
        a1_prop=getattr(atom1,attr,None)
        a2_prop=getattr(atom2,attr,None)
        if (attr!='chainID' and (not a1_prop or not a2_prop)) or a1_prop!=a2_prop:
            return None
    return 1

#
# ----
#

TITRATIONFILE = os.path.join(scriptpath,"TITRATION.DAT")

class pKaRoutines:
    """
        Class for running all pKa related functions
    """
    def __init__(self, protein, routines, forcefield,apbs_setup):
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
        self.apbs_setup=apbs_setup
        self.pKagroups = self.readTitrationDefinition()
        
        self.pKas = []

        myHydrogenRoutines = hydrogenRoutines(routines)
        self.hydrogenRoutines = myHydrogenRoutines
        #
        # Not sure this is the best place for the interaction energies...
        #
        self.matrix={}
        #
        # Set the verbosity level
        #
        self.verbose=None
        
        return
    
    #
    # -----------------------------------------
    #
    

    def insert_new_titratable_group(self,ligand_titratable_groups):
        """Insert the new titratable groups in to self.pkagroups"""
        group_type=ligand_titratable_groups['type']
        if self.pKagroups.has_key(group_type):
            #
            # Now modify the group so that it will correspond to the group
            # we have in the ligand
            #
            ligand_name='LIG' #Note: we have to implement automatic determination of ligand name
            import copy
            new_group=copy.deepcopy(self.pKagroups[group_type])
            new_group.DefTitrations[0].modelpKa=ligand_titratable_groups['modelpka']
            new_group.name='LIG'
            new_group.resname='LIG'
            self.pKagroups['LIG']=copy.deepcopy(new_group)
            atom_map=ligand_titratable_groups['matching_atoms']
            #
            # Insert definition into HYDROGEN arrays
            #
            for hdef in self.hydrogenRoutines.hydrodefs:
                if hdef.name==group_type:
                    newdef=copy.deepcopy(hdef)
                    newdef.name=ligand_name
                   
                    #
                    # Change the names in each of the conformatinos
                    #
                    # The name of the H is not changed!
                    #
                    for conformation in newdef.conformations:
                        #
                        # Change the name of the atom that the H is bound to
                        #
                        if atom_map.has_key(conformation.boundatom):
                            conformation.boundatom=atom_map[conformation.boundatom]
                        #
                        # Change the name of the hydrogen
                        #
                        oldhname=conformation.hname
                        conformation.hname='H'+conformation.boundatom
                        #
                        # And then for the individual atom names
                        #
                        for atom in conformation.atoms:
                            if atom_map.has_key(atom.name):
                                atom.name=atom_map[atom.name]
                            elif atom.name==oldhname:
                                atom.name=conformation.hname
                    self.hydrogenRoutines.hydrodefs.append(copy.deepcopy(newdef))
        return

    #
    # -----------------------------------------
    #
        
    def runpKa(self,ghost=None):
        #
        #    Main driver for running pKa calculations
        #
        self.findTitratableGroups()
        
        if ghost:
            """ Calculate Pairwise Interactions """
            potentialDifference=self.calculatePotentialDifferences()
            return potentialDifference
                        
        else:
            self.calculateIntrinsicpKa()
    
            """ Calculate Pairwise Interactions """
            self.calculatePairwiseInteractions()
    
            """ Calculate Full pKa Value """
            self.calculatepKaValue()
            return

    #
    # -----------------------------------
    #
    def calculatePotentialDifferences(self):
        #
        # calculate potential difference of backbone atoms when each titratable group is set to charged and neutral states
        #
        import cPickle
        potentialDifference={}
        for pKa in self.pKas:
            self.get_interaction_energies_setup(pKa,mode='pKD')
            all_potentials=self.all_potentials[pKa].copy()
            residue = pKa.residue
            pKaGroup = pKa.pKaGroup
            ambiguity = pKa.amb
            
            titgroup='%s:%s:%s' %(residue.chainID, string.zfill(residue.resSeq,4),pKaGroup.name)
            if not potentialDifference.has_key(titgroup):
                potentialDifference[titgroup]={}
                for atom in self.protein.getAtoms():
                    if atom.name in ['N','H','C']:
                        atom_uniqueid=atom.chainID+':'+string.zfill(atom.resSeq,4)+':'+atom.name
                        potentialDifference[titgroup][atom_uniqueid]=0.00
                
            #
            #loop over each titration
            #
            for titration in pKaGroup.DefTitrations:
                startstates = titration.startstates
                endstates = titration.endstates
                possiblestates = titration.allstates
                #
                #loop over each state
                #                           
                for state in possiblestates:
                    for atom in self.protein.getAtoms():
                        if atom.name in ['N','H','C']:
                            atom_uniqueid=atom.chainID+':'+string.zfill(atom.resSeq,4)+':'+atom.name
                            if state in startstates:
                                potentialDifference[titgroup][atom_uniqueid]-=all_potentials[titration][state][pKa][titration][state][atom_uniqueid]/len(startstates)
                            if state in endstates:
                                potentialDifference[titgroup][atom_uniqueid]+=all_potentials[titration][state][pKa][titration][state][atom_uniqueid]
        
        return potentialDifference
    #
    # -----------------------------------
    #
    def calculatePairwiseInteractions(self):
        #
        # Calculate the pairwise interaction energies
        #
        for pKa in self.pKas:
            self.get_interaction_energies_setup(pKa)
        return

    #
    # -----
    #

    def get_interaction_energies_setup(self,pKa,mode='pkacalc'):
        """Perform the setup for the energy calculation"""
        residue = pKa.residue
        pKaGroup = pKa.pKaGroup
        ambiguity = pKa.amb

        #
        # Loop over each titration
        #
        self.all_potentials={}
        if not self.matrix.has_key(pKa):
            self.matrix[pKa]={}
            self.all_potentials[pKa]={}
        #
        for titration in pKaGroup.DefTitrations:
            if not self.matrix[pKa].has_key(titration):
                self.matrix[pKa][titration]={}
                self.all_potentials[pKa][titration]={}
            #
            # Get the atomnames
            #
            atomnames = self.getAtomsForPotential(pKa,titration)
            #
            # Get all states
            #
            possiblestates = titration.allstates

            for state in possiblestates:
                #
                # Here we switch the center group to a particular state
                #
                self.hydrogenRoutines.switchstate('pKa', ambiguity, state) 
                intenename=pdbfile_name+'.intene_%s_%s_%s_%s' %(titration.name,pKa.residue.chainID,pKa.residue.resSeq,state)
                import os
                if not os.path.isfile(intenename):

                    myRoutines = Routines(self.protein, 0)
                    myRoutines.updateResidueTypes()
                    myRoutines.updateSSbridges()
                    myRoutines.updateBonds()
                    myRoutines.updateInternalBonds()

                    pKa.residue.fixed = 2

                    self.hydrogenRoutines.initializeFullOptimization()

                    self.hydrogenRoutines.optimizeHydrogens()

                    myRoutines.debumpProtein()

                self.zeroAllRadiiCharges()
                self.setAllRadii()
                self.setCharges(residue, atomnames)
                #
                # get_interaction_energies get the potential at all titratable groups due the charges
                # this group
                #
                self.matrix[pKa][titration][state],self.all_potentials[pKa][titration][state]=self.get_interaction_energies(pKa,titration,state,mode)

        return

    #
    # ----
    #

    def get_interaction_energies(self,pKa_center,titration_center,state_center,mode):
        """Get the potentials and charges at all titratable groups"""
        print '------------>Charge - charge interactions for group: %s, state: %s' %(pKa_center.residue.resSeq,state_center)
        intenename=pdbfile_name+'.intene_%s_%s_%s_%s' %(titration_center.name,pKa_center.residue.chainID,pKa_center.residue.resSeq,state_center)
        allpotsname=pdbfile_name+'.intene_%s_%s_%s_%s_allpots' %(titration_center.name,pKa_center.residue.chainID,pKa_center.residue.resSeq,state_center)
        read_allpots=None
        #

        import os
        if os.path.isfile(intenename):
            fd=open(intenename)
            import cPickle
            savedict=cPickle.load(fd)
            fd.close()
            savedict_loaded=1
            #
            # pKD?
            #
            if mode=='pKD' and os.path.isfile(allpotsname):
                try:
                    import sys
                    sys.stdout.flush()
                    fd=open(allpotsname)
                    allsavedict=cPickle.load(fd)
                    fd.close()
                    read_allpots=1
                except EOFError:
                    print
                    print 'File %s is corrupt.\nDeleting and continuing...' %allpotsname
                    import os
                    os.unlink(allpotsname)
                    allsavedict={}
            else:
                allsavedict={}
        else:
            print 'Not found',intenename
            savedict={}
            savedict_loaded=None
            allsavedict={}
        #
        # Set the calc type and center
        #
        self.apbs_setup.set_type('intene')
        #
        # Run APBS and get the interaction with all other states
        #
        if debug:
            CM.set_calc('IE %s %s' %(pKa_center.residue.resSeq,state_center))
        if savedict=={} or (allsavedict=={} and mode=='pKD'):
            potentials=self.getAPBSPotentials(pKa_center,titration_center,state_center,cleanup=None)

        #
        # construct this side
        #
        energies={}
        all_potentials={}
        #
        # Loop over all groups
        #
        for pKa in self.pKas:
            residue = pKa.residue
            pKaGroup = pKa.pKaGroup
            ambiguity = pKa.amb
            #
            # Loop over each titration
            #
            if not energies.has_key(pKa):
                energies[pKa]={}
                all_potentials[pKa]={}
            #

            for titration in pKaGroup.DefTitrations:
                if not energies[pKa].has_key(titration):
                    energies[pKa][titration]={}
                    all_potentials[pKa][titration]={}
                #
                # Get all states
                #
                possiblestates = titration.allstates
                atomnames = self.getAtomsForPotential(pKa,titration)

                for state in possiblestates:
                    all_potentials[pKa][titration][state]={}
                    name='%s_%s_%s_%s' %(titration.name,pKa.residue.chainID,pKa.residue.resSeq,state)

                    if savedict.has_key(name):
                        energies[pKa][titration][state]= savedict[name]
                        if mode=='pKD':
                            if allsavedict.has_key(name):
                                all_potentials[pKa][titration][state]=allsavedict[name]
                                continue
                        else:
                            continue
                    elif savedict_loaded:
                        print 'Not found',name,pKa.residue.resSeq
                        print savedict.keys()
                        print 'Interaction energy file corrupt - delete all intene files and restart'
                        stop

                    #
                    # Switch the states of all other titratable groups to a neutral state
                    #
                    for other_pKa in self.pKas:
                        if pKa==other_pKa:
                            continue
                        for other_titration in other_pKa.pKaGroup.DefTitrations:
                            self.hydrogenRoutines.switchstate('pKa',other_pKa.amb,self.neutral_ref_state[other_pKa][other_titration])

                    #
                    # Switch to the particular state
                    #

                    self.hydrogenRoutines.switchstate('pKa', ambiguity, state)

                    potentials=self.getmoreAPBSPotentials()

                    self.zeroAllRadiiCharges()
                    self.setAllRadii()
                    self.setCharges(residue, atomnames)

                    #
                    # Get atoms for potential
                    #
                    #print 'Atoms for measuring potential',atomnames
                    atomlist=[]
                    for atomname in atomnames:
                        if residue.getAtom(atomname) in atomlist: pass
                        else:
                            atomlist.append(residue.getAtom(atomname))

                    energy=0.0
                    count=0               
                    for atom in self.protein.getAtoms():
                        for atom2 in atomlist:
                            if is_sameatom(atom,atom2):
                                energy=energy+potentials[count]*atom.get("ffcharge")
                        count=count+1

                    #
                    # We set all energies with self to zero
                    #
                    if pKa==pKa_center:
                        energies[pKa][titration][state]=0.0
                    else:
                        energies[pKa][titration][state]=energy
                    #
                    # Save in dict
                    #
                    savedict[name]=energies[pKa][titration][state]
                    #
                    # If running in pKD mode then we also return all potentials
                    #
                    if mode=='pKD':
                        count=0
                        for atom in self.protein.getAtoms():
                            atom_uniqueid=atom.chainID+':'+string.zfill(atom.resSeq,4)+':'+atom.name
                            all_potentials[pKa][titration][state][atom_uniqueid]=potentials[count]
                            count=count+1
                        import copy
                        allsavedict[name]=all_potentials[pKa][titration][state]

        #
        # Get rid of APBS instance
        #
        if getattr(self,'APBS',None):
            try:
                self.APBS.cleanup()
                self.APBS=None
            except:
                pass
        #
        # Dump a pickle file
        #
        fd=open(intenename,'w')
        import cPickle
        cPickle.dump(savedict,fd)
        fd.close()
        if mode=='pKD' and not read_allpots:
            fd=open(allpotsname,'w')
            import cPickle
            cPickle.dump(allsavedict,fd)
            fd.close()
        return energies,all_potentials

    #
    # ----------------------------------
    #

    def calculatepKaValue(self):
        """
        #  Calculate the pKa Value
        #
        # We use a c++ class for the MC steps since it's a lot faster..
        """
        #
        # First we need to correct all the interaction energies to make sure that terms aren't
        # counted twice
        #
        correct_matrix=self.correct_matrix()
        #
        # Now calculate the titration curves
        #
        import pMC_mult
        #
        # Matrix needs to be linear for transport to c++
        #
        linear=[]
        intpkas=[]
        acidbase=[]
        state_counter=[]
        is_charged_state=[]
        #print
        #print
        #print 'Residue\tIntrinsic pKa'
        import math
        ln10=math.log(10)
        for pKa in self.pKas:
            pKaGroup = pKa.pKaGroup
            for titration in pKaGroup.DefTitrations:
                #
                # Acid/Base
                #
                if pKaGroup.type=='acid':
                    acidbase.append(-1)
                else:
                    acidbase.append(1)

                possiblestates = titration.allstates
                #
                # Store all the intrinsic pKa values
                #
                pos_states=possiblestates
                pos_states.sort()
                for state1 in pos_states:
                    #
                    # Is this a charged state?
                    #
                    crg=self.is_charged(pKa,titration,state1)
                    is_charged_state.append(crg)
                    intpkas.append(pKa.intrinsic_pKa[state1])
                #
                # Record the number of states for each titratable group
                #
                state_counter.append(len(possiblestates))
                #
                #
                for pKa2 in self.pKas:
                    pKaGroup2 = pKa2.pKaGroup
                    for titration2 in pKaGroup2.DefTitrations:
                        #
                        # Loop over all states for group 1
                        #
                        pos_states=possiblestates
                        pos_states.sort()
                        for state1 in pos_states:
                            #
                            # Loop over all the states that this state interacts with
                            #
                            states2=titration2.allstates
                            states2.sort()
                            for state2 in states2:
                                linear.append(correct_matrix[pKa][titration][state1][pKa2][titration2][state2])
        mcsteps=50000 #This should be adjusted depending on the size of the protein
        phstart=0.1
        phend=20.1
        phstep=0.1
        #
        # Call our little C++ module
        #
        print "acidbase: %s" % (acidbase)
        FAST=pMC_mult.MC(intpkas,linear,acidbase,state_counter,is_charged_state)
        FAST.set_MCsteps(int(mcsteps))
        print 'Calculating titration curves'
        pKavals=FAST.calc_pKas(phstart,phend,phstep)
        count=0
        #
        # done
        #
        print
        print
        print 'Final pKa values'
        print
        pkas={}
        for pKa in self.pKas:
            pKaGroup = pKa.pKaGroup
            Gtype=pKa.pKaGroup.type
            for titration in pKaGroup.DefTitrations:
                name=pKa.uniqueid
                pkas[name]={'pKa':pKavals[count]}
                print "name: %s, PKAS[name]: %s" % (name, pkas[name])
                pkas[name]['modelpK']=titration.modelpKa
                #
                # Find an uncharged reference state
                #
                ref_state=self.neutral_ref_state[pKa][titration]
                #
                all_states=titration.allstates
                all_states.sort()
                for state in all_states:
                    if self.is_charged(pKa,titration,state)==1:
                        dpKa_desolv=(pKa.desolvation[state]-pKa.desolvation[ref_state])/ln10
                        dpKa_backgr=(pKa.background[state]-pKa.background[ref_state])/ln10
                        #
                        # Make acid and base modifications
                        #
                        if Gtype=='base':
                            dpKa_desolv=-dpKa_desolv
                            dpKa_backgr=-dpKa_backgr

                pkas[name]['desolv']=dpKa_desolv
                pkas[name]['backgr']=dpKa_backgr

                print 'Desolvation',pKa.desolvation
                print 'Background',pKa.background

                possiblestates = titration.allstates
                #
                # Record the number of states for each titratable group
                #
                #state_counter.append(len(possiblestates))

                pos_states=possiblestates
                chg_intpkas=[]
                neut_intpkas=[]
                pos_states.sort()
                print 'States',pos_states
                for state in pos_states:
                    crg=self.is_charged(pKa,titration,state)
                    if abs(crg)>0.1:
                        chg_intpkas.append(pKa.intrinsic_pKa[state])
                    else:
                        neut_intpkas.append(pKa.intrinsic_pKa[state])
                    print 'State: %s, charge: %5.2f, intpka: %5.3f' %(state,crg,pKa.intrinsic_pKa[state])
              
                pkas[name]['intpka']=pKa.simulated_intrinsic_pKa
                pkas[name]['delec']=pKavals[count]-pkas[name]['intpka']
                #
                # Print
                #
                pKa.pKa=pKavals[count]
                count=count+1
                print 'Simulated intrinsic pKa: %5.3f, delec: %5.3f' %(pkas[name]['intpka'],pkas[name]['delec'])
                print '%s final pKa: %5.2f' %(pKa.uniqueid,pKa.pKa)
                print '============================================='
                print
        #raw_input('Continue?')
        #
        # Write a WHAT IF -style pKa file
        import pKaIO_compat
        X=pKaIO_compat.pKaIO()
        X.write_pka('%s.PKA.DAT' %pdbfile_name,pkas,format='pdb2pka')
        #
        # Desolv file % Backgr file
        #
        X.desolv={}
        X.backgr={}
        for name in pkas.keys():
            X.desolv[name]=pkas[name]['desolv']
            X.backgr[name]=pkas[name]['backgr']
        X.write_desolv('%s.DESOLV.DAT' %pdbfile_name,format='pdb2pka')
        X.write_backgr('%s.BACKGR.DAT' %pdbfile_name,format='pdb2pka')
        #
        # Get the charges
        #
        charges={}
        pH_start=pKavals[count]
        pH_step=pKavals[count+1]
        num_pHs=pKavals[count+2]
        count=count+2
        for pKa in self.pKas:
            pKaGroup = pKa.pKaGroup
            for titration in pKaGroup.DefTitrations:
                name=pKa.uniqueid
                charges[name]={}
                for x in range(int(num_pHs)):
                    count=count+1
                    pH=pKavals[count]
                    count=count+1
                    charges[name][pH]=pKavals[count]
                    pH=pH+pH_step
                if pKavals[count+1]==999.0 and pKavals[count+2]==-999.0:
                    count=count+2
                else:
                    print 'Something is wrong'
                    print pKavals[count:count+30]
                    stop
        #
        # Write a WHAT IF titration curve file
        #
        X.write_titration_curve('%s.TITCURV.DAT' %pdbfile_name,charges,format='pdb2pka')
        #
        # Write the charge matrix
        #
        X.write_pdb2pka_matrix('%s.MATRIX.DAT' %pdbfile_name, correct_matrix)
        return

    #
    # ----
    #

    def correct_matrix(self):
        """Correct the matrix so that all energies are correct, i.e.
        that we subtract the neutral-charged and charged-neutral interaction energies
        in the right places.
        After having fixed this we make the matrix symmetric"""
        corrected_matrix={}
        #
        # Loop over the original matrix
        #
        for pKa1 in self.pKas:
            corrected_matrix[pKa1]={}
            pKaGroup = pKa1.pKaGroup
            for titration1 in pKaGroup.DefTitrations:
                corrected_matrix[pKa1][titration1]={}
                #
                # Get the reference state for first state
                #
                ref_state1=self.neutral_ref_state[pKa1][titration1]
                
                possible_states = titration1.allstates
                possible_states.sort()
                for state1 in possible_states:
                    corrected_matrix[pKa1][titration1][state1]={}
                    for pKa2 in self.pKas:
                        corrected_matrix[pKa1][titration1][state1][pKa2]={}
                        pKaGroup2 = pKa2.pKaGroup
                        for titration2 in pKaGroup2.DefTitrations:
                            corrected_matrix[pKa1][titration1][state1][pKa2][titration2]={}
                            #
                            # Get the reference state for the second state
                            #
                            ref_state2=self.neutral_ref_state[pKa2][titration2]
                            states2=titration2.allstates
                            states2.sort()
                            for state2 in states2:
                                #
                                # Now figure out what the correct energy for the interaction is
                                #
                                if state1==ref_state1 and state2==ref_state2:
                                    # No interaction for ref-ref
                                    value=0.0
                                else:
                                    value=self.matrix[pKa1][titration1][state1][pKa2][titration2][state2]
                                    #
                                    # Subtract Eint(state1,ref2)+Eint(state2,ref1) and add Eint(ref1,ref2)
                                    #
                                    state1_ref2=self.matrix[pKa1][titration1][state1][pKa2][titration2][ref_state2]
                                    state2_ref1=self.matrix[pKa1][titration1][ref_state1][pKa2][titration2][state2]
                                    ref1_ref2=self.matrix[pKa1][titration1][ref_state1][pKa2][titration2][ref_state2]
                                    value=value-state1_ref2-state2_ref1+ref1_ref2
                                #
                                # Insert this value in the corrected matrix
                                #
                                #print pKa1.uniqueid,titration1.name,state1,pKa2.uniqueid,titration2.name,state2,value
                                corrected_matrix[pKa1][titration1][state1][pKa2][titration2][state2]=value
        #
        # Make matrix symmetric
        #
        #print
        #print 'Interaction energy matrix'
        #print '%25s %25s %10s %10s %16s' %('Group1','Group 2','State G1','State G2','Interaction energy (kT)')
        symmetric_matrix={}
        for pKa1 in self.pKas:
            symmetric_matrix[pKa1]={}
            pKaGroup = pKa1.pKaGroup
            for titration1 in pKaGroup.DefTitrations:
                symmetric_matrix[pKa1][titration1]={}
                possible_states = titration1.allstates
                possible_states.sort()
                for state1 in possible_states:
                    symmetric_matrix[pKa1][titration1][state1]={}
                    for pKa2 in self.pKas:
                        symmetric_matrix[pKa1][titration1][state1][pKa2]={}
                        pKaGroup2 = pKa2.pKaGroup
                        for titration2 in pKaGroup2.DefTitrations:
                            symmetric_matrix[pKa1][titration1][state1][pKa2][titration2]={}
                            #
                            # Get the reference state for the second state
                            #
                            states2=titration2.allstates
                            states2.sort()
                            for state2 in states2:
                                #
                                # Now figure out what the correct energy for the interaction is
                                #
                                value1=corrected_matrix[pKa1][titration1][state1][pKa2][titration2][state2]
                                value2=corrected_matrix[pKa2][titration2][state2][pKa1][titration1][state1]
                                #
                                # Insert the average value in the symetric matrix
                                #
                                #print '%25s %25s %10s %10s %6.3f' %(pKa1.uniqueid,pKa2.uniqueid,state1,state2,(value1+value2)/2.0)
                                symmetric_matrix[pKa1][titration1][state1][pKa2][titration2][state2]=(value1+value2)/2.0
        return symmetric_matrix

    #
    # ----------------------------------
    # 

    def is_charged(self,pKa,titration,state):
        """
        # Check if this state is a charged state
        """
        #
        # Check if there are any other titratable groups in this residue
        #
        for other_pKa in self.pKas:
            if pKa==other_pKa:
                continue
            for other_titration in other_pKa.pKaGroup.DefTitrations:
                self.hydrogenRoutines.switchstate('pKa',other_pKa.amb,self.neutral_ref_state[other_pKa][other_titration])
        #
        # Get charge
        #
        ambiguity = pKa.amb
        self.hydrogenRoutines.switchstate('pKa', ambiguity, state)
        residue = pKa.residue
        pKaGroup = pKa.pKaGroup
        sum=0.0
        for atom in residue.getAtoms():
            atomname = atom.get("name")
            if atomname.find('FLIP')!=-1:
                continue
            charge, radius = self.forcefield.getParams1(residue, atomname)
            sum=sum+charge
        if abs(sum)>0.05:
            return 1
        return 0

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
            print "======== Residue: %s ========" % (pKa.residue)
            print 'State\tModel pKa\tDesolvation\tBackground'
            for titration in pKa.pKaGroup.DefTitrations:
                for state in titration.allstates:
                    print '%10s\t%5.3f\t\t%5.3f\t\t%5.3f' %(state,titration.modelpKa,pKa.desolvation[state],pKa.background[state])
        print
        print
        #
        # We calculate an intrinsic pKa for every possible <startstate> -> <endstate> transition
        #
        import math
        ln10=math.log(10)
        for pKa in self.pKas:
            pKaGroup=pKa.pKaGroup
            Gtype=pKa.pKaGroup.type
            #
            # We measure intrinsic pKa values against a single reference state
            #
            for titration in pKaGroup.DefTitrations:
                #
                # Find an uncharged reference state
                #
                ref_state=self.neutral_ref_state[pKa][titration]
                #
                #
                #
                all_states=titration.allstates
                all_states.sort()
                for state in all_states:
                    if self.is_charged(pKa,titration,state)==1:
                        dpKa_desolv=(pKa.desolvation[state]-pKa.desolvation[ref_state])/ln10
                        dpKa_backgr=(pKa.background[state]-pKa.background[ref_state])/ln10
                        #
                        # Make acid and base modifications
                        #
                        if Gtype=='base':
                            dpKa_desolv=-dpKa_desolv
                            dpKa_backgr=-dpKa_backgr
                        #
                        # Now calculate intrinsic pKa
                        #
                        intpKa=titration.modelpKa+dpKa_desolv+dpKa_backgr   
                        print 'Energy difference for %s -> %s [reference state] is %5.2f pKa units' %(state,ref_state,intpKa)
                        pKa.intrinsic_pKa[state]=intpKa
                    else:
                        #
                        # Neutral states 
                        #
                        dpKa_desolv=(pKa.desolvation[state]-pKa.desolvation[ref_state])/ln10
                        dpKa_backgr=(pKa.background[state]-pKa.background[ref_state])/ln10
                        #
                        # Make acid and base modifications
                        #
                        if Gtype=='base':
                            dpKa_desolv=-dpKa_desolv
                            dpKa_backgr=-dpKa_backgr
                        dpKa=dpKa_desolv+dpKa_backgr
                        print 'Energy difference for %s -> %s [reference state] is %5.2f pKa units' %(state,ref_state,dpKa)
                        pKa.intrinsic_pKa[state]=dpKa
            # -----------------------------------------------------------------
            # Get the intrinsic pKa with a small MC calculation
            #
            acidbase=[]
            is_charged=[]
            intpKas=[]
            for titration in pKaGroup.DefTitrations:
                #
                # Acid/Base
                #
                if pKaGroup.type=='acid':
                    acidbase.append(-1)
                else:
                    acidbase.append(1)
                #
                # 
                #
                possiblestates = titration.allstates
                #
                # Record the number of states for each titratable group
                #
                pos_states=possiblestates
                pos_states.sort()
                for state in pos_states:
                    #
                    # Is this a charged state?
                    #
                    crg=self.is_charged(pKa,titration,state)
                    is_charged.append(crg)
                    intpKas.append(pKa.intrinsic_pKa[state])
            intpka=titrate_one_group(name='%s' %(pKa.residue),intpkas=intpKas,is_charged=is_charged,acidbase=acidbase)
            pKa.simulated_intrinsic_pKa=intpka
        return

    #
    # --------------------------------
    #

    def calculateBackground(self,onlypKa=None):
        #
        #    Calculate background interaction energies
        #
        backgroundname=pdbfile_name+'.background'
        if os.path.isfile(backgroundname):
            fd=open(backgroundname)
            import pickle
            savedict=pickle.load(fd)
            fd.close()
        else:
            savedict={}

        for pKa in self.pKas:
            if onlypKa:
                if not pKa==onlypKa:
                    continue
            residue = pKa.residue
            pKaGroup = pKa.pKaGroup
            ambiguity = pKa.amb

            print "-----> Finding Background Interaction Energy for %s %s" %(residue.name, residue.resSeq)
            #
            # Loop over all titrations in this group
            #
            for titration in pKaGroup.DefTitrations:
                #
                # Get all the states for this titration
                #
                possiblestates = titration.startstates + titration.endstates
                #
                # Loop over all states and calculate the Background Interaction energy for each
                #
                for state in possiblestates:
                    #
                    # Set the name for this energy
                    #
                    name='%s_%s_%s_%s' %(titration.name,pKa.residue.chainID,pKa.residue.resSeq,state)
                    if savedict.has_key(name):
                        pKa.background[state] = savedict[name]
                        continue
                    #
                    # Get the atoms where we will measure the potential
                    #
                    firststate = possiblestates[0]
                    atomnames = self.getAtomsForPotential(pKa,titration)
                    atomlist=[]
                    for atomname in atomnames:
                        atomlist.append(residue.getAtom(atomname))
                    #
                    # Switch the states of all other titratable groups to the neutral reference state
                    #
                    for other_pKa in self.pKas:
                        if pKa==other_pKa:
                            continue
                        for other_titration in other_pKa.pKaGroup.DefTitrations:
                            self.hydrogenRoutines.switchstate('pKa',other_pKa.amb,self.neutral_ref_state[other_pKa][other_titration])

                    #
                    # Switch the state for the group in question
                    #
                    print "----------> Calculating Background for state %s" % state
                    self.hydrogenRoutines.switchstate('pKa', ambiguity, state) 

                    myRoutines = Routines(self.protein, 0)
                    myRoutines.updateResidueTypes()
                    myRoutines.updateSSbridges()
                    myRoutines.updateBonds()
                    myRoutines.updateInternalBonds()

                    pKa.residue.fixed = 2

                    self.hydrogenRoutines.initializeFullOptimization()

                    self.hydrogenRoutines.optimizeHydrogens()

                    myRoutines.debumpProtein()

                    self.zeroAllRadiiCharges()
                    self.setAllRadii()
                    #
                    # Set charges on all other residues
                    #
                    for chain in self.protein.getChains():
                        for otherresidue in chain.get("residues"):
                            if residue == otherresidue:
                                continue
                            otherlist = []
                            for atom in otherresidue.atoms:
                                if not atom:
                                    continue
                                if atom.get('name').find('FLIP')==-1:
                                    otherlist.append(atom.get("name"))
                            self.setCharges(otherresidue, otherlist)

                    #
                    # Center the map on our residue
                    #
                    center=self.get_atoms_center(atomlist)

                    all_center,extent=self.apbs_setup.getCenter()

                    #
                    # For small proteins we set the center to the center of the molecule
                    #
                    if extent[0]<20.0 or extent[1]<20.0 or extent[2]<20.0:
                        self.apbs_setup.setfineCenter(all_center)
                    else:
                        self.apbs_setup.setfineCenter(center)
                    self.apbs_setup.set_type('background')
                    #
                    # Run APBS
                    #
                    if debug:
                        CM.set_calc('background %s %s' %(pKa.residue.resSeq,state))
                    potentials=self.getAPBSPotentials(pKa,titration,state)
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
                    #
                    # Save it under a unique name
                    #
                    print 'Saving energy as',name
                    savedict[name]=energy

        #
        # Dump a pickle file
        #
        fd=open(backgroundname,'w')
        import pickle
        pickle.dump(savedict,fd)
        fd.close()
        return


    #
    # --------------------------------
    #
                    
    def calculateDesolvation(self,onlypKa=None):
        """
        #
        #   Calculate the Desolvation Energies
        #
        """
        desolvname=pdbfile_name+'.desolv'
        if os.path.isfile(desolvname):
            fd=open(desolvname)
            import pickle
            savedict=pickle.load(fd)
            fd.close()
        else:
            savedict={}
        #
        #
        #
        for pKa in self.pKas:
            if onlypKa:
                if not pKa==onlypKa:
                    continue
            residue = pKa.residue
            pKaGroup = pKa.pKaGroup
            ambiguity = pKa.amb

            print "-----> Calculating Desolvation Energy for %s %s" %(residue.name, residue.resSeq)
            for titration in pKaGroup.DefTitrations:
                possiblestates = titration.allstates
                #
                # Get atoms for potential
                #
                atomnames = self.getAtomsForPotential(pKa,titration)
                atomlist=[]
                for atomname in atomnames:
                    atomlist.append(residue.getAtom(atomname))
                #
                # Switch all the other groups to the neutral reference state
                #
                for other_pKa in self.pKas:
                    if pKa==other_pKa:
                        continue
                    for other_titration in other_pKa.pKaGroup.DefTitrations:
                        self.hydrogenRoutines.switchstate('pKa',other_pKa.amb,self.neutral_ref_state[other_pKa][other_titration])
                #
                # Calculate the self energy for each state
                #        
                #fixedstates = self.hydrogenRoutines.getstates(ambiguity)
                for state in possiblestates:
                    name='%s_%s_%s_%s' %(titration.name,pKa.residue.chainID,pKa.residue.resSeq,state)
                    if savedict.has_key(name):
                        pKa.desolvation[state] = savedict[name]
                        continue
                    print "---------> Calculating desolvation energy for residue %s state %s in solvent" %(residue.name,state)
                    #
                    # Center the map on our set of atoms
                    #
                    center=self.get_atoms_center(atomlist)
                    self.apbs_setup.setfineCenter(center)
                    self.apbs_setup.set_type('desolv')
                    #
                    # Switch to the state
                    # Assign, radii, charges

                    self.hydrogenRoutines.switchstate('pKa', ambiguity, state) 

                    myRoutines = Routines(self.protein, 0)
                    myRoutines.updateResidueTypes()
                    myRoutines.updateSSbridges()
                    myRoutines.updateBonds()
                    myRoutines.updateInternalBonds()

                    pKa.residue.fixed = 2

                    self.hydrogenRoutines.initializeFullOptimization()

                    self.hydrogenRoutines.optimizeHydrogens()

                    myRoutines.debumpProtein()

                    self.zeroAllRadiiCharges()
                    self.setCharges(residue, atomnames)
                    self.setRadii(residue, atomnames)

                    #
                    # Run APBS first time for the state in solvent
                    #

                    if debug:
                        CM.set_calc('Desolv solv %s %s' %(pKa.residue.resSeq,state))

                    solutionEnergy=self.get_elec_energy(self.getAPBSPotentials(pKa,titration,state),atomlist)
                    #
                    # Now we set all radii (= in protein)
                    #
                    self.setAllRadii()
                    #
                    # Run APBS again, - this time for the state in the protein
                    #

                    print '--------> Calculating self energy for residue %s %d state %s in the protein' %(residue.name,residue.resSeq,state)

                    if debug:
                        CM.set_calc('Desolv prot %s %s' %(pKa.residue.resSeq,state))

                    proteinEnergy = self.get_elec_energy(self.getAPBSPotentials(pKa,titration,state),atomlist)
                    #
                    # Calculate the difference in self energy for this state
                    #
                    desolvation = (proteinEnergy - solutionEnergy)/2.0 # Reaction field energy
                    print 'Desolvation for %s %d in state %s is %5.3f'  \
                          %(residue.name,residue.resSeq,state,desolvation)
                    print
                    print '======================================='
                    pKa.desolvation[state] = desolvation
                    print 'Saving energy as',name
                    #stop
                    savedict[name]=desolvation
                    #
                    # Check for bumps and modify energy if needed
                    #
                    pass
        #
        # Dump a pickle file
        #
        fd=open(desolvname,'w')
        import pickle
        pickle.dump(savedict,fd)
        fd.close()
        return

    #
    # ---------
    #

    def get_atoms_center(self,atomlist):
        #
        # Get the centre of a list of atoms
        #
        minmax={'x':[999.9,-999.9],'y':[999.9,-999.9],'z':[999.9,-999.9]}
        for atom in atomlist:
            if atom:
                for axis in ['x','y','z']:
                    coord=getattr(atom,axis)
                    if coord<minmax[axis][0]:
                        minmax[axis][0]=coord
                    if coord>minmax[axis][1]:
                        minmax[axis][1]=coord
        #
        # Calc the geometric center and extent
        #
        center={}
        extent={}
        for axis in minmax.keys():
            extent[axis]=minmax[axis][1]-minmax[axis][0]
            center[axis]=extent[axis]/2.0+minmax[axis][0]
        return [center['x'],center['y'],center['z']]
        

    #
    # -------------------------------
    #

    def get_elec_energy(self,potentials,atomlist):
        #
        # Given the electrostatic potential from getAPBSPotentials and a list
        # of atoms, this routine returns the energy in kT
        #
        # This function could be made a lot smarter!! (JN)
        # 
        energy=0.0
        count=0
        totphi=0.0
        totcrg=0.0
        found=0
        for atom in self.protein.getAtoms():
            if not atom:
                continue
            for atom_2 in atomlist:
                if not atom_2:
                    continue
                if is_sameatom(atom,atom_2):
                    totcrg=totcrg+abs(atom.get("ffcharge"))
                    totphi=totphi+abs(potentials[-1][count])
                    energy=energy+(potentials[-1][count])*atom.get("ffcharge")
                    found=found+1
                    break
            #
            # This counter is outside the atom_2 loop!!
            #
            count=count+1
            if found==len(atomlist):
                break
        if abs(totphi)<0.01 or abs(totcrg)<0.01:
            print 'total abs phi',totphi
            print 'total abs crg',totcrg
            raise 'Something is rotten'
        return energy

    #
    # ----------------------------------
    #

    def getAPBSPotentials(self,group,titration,state,cleanup=1):
        """
        #    Run APBS and get the potentials 
        #
        #    Returns
        #        list of potentials (list of floats)
        #"""
        #
        # Do we have results for this calculation?
        #
        result_file=os.path.join(phidir,'%s_%s_%s.potentials' %(self.apbs_setup.type,group.uniqueid,state))
        import cPickle
        loaded=None
        if os.path.isfile(result_file):
            #
            # Yes!
            #
            fd=open(result_file,'rb')
            try:
                potentials=cPickle.load(fd)
                fd.close()
                loaded=1
            except:
                fd.close()
        loaded=None
        #
        # Run calc again if needed
        #
        if not loaded:
            apbs_inputfile=self.apbs_setup.printInput()
            self.APBS=runAPBS()
            potentials = self.APBS.runAPBS(self.protein, apbs_inputfile, CM)
            if cleanup:
                self.APBS.cleanup()
                self.APBS=None
            fd=open(result_file,'wb')
            cPickle.dump(potentials,fd)
            fd.close()
        return potentials

    #
    # -----
    #

    def getmoreAPBSPotentials(self):
        if not self.APBS:
            raise 'APBS instance killed'
        return self.APBS.get_potentials(self.protein)

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
        for atom in residue.getAtoms():
            atomname = atom.get("name")
            if atomname not in atomlist: continue
            charge, radius = self.forcefield.getParams1(residue, atomname)
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
        for atom in residue.getAtoms():
            atomname = atom.get("name")
            if atomname not in atomlist:
                continue
            charge, radius = self.forcefield.getParams1(residue, atomname)
            if charge != None:
                atom.set("ffcharge", charge)
            else:
                text = "Could not find charge for atom %s" % atomname
                text += " in residue %s %i" % (residue.name, residue.resSeq)
                text += " while attempting to set charge!"
                raise ValueError, text
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
                    if atomname.find('FLIP')!=-1:
                        continue
                    else:
                        charge, radius = self.forcefield.getParams1(residue, atomname)
                    ###PC
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
                
    def getAtomsForPotential(self, pKa,titration, get_neutral_state=None):
        #
        #    Find the atoms that are needed for measuring the potential,
        #    only selecting atoms where the charge changes.
        #    Parameters
        #        pKa:  The pKa object (pKa)
        #    Returns:
        #        atomnames:  A list of atomnames to measure (list)
        neutral_state=None
        atomnames = []
        newatomnames = []
        initialmap = {}
        residue = pKa.residue
        pKaGroup = pKa.pKaGroup
        ambiguity = pKa.amb
        #states = self.hydrogenRoutines.getstates(ambiguity)
        #
        # Change to the start state
        #
        start_state=titration.startstates[0]
        self.hydrogenRoutines.switchstate('pKa', ambiguity, start_state)
        sum=0.0
        for atom in residue.getAtoms():
            atomname = atom.get("name")
            if atomname.find('FLIP')!=-1:
                continue

            charge, radius = self.forcefield.getParams1(residue, atomname)
            initialmap[atomname] = charge
            if charge is None:
                print atomname,charge
                print residue.isCterm
                raise 'Charge on atom is None'
            sum=sum+charge
        if abs(sum)<0.001:
            neutral_state=start_state
        #
        # Check if charges change in all other states
        #

        for state in titration.endstates+titration.startstates[1:]:
            self.hydrogenRoutines.switchstate('pKa', ambiguity, state)
            #
            # Check that no charges changed and that no atoms were added
            #
            sum=0.0
            for atom in residue.getAtoms():
                atomname = atom.get("name")
                if atomname.find('FLIP')!=-1:
                    continue

                charge, radius = self.forcefield.getParams1(residue, atomname)
                sum=sum+charge
                if initialmap.has_key(atomname):
                    initcharge = initialmap[atomname]
                    if charge != initcharge:
                        if not atomname in atomnames:
                            atomnames.append(atomname)
                else:
                    if not atomname in atomnames:
                        atomnames.append(atomname)
            #
            # Check that no atoms were removed
            #
            for atom in initialmap.keys():
                if not atom in residue.get('map'):
                    atomnames.append(atom)
        #
        # Make sure that the charges add up to integers by adding extra atoms
        #
        sum=0.01
        while sum>0.001:
            sum=0.0
            added=None
            neutral_state=None
            #
            # Loop over all states to find atoms to add
            #
            for state in titration.endstates+titration.startstates:
                self.hydrogenRoutines.switchstate('pKa', ambiguity, state)
                #
                # Sum this state
                #
                this_sum=0.0
                for atom in residue.atoms:
                    atomname = atom.get("name")
                    if atomname.find('FLIP')!=-1:
                        continue
                    if not atomname in atomnames:
                        continue
                    charge, radius = self.forcefield.getParams1(residue, atomname)
                    this_sum=this_sum+charge
                #
                # Is this the first neutral state?
                #
                if abs(this_sum)<0.0001:
                    if not neutral_state:
                        neutral_state=state
                #
                # Is this an integer charge?
                #
                diff=float(abs(1000.0*this_sum)-abs(1000.0*int(this_sum)))/1000.0
                sum=sum+diff
                if diff>0.001:
                    #
                    # Find all atoms one bond away
                    #
                    add_atoms=[]
                    for atom in residue.atoms:
                        atomname=atom.get('name')
                        if atomname.find('FLIP')!=-1:
                            continue
                        if not atomname in atomnames:
                            continue
                        #
                        # Add all atoms that are not already added
                        #
                        for bound_atom in atom.bonds:
                            if type(bound_atom) is str:
                                if not bound_atom in atomnames:
                                    add_atoms.append(bound_atom)
                                    added=1
                            else:
                                if not bound_atom.name in atomnames:
                                    add_atoms.append(bound_atom.name)
                                    added=1
                    #
                    # Update atomnames
                    #
                    for addatom in add_atoms:
                        if not addatom in atomnames:
                            atomnames.append(addatom)
                #
                # Next state
                #
                pass
            #
            # Did we add anything?
            #
            if added is None and sum>0.001:
                print sum
                print atomnames
                raise 'Could not find integer charge state'
        #
        # Did we just want a neutral state identification?
        #
        if get_neutral_state:
            if not neutral_state:
                raise "no neutral state for",residue.resSeq
            return neutral_state
        #
        # No, we wanted the atomnames
        #
        if atomnames==[]:
            print 'Did not find any atoms for ',residue.resSeq
            raise 'Something wrong with charges'

        for atomname in atomnames:
            if not atomname in newatomnames: 
                newatomnames.append(atomname)
        return newatomnames

    #
    # -------------------------------
    #


    def findTitratableGroups(self):
        """
            Find all titratable groups in the protein based on the definition

            We do a simple name-matching on residue names and titration group names, and
            also a Cterm/Nterm matching.

            We need to build in checks for post-translational modifications
            
            Returns
                pKalist:  A list of pKa objects (list)
        """
        pKalist = []
        
        print "Finding Titratable groups....",
        import sys
        sys.stdout.flush()
        if self.verbose:
            print
        for chain in self.protein.getChains():
            for residue in chain.get("residues"):
                resname = residue.get("name")
                for group in self.pKagroups:
                    if resname == group:
                        amb=self.find_hydrogen_amb_for_titgroup(residue,group)
                        thispKa = pKa(residue, self.pKagroups[group], amb)
                        pKalist.append(thispKa)
                        if self.verbose:
                            print "\t%s %s" % (resname, residue.resSeq)
                    if residue.isNterm and group=='NTR':
                        #
                        # N-terminus
                        #
                        amb=self.find_hydrogen_amb_for_titgroup(residue,group)
                        thispKa=pKa(residue,self.pKagroups[group],amb)
                        pKalist.append(thispKa)
                        if self.verbose:
                            print "\t%s %s" % (resname, residue.resSeq)
                    if residue.isCterm and group=='CTR':
                        #
                        # C-terminus
                        #
                        amb=self.find_hydrogen_amb_for_titgroup(residue,group)
                        thispKa=pKa(residue,self.pKagroups[group],amb)
                        pKalist.append(thispKa)
                        if self.verbose:
                            print "\t%s %s" % (resname, residue.resSeq)
        #
        # Find a neutral state for each group
        #
        self.neutral_ref_state={}
        for this_pka in pKalist:
            residue = this_pka.residue
            pKaGroup = this_pka.pKaGroup
            ambiguity = this_pka.amb
            self.neutral_ref_state[this_pka]={}
            for titration in pKaGroup.DefTitrations:
                neutral_state = self.getAtomsForPotential(this_pka,titration,get_neutral_state=1)
                self.neutral_ref_state[this_pka][titration]=neutral_state
        #
        # Store pKa groups in self.pKas
        #
        self.pKas=pKalist
        print 'done'
        return 

    #
    # ----------------------------------
    #

    def find_hydrogen_amb_for_titgroup(self,residue,group):
        """Find the hydrogen ambiguity that controls the protonation state for the
        titratable group within the given residue"""
        amb = None
        self.hydrogenRoutines.readHydrogenDefinition()
        for hydrodef in self.hydrogenRoutines.hydrodefs:
            hydname = hydrodef.name
            if hydname == group: # or group in residue.patches: (for ASP/ASH, GLU/GLH)
                amb = hydrogenAmbiguity(residue, hydrodef,self.routines)
            elif group == 'ASP':
                if hydname == 'ASH':
                    amb = hydrogenAmbiguity(residue, hydrodef,self.routines)
                    self.routines.applyPatch('ASH', residue)
            elif group == 'GLU':
                if hydname == 'GLH':
                    amb = hydrogenAmbiguity(residue, hydrodef,self.routines)
                    self.routines.applyPatch('GLH', residue)
        if amb == None:
            text = "Could not find hydrogen ambiguity "
            text += "for titratable group %s!" % group
            raise ValueError, text
        return amb

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
                    #
                    # Skip comments
                    #
                    while line[:2]=='//':
                        line=file.readline()
                    
                    startstates = []
                    endstates = []
                    modelpKa = None

                    if line[:11] != 'Transition:':
                        text = "Wrong line found when looking for 'Transition:'"
                        raise ValueError, "%s: %s" % (text, line)
                    
                    split=string.split(line[11:],'->')
                    for number in string.split(split[0], ','):
                        startstates.append(string.strip(number))                        
                    for number in string.split(split[1], ','):
                        endstates.append(string.strip(number))

                    line = file.readline()
                    #
                    # Skip comments
                    #
                    while line[:2]=='//':
                        line=file.readline()
                    #
                    # Must be the model pKa line
                    #
                    if line[:10]!='Model_pKa:':
                        text = "Wrong line found when looking for 'Model_pKa'"
                        raise ValueError, "%s: %s" % (text, line)
                       
                    modelpKa = float(string.split(line)[1])

                    thisTitration = DefTitration(startstates, endstates,modelpKa,name)
                    titrations.append(thisTitration)
                    
                    line = file.readline()
                    if string.strip(line) == 'END': break

                thisGroup = pKaGroup(name, resname, type, titrations)
                mygroups[name] = thisGroup

                line = file.readline()
                if string.strip(line) == 'END OF FILE': break
        
        return mygroups


#
# -------------
#

def usage(x):
    #
    # Print usage guidelines
    #
    print 'Usage: pka.py --ff <forcefield> --lig <ligand in MOL2> --pdie <protein diel cons> --maps <1 for using 3D maps>'
    print '--xdiel <xdiel maps> --ydiel <ydiel maps> --zdiel <zdiel maps> --kappa <ion-accessibility map> <pdbfile>'
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
    longOptlist = ["help","verbose","ff=",'lig=',"pdie=","maps=","xdiel=","ydiel=","zdiel=","kappa=",]

    try:
        opts, args = getopt.getopt(sys.argv[1:], shortOptlist, longOptlist)
    except getopt.GetoptError, details:
        sys.stderr.write("GetoptError:  %s\n" % details)
        usage(2)
        sys.exit(0)
    #
    #
    #
    if len(args) < 1 or len(args) > 9:
        sys.stderr.write("Incorrect number (%d) of arguments!\n" % len(args))
        usage(2)
        sys.exit(0)

    verbose = 0
    ligfilename=None
    ff = None
    pdie = None
    maps=None
    xdiel=None
    ydiel=None
    zdiel=None
    kappa=None
    
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
        elif o == "--pdie":
            pdie = int(a)
        elif o == "--maps":
            maps = int(a)
        elif o == "--xdiel":
            xdiel = string.lower(a)
        elif o == "--ydiel":
            ydiel = string.lower(a)
        elif o == "--zdiel":
            zdiel = string.lower(a)
        elif o == "--kappa":
            kappa = string.lower(a)
        elif o == "--lig":
            ligfilename=a


    #
    # No forcefield? Set default forcefield to parse
    #
    if ff == None:
        ff = "parse"

    #
    # No dielectric constant
    #
    if pdie == None:
        pdie = 8
    #
    # Find the PDB file
    #
    path = args[0]
    #
    # Call the pre_init function
    #
    return pre_init(pdbfilename=path,ff=ff,ligand=ligfilename,pdie=pdie,maps=maps,xdiel=xdiel,ydiel=ydiel,zdiel=zdiel,kappa=kappa,)

#
# ----
#

def remove_hydrogens(pdb):

    fd = open(pdb,'r')
    l_lines_i = fd.readlines()
    fd.close()

    l_lines_o = []
    for s_line in l_lines_i:
        record = s_line[:6].strip()
        if record in ['HETATM']:
            continue
        if record in ['ATOM','ANISOU','SIGUIJ','SIGATM',]:
            element = s_line[76:78].strip()
            if element == 'H':
                continue
        l_lines_o += [s_line]

    fd = open(pdb,'w')
    fd.writelines(l_lines_o)
    fd.close()

    return


def pre_init(pdbfilename=None,ff=None,ligand=None,verbose=1,pdie=None,maps=None,xdiel=None,ydiel=None,zdiel=None,kappa=None,):
    #
    # remove hydrogen atoms
    #
    remove_hydrogens(pdbfilename)
   
    #
    # Get the PDBfile
    #
    global pdbfile_name
    pdbfile_name=pdbfilename
    pdbfile = getPDBFile(pdbfilename)
    pdblist, errlist = readPDB(pdbfile)
    
    if len(pdblist) == 0 and len(errlist) == 0:
        print "Unable to find file %s!\n" % path
        os.remove(path)
        sys.exit(2)

    if len(errlist) != 0 and verbose:
        print "Warning: %s is a non-standard PDB file.\n" %pdbfilename
        print errlist

    if verbose:
        print "Beginning PDB2PQR...\n"
    #
    # Read the definition file
    #
    myDefinition = Definition()
    ligand_titratable_groups=None
    #
    #
    # Choose whether to include the ligand or not
    #
    # Add the ligand to the pdb2pqr arrays
    #
    Lig=None
    MOL2FLAG = False 
    if not ligand:
        dummydef = Definition()
        myProtein = Protein(pdblist, dummydef)
    else:
        #
        # Read the ligand into Paul's code
        #
        Lig = ligandclean.ligff.ligand_charge_handler()
        Lig.read(ligand)
        #
        # Create the ligand definition from the mol2 data
        #
        import NEWligand_topology
        MOL2FLAG = True
        #
        X=NEWligand_topology.get_ligand_topology(Lig.lAtoms,MOL2FLAG)
        #
        # Add it to the 'official' definition
        #
        ligresidue=myDefinition.parseDefinition(X.lines, 'LIG', 2)
        myDefinition.AAdef.addResidue(ligresidue)
        #
        # Look for titratable groups in the ligand
        #
        print '==============================\n================================\n======================='
        ligand_titratable_groups=X.find_titratable_groups()
        print '==============================\n================================\n======================='
        print "ligand_titratable_groups", ligand_titratable_groups
        #
        # ------------------------------------------------------
        # Creation of ligand definition and identification of ligand titgrps done
        # Start loading everything into PDB2PQR
        #
        # Append the ligand data to the end of the PDB data
        #
        newpdblist=[]
        # First the protein
        for line in pdblist:
            if isinstance(line, END) or isinstance(line,MASTER):
                continue
            newpdblist.append(line)
        # Now the ligand
        for e in Lig.lAtoms:
            newpdblist.append(e)
        #
        # Add a TER and an END record for good measure
        #
        newpdblist.append(TER)
        newpdblist.append(END)
        #
        # Let PDB2PQR parse the entire file
        #
        myProtein = Protein(newpdblist)

        #
        # Post-Processing for adding sybylTypes to lig-atoms in myProtein
        # Jens: that's the quick and easy solution
        for rrres in  myProtein.chainmap['L'].residues:
            for aaat in rrres.atoms:
                for ligatoms in Lig.lAtoms:
                    if ligatoms.name == aaat.name:
                        aaat.sybylType = ligatoms.sybylType
                        #
                        # setting the formal charges
                        if ligatoms.sybylType == "O.co2":
                            aaat.formalcharge = -0.5
                        else: aaat.formalcharge = 0.0
                        xxxlll = []
                        for xxx in ligatoms.lBondedAtoms:
                            xxxlll.append(xxx.name)
                        aaat.intrabonds = xxxlll
                        #
                        # charge initialisation must happen somewhere else
                        # but i don't know where...
                        aaat.charge = 0.0
        #
    #
    # =======================================================================
    #
    # We have identified the structural elements, now contiue with the setup
    #
    # Print something for some reason?
    #
    if verbose:
        print "Created protein object -"
        print "\tNumber of residues in protein: %s" % myProtein.numResidues()
        print "\tNumber of atoms in protein   : %s" % myProtein.numAtoms()
    #
    # Set up all other routines
    #
    myRoutines = Routines(myProtein, verbose) #myDefinition)
    myRoutines.updateResidueTypes()
    myRoutines.updateSSbridges()
    myRoutines.updateBonds()
    myRoutines.setTermini()
    myRoutines.updateInternalBonds()

    myRoutines.applyNameScheme(Forcefield(ff, myDefinition, None))
    myRoutines.findMissingHeavy()
    myRoutines.addHydrogens()

    #myRoutines.randomizeWaters()
    myProtein.reSerialize()
    #
    # Inject the information on hydrogen conformations in the HYDROGENS.DAT arrays
    # We get this information from ligand_titratable_groups
    #
    from src.hydrogens import hydrogenRoutines
    myRoutines.updateInternalBonds()
    myRoutines.calculateDihedralAngles()
    myhydRoutines = hydrogenRoutines(myRoutines)
    #
    # Here we should inject the info!!
    #
    myhydRoutines.initializeFullOptimization()
    myhydRoutines.optimizeHydrogens()
    #
    # Choose the correct forcefield
    #
    if not Lig:
        myForcefield = Forcefield(ff, myDefinition, None)
    else:
        # Here I am now passing the instance of Lig that's already been create - I have not
        # yet changed it in ligforcefield.... It will give an error for now
        #
        myForcefield = ligandclean.ligff.ligforcefield(ff,Lig)
        print "That's Lig: ", Lig
        print "That's X:   " , X
    
    if verbose:
        print "Created protein object (after processing myRoutines) -"
        print "\tNumber of residues in protein: %s" % myProtein.numResidues()
        print "\tNumber of atoms in protein   : %s" % myProtein.numAtoms()
    #
    # Create the APBS input file
    #
    import src.psize
    size=src.psize.Psize()

    method=""
    async=0
    split=0
    import inputgen_pKa
    igen = inputgen_pKa.inputGen(pdbfilename)
    #
    # For convenience
    #
    igen.pdie = pdie
    print 'Setting protein dielectric constant to ',igen.pdie
    igen.sdie=80
    igen.maps=maps
    if maps==1:
        print "Using dielectric and mobile ion-accessibility function maps in PBE"
        if xdiel:
            igen.xdiel = xdiel
        else:
            sys.stderr.write ("X dielectric map is missing\n")
            usage(2)
            sys.exit(0)
        if ydiel:
            igen.ydiel = ydiel
        else:
            sys.stderr.write ("Y dielectric map is missing\n")
            usage(2)
            sys.exit(0)
        if zdiel:
            igen.zdiel = zdiel
        else:
            sys.stderr.write ("Z dielectric map is missing\n")
            usage(2)
            sys.exit(0)
        
        print 'Setting dielectric function maps: %s, %s, %s'%(igen.xdiel,igen.ydiel,igen.zdiel)
        
        if kappa:
            igen.kappa = kappa
        else:
            sys.stderr.write ("Mobile ion-accessibility map is missing\n")
            usage(2)
            sys.exit(0)
            
        print 'Setting mobile ion-accessibility function map to: ',igen.kappa
    #
    # Return all we need
    #
    return myProtein, myRoutines, myForcefield,igen, ligand_titratable_groups

#
# -----------------------------------------------
#

def test_interface():
    """Test the interface with pKaTool"""
    import pKaTool.pKa_calc
    X=pKaTool.pKa_calc.Monte_Carlo_Mult_CPP()
    
    X.intrinsic_pKa={':0001:ASP':[0.0,4.0,5.0]}
    X.charged_state={':0001:ASP':[0,1,1]}
    X.acid_base={':0001:ASP':-1}
    X.intene_mult={':0001:ASP':{':0001:ASP':[[0,0,0],[0,0,0],[0,0,0]]}}
    X._calc_pKas(0.0,10.0,0.5)
    return

def titrate_one_group(name,intpkas,is_charged,acidbase):
    """Titrate a single group and return the pKa value for it"""
    names=[name]
    num_states=len(intpkas)
    state_counter=[num_states]
    linear=[] # The linear matrix
    for group2 in range(1):
        for state1 in range(num_states):
            for state2 in range(num_states):
                linear.append(0.0)
    #
    # Set the MC parameters
    #
    mcsteps=5000
    phstart=0.1
    phend=20.0
    phstep=0.1
    #
    # Call our little C++ module
    #
    import pMC_mult
    FAST=pMC_mult.MC(intpkas,linear,acidbase,state_counter,is_charged)
    FAST.set_MCsteps(int(mcsteps))
    print 'Calculating intrinsic pKa value'
    pKavals=FAST.calc_pKas(phstart,phend,phstep)
    count=0
    intpka=pKavals[0]
    print 'Simulated intrinsic pKa value: %5.2f' %intpka
    count=1
    #
    # Get the charges
    #
    charges={}
    pH_start=pKavals[count]
    pH_step=pKavals[count+1]
    num_pHs=pKavals[count+2]
    count=count+2
    pHs=[]
    charges=[]
    pH=pH_start
    for x in range(int(num_pHs)):
        count=count+1
        pHs.append(pKavals[count])
        count=count+1
        charges.append(pKavals[count])
        pH=pH+pH_step
    if pKavals[count+1]==999.0 and pKavals[count+2]==-999.0:
        count=count+2
    else:
        print 'Something is wrong'
        print pKavals[count:count+30]
        raise Exception('Incorrect data format from pMC_mult')
    return intpka
            
if __name__ == "__main__":
    state=1
    if state==0:
        import pMC_mult
        #
        # System definition
        #
        groups=2
        acidbase=[-1,1] # 1 if group is a base, -1 if group is an acid
        intpkas=[3.4,0.0,0.0,0.0,0.0,
                 9.6,0.0,0.0,0.0,0.0] 
        is_charged_state=[1,0,0,0,0,
                          1,0,0,0,0]
        #
        # Automatic configuration from here on
        #
        states=len(intpkas)/groups #States per group
        state_counter=[]
        linear=[]
        names=[]
        for group in range(groups):
            #
            # Names
            #
            names.append('Group %d' %group)
            #
            # Add state_counter
            #
            state_counter.append(states)
            #
            # Matrix
            #
            for group2 in range(groups):
                for state1 in range(states):
                    for state2 in range(states):
                        if state1==0 and state2==0 and group!=group2:
                            linear.append(-2.3)
                        else:
                            linear.append(0.0)
        mcsteps=50000
        phstart=2.0
        phend=14.0
        phstep=0.1
        #
        # Call our little C++ module
        #
        FAST=pMC_mult.MC(intpkas,linear,acidbase,state_counter,is_charged_state)
        FAST.set_MCsteps(int(mcsteps))
        print 'Starting to calculate pKa values'
        pKavals=FAST.calc_pKas(phstart,phend,phstep)
        count=0
        print '\n****************************'
        print 'Final pKa values'
        pkas={}
        for name in names:
            pkas[name]={'pKa':pKavals[count]}
            print '%s pKa: %5.2f ' %(name,pKavals[count])
            count=count+1
        #
        # Write the WHAT IF pKa file
        #
        for name in names:
            pkas[name]['modelpK']=0.0
            pkas[name]['desolv']=0.0
            pkas[name]['backgr']=0.0
            pkas[name]['delec']=0.0
        import pKaTool.pKaIO
        X=pKaTool.pKaIO.pKaIO()
        X.write_pka('test.pdb.PKA.DAT',pkas)
        #
        # Get the charges
        #
        charges={}
        pH_start=pKavals[count]
        pH_step=pKavals[count+1]
        num_pHs=pKavals[count+2]
        count=count+2
        for name in names:
            charges[name]={}
            for x in range(int(num_pHs)):
                count=count+1
                pH=pKavals[count]
                count=count+1
                charges[name][pH]=pKavals[count]
                pH=pH+pH_step
            if pKavals[count+1]==999.0 and pKavals[count+2]==-999.0:
                count=count+2
            else:
                print 'Something is wrong'
                print pKavals[count:count+30]
                raise Exception('Incorrect data format from pMC_mult')
    elif state==1:
        #
        # Do a real pKa calculation
        #
        protein, routines, forcefield,apbs_setup, ligand_titratable_groups = startpKa()
        mypkaRoutines = pKaRoutines(protein, routines, forcefield,apbs_setup)
        #
        # Debugging
        #
        if debug:
            CM.init_protein(mypkaRoutines)
        #
        # Deal with ligands
        #
        if ligand_titratable_groups:
            print "lig (before mypKaRoutines) ", ligand_titratable_groups['titratableatoms']
            mypkaRoutines.insert_new_titratable_group(ligand_titratable_groups)
        mypkaRoutines.runpKa()
    elif state==2:
        #
        # Just assign charges
        #
        protein, routines, forcefield,apbs_setup, ligand_titratable_groups = startpKa()
        for chain in protein.getChains():
            for residue in chain.get("residues"):
                for atom in residue.get("atoms"):
                    atomname = atom.get("name")
                    charge, radius = forcefield.getParams(residue, atomname)
                    print '%2s %4s %3d %4s %5.2f %5.2f' %(chain.chainID,residue.name,residue.resSeq,atomname,charge,radius)
