"""
    Hydrogen optimization routines for PDB2PQR

    This module contains the hydrogen optimization routines and classes for
    PDB2PQR. It uses a deterministic algorithm to perform flips, rotate
    alcoholic hydrogens, and optimize waters.

    Based on C code from Jens Erik Nielsen
    UCSD/HHMI
    
    Ported to Python by Todd Dolinsky (todd@ccb.wustl.edu)
    Washington University in St. Louis
"""

__date__ = "11 May 2005"
__author__ = "Jens Erik Nielsen, Todd Dolinsky"

HYDROGENFILE = "dat/HYDROGENS.DAT"

import os
import string

from definitions import *
from utilities import *
from math import *
from quatfit import *
from random import *
from time import *

def sortdict(dict):
    """
        Sort a dictionary by its values

        Parameters
            dict:  The dictionary to sort (dict)
        Returns
            items:  A list of the items in tuples (list)
    """
    items = [(v, k) for k, v in dict.items()]
    items.sort()
    items.reverse()             
    items = [(k, v) for v, k in items]
    return items

class hbond:
    """
        A small class containing the hbond structure
    """
    def __init__(self, atom1, atom2, dist):
        """
            Initialize the class

            Parameters
                atom1:  The first atom in the potential bond (Atom)
                atom2:  The second atom in the potential bond (Atom)
                dist:  The distance between the two atoms (float)
        """
        self.atom1 = atom1
        self.atom2 = atom2
        self.dist = dist
      
class hydrogenAmbiguity:
    """
        A class containing information about the ambiguity
    """
    def __init__(self, residue, hdef, routines):
        """
            Initialize the class - if the residue has a rotateable hydrogen,
            remove it.  If it can be flipped, pre-flip the residue by creating
            all additional atoms.

            Parameters
                residue:  The residue in question (residue)
                hdef:     The hydrogen definition matching the residue
                routines: Pointer to the general routines object
        """
        self.residue = residue
        self.hdef = hdef
        self.nearatoms = []
        self.routines = routines
        self.fixed = 0
        self.proth = None

        if residue.name in ["HIS","HS2N","HSE","HSD","HSP","HIE","HID","HIP"]:
            nd1atom = residue.getAtom("ND1").getCoords()
            cd2atom = residue.getAtom("CD2").getCoords()
            ce1atom = residue.getAtom("CE1").getCoords()
            ne2atom = residue.getAtom("NE2").getCoords()
            hd2atom = residue.getAtom("HD2").getCoords()
            he1atom = residue.getAtom("HE1").getCoords()
            if "HE2" in residue.map:
                he2atom = residue.getAtom("HE2").getCoords()
            if "HD1" in residue.map:
                hd1atom = residue.getAtom("HD1").getCoords()
            self.makeFlip()
            residue.createAtom("ND1FLIP",nd1atom,"ATOM")
            residue.createAtom("CD2FLIP",cd2atom,"ATOM")
            residue.createAtom("CE1FLIP",ce1atom,"ATOM")
            residue.createAtom("NE2FLIP",ne2atom,"ATOM")
            residue.createAtom("HD2FLIP",hd2atom,"ATOM")
            residue.createAtom("HE1FLIP",he1atom,"ATOM")
            if "HE2" in residue.map:
                residue.createAtom("HE2FLIP",he2atom,"ATOM")
                residue.getAtom("NE2FLIP").hdonor = 1
                residue.getAtom("NE2FLIP").intrabonds = ["HE2FLIP", "CE1FLIP", "CD2FLIP"]
                residue.getAtom("HE2FLIP").intrabonds = ["NE2FLIP"]
            else:
                residue.getAtom("NE2FLIP").hacceptor = 1
                residue.getAtom("NE2FLIP").intrabonds = ["CE1FLIP", "CD2FLIP"]
            if "HD1" in residue.map:
                residue.createAtom("HD1FLIP",hd1atom,"ATOM")
                residue.getAtom("ND1FLIP").hdonor = 1
                residue.getAtom("ND1FLIP").intrabonds = ["HD1FLIP", "CE1FLIP", "CG"]
                residue.getAtom("HD1FLIP").intrabonds = ["ND1FLIP"]
            else:
                residue.getAtom("ND1FLIP").hacceptor = 1
                residue.getAtom("ND1FLIP").intrabonds = ["CE1FLIP", "CG"]

            # Set intrabonds
            residue.getAtom("CD2FLIP").intrabonds = ["CG","HD2FLIP", "NE2FLIP"]
            residue.getAtom("CE1FLIP").intrabonds = ["HE1FLIP","ND1FLIP","NE2FLIP"]
            residue.getAtom("HD2FLIP").intrabonds = ["CD2FLIP"]
            residue.getAtom("HE1FLIP").intrabonds = ["CE1FLIP"]

            if residue.name == "HS2N": # Set to acceptor and donor
                residue.getAtom("ND1FLIP").hacceptor = 1
                residue.getAtom("ND1FLIP").hdonor = 1
                residue.getAtom("NE2FLIP").hdonor = 1
                residue.getAtom("NE2FLIP").hacceptor = 1

        elif residue.name == "GLN":
            oatom = residue.getAtom("OE1").getCoords()
            natom = residue.getAtom("NE2").getCoords()
            hatom1 = residue.getAtom("HE21").getCoords()
            hatom2 = residue.getAtom("HE22").getCoords()
            self.makeFlip()
            residue.createAtom("OE1FLIP", oatom, "ATOM")
            residue.createAtom("NE2FLIP", natom, "ATOM")
            residue.createAtom("HE21FLIP", hatom1, "ATOM")
            residue.createAtom("HE22FLIP", hatom2, "ATOM")
            residue.getAtom("OE1FLIP").hacceptor = 1
            natom = residue.getAtom("NE2FLIP")
            natom.hdonor = 1
            natom.intrabonds = ["HE21FLIP", "HE22FLIP"]
            residue.getAtom("HE21FLIP").intrabonds = ["NE2FLIP"]
            residue.getAtom("HE22FLIP").intrabonds = ["NE2FLIP"]
        elif residue.name == "ASN":
            oatom = residue.getAtom("OD1").getCoords()
            natom = residue.getAtom("ND2").getCoords()
            hatom1 = residue.getAtom("HD21").getCoords()
            hatom2 = residue.getAtom("HD22").getCoords()
            self.makeFlip()
            residue.createAtom("OD1FLIP", oatom, "ATOM")
            residue.createAtom("ND2FLIP", natom, "ATOM")
            residue.createAtom("HD21FLIP", hatom1, "ATOM")
            residue.createAtom("HD22FLIP", hatom2, "ATOM")
            residue.getAtom("OD1FLIP").hacceptor = 1
            natom = residue.getAtom("ND2FLIP")
            natom.hdonor = 1
            natom.intrabonds = ["HD21FLIP", "HD22FLIP"]
            residue.getAtom("HD21FLIP").intrabonds = ["ND2FLIP"]
            residue.getAtom("HD22FLIP").intrabonds = ["ND2FLIP"]

        elif hdef.type == 2:  # Remove the hydrogen atom - it will be readded
            if residue.name == "SER":
                residue.removeAtom("HG")
                bonds = residue.getAtom("OG").intrabonds
                bonds.pop(bonds.index("HG"))
            elif residue.name == "THR":
                residue.removeAtom("HG1")
                bonds = residue.getAtom("OG1").intrabonds
                bonds.pop(bonds.index("HG1"))
            elif residue.name == "TYR":
                residue.removeAtom("HH")
                bonds = residue.getAtom("OH").intrabonds
                bonds.pop(bonds.index("HH"))

        elif hdef.type == 13:  # Remove the current hydrogen atom and pre-add
            if residue.name == "GLH":
                bond1 = "OE1"
                bond2 = "OE2"
                ref = "CD"
                h1 = "HE1"
                h2 = "HE2"
            elif residue.name == "ASH":
                bond1 = "OD1"
                bond2 = "OD2"
                ref = "CG"
                h1 = "HD1"
                h2 = "HD2"
            elif residue.isCterm:
                bond1 = "O"
                bond2 = "OXT"
                ref = "C"
                h1 = "HO1"
                h2 = "HO2"
            else:
                txt = "Unsupported type 13 residue for optimization!"
                raise ValueError, txt
                
            if h1 in residue.map:
                residue.removeAtom(h1)
                residue.getAtom(bond1).intrabonds.pop(residue.getAtom(bond1).intrabonds.index(h1))
            if h2 in residue.map:
                residue.removeAtom(h2)
                residue.getAtom(bond2).intrabonds.pop(residue.getAtom(bond2).intrabonds.index(h2))
            if "HO" in residue.map and residue.isCterm:
                residue.removeAtom("HO")
               
            # Set up the coordinate frames
            refcoords = []
            defcoords = []
            refcoords.append(residue.getAtom(bond1).getCoords())
            refcoords.append(residue.getAtom(bond2).getCoords())
            refcoords.append(residue.getAtom(ref).getCoords())
            defcoords = [[14.109, 14.303, 18.212], [12.267, 14.963, 19.265], \
                         [13.140, 14.094, 18.958]]
            
            # Add the atoms
            defatomcoords = [14.246,  15.203,  17.799]
            newcoords = findCoordinates(3, refcoords, defcoords, defatomcoords)
            hname = "%s1" % h1
            residue.createAtom(hname, newcoords, "ATOM")
            residue.getAtom(hname).intrabonds = [bond1]
            residue.getAtom(bond1).intrabonds.append(hname)

            defatomcoords = [14.797, 13.618, 17.970]
            newcoords = findCoordinates(3, refcoords, defcoords, defatomcoords)
            hname = "%s2" % h1
            residue.createAtom(hname, newcoords, "ATOM")
            residue.getAtom(hname).intrabonds = [bond1]
            residue.getAtom(bond1).intrabonds.append(hname)
            
            defatomcoords = [12.404, 15.863, 18.852]
            newcoords = findCoordinates(3, refcoords, defcoords, defatomcoords)
            hname = "%s1" % h2
            residue.createAtom(hname, newcoords, "ATOM")
            residue.getAtom(hname).intrabonds = [bond2]
            residue.getAtom(bond2).intrabonds.append(hname)

            defatomcoords = [11.486, 14.795, 19.866]
            newcoords = findCoordinates(3, refcoords, defcoords, defatomcoords)
            hname = "%s2" % h2
            residue.createAtom(hname, newcoords, "ATOM")
            residue.getAtom(hname).intrabonds = [bond2]
            residue.getAtom(bond2).intrabonds.append(hname)

            residue.getAtom(bond1).hdonor = 1
            residue.getAtom(bond1).hacceptor = 1
            residue.getAtom(bond2).hdonor = 1
            residue.getAtom(bond2).hacceptor = 1

    def __str__(self):
        """
            Print the ambiguity for debugging purposes
        """
        text = "%s %i %s (%i)" % (self.residue.name, self.residue.resSeq, \
                                  self.residue.chainID, self.hdef.type)
        return text
          
    def makeFlip(self):
        """
            Flip the residue (rotate around the chiangle by 180)
        """
        residue = self.residue
        name = residue.get("name")
        defresidue = self.routines.aadef.getResidue(name)
        chinum = self.hdef.chiangle - 1
        oldangle = residue.get("chiangles")[chinum]
        angle = 180.0 + oldangle
        self.routines.setChiangle(residue, chinum, angle, defresidue)

    def fixFlip(self, boundatom, donorflag=None):
        """
            Fix the flippable residue.  The residue will no loner be
            rotated.  This deletes the worst-state atoms, but does not
            rename them.

            Parameters:
                boundatom:  The atom that made the best bond (Atom)
                donorflag:  Optional flag used for HS2N to determine whether
                            the boundatom is donor/acceptor (int)
        """
        if self.fixed: return
        residue = self.residue
        self.fixed = 1
        atomlist = []
        intrabondlist = []
        flag = 0
        if boundatom == None or boundatom.name.endswith("FLIP"): flag = 1
        for atom in residue.atoms: atomlist.append(atom)
        
        # Handle the special case for HS2N!
       
        if residue.name == "HS2N":
            if donorflag == None:
                txt = "Donor flag must be set for HS2N!"
                raise ValueError, txt

            if flag: flip = "FLIP"
            else: flip = ""
            residue.renameResidue("HIS")
            if boundatom == None: # Set to HSD
                residue.removeAtom("HE2%s" % flip)
                ne2atom = residue.getAtom("NE2%s" % flip)
                ne2atom.intrabonds.pop(ne2atom.intrabonds.index("HE2%s" % flip))
                ne2atom.hdonor = 0
                residue.getAtom("ND1%s" % flip).hacceptor = 0   
            elif (boundatom.name.startswith("ND1") and donorflag) or \
                 (boundatom.name.startswith("NE2") and not donorflag): # Set to HSD
                residue.removeAtom("HE2%s" % flip)
                ne2atom = residue.getAtom("NE2%s" % flip)
                ne2atom.intrabonds.pop(ne2atom.intrabonds.index("HE2%s" % flip))
                ne2atom.hdonor = 0
                residue.getAtom("ND1%s" % flip).hacceptor = 0
            elif (boundatom.name.startswith("NE2") and donorflag) or \
                 (boundatom.name.startswith("ND1") and not donorflag): # Set to HSE
                residue.removeAtom("HD1%s" % flip)
                nd1atom = residue.getAtom("ND1%s" % flip)
                nd1atom.intrabonds.pop(nd1atom.intrabonds.index("HD1%s" % flip))
                nd1atom.hdonor = 0
                residue.getAtom("NE2%s" % flip).hacceptor = 0
            else:
                text = "Invalid bound atom %s for HS2N!" % boundatom.name
                raise ValueError, text

        for atom in atomlist:
            atomname = atom.name
            if atomname.endswith("FLIP") and flag: # Delete the other list
                residue.removeAtom(atomname[:-4])
                intrabondlist.append(atomname[:-4])
            elif atomname.endswith("FLIP"):  # Delete the flip
                residue.removeAtom(atomname)
                intrabondlist.append(atomname)
            else: continue

        # Since the atoms are gone, remove them from the intrabond list

        for atom in atomlist:
            if atom not in residue.atoms: continue
            for intrabond in atom.intrabonds:
                if intrabond in intrabondlist:
                    atom.intrabonds.pop(atom.intrabonds.index(intrabond))

    def fixProt(self, hydatom):
        """
            Fix a type 13 residue by keeping the hydatom given.
            All others can be deleted
        """
        if self.fixed: return
        self.fixed = 1
        residue = self.residue

       
        if residue.name == "GLH":
            bond1 = "OE1"
            bond2 = "OE2"
            ref = "CD"
            h1 = "HE1"
            h2 = "HE2"
        elif residue.name == "ASH":
            bond1 = "OD1"
            bond2 = "OD2"
            ref = "CG"
            h1 = "HD1"
            h2 = "HD2"
        elif residue.isCterm:
            bond1 = "O"
            bond2 = "OXT"
            ref = "C"
            h1 = "HO1"
            h2 = "HO2"

        if hydatom == None: # Choose h1-1
            hydatom = residue.getAtom("%s1" % h1)

        # Remove all other hydrogens
        for i in range(1,3):
            hname = "%s%i" % (h1, i)
            if hname != hydatom.name:
                residue.removeAtom(hname)
                residue.getAtom(bond1).intrabonds.pop(residue.getAtom(bond1).intrabonds.index(hname))
            hname = "%s%i" % (h2, i)
            if hname != hydatom.name:
                residue.removeAtom(hname)
                residue.getAtom(bond2).intrabonds.pop(residue.getAtom(bond2).intrabonds.index(hname))

        self.proth = hydatom.name
  
    def setNearatoms(self, allatoms, ambmap):
        """
            Set the nearby atoms to this residue.  The only donors/acceptors
            that will be changing positions are the flips.

            Parameters
                allatoms:  A list of all donors/acceptors (list)
                ambmap:  A mapping of all residues with ambiguities
        """
        outlist = []
        nearatoms = []
        resname = self.residue.get("name")
        confs = []
        if self.hdef.type == 11:
            for atom in self.residue.atoms:
                atomname = atom.name
                if atomname.endswith("FLIP") and not atomname[0] == "H": confs.append(atomname)
        for conf in self.hdef.conformations:
            atomname = conf.boundatom
            if atomname not in confs: confs.append(atomname)
        map = []
        for atomname in confs:
            atom = self.residue.getAtom(atomname)
            for nearatom in allatoms:
                nearres = nearatom.get("residue")
                nearname = nearres.get("name")
                if nearres == self.residue: continue
                if atom.hdonor and not atom.hacceptor and not nearatom.hacceptor: continue
                elif atom.hacceptor and not atom.hdonor and not nearatom.hdonor: continue
                elif not (atom.hdonor or atom.hacceptor): continue
                dist = distance(atom.getCoords(), nearatom.getCoords())
                item = (atom, nearatom, dist)
                if dist < 4.3 and item not in map:
                    myhbond = hbond(atom, nearatom, dist)
                    nearatoms.append(myhbond)
                    map.append(item)
                    if nearatom.residue in ambmap and (nearatom.residue.name == "WAT" or nearatom.name not in ["O","N"]):
                        outlist.append(ambmap[nearatom.residue])
                                                                                                              
        self.nearatoms = nearatoms
        return outlist
        
class hydrogenRoutines:
    """
        The main class of routines in the hydrogen optimization process.
    """
    
    def __init__(self, routines):
        """
            Initialize the routines and run the hydrogen optimization

            Parameters
                routines: The parent routines object (Routines)             
        """
        self.hdebug = 1
        self.routines = routines
        self.protein = routines.protein
        self.hydrodefs = []
        self.groups = []
   
    def debug(self, text):
        """
            Print text to stdout for debugging purposes.

            Parameters
                text:  The text to output (string)
        """
        if self.hdebug:
            print text

    def getstates(self, amb):
        """
            Get all possible states for a conformation/protonation
            ambiguity and store them in a list. Each.
            hydrogen type must be explicitly defined.

            Parameters
                amb   : The ambiguity to get the states of (tuple)
            Returns
                states: A list of states, where each state
                        is a list of conformations of the atom. (list)
        """
        states = []
        residue = getattr(amb,"residue")
        hdef = getattr(amb,"hdef")
        type = hdef.type

        confs = hdef.conformations

        # Now make the states

        if type == 1: # CTR/ASP/GLU
            states.append([()])
            states.append([confs[0]]) # H on O1 cis
            states.append([confs[1]]) # H on O1 trans
            states.append([confs[2]]) # H on O2 cis
            states.append([confs[3]]) # H on O2 trans
        elif type == 3: # NTR
            states.append([()]) # No Hs
            states.append([confs[0]]) # H3 only
            states.append([confs[1]]) # H2 only
            states.append([confs[0], confs[1]]) # H3 and H2
        elif type == 4: # HIS 
            states.append([confs[0]]) # HSD
            states.append([confs[1]]) # HSE
            states.append([()]) # Negative HIS
            states.append([confs[0], confs[1]]) #HSP
        elif type == 10: #PNTR
            states.append([()])
            states.append([confs[0]])
            states.append([confs[1]])
            states.append([confs[0], confs[1]])
        elif type == 13: #ASH/GLH/CTR
            states.append([confs[0]]) # H on O1 cis
            states.append([confs[1]]) # H on O1 trans
            states.append([confs[2]]) # H on O2 cis
            states.append([confs[3]]) # H on O2 trans
        elif type == 14: #HSD/HSE ONLY!
            states.append([confs[0]]) # HSD
            states.append([confs[1]]) # HSE
        return states          
                
    def switchstate(self, states, amb, id):
        """
            Switch a residue to a new state by first removing all
            hydrogens.

            Parameters
                states: The list of states (list)
                amb   : The amibiguity to switch (tuple)
                id    : The state id to switch to (int)
        """
        if id > len(states):
            raise ValueError, "Invalid State ID!"
        
        # First Remove all Hs
        residue = getattr(amb,"residue")
        hdef = getattr(amb,"hdef")
        type = hdef.type
        for conf in hdef.conformations:
            hname = conf.hname
            boundname = conf.boundatom
            if residue.getAtom(hname) != None:
                residue.removeAtom(hname)
            residue.getAtom(boundname).hacceptor = 1
            residue.getAtom(boundname).hdonor = 0

        # Update the IntraBonds
        name = residue.get("name")
        defresidue = self.routines.aadef.getResidue(name)
        residue.updateIntraBonds(defresidue)
                              
        # Now build appropriate atoms
        state = states[id]
        for conf in state:
            refcoords = []
            defcoords = []
            defatomcoords = []
            if conf == (): continue # Nothing to add
            hname = conf.hname
            for atom in conf.atoms:
                atomname = atom.get("name")
                resatom = residue.getAtom(atomname)
                if atomname == hname:
                    defatomcoords = atom.getCoords()
                elif resatom != None:
                    refcoords.append(resatom.getCoords())
                    defcoords.append(atom.getCoords())
                else:
                    raise ValueError, "Could not find necessary atom!"
        
            newcoords = findCoordinates(3, refcoords, defcoords, defatomcoords)
            boundname = conf.boundatom
            residue.createAtom(hname, newcoords, "ATOM")
            residue.addDebumpAtom(residue.getAtom(hname))
            residue.getAtom(boundname).addIntraBond(hname)    
            residue.getAtom(boundname).hacceptor = 0
            residue.getAtom(boundname).hdonor = 1

 
    def printNetwork(self, network):
        """
            Print the network of ambiguities
        """
        text = ""
        for amb in network:
            text += "%s, " % amb
        text = text[:-2]
        return text

    def optimizeHydrogens(self, pkaflag):
        """
            Optimize hydrogens according to HYDROGENS.DAT.  This
            function serves as the main driver for the optimizing
            script.

            Parameters:
                pkaflag: 1 if pka calculations have already been done.  This
                         will ignore all protonation ambiguities (but not
                         necessarily the placement of hydrogens).  Otherwise
                         set to 0.
        """
        starttime = time()
        allatoms = self.findAmbiguities(0, pkaflag)
        self.printAmbiguities()

        ambmap = {}

        # Make the network map
        for amb in self.groups:   
            residue = getattr(amb, "residue")
            ambmap[residue] = amb

        nearmap = {}
        for amb in self.groups:
            nearlist = amb.setNearatoms(allatoms, ambmap)
            nearmap[amb] = nearlist

        done = []
        networks = []
        for amb in self.groups:
            if amb in nearmap and amb not in done:
                list = analyzeMap(nearmap,amb,[])
                for item in list:
                    done.append(item)
                networks.append(list)
            elif amb not in nearmap and amb not in done:
                networks.append([i])
                
        for network in networks:
            self.debug(self.printNetwork(network))

        # Start the optimization        
        for cluster in networks:
            self.debug("*** STARTING NETWORK %s ***\n" % self.printNetwork(cluster))
            # STEP 1 - Set up the clustermap
            clustermap = {}
            for amb in cluster:
                residue = getattr(amb, "residue")
                clustermap[residue] = amb
                
            # STEP 2 - Order the distance of the hbonds
            # Do potential bonds to backbone first

            seen = []
            backbonemap = {}
            for amb in cluster:
                residue = getattr(amb,"residue")
                nearatoms = getattr(amb,"nearatoms")
                for hbond in nearatoms:
                    if not hbond.atom2.residue in clustermap or (hbond.atom2.residue in clustermap and hbond.atom2.residue.name != "WAT" and hbond.atom2.name in ["O","N"]): # Only non-dual ambs first
                        backbonemap[hbond] = hbond.dist
                        seen.append(hbond)
                          
            backbonelist = sortdict(backbonemap)
            backbonelist.reverse()
              
            for item in backbonelist:
                hbond = item[0]
                atom1 = hbond.atom1
                atom2 = hbond.atom2
                dist = item[1]
                
                if atom1 == None:
                    self.debug("already fixed")
                    continue

                residue = atom1.residue
                amb = clustermap[residue]
                type = amb.hdef.type
                
                if amb.fixed: continue
                
                self.debug("Working on %s %i %s %s %s %i %s %s %.2f" % (atom1.residue.name, atom1.residue.resSeq, atom1.residue.chainID, atom1.name, atom2.residue.name, atom2.residue.resSeq, atom2.residue.chainID, atom2.name, dist))
                if atom2.hdonor and atom2.hacceptor: # We have to try both!
                    self.debug("IN OPTION 0")

                    # First try atom2 as donor:
                    result = self.tryAcceptor(atom1, atom2, amb)
                    if result: continue
                    
                    # Next try atom2 as acceptor
                    self.tryDonor(atom1, atom2, amb)

                elif atom2.hacceptor: # Treat Atom 1 as donor
                    self.debug("IN OPTION 1!")
                    self.tryDonor(atom1, atom2, amb)
                elif atom2.hdonor:  # Treat Atom 1 as acceptor
                    self.debug("IN OPTION 2!")
                    self.tryAcceptor(atom1, atom2, amb)
                else: 
                    self.debug("Resolved in flipping")
                    

            # STEP 3: All rotateable to other rotateable bonds
            self.debug("\nNOW MOVING ON TO ALL ROTATE to ROTATE bonds")
            distmap = {}
            for amb in cluster:
                if amb.fixed: continue
                residue = getattr(amb,"residue")
                nearatoms = getattr(amb,"nearatoms")
                for hbond in nearatoms:
                    if hbond.atom2.residue in clustermap and (clustermap[hbond.atom2.residue].hdef.type != 12 or amb.hdef.type != 12) \
                           and hbond not in seen:
                        distmap[hbond] = hbond.dist
                            
            distlist = sortdict(distmap)
            distlist.reverse()
          
            for item in distlist:
                hbond = item[0]
                atom1 = hbond.atom1
                atom2 = hbond.atom2
                dist = item[1]
                residue = atom1.residue
                amb = clustermap[residue]
                type = amb.hdef.type
                amb2 = clustermap[atom2.residue]
                type2 = amb2.hdef.type
                
                self.debug("Working on %s %i %s %s %s %i %s %s %.2f" % (atom1.residue.name, atom1.residue.resSeq, atom1.residue.chainID, atom1.name, atom2.residue.name, atom2.residue.resSeq, atom2.residue.chainID, atom2.name, dist))

                if amb.fixed: continue

                if atom2 not in atom2.residue.atoms:
                    self.debug("Was already fixed!")
                    continue
                
                # We would rather do WAT-X than X-WAT and X-11 than 11-X
                # We would also rather do X-13 than 13-X
                if type2 == 12 or (type in [11,14] and type2 not in [11,14]):
                    self.debug("Ignoring this one - will get the other")
                    continue
                if type == 13 and type2 != 13:
                    self.debug("Ignoring this one - will get the other")
                    continue
                
                if atom2.hdonor and atom2.hacceptor: # We have to try both!
                    self.debug("IN OPTION 0")
                    if atom1.hdonor: # Try atom1 as donor first
                        setflag = self.tryDonor(atom1, atom2, amb)
                        if setflag:
                            accflag = self.tryAcceptor(atom2, atom1,  amb2)
                            if not accflag:
                                self.unAdd(amb, atom1)
                            else:
                                continue
                                
                    setflag = self.tryAcceptor(atom1, atom2, amb)
                    if setflag:
                        accflag = self.tryDonor(atom2, atom1, amb2)
                        if not accflag:
                            self.unAdd(amb, atom1)                       
               
                elif atom2.hacceptor: # Treat Atom 1 as donor
                    self.debug("IN OPTION 1")
                    setflag = self.tryDonor(atom1, atom2, amb)
                    if setflag:
                        accflag = self.tryAcceptor(atom2, atom1, amb2)
                        if not accflag:
                            self.unAdd(amb, atom1)
                                     
                elif atom2.hdonor:  # Treat Atom 1 as acceptor
                    self.debug("IN OPTION 2")
                    setflag = self.tryAcceptor(atom1, atom2, amb)
                    if setflag:
                        accflag = self.tryDonor(atom2, atom1, amb2)
                        if not accflag:
                            self.unAdd(amb, atom1)
                         
                else:
                    self.debug("Resolved in flipping")
                   
               
            # STEP 4: All Remaining water to water bonds
            self.debug("\nNOW MOVING ON TO ALL WAT to WAT bonds")
            watermap = {}
            bondlist = []
            for amb in cluster:
                if amb.fixed: continue
                residue = getattr(amb,"residue")
                nearatoms = getattr(amb,"nearatoms")
                if nearatoms == [] and amb.hdef.type == 12:
                    watermap[residue.getAtom("O")] = []
                for hbond in nearatoms:
                    if amb.hdef.type == 12:
                        if hbond.atom1 not in watermap: watermap[hbond.atom1] = []
                    if hbond.atom2.residue in clustermap and (clustermap[hbond.atom2.residue].hdef.type == 12 and amb.hdef.type == 12):
                        atom1 = hbond.atom1
                        atom2 = hbond.atom2
                        self.debug("%s %i %s %s %i %s %.3f" % (atom1.residue.name, atom1.residue.resSeq, atom1.name, atom2.residue.name, atom2.residue.resSeq, atom2.name, hbond.dist))
                        watermap[atom1].append(hbond)
                      

            while watermap != {}:
                maxnum = -1
                maxatom = None
                
                # Pick the guy with the most intrabonds placed
                for atom in watermap:
                    if len(atom.intrabonds) > maxnum:
                        maxnum = len(atom.intrabonds)
                        maxatom = atom

                # Make the list of nearby waters, put the best one first
                
                if "H2" in maxatom.intrabonds: # We're already done!
                    del watermap[maxatom]
                    continue
                    
                nearwaters = []
                bestdist = 999.9
                for i in range(len(watermap[maxatom])):
                    hbond = watermap[maxatom][i]  
                    atom1 = hbond.atom1
                    atom2 = hbond.atom2
                    dist = hbond.dist
                    
                    if dist < bestdist:
                        bestdist = dist
                        nearwaters = [atom2] + nearwaters
                    else:
                        nearwaters.append(atom2)

                # Finish this water
                self.debug("Working on %s %i %s" % (maxatom.residue.name, maxatom.residue.resSeq, maxatom.name))                
                    
                self.fixWater(maxatom, nearwaters) # Try adding hydrogen to atom1, LP to atom2
                del watermap[maxatom]

            # STEP 5:  Fix all remaining flips/ alcoholics
                
            for residue in clustermap:
                amb = clustermap[residue]
                if amb.hdef.type == 2:
                    self.fixAlcoholic(amb)
                if amb.hdef.type in [11,14] and not amb.fixed:
                    self.setFlip(amb)
                if amb.hdef.type == 13 and not amb.fixed:
                    self.setProt(amb)
                        
            # STEP 6:  Remove all LPs from waters/alcoholics, rename all FLIPS back to standard

            for residue in clustermap:
                amb = clustermap[residue]
                if amb.hdef.type == 12 or amb.hdef.type == 2: # Remove LPS
                    atomlist = []
                    for atom in residue.atoms: atomlist.append(atom)
                    for atom in atomlist:
                        if atom.name.startswith("LP"): residue.removeAtom(atom.name)
                elif amb.hdef.type in [11,14]: # Flip names back
                    for atom in residue.atoms:
                        atomname = atom.name
                        if atomname.endswith("FLIP"): residue.renameAtom(atomname, atomname[:-4])
                elif amb.hdef.type == 13:  # Remove the extra digit from the extra hyd
                    residue.renameAtom(amb.proth, amb.proth[:-1])

            self.debug("***** NETWORK COMPLETE *****\n")

      
     
    def unAdd(self, amb, atom):
        """
            Remove a recently added atom; Necessary for case
            where we are dealing with two ambiguities and the first is
            able to place a bond.

            Parameters
                amb:  The ambiguity in question
                atom: The atom that is bonded to the removed atom
        """
        residue = amb.residue
        type = amb.hdef.type
        if type == 11:
            amb.fixed = 0
            self.debug("Unfixed residue %i!" % residue.resSeq)
        else:
            rmname = atom.intrabonds.pop(-1)
            residue.removeAtom(rmname)
            self.debug("Removed %s!" % rmname)

    def setProt(self, amb):
        """
            Since the protontated residue cannot make any bonds, set
            to the best position as defined by the energy function
        """
        hbonds = amb.nearatoms
        bestenergy = 0
        bestatom = None
        bestnearatom = None
        for hbond in hbonds:
            atom = hbond.atom1
            nearatom = hbond.atom2
            if nearatom not in nearatom.residue.atoms: continue # Was deleted
            energy = self.getPairEnergy(atom, nearatom)
            if energy < bestenergy:
                bestatom = atom
                bestnearatom = nearatom
                bestenergy = energy

        # Now figure out which H it was - closest dist

        if bestatom == None:
            amb.fixProt(None)
            return

        besth = None
        bestdist = 99.9
        for bond in bestatom.intrabonds:
            if not bond.startswith("H"): continue
            dist = distance(bestatom.residue.getAtom(bond).getCoords(), bestnearatom.getCoords())
            if dist < bestdist:
                besth = bestatom.residue.getAtom(bond)
                bestdist = dist
        amb.fixProt(besth)
      
    def setFlip(self, amb):
        """
            Since the flip cannot make any bonds, set to the best position
            as defined by the energy function
        """
        hbonds = amb.nearatoms
        bestenergy = 0
        bestatom = None
        donorflag = 0
        for hbond in hbonds:
            atom = hbond.atom1
            nearatom = hbond.atom2
            if nearatom not in nearatom.residue.atoms: continue # Was deleted
            energy = self.getPairEnergy(atom, nearatom)
            if energy < bestenergy:
                bestatom = atom
                bestenergy = energy
                donorflag = 1   
            energy = self.getPairEnergy(nearatom, atom)
            if energy < bestenergy:
                bestatom = atom
                bestenergy = energy  
        amb.fixFlip(bestatom, donorflag)

    def printAtoms(self, clustermap):
        """
            For debugging, print all residues examined in this
            cluster to stdout
        """
        seen = []
        for residue in clustermap:
            if residue not in seen:
                seen.append(residue)
                for atom in residue.atoms:
                    print atom
            for hbond in clustermap[residue].nearatoms:
                if hbond.atom2.residue not in seen:
                    seen.append(hbond.atom2.residue)
                    for atom in hbond.atom2.residue.atoms:
                        print atom
         
    def tryDonor(self, atom1, atom2, amb):
        """
            Try setting atom1 as a donor, atom2 as an acceptor

            Parameters
                atom1:  The donor atom
                atom2:  The acceptor atom
                amb:    The ambiguity for the donor atom

            Returns
                result : 1 if an hbond is made
        """
        result = 0
        type = amb.hdef.type
        if not atom1.hdonor or not atom2.hacceptor: return result
        if type == 12:
            if self.optWater(atom1, atom2, 1): result = 1
        elif type == 2:
            if self.addAlcoholic(atom1, atom2, 1): result = 1
        elif type in [11,13,14]:
            hlist = self.isHbond(atom1, atom2)
            if hlist != []:
                result = 1
                if type == 13:
                    amb.fixProt(hlist[0])
                else:
                    amb.fixFlip(atom1, 1)
                
        else:
            print "Unable to process type %i!" % type
            sys.exit()
        return result
                        
    def tryAcceptor(self, atom1, atom2, amb):
        """
            Try setting atom1 as an acceptor, atom2 as a donor

            NOTE: For now we can get away with assuming that atom2's
                  hydrogens are already present (in rotate-to-backbone this
                  is always the case, and since we flip all rotate-to-waters
                  to water-to-rotate this will also work.

            Parameters
                atom1:  The acceptor atom
                atom2:  The donor atom
                amb:    The ambiguity for the acceptor atom
    
            Returns
                result : 1 if an hbond is made
        """
        result = 0
        if not atom1.hacceptor or not atom2.hdonor: return result
        type = amb.hdef.type
    
        if type == 2:#Try adding the lone pair
            for hyd in atom2.intrabonds:
                if hyd.startswith("H"):
                    if self.addAlcoholic(atom1, atom2.residue.getAtom(hyd), 0):
                        result = 1
                        return result
       
        hlist = self.isHbond(atom2, atom1)
        if hlist != []:
            if type in [11,14]:
                result = 1
                # Fix flip and remove extra atoms
                amb.fixFlip(atom1, 0)
            elif type == 13:
                result = 1
            elif type == 12: #Try adding the lone pair
                for hyd in atom2.intrabonds:
                    if hyd.startswith("H"): 
                        if self.optWater(atom1, atom2.residue.getAtom(hyd), 0):
                            result = 1
                            break
            elif type != 2:
                print "Unable to process type %i!" % type
                sys.exit()
        return result


    def makeAtomWithOneBond(self, oxygen, add, nearatoms, watflag):
        """
            Make an atom by rotating about the one existing bonded atom
            to the oxygen.

            Parameters
               oxygen:    The oxygen atom in question
               add:       The name of the atom to add
               nearatoms: A list of nearby atoms to evaluate
               watflag:   A flag; 1 if oxygen is water, 0 otherwise
        """
        residue = oxygen.residue
        bonds = oxygen.intrabonds
        refcoords = []
        defcoords = []

        defatomcoords = [0.9428,0,0]
        refcoords.append(oxygen.getCoords())
        refcoords.append(residue.getAtom(bonds[0]).getCoords())
        defcoords.append([0,0,.3333]) # Oxygen
        if watflag:
            defcoords.append([0,0,1.3333])
            atomtype = "HETATM"
        else:
            defcoords.append([0,0,1.7093]) # Not quite tetrahedral
            atomtype = "ATOM"

        newcoords = findCoordinates(2, refcoords, defcoords, defatomcoords)
        residue.createAtom(add,newcoords,atomtype) 
        oxygen.intrabonds.append(add)
        addatom = residue.getAtom(add)
        addatom.intrabonds.append("O")
          
        # Check to make sure this makes an H bond

        bestcoords = []
        besten = 999.9
        fixed = residue.getAtom(bonds[0])

        for i in range(36):
            self.rotateResidue(residue, fixed, oxygen, 10.0)          
            energy = 0
            for nearatom in nearatoms:
                if nearatom not in nearatom.residue.atoms: continue # Could happen due to flips
                if add.startswith("H"): # Use the energy function
                    energy += self.getPairEnergy(oxygen, nearatom) + self.getPairEnergy(nearatom, oxygen)
                else: # Adding LP, so put at minimum distance
                    energy += distance(addatom.getCoords(), nearatom.getCoords())
            if energy < besten:
                bestcoords = addatom.getCoords()
                besten = energy

        # Switch to bestcoords
        addatom.x = bestcoords[0]
        addatom.y = bestcoords[1]
        addatom.z = bestcoords[2]

    def getPositionsWithTwoBonds(self, oxygen):
        """
            Return the two positions that are possible when using the two
            existing bonded atoms

            Parameters
                oxygen:  The oxygen for the residue in question
            Returns:
                loc1:  Possible location 1 for the new atom
                loc2:  Possible location 2 for the new atom
        """
        # There are only two possible cases remaining - we can find them
        #   by rotating one of the existing bonds about the other
        residue = oxygen.residue
        bonds = oxygen.intrabonds
        
        fixed = residue.getAtom(bonds[0])
        rotate = residue.getAtom(bonds[1])
        origcoords = rotate.getCoords()
        
        # Rotate by 120 degrees twice
      
        self.rotateResidue(residue, fixed, oxygen, 120)
        loc1 = rotate.getCoords()  
        self.rotateResidue(residue, fixed, oxygen, 120)
        loc2 = rotate.getCoords()
    
        # Set rotate back to original
        rotate.x = origcoords[0]
        rotate.y = origcoords[1]
        rotate.z = origcoords[2]

        return loc1, loc2

    def fixPositionsWithTwoBonds(self, oxygen, add, loc1, loc2, nearatoms, watflag):
        """
            Fix the new atom position using either loc1 or loc2

            Paramters
                oxygen:  The oxygen for the residue in question
                add:     The name of the atom to add
                loc1:  Possible location 1 for the new atom
                loc2:  Possible location 2 for the new atom
                nearatoms:  A list of nearby atoms
                watflag: A flag; 1 if oxygen is water, 0 otherwise
        """
        # Try placing the atom to be added at each of the spaces
        if watflag == 0:
            atomtxt = "ATOM"
        else:
            atomtxt = "HETATM"
            
        residue = oxygen.residue
        residue.createAtom(add, loc1,atomtxt)
        oxygen.intrabonds.append(add)
        residue.getAtom(add).intrabonds.append("O")

        energy1 = 0
        for nearatom in nearatoms:
            if nearatom not in nearatom.residue.atoms: continue # Could happen due to flips
            energy1 += self.getPairEnergy(oxygen, nearatom) + self.getPairEnergy(nearatom, oxygen)
      
        # Try the other location
        residue.removeAtom(add)
        residue.createAtom(add, loc2, atomtxt)
        residue.getAtom(add).intrabonds.append("O")
        
        energy2 = 0
        for nearatom in nearatoms:
            if nearatom not in nearatom.residue.atoms: continue # Could happen due to flips
            energy2 += self.getPairEnergy(oxygen, nearatom) + self.getPairEnergy(nearatom, oxygen)
           
        # If the other one was better use that
        if energy1 < energy2:
            residue.removeAtom(add)
            residue.createAtom(add, loc1,atomtxt)
            residue.getAtom(add).intrabonds.append("O")

    def optPositionsWithTwoBonds(self, oxygen, add, loc1, loc2, nearatom, watflag):
        """
            Try adding an atom to either position loc1 or loc2

            Paramters
                oxygen:  The oxygen for the residue in question
                add:     The name of the atom to add
                loc1:  Possible location 1 for the new atom
                loc2:  Possible location 2 for the new atom
                nearatoms:  A list of nearby atoms
                watflag: A flag; 1 if oxygen is water, 0 otherwise

            Returns 1 if atom is added, 0 otherwise
        """
        # Try placing the atom to be added at each of the spaces
        if watflag == 0:
            atomtxt = "ATOM"
        else:
            atomtxt = "HETATM"
            
        residue = oxygen.residue
        residue.createAtom(add, loc1,atomtxt)
        oxygen.intrabonds.append(add)
        residue.getAtom(add).intrabonds.append("O")
        if nearatom.name.startswith("H"): hlist = self.isHbond(nearatom.residue.getAtom(nearatom.intrabonds[0]), oxygen)
        else: hlist = self.isHbond(oxygen, nearatom)
        
        if hlist != []: 
            self.debug("ADDED %s!" % add)
            return 1
        else: # Try the other location
            residue.removeAtom(add)
            residue.createAtom(add, loc2, atomtxt)
            residue.getAtom(add).intrabonds.append("O")
            if nearatom.name.startswith("H"): hlist2 = self.isHbond(nearatom.residue.getAtom(nearatom.intrabonds[0]), oxygen)
            else: hlist2 = self.isHbond(oxygen, nearatom)
            if hlist2 != []:
                self.debug("ADDED %s!" % add)
                return 1
            else: # No hbond here
                oxygen.intrabonds.pop(oxygen.intrabonds.index(add))
                residue.removeAtom(add)
                return 0

    def makeAtomWithThreeBonds(self, oxygen, add):
        """
            Since the atom already has three bonds, place in
            the lone remaining available spot.

            Parameters
                oxygen:  The oxygen to be added to
                add:     The name of the atom to be added
        """
        bonds = oxygen.intrabonds
        residue = oxygen.residue
        fixed = residue.getAtom(bonds[0])
        rotate1 = residue.getAtom(bonds[1])
        rotate2 = residue.getAtom(bonds[2])
        origcoords1 = rotate1.getCoords()
        origcoords2 = rotate2.getCoords()
        
        # Rotate by 120 degrees twice
        self.rotateResidue(residue, fixed, oxygen, 120)
        loc1 = rotate1.getCoords()
        self.rotateResidue(residue, fixed, oxygen, 120)
        loc2 = rotate1.getCoords()
        
        # Set rotates back to original
        rotate1.x = origcoords1[0]
        rotate1.y = origcoords1[1]
        rotate1.z = origcoords1[2]
        rotate2.x = origcoords2[0]
        rotate2.y = origcoords2[1]
        rotate2.z = origcoords2[2]
        
        # Check to see which position is already placed
        newloc = loc1
        if distance(origcoords2, loc2) > 0.5: newloc = loc2
        
        residue.createAtom(add, newloc,"ATOM") 
        oxygen.intrabonds.append(add)
        residue.getAtom(add).intrabonds.append("O")
        
    def fixAlcoholic(self, amb):
        """
            If H is not present, add it at the best energy location
        """
        residue = amb.residue
        if residue.name == "SER":
            oxygen = residue.getAtom("OG")
            add = "HG"
        elif residue.name == "THR":
            oxygen = residue.getAtom("OG1")
            add = "HG1"           
        elif residue.name == "TYR":
            oxygen = residue.getAtom("OH")
            add = "HH"
        if add in residue.map: return

        hbonds = amb.nearatoms
        nearatoms = []
        for hbond in hbonds:
            nearatoms.append(hbond.atom2)
           
        bonds = oxygen.intrabonds
        if len(bonds) == 1: # No H or LP attached
            self.makeAtomWithOneBond(oxygen, add, nearatoms, 0)
     
        elif len(bonds) == 2:  # One LP present, add H
            loc1, loc2 = self.getPositionsWithTwoBonds(oxygen)
            self.fixPositionsWithTwoBonds(oxygen, add, loc1, loc2, nearatoms, 0)

        elif len(bonds) == 3:  # Both LPs present, add H to remaining spot
            self.makeAtomWithThreeBonds(oxygen, add)

    def addAlcoholic(self, oxygen, nearatom, flag):
        """
            oxygen is the residue's oxygen
            nearatom is either an acceptor or a hydrogen that is
               to be donated
            if flag is 1 treat oxygen as donor
            else as acceptor
        """
        bonds = oxygen.intrabonds
        oxcoords = oxygen.getCoords()
        residue = oxygen.residue
        refcoords = []
        defcoords = []
        if residue.name == "SER": add = "HG"
        elif residue.name == "THR": add = "HG1"           
        elif residue.name == "TYR": add = "HH"

        if flag == 0 and ("LP2" in bonds or "LP1" in bonds): add = "LP2"
        elif flag == 0: add = "LP1"

        if add in bonds:
            self.debug("%s was already placed!" % add)
            return 0
        
        if len(bonds) == 1: # No H or LP attached
            self.makeAtomWithOneBond(oxygen, add, [nearatom], 0)
   
            # If this still isn't a bond, give up
            if nearatom.name.startswith("H"): hlist = self.isHbond(nearatom.residue.getAtom(nearatom.intrabonds[0]), oxygen)
            else: hlist = self.isHbond(oxygen, nearatom)
            
            if hlist == []:
                self.debug("Couldn't fix!")
                residue.removeAtom(add)
                oxygen.intrabonds.pop(oxygen.intrabonds.index(add))
                return 0
            else:
                self.debug("ADDED %s!" % add)
                return 1
            
        elif len(bonds) == 2:  # LP present, add H

            loc1, loc2 = self.getPositionsWithTwoBonds(oxygen)
            return self.optPositionsWithTwoBonds(oxygen, add, loc1, loc2, nearatom, 0)
    
        elif len(bonds) == 3:  # Both LPs present, add H to remaining spot
            self.makeAtomWithThreeBonds(oxygen, add)
           
            hlist = self.isHbond(oxygen, nearatom)
            if hlist != []:
                self.debug("ADDED %s!" % add)
                return 1
            else:
                residue.removeAtom(add)
                oxygen.intrabonds.pop(oxygen.intrabonds.index(add))
                return 0
                
    def optWater(self, oxygen, nearatom, flag):
        """
            oxygen is the water's oxygen
            nearatom is either an acceptor or a hydrogen that is
               to be donated
            if flag is 1 treat oxygen as donor
            else as acceptor

            The 4 coords for water are:
               [0.9428,0,0]
               [-0.4714, 0.8165, 0.000]
               [-0.4714,-0.8165, 0.000]
               [0.000, 0.000, 1.333]
               with oxygen at [0,0, 0.3333]
        """
        bonds = oxygen.intrabonds   
        residue = oxygen.residue
     
        if flag == 1 and ("H2" in bonds or "H1" in bonds): add = "H2"
        elif flag == 1: add = "H1"
        elif flag == 0 and ("LP2" in bonds or "LP1" in bonds): add = "LP2"
        elif flag == 0: add = "LP1"

        if add in bonds:
            self.debug("%s was already placed!" % add)
            return 0

        if len(bonds) == 0: # Use the near atom
            refcoords = []
            defcoords = []
            oxcoords = oxygen.getCoords()
            nearcoords = nearatom.getCoords()
            dist = distance(oxcoords, nearcoords)
            defcoords.append([0,0,0.3333+dist])
            defcoords.append([0,0,0.3333]) # Oxygen
            refcoords.append(nearcoords)
            refcoords.append(oxcoords)
            defatomcoords = [0,0,1.3333] # Location to be placed
            
            newcoords = findCoordinates(2, refcoords, defcoords, defatomcoords)
            residue.createAtom(add,newcoords,"HETATM") 
            oxygen.intrabonds.append(add)
            residue.getAtom(add).intrabonds.append("O")
            
            if nearatom.name.startswith("H"): hlist = self.isHbond(nearatom.residue.getAtom(nearatom.intrabonds[0]), oxygen)
            else: hlist = self.isHbond(oxygen, nearatom)
            
            if hlist == []:
                self.debug("Couldn't fix!")
                residue.removeAtom(add)
                oxygen.intrabonds.pop(oxygen.intrabonds.index(add))
                return 0
            else:
                self.debug("ADDED %s!" % add)
                return 1
            
        elif len(bonds) == 1:  # Use the other bond
            self.makeAtomWithOneBond(oxygen, add, [nearatom], 1)
            
            # If this still isn't a bond, give up
            if nearatom.name.startswith("H"): hlist = self.isHbond(nearatom.residue.getAtom(nearatom.intrabonds[0]), oxygen)
            else: hlist = self.isHbond(oxygen, nearatom)
            
            if hlist == []:
                self.debug("Couldn't fix!")
                residue.removeAtom(add)
                oxygen.intrabonds.pop(oxygen.intrabonds.index(add))
                return 0
            else:
                self.debug("ADDED %s!" % add)
                return 1

        elif len(bonds) == 2:  # We have LP1 and H1 already

            loc1, loc2 = self.getPositionsWithTwoBonds(oxygen)
            return self.optPositionsWithTwoBonds(oxygen, add, loc1, loc2, nearatom, 0)

    def fixWater(self, oxygen, nearatoms):
        """
            oxygen is the water's oxygen
            nearatom is either an acceptor or a hydrogen that is
               to be donated
            if flag is 1 treat oxygen as donor
            else as acceptor

          The 4 coords for water are:
               [0.9428,0,0]
               [-0.4714, 0.8165, 0.000]
               [-0.4714,-0.8165, 0.000]
               [0.000, 0.000, 1.333]
               with oxygen at [0,0, 0.3333]
        """
        bonds = oxygen.intrabonds
        residue = oxygen.residue
      
        if "H1" in bonds: add = "H2"
        else: add = "H1"
       
        if len(bonds) == 0: # Use the nearest atom

            refcoords = []
            defcoords = []
            oxcoords = oxygen.getCoords()
            defcoords.append([0,0,.3333]) # Oxygen
            refcoords.append(oxcoords)
            defatomcoords = [0,0,1.3333] # Location to be placed
            
            # Unless that nearest atom is a water donating an H to this atom!
            if nearatoms != []:
                hlist = self.isHbond(nearatoms[0], oxygen)
                nearcoords = nearatoms[0].getCoords()
                if hlist != []:
                    # Point away from this bond!
                    add = "LP1"
            else: # It doesn't matter where it points
                nearcoords = [1.0,1.0,1.0]
            dist = distance(oxcoords, nearcoords)
            defcoords.append([0,0,.3333+dist]) 
            refcoords.append(nearcoords)

            newcoords = findCoordinates(2, refcoords, defcoords, defatomcoords)
            residue.createAtom(add,newcoords,"HETATM") 
            oxygen.intrabonds.append(add)
            residue.getAtom(add).intrabonds.append("O")

            # Now add H2
            self.fixWater(oxygen, nearatoms)

        elif len(bonds) == 1:  # Use the other bond
            self.makeAtomWithOneBond(oxygen, add, nearatoms, 1)
            if add == "H1": self.fixWater(oxygen, nearatoms)

        elif len(bonds) == 2:  # We have LP1 and H1 already

            loc1, loc2 = self.getPositionsWithTwoBonds(oxygen)
            self.fixPositionsWithTwoBonds(oxygen, add, loc1, loc2, nearatoms, 0)


    def rotateResidue(self, residue, fixed, oxygen, newangle):
        """
            Rotate a residue about a bond to a certain angle

            Parameters
                residue:  The residue to be rotated
                fixed:    The name of the atom that is fixed
                oxygen:   The name of the oxygen; all other atoms
                          will be rotated about the fixed-oxygen bond
                newangle: The increment to change the angle
        """
        movenames = []
        movecoords = []
        
        initcoords = subtract(oxygen.getCoords(), fixed.getCoords())

        # Determine which atoms to rotate
        for atomname in oxygen.intrabonds:
            if atomname != fixed.name: movenames.append(atomname)

        for name in movenames:
            atom = residue.getAtom(name)
            movecoords.append(subtract(atom.getCoords(), fixed.getCoords()))

        newcoords = qchichange(initcoords, movecoords, newangle)
        for i in range(len(movenames)):
            name = movenames[i]
            atom = residue.getAtom(name)
            self.routines.removeCell(atom)
            x = (newcoords[i][0] + fixed.get("x"))
            y = (newcoords[i][1] + fixed.get("y"))
            z = (newcoords[i][2] + fixed.get("z"))
            atom.set("x", x)
            atom.set("y", y)
            atom.set("z", z)
            self.routines.addCell(atom)
        
    def isHbond(self, donor, acc):
        """
            Determine whether this donor acceptor pair is a
            hydrogen bond
        """
        hbonds = []
        donorhs = []
        acchs = []
        for bond in donor.get("intrabonds"):
            if bond[0] == "H": donorhs.append(bond)
        for bond in acc.get("intrabonds"):
            if bond[0] == "H": acchs.append(bond)
        for donorh in donorhs:
            donorhatom = donor.get("residue").getAtom(donorh)
            dist = distance(donorhatom.getCoords(), acc.getCoords())
            if dist > 3.3: continue
            if len(acchs) != 0:
                for acch in acchs:
                    acchatom = acc.get("residue").getAtom(acch)
                    hdist = distance(donorhatom.getCoords(), acchatom.getCoords())
                    if hdist < 1.5 and acc.residue.name not in ["HS2N", "GLH", "ASH"] \
                       and not acc.residue.isCterm: continue
                    angle = self.getHbondangle(acc, donor, donorhatom)
                    if angle > 20.0: continue
                    angle2 = self.getHbondangle(donorhatom, acchatom, acc)
                    if angle2 > 110.0: continue
                    hbonds.append(donorhatom)
                    self.debug("Found HBOND! %.4f %.4f" % (dist, angle))
            else:
                angle = self.getHbondangle(acc, donor, donorhatom)
                #print dist, angle
                if angle <= 21.0:
                    hbonds.append(donorhatom)
                    self.debug("Found HBOND! %.4f %.4f" % (dist, angle))
        return hbonds               

    def getPairEnergy(self, donor, acceptor):
        """
            Get the energy between two atoms

            Parameters
                donor:    The first atom in the pair (Atom)
                acceptor: The second atom in the pair (Atom)
            Returns
                energy:   The energy of the pair (float)
        """
        max_hbond_energy = -10.0
        max_ele_energy = -1.0
        maxangle = 20.0
        max_dha_dist = 3.3
        max_ele_dist = 5.0
        energy = 0.0
        donorh = None
        acceptorh = None
        donorhs = []
        acceptorhs = []
        
        if not (donor.get("hdonor") and acceptor.get("hacceptor")):
            return energy

        # See if hydrogens are presently bonded to the acceptor and donor

        for bond in donor.get("intrabonds"):
            if bond[0] == "H":
                donorhs.append(bond)
        for bond in acceptor.get("intrabonds"):
            if bond[0] == "H":
                acceptorhs.append(bond)
            
        if donorhs == []: return energy

        # Case 1: Both donor and acceptor hydrogens are present

        if acceptorhs != []:
            for donorh in donorhs:
                for acceptorh in acceptorhs:
                    donorhatom = donor.get("residue").getAtom(donorh)
                    acceptorhatom = acceptor.get("residue").getAtom(acceptorh)
                    if donorhatom == None or acceptorhatom == None:
                        text = "Couldn't find bonded hydrogen even though "
                        text = text + "it is present in intrabonds!"
                        raise ValueError, text
            
                    dist = distance(donorhatom.getCoords(), acceptor.getCoords())
                    if dist > max_dha_dist and dist < max_ele_dist: # Are the Hs too far?
                        energy += max_ele_energy/(dist*dist)
                        continue
                    
                    hdist = distance(donorhatom.getCoords(), acceptorhatom.getCoords())
                    if hdist < 1.5 and acceptor.residue.name not in ["HS2N", "GLH", "ASH"] \
                       and not acceptor.residue.isCterm:  # Are the Hs too close?
                        energy += -1 * max_hbond_energy
                        continue
                    
                    angle1 = self.getHbondangle(acceptor, donor, donorhatom)
                    if angle1 <= maxangle:
                        angleterm = (maxangle - angle1)/maxangle
                        angle2 = self.getHbondangle(donorhatom, acceptorhatom, acceptor)
                        if angle2 < 110.0: angle2 = 1.0
                        else: angle2=-1.0/110.0*(angle2-110.0)
                        energy+=max_hbond_energy/pow(dist,3)*angleterm*angle2
                    
            return energy

        # Case 2: Only the donor hydrogen is present

        elif acceptorhs == []:
            for donorh in donorhs:
                donorhatom = donor.get("residue").getAtom(donorh)
                if donorhatom == None:
                    text = "Couldn't find bonded hydrogen even though "
                    text = text + "it is present in intrabonds!"
                    raise ValueError, text
                
                dist = distance(donorhatom.getCoords(), acceptor.getCoords())
                if dist > max_dha_dist and dist < max_ele_dist: # Or too far?
                    energy += max_ele_energy/(dist*dist)               
                    continue
              
                angle1 = self.getHbondangle(acceptor, donor, donorhatom)
                if angle1 <= maxangle:
                    angleterm = (maxangle - angle1)/maxangle
                    energy += max_hbond_energy/pow(dist,2)*angleterm
                  
            return energy

    def getHbondangle(self, atom1, atom2, atom3):
        """
            Get the angle between three atoms

            Parameters
                atom1:  The first atom (atom)
                atom2:  The second (vertex) atom (atom)
                atom3:  The third atom (atom)
            Returns
                angle:  The angle between the atoms (float)
        """
        angle = 0.0
        atom2Coords = atom2.getCoords()
        coords1 = subtract(atom3.getCoords(), atom2Coords)
        coords2 = subtract(atom1.getCoords(), atom2Coords)
        norm1 = normalize(coords1)
        norm2 = normalize(coords2)
        dotted = dot(norm1, norm2)
        if dotted > 1.0: # If normalized, this is due to rounding error
            dotted = 1.0
        rad = abs(acos(dotted))
        angle = rad*180.0/pi
        if angle > 180.0:
            angle = 360.0 - angle
        return angle

    def getPenalty(self, residue):
        """
            Add penalties for unusual protonation states.

            Parameters
                atom:    The residue to examine (Atom)
            Returns
                penalty: The amount of the penalty (float)
        """
        acidpenalty = 25.0
        hispos = 0.1
        hisminus = 10.0
        nterm = 5.0
        penalty = 0.0

        resname = residue.get("name")

        if residue.get("isNterm"):
            charge = 1
            if residue.getAtom("H2") == None: charge = charge - 1
            if residue.getAtom("H3") == None: charge = charge - 1
            penalty = penalty + (1- charge)*nterm

        if resname == "HIS":
            hd1 = residue.getAtom("HD1")
            he2 = residue.getAtom("HE2")
            if hd1 == None and he2 == None: 
                penalty = penalty + hisminus
            elif hd1 != None and he2 != None:
                penalty = penalty + hispos

        if resname != "WAT":
            for atom in residue.get("atoms"):
                atomname = atom.get("name")
                if atomname in ["OD1","OD2","OE1","OE2","O","OXT"] and atom.get("hdonor"):
                    penalty = penalty + acidpenalty
                    break

        return penalty
                

    def findAmbiguities(self, water, pkaflag):
        """
            Find the amibiguities within a protein according to the
            DAT file, and set all boundatoms to their hydrogen donor/
            acceptor state.  Store the ambiguities as (residue, hydrodef)
            tuples in self.groups.

            Returns
                allatoms:  A list of all donors and acceptors in the
                           protein (list)
                water:     If 1, only put waters in groups, but fill allatoms
                           appropriately
                pkaflag:   If 1, ignore all protonation ambiguities.
        """
        self.routines.setDonorsAndAcceptors()
        allatoms = []
        hydrodefs = self.hydrodefs
        for chain in self.protein.getChains():
            for residue in chain.get("residues"):
                resname = residue.get("name")
                nter = residue.get("isNterm")
                cter = residue.get("isCterm")
                type = residue.get("type")
                if type == 2: continue
                for group in hydrodefs:
                    groupname = group.name
                    htype = group.type
                    if htype in [1,3,4] and pkaflag: continue

                    # FOR NOW REMOVE the other ones as well
                    if htype not in [2,11,14,13] and not pkaflag: continue
                    
                    if resname == groupname or \
                       (groupname == "APR") or \
                       (groupname == "APP" and resname != "PRO") or \
                       (groupname.endswith("FLIP") and resname == groupname[0:3]) or \
                       (groupname == "HISFLIP" and resname in ["HIP","HID","HIE","HSP","HSE","HSD"]) or \
                       (groupname == "NTR" and nter and resname != "PRO") or \
                       (groupname == "PNTR" and nter and resname == "PRO") or \
                       (groupname == "CTR" and cter) or \
                       (groupname == "CTN" and cter == 2):

                        if group.method != 0:
                            amb = hydrogenAmbiguity(residue, group, self.routines)
                            self.groups.append(amb)
        for atom in self.protein.getAtoms():
            if atom.residue.type not in [1,3]: continue
            if atom.name == "O" and atom.residue.name == "WAT": atom.hdonor = 1
            if atom.get("hacceptor") or atom.get("hdonor"): allatoms.append(atom)
            
        return allatoms
                                
    def printAmbiguities(self):
        """
            Print the list of ambiguities to stdout
        """
        if self.hdebug == 0: return
        i = 0
        for amb in self.groups:
            residue = getattr(amb,"residue")
            hydrodef = getattr(amb,"hdef")
            conf = hydrodef.conformations[0]
            self.routines.write("Ambiguity #: %i, chain: %s, residue: %i %s, hyd_type: %i state: %i Grp_name: %s Hname: %s Boundatom: %s\n" % (i, residue.chainID, residue.resSeq, residue.name, hydrodef.type, hydrodef.standardconf, hydrodef.name, conf.hname, conf.boundatom))
            i += 1
            
    def parseHydrogen(self, lines):
        """
            Parse a list of lines in order to make a hydrogen
            definition

            Parameters
                lines:  The lines to parse (list)
            Returns
                mydef:  The hydrogen definition object (HydrogenDefinition)
        """
        maininfo = string.split(lines[0])
        name = maininfo[0]
        group = maininfo[1]
        numhydrogens = int(maininfo[2])
        standardconf = int(maininfo[4])
        type = int(maininfo[5])
        chiangle = int(maininfo[6])
        method = int(maininfo[7])

        mydef = HydrogenDefinition(name, group, numhydrogens, standardconf, \
                                   type, chiangle, method)

        conf = []
        for newline in lines[1:]:
            if newline.startswith(">"):
                if conf == []: continue
                confinfo = string.split(conf[0])
                hname = confinfo[0]
                boundatom = confinfo[1]
                bondlength = float(confinfo[2])
                myconf = HydrogenConformation(hname, boundatom, bondlength)
                count = 0
                for line in conf[1:]:
                    textatom = string.split(line)
                    name = textatom[0]
                    x = float(textatom[1])
                    y = float(textatom[2])
                    z = float(textatom[3])
                    atom = DefinitionAtom(count, name, "", x,y,z)
                    myconf.addAtom(atom)

                mydef.addConf(myconf)
                conf = []
            else:
                conf.append(newline)
                
        if conf != []:
            confinfo = string.split(conf[0])
            hname = confinfo[0]
            boundatom = confinfo[1]
            bondlength = float(confinfo[2])
            myconf = HydrogenConformation(hname, boundatom, bondlength)
            count = 0
            for line in conf[1:]:
                textatom = string.split(line)
                name = textatom[0]
                x = float(textatom[1])
                y = float(textatom[2])
                z = float(textatom[3])
                atom = DefinitionAtom(count, name, "", x,y,z)
                myconf.addAtom(atom)
            mydef.addConf(myconf)
         
        if lines[1:] == []:  # FLIPS
            # ASN: OD1, ND2
            # GLN: OE1, NE2
            # HIS: ND1, CD2, CE1, NE2
            #boundatom = "CG"
            #if name == "GLNFLIP": boundatom = "CD"
            #myconf = HydrogenConformation(None, boundatom, 0.0)
            #mydef.addConf(myconf)
            if name == "ASNFLIP":
                for boundatom in ["OD1","ND2"]:
                    myconf = HydrogenConformation(None, boundatom, 0.0)
                    mydef.addConf(myconf)
            elif name == "GLNFLIP":
                for boundatom in ["OE1","NE2"]:
                    myconf = HydrogenConformation(None, boundatom, 0.0)
                    mydef.addConf(myconf)
            elif name == "HISFLIP":
                for boundatom in ["ND1","CD2","CE1","NE2"]:
                    myconf = HydrogenConformation(None, boundatom, 0.0)
                    mydef.addConf(myconf)      
        return mydef

    def readHydrogenDefinition(self):
        """
            Read the Hydrogen Definition file

            Returns
                hydrodef:  The hydrogen definition ()
        """
        defpath = HYDROGENFILE
        if not os.path.isfile(defpath):
            for path in sys.path:
                testpath = "%s/%s" % (path, defpath)
                if os.path.isfile(testpath):
                    defpath = testpath
                    break
        if not os.path.isfile(defpath):
            raise ValueError, "%s not found!" % defpath

       
        file = open(defpath)
        lines = file.readlines()
        file.close()
        info = []

        for line in lines:
            if line.startswith("//"): pass
            elif line.startswith("*") or line.startswith("!"):
                if info == []: continue
                mydef = self.parseHydrogen(info)
                self.hydrodefs.append(mydef)
                info = []
            else:
                info.append(string.strip(line))

class HydrogenDefinition:
    """
        HydrogenDefinition class

        The HydrogenDefinition class provides information on possible
        ambiguities in amino acid hydrogens.  It is essentially the hydrogen
        definition file in object form.
    """
    
    def __init__(self, name, group, numhydrogens, standardconf, type, \
                 chiangle, method):
        """
            Initialize the object with information from the definition file

            Parameters:
                name:          The name of the grouping (string)
                group:         The group of the definition
                               (acid/base/none, string)
                numhydrogens:  The number of hydrogens that can be added (int)
                standardconf:  The number of standard conformations (int)
                type        :  Type of Hydrogen (int)
                chiangle    :  The chiangle to be changed (int)
                method      :  The standard optimization method (int)

                See HYDROGENS.DAT for more information
        """
        self.name = name
        self.group = group
        self.numhydrogens = numhydrogens
        self.standardconf = standardconf
        self.type = type
        self.chiangle = chiangle
        self.method = method
        self.conformations = []

    def __str__(self):
        """
            Used for debugging purposes

            Returns
                output:  The information about this definition (string)
        """
        output =  "Name:                  %s\n" % self.name
        output += "Group:                 %s\n" % self.group
        output += "# of Hydrogens:        %i\n" % self.numhydrogens
        output += "# of Conformations:    %i\n" % len(self.conformations)
        output += "Standard Conformation: %i\n" % self.standardconf
        output += "Type of Hydrogen:      %i\n" % self.type
        output += "Chiangle to change:    %i\n" % self.chiangle
        output += "Optimization method:   %i\n" % self.method
        output += "Conformations:\n"
        for conf in self.conformations:
            output += "\n%s" % conf
        output += "*****************************************\n"
        return output

    def addConf(self, conf):
        """
            Add a HydrogenConformation to the list of conformations

            Parameters
                conf:  The conformation to be added (HydrogenConformation)
        """
        self.conformations.append(conf)

class HydrogenConformation:
    """
        HydrogenConformation class

        The HydrogenConformation class contains data about possible
        hydrogen conformations as specified in the hydrogen data file.
    """

    def __init__(self, hname, boundatom, bondlength):
        """
           Initialize the object

           Parameters
               hname      : The hydrogen name (string)
               boundatom  : The atom the hydrogen is bound to (string)
               bondlength : The bond length (float)
        """
        self.hname = hname
        self.boundatom = boundatom
        self.bondlength = bondlength
        self.atoms = []

    def __str__(self):
        """
            Used for debugging purposes

            Returns
                output:  Information about this conformation (string)
        """
        output  = "Hydrogen Name: %s\n" % self.hname
        output += "Bound Atom:    %s\n" % self.boundatom
        output += "Bond Length:   %.2f\n" % self.bondlength
        for atom in self.atoms:
            output += "\t%s\n" % atom
        return output
    
    def addAtom(self, atom):
        """
            Add an atom to the list of atoms

            Parameters
                atom: The atom to be added (DefinitionAtom)
        """
        self.atoms.append(atom)
        
