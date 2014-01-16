#!/usr/bin/python
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

#propka3.0, revision 182                                                                      2011-08-09
#-------------------------------------------------------------------------------------------------------
#--                                                                                                   --
#--                                   PROPKA: A PROTEIN PKA PREDICTOR                                 --
#--                                                                                                   --
#--                              VERSION 3.0,  01/01/2011, COPENHAGEN                                 --
#--                              BY MATS H.M. OLSSON AND CHRESTEN R. SONDERGARD                       --
#--                                                                                                   --
#-------------------------------------------------------------------------------------------------------
#
#
#-------------------------------------------------------------------------------------------------------
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
#-------------------------------------------------------------------------------------------------------


import sys, pdb, protonate, lib, bonds
from vector_algebra import *
pka_print = lib.pka_print

all_sybyl_types = [
    'C.3',   #  carbon sp3
    'H',     #  hydrogen
    'C.2',   #  carbon sp2
    'H.spc', #  hydrogen in Single Point Charge (SPC) water model
    'C.1',   #  carbon sp
    'H.t3p', #  hydrogen in Transferable intermolecular Potential (TIP3P) water model
    'C.ar',  #  carbon aromatic
    'LP',    #  lone pair
    'C.cat', #  carbocation (C+) used only in a guadinium group
    'Du',    #  dummy atom
    'N.3',   #  nitrogen sp3
    'Du.C',  #  dummy carbon
    'N.2',   #  nitrogen sp2
    'Any',   #  any atom
    'N.1',   #  nitrogen sp
    'Hal',   #  halogen
    'N.ar',  #  nitrogen aromatic
    'Het',   #  heteroatom = N, O, S, P
    'N.am',  #  nitrogen amide
    'Hev',   #  heavy atom (non hydrogen)
    'N.pl3', #  nitrogen trigonal planar
    'Li',    #  lithium
    'N.4',   #  nitrogen sp3 positively charged
    'Na',    #  sodium
    'O.3',   #  oxygen sp3
    'Mg',    #  magnesium
    'O.2',   #  oxygen sp2
    'Al',    #  aluminum
    'O.co2', #  oxygen in carboxylate and phosphate groups
    'Si',    #  silicon
    'O.spc', #  oxygen in Single Point Charge (SPC) water model
    'K',     #  potassium
    'O.t3p', #  oxygen in Transferable Intermolecular Potential (TIP3P) water model
    'Ca',    #  calcium
    'S.3',   #  sulfur sp3
    'Cr.th', #  chromium (tetrahedral)
    'S.2',   #  sulfur sp2
    'Cr.oh', #  chromium (octahedral)
    'S.O',   #  sulfoxide sulfur
    'Mn',    #  manganese
    'S.O2',  #  sulfone sulfur
    'Fe',    #  iron
    'P.3',   #  phosphorous sp3
    'Co.oh', #  cobalt (octahedral)
    'F',     #  fluorine
    'Cu',    #  copper
    'Cl',    #  chlorine
    'Zn',    #  zinc
    'Br',    #  bromine
    'Se',    #  selenium
    'I',     #  iodine
    'Mo',    #  molybdenum
    'Sn']    #  tin


#propka_input_types = ['1P','1N','2P','2N']
#for type in all_sybyl_types:
#    temp = type.replace('.','')
#    if len(temp)>3:
#        temp = temp[0:3]
#    propka_input_types.append(temp)
#
#for t in propka_input_types:
#    print (t)


propka_input_types = [
    '1P',
    '1N',
    '2P',
    '2N',
    'C3',
    'H',
    'C2',
    'Hsp',
    'C1',
    'Ht3',
    'Car',
    'LP',
    'Cca',
    'Du',
    'N3',
    'DuC',
    'N2',
    'Any',
    'N1',
    'Hal',
    'Nar',
    'Het',
    'Nam',
    'Hev',
    'Npl',
    'Li',
    'N4',
    'Na',
    'O3',
    'Mg',
    'O2',
    'Al',
    'Oco',
    'Si',
    'Osp',
    'K',
    'Ot3',
    'Ca',
    'S3',
    'Crt',
    'S2',
    'Cro',
    'SO',
    'Mn',
    'SO2',
    'Fe',
    'P3',
    'Coo',
    'F',
    'Cu',
    'Cl',
    'Zn',
    'Br',
    'Se',
    'I',
    'Mo',
    'Sn']


ions = ['CA','NA']

max_C_double_bond = 1.3
max_C_triple_bond = 1.2

max_C_double_bond_squared = max_C_double_bond*max_C_double_bond
max_C_triple_bond_squared = max_C_triple_bond*max_C_triple_bond

class ligand:
    def __init__(self, atoms):
        self.atoms = atoms
        for atom in self.atoms:
            atom.residue = self

        #self.remove_ions()
        self.configuration_keys = atoms[0].configurations.keys()
        
        # create ligand residue objects
        self.ligand_residues = []
        self.split_into_residues()

        return

    def __str__(self):
        res  = '----Ligand----\n'
        for atom in self.atoms:
            res += '%s\n'%atom
        res+='--------------'
        return res


    def split_into_residues(self):
        
        residue = []
        if len(self.atoms)>0:
            current_residue_number = self.atoms[0].resNumb

        for atom in self.atoms:
            if atom.resNumb != current_residue_number:
                self.ligand_residues.append(ligand_residue(residue))
                residue = []
                current_residue_number = atom.resNumb

            residue.append(atom)

        # remember to include the last ligand residue
        self.ligand_residues.append(ligand_residue(residue))


        return


    def remove_ions(self):
        self.atoms = [atom for atom in self.atoms if atom.get_element() not in ions]
        return

    def writePDB(self, pdbname):
        pka_print("writing pdbfile %s" % (pdbname))
        file = open(pdbname, 'w')

        configurations = lib.get_sorted_configurations(self.configuration_keys)
        if len(configurations)==1:        
            self.write_atoms(file)
        else:
            for configuration in configurations:
                self.setConfiguration(configuration)
                file.write('MODEL%9d\n'%int(configuration[1]))
                self.write_atoms(file)
                file.write('ENDMDL\n')
        file.close()

        return
    
    def write_atoms(self,file):
        atom_number=1
        for atom in my_ligand.atoms:
            atom.writePDB(file,atom_number)
            atom_number+=1
        return


    def assign_atom_names(self):
        """
        Assigns sybyl names to ligand atoms based on elements and coordinates
        copied from propka/ligand.py
        modified:
        P -> P3
        O.co2+ - > O.co2
        """
        # find bonding atoms
        self.my_bond_maker = bonds.bondmaker()
        self.my_bond_maker.find_bonds_for_atoms(self.atoms)
        for atom in self.atoms:
            # check if we already have assigned a name to this atom
            if hasattr(atom, 'sybyl_assigned'):
                print(atom.resName, atom.numb, atom.name, 'alreadyassigned')
                continue

            # find some properties of the atom
            ring_atoms = self.is_ring_member(atom)
            aromatic = self.is_aromatic_ring(ring_atoms)
            planar = self.is_planar(atom)
            bonded_elements = {}
            for i in range(len(atom.bonded_atoms)):
                bonded_elements[i] = atom.bonded_atoms[i].get_element()

            # Aromatic carbon/nitrogen
            if aromatic:
                #print "--- if aromatic",atom.resName, atom.numb, atom.name, atom.get_element()
                #for ra in ring_atoms:
                #    if ra.get_element() in ['C', 'N']:
                #        self.set_type(ra, ra.get_element() + '.ar')
                #continue
                '''SH: In the original version, if a ring is planar (eg., 4FXF_D_FBP_606)
                but contains other atoms than C or N the loop breaks without assigning these atoms
                '''
                if atom.get_element() in ['C', 'N']:
                    self.set_type(atom, atom.get_element() + '.ar')
                    continue

            # check for amide
            if atom.get_element() in ['O', 'N']:
                amide = 0
                for b in atom.bonded_atoms:
                    if b.element == 'C':
                        for bb in b.bonded_atoms:
                            if (bb.get_element() == 'N' and atom.get_element() == 'O'):
                                self.set_type(bb, 'N.am')
                                self.set_type(b, 'C.2')
                                self.set_type(atom, 'O.2')
                                amide = 1
                            if (bb.get_element() == 'O' and atom.get_element() == 'N'):
                                self.set_type(atom, 'N.am')
                                self.set_type(b, 'C.2')
                                self.set_type(bb, 'O.2')
                                amide = 1
                if amide == 1:
                    continue


            if atom.get_element() == 'C':
                # check for amide
                if 'O' in bonded_elements.values() and 'N' in bonded_elements.values():
                    self.set_type(atom, 'C.2')
                    for b in atom.bonded_atoms:
                        if b.get_element() == 'N':
                            self.set_type(b, 'N.am')
                        if b.get_element() == 'O':
                            self.set_type(b, 'O.2')
                    continue

                # check for carboxyl
                if len(atom.bonded_atoms) == 3 and list(bonded_elements.values()).count('O') == 2:
                    i1 = list(bonded_elements.values()).index('O')
                    i2 = list(bonded_elements.values()).index('O', i1 + 1)
                    if len(atom.bonded_atoms[i1].bonded_atoms) == 1 and len(atom.bonded_atoms[i2].bonded_atoms) == 1:
                        #self.set_type(atom.bonded_atoms[i1], 'O.co2+')
                        self.set_type(atom.bonded_atoms[i1], 'O.co2') #SH: problems in 1a3w
                        self.set_type(atom.bonded_atoms[i2], 'O.co2')
                        self.set_type(atom, 'C.2')
                        continue

                # sp carbon
                if len(atom.bonded_atoms) <= 2:
                    for b in atom.bonded_atoms:
                        if self.my_bond_maker.squared_distance(atom, b) < max_C_triple_bond_squared:
                            self.set_type(atom, 'C.1')
                            self.set_type(b, b.get_element() + '.1')
                    if hasattr(atom, 'sybyl_assigned'):
                        continue

                # sp2 carbon
                if planar:
                    self.set_type(atom, 'C.2')
                    # check for N.pl3
                    for b in atom.bonded_atoms:
                        if b.get_element() == 'N':
                            if len(b.bonded_atoms) < 3 or self.is_planar(b):
                                self.set_type(b, 'N.pl3')
                    continue

                # sp3 carbon
                self.set_type(atom, 'C.3')
                continue

            # Nitrogen
            if atom.get_element() == 'N':
                # check for planar N
                if len(atom.bonded_atoms) == 1:
                    if self.is_planar(atom.bonded_atoms[0]):
                        self.set_type(atom, 'N.pl3')
                        continue

                if planar:
                    self.set_type(atom, 'N.pl3')
                    continue

                self.set_type(atom, 'N.3')
                continue

            # Oxygen
            if atom.get_element() == 'O':
                self.set_type(atom, 'O.3')
                # check for X=O
                if len(atom.bonded_atoms) == 1:
                    if self.my_bond_maker.squared_distance(atom, atom.bonded_atoms[0]) < max_C_double_bond_squared:
                        self.set_type(atom, 'O.2')
                        if atom.bonded_atoms[0].get_element() == 'C':
                            self.set_type(atom.bonded_atoms[0], 'C.2')
                continue


            # Sulphur
            if atom.get_element() == 'S':
                #check for SO2
                if list(bonded_elements.values()).count('O') == 2:
                        i1 = list(bonded_elements.values()).index('O')
                        i2 = list(bonded_elements.values()).index('O', i1 + 1)
                        self.set_type(atom.bonded_atoms[i1], 'O.2')
                        self.set_type(atom.bonded_atoms[i2], 'O.2')
                        self.set_type(atom, 'S.o2')
                        continue

                # check for SO4
                if list(bonded_elements.values()).count('O') == 4:
                    no_O2 = 0
                    for i in range(len(atom.bonded_atoms)):
                        if len(atom.bonded_atoms[i].bonded_atoms) == 1 and no_O2 < 2:
                            self.set_type(atom.bonded_atoms[i], 'O.2')
                            no_O2 += 1
                        else:
                            self.set_type(atom.bonded_atoms[i], 'O.3')

                self.set_type(atom, 'S.3')
                continue

            # Phosphorous (phosphorous sp3)
            #@attention: This was added and may not consider all types of phosphorous
            if atom.get_element() == 'P':
                self.set_type(atom, 'P.3')
                continue
            

            element = atom.get_element().capitalize()
            self.set_type(atom, element)
            #pka_print('Using element as type for %s'%atom.get_element())
        return

    def set_type(self,atom,type):
        pka_print(atom,'->',type)
        atom.name = type
        atom.sybyl_assigned=1
        return


    def is_ring_member(self, atom):
        return self.identify_ring(atom,atom,0,[])

    def identify_ring(self, this_atom, original_atom, number, past_atoms):
        number+=1
        past_atoms=past_atoms+[this_atom]
        return_atoms = []

        for atom in this_atom.bonded_atoms:
            if atom == original_atom and number>2:
                return past_atoms

            if atom not in past_atoms:
                these_return_atoms = self.identify_ring(atom, original_atom, number, past_atoms)
                if len(these_return_atoms) > 0:
                    if len(return_atoms)>len(these_return_atoms) or len(return_atoms)==0:
                        return_atoms = these_return_atoms

        return return_atoms


    def is_planar(self, atom):
        """ Finds out if atom forms a plane together with its bonded atoms"""
        atoms = [atom]+atom.bonded_atoms
        return self.are_atoms_planar(atoms)

    def are_atoms_planar(self, atoms):
        if len(atoms)==0:
            return False
        if len(atoms)<4:
            return False
        v1 = vector(atom1=atoms[0], atom2=atoms[1])
        v2 = vector(atom1=atoms[0], atom2=atoms[2])
        n = (v1**v2).rescale(1.0)

        margin = 0.20
        for b in atoms[3:]:
            v = vector(atom1=atoms[0], atom2=b).rescale(1.0)
            #pka_print(atoms[0],abs(v*n) )
            if abs(v*n)>margin:
                return False
            
        return True

    def is_aromatic_ring(self, atoms):
        if len(atoms)<5:
            return False
        
        for i in range(len(atoms)):
            if not self.are_atoms_planar(atoms[i:]+atoms[:i]):
                return False

        return True



class ligand_residue:
    def __init__(self, atoms):
        self.resNumb = -1
        self.resName = ''

        self.atoms = atoms

        if len(self.atoms)>0:
            self.resNumb = self.atoms[0].resNumb
            self.resName = self.atoms[0].resName

        pka_print('Created ligand residue %s with %2d atoms'%(self, len(self.atoms)))

        return

    def __str__(self):
        return '%s-%4d'%(self.resName,self.resNumb)


if __name__ == '__main__':
    if len(sys.argv)<2:
        sys.exit(0)

    protonator = protonate.Protonate()

    pka_print(sys.argv[1])
    atoms = pdb.readPDB(sys.argv[1],tags=["ATOM","HETATM"])
    
    my_ligand = ligand(atoms)
    
    

    #assign sybyl names
    protonator.remove_all_hydrogen_atoms_from_ligand(my_ligand)
    my_ligand.assign_atom_names()

    my_ligand.writePDB('before_ligand_protonation.pdb')

    #convert to propka input names
#    for atom in my_ligand.atoms:
#        temp = atom.name
#        temp = temp.replace('.','')
#        if len(temp)>3:
#            temp = temp[0:3]
#        atom.name = temp

    # protonate
    protonator.protonate_ligand(my_ligand)
    my_ligand.writePDB('after_ligand_protonation.pdb')

