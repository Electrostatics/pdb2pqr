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


import  sys, math
from vector_algebra import *
import bonds as bonds
import pdb as pdb


def makeProtonator(scheme=None):
    """
    selects a protonator class based on 'scheme'
    """
    if scheme in ["old-school", "old"]:
      return  old_scheme()
    else:
      return  new_scheme()




class old_scheme:
    """ 
    Protonates a protein using the old propka2-scheme
    """

    def __init__(self, atoms=None, pdbfile=None, name=None, options=None):
        """
        constructer of the protein object.
        """
        self.name               = "old-school protonator"
        self.protonate_residues = ['ASN', 'GLN', 'TRP', 'HIS', 'ARG']

        return


    def protonate(self, protein=None):
        """
        protonate a given protein
        """
        C = None; O = None

        for chain in protein.chains:
          C = None; O = None
          for residue in chain.residues:
            if residue.type == "amino-acid":
              C, O = self.protonateBackBone(residue, C, O)
              if   residue.resType == "AMD":
                self.protonateAMD(residue)
              elif residue.resType == "TRP":
                self.protonateTRP(residue)
              elif residue.resType == "HIS":
                self.protonateHIS(residue)
              elif residue.resType == "ARG":
                self.protonateARG(residue)
              elif residue.resName in self.protonate_residues:
                print("no protocol to protonate '%s' in 'old-scheme'" % (residue.label))
                sys.exit(8)

        return


    def protonateBackBone(self, residue, C, O):
        """
        Protonates an atom, X1, given a direction (X2 -> X3) [X1, X2, X3]
        """
        N = residue.getAtom(name='N')

        if C == None and O == None:
          """ do nothing, first residue """
        elif N == None:
          print( "could not find N atom in '%s' (protonateBackBone())" % (residue.label) ); sys.exit(9)
        elif residue.resName == "PRO":
          """ do nothing, proline doesn't have a proton """
        else:
          H = self.protonateDirection(atoms=[N, O, C])
          residue.atoms.append(H)

        return  residue.getAtom(name='C'), residue.getAtom(name='O')


    def protonateAMD(self, residue):
        """
        Protonates the side-chain of 'ASN' and 'GLN'
        """
        if residue.resName == "ASN":
          C = residue.getAtom(name='CG')
          O = residue.getAtom(name='OD1')
          N = residue.getAtom(name='ND2')
          H1 = self.protonateDirection(name='1HD2', atoms=[N, O, C])
          H2 = self.protonateAverageDirection(name='2HD2', atoms=[N, C, O])
        elif residue.resName == "GLN":
          C = residue.getAtom(name='CD')
          O = residue.getAtom(name='OE1')
          N = residue.getAtom(name='NE2')
          H1 = self.protonateDirection(name='1HE2', atoms=[N, O, C])
          H2 = self.protonateAverageDirection(name='2HE2', atoms=[N, C, O])

        residue.atoms.extend([H1, H2])

        return


    def protonateTRP(self, residue):
        """
        Protonates the side-chain of 'TRP'
        """
        CD = residue.getAtom(name='CD1')
        NE = residue.getAtom(name='NE1')
        CE = residue.getAtom(name='CE2')

        HE = self.protonateSP2(name='HE1', atoms=[CD, NE, CE])
        residue.atoms.append(HE)

        return


    def protonateHIS(self, residue):
        """
        Protonates the side-chain of 'HIS'
        """
        CG = residue.getAtom(name='CG')
        ND = residue.getAtom(name='ND1')
        CD = residue.getAtom(name='CD2')
        CE = residue.getAtom(name='CE1')
        NE = residue.getAtom(name='NE2')

        HD = self.protonateSP2(name='HD1', atoms=[CG, ND, CE])
        HE = self.protonateSP2(name='HE2', atoms=[CD, NE, CE])

        residue.atoms.extend([HD, HE])

        return


    def protonateARG(self, residue):
        """
        Protonates the side-chain of 'ARG'
        """
        CD  = residue.getAtom(name='CD')
        CZ  = residue.getAtom(name='CZ')
        NE  = residue.getAtom(name='NE')
        NH1 = residue.getAtom(name='NH1')
        NH2 = residue.getAtom(name='NH2')

        H1 = self.protonateSP2(name='HE', atoms=[CD, NE, CZ])
        H2 = self.protonateDirection(name='1HH1', atoms=[NH1, NE, CZ])
        H3 = self.protonateDirection(name='2HH1', atoms=[NH1, NE, CD])
        H4 = self.protonateDirection(name='1HH2', atoms=[NH2, NE, CZ])
        H5 = self.protonateDirection(name='2HH2', atoms=[NH2, NE, H1])

        residue.atoms.extend([H1, H2, H3, H4, H5])

        return


    def protonateDirection(self, name="H", atoms=None):
        """
        Protonates an atom, X1, given a direction (X2 -> X3) [X1, X2, X3]
        """
        X1, X2, X3 = atoms
        H  = X1.makeCopy(name=name, element='H')

        for key in H.configurations.keys():
          for atom in [H, X1, X2, X3]:
            atom.setConfiguration(key)

          dX = (X3.x - X2.x)
          dY = (X3.y - X2.y)
          dZ = (X3.z - X2.z)
          length = math.sqrt( dX*dX + dY*dY + dZ*dZ )
          H.x += dX/length
          H.y += dY/length
          H.z += dZ/length

          H.setConfigurationPosition(key)

        return  H


    def protonateAverageDirection(self, name="H", atoms=None):
        """
        Protonates an atom, X1, given a direction (X1/X2 -> X3) [X1, X2, X3]
        Note, this one uses the average of X1 & X2 (N & O) unlike the previous
        N - C = O
        """
        X1, X2, X3 = atoms
        H  = X1.makeCopy(name=name, element='H')

        for key in H.configurations.keys():
          for atom in [H, X1, X2, X3]:
            atom.setConfiguration(key)

          dX = (X3.x + X1.x)*0.5 - X2.x
          dY = (X3.y + X1.y)*0.5 - X2.y
          dZ = (X3.z + X1.z)*0.5 - X2.z

          length = math.sqrt( dX*dX + dY*dY + dZ*dZ )
          H.x += dX/length
          H.y += dY/length
          H.z += dZ/length

          H.setConfigurationPosition(key)

        return H


    def protonateSP2(self, name="H", atoms=None):
        """
        Protonates a SP2 atom, H2, given a list of [X1, X2, X3]
        X1-X2-X3
        """
        X1, X2, X3 = atoms
        H  = X2.makeCopy(name=name, element='H')

        for key in H.configurations.keys():
          for atom in [H, X1, X2, X3]:
            atom.setConfiguration(key)

          dX = (X1.x + X3.x)*0.5 - X2.x
          dY = (X1.y + X3.y)*0.5 - X2.y
          dZ = (X1.z + X3.z)*0.5 - X2.z
          length = math.sqrt( dX*dX + dY*dY + dZ*dZ )
          H.x -= dX/length
          H.y -= dY/length
          H.z -= dZ/length

          H.setConfigurationPosition(key)

        return H






class new_scheme:
    """
    Protonate atoms according to VSEPR theory 
    """
    
    def __init__(self):

        self.name              = "new-school protonator"
        self.valence_electrons = {'H':1,
                                  'C':4,
                                  'N':5,
                                  'O':6,
                                  'F':7,
                                  'P':5,
                                  'S':6,
                                  'CL':7}

        self.standard_charges= {'ARG-NH1':1.0,
                                'ASP-OD2':-1.0,
                                'GLU-OE2':-1.0,
                                'HIS-NE2':0.0,
                                'LYS-NZ':1.0,
                                'NTERM':1.0,
                                'CTERM':-1.0} 
        
        self.sybyl_charges = {'N.pl3':+1,
                              'N.3':+1,
                              'N.4':+1,
                              'N.ar':+1,
                              'O.co2+':-1}



#        self.standard_conjugate_charges= {'ARG-NH1':1.0}


        self.bond_lengths = {'C':1.09,
                             'N':1.01,
                             'O':0.96,
                             'F':0.92,
                             'Cl':1.27,
                             'Br':1.41,
                             'I':1.61}

        self.ions = ['NA','CA']


        # protonation_methods[steric_number] = method
        self.protonation_methods = {4:self.tetrahedral,
                                    3:self.trigonal}




        self.my_bond_maker = bonds.bondmaker()
        return



    def protonate(self, protein=None):
        """
        protonate a given protein
        """
        self.protonate_protein(protein)

        return


    def protonate_protein(self, protein):
        """ Will protonate all atoms in the protein """

        #print('----- Protontion started -----')
        # Remove all currently present hydrogen atoms
        #self.remove_all_hydrogen_atoms_from_protein(protein)
     
        # make bonds
        self.my_bond_maker.find_bonds_for_protein(protein)

        # set charges
        self.set_charges(protein)
        
        # protonate all atom
        non_H_atoms = []
        for chain in protein.chains:
            for residue in chain.residues:
                 if   residue.type == 'amino-acid':
                     for atom in residue.atoms:
                         non_H_atoms.append(atom)

        for atom in non_H_atoms:
          self.protonate_atom(atom)
            
        # set correct hydrogen names
        self.set_proton_names(non_H_atoms)

        return


    def protonate_ligand(self, ligand):
        """ Will protonate all atoms in the ligand """

        #print('----- Protonation started -----')
        # Remove all currently present hydrogen atoms
        self.remove_all_hydrogen_atoms_from_ligand(ligand)

        print(ligand)

        # make bonds
        self.my_bond_maker.find_bonds_for_ligand(ligand)

        # set charges
        self.set_ligand_charges(ligand)

        # protonate all atoms
        atoms = []
        for atom in ligand.atoms:
            if atom.type == 'atom':
                atoms.append(atom)
        for atom in atoms:
            self.protonate_atom(atom)
  
        # fix hydrogen names
        self.set_proton_names(ligand.atoms)
                    
        return


    def remove_all_hydrogen_atoms_from_protein(self, protein):
        for chain in protein.chains:
            for residue in chain.residues:
                residue.atoms = [atom for atom in residue.atoms if atom.get_element() != 'H']

        return


    def remove_all_hydrogen_atoms_from_ligand(self, ligand):
        ligand.atoms = [atom for atom in ligand.atoms if atom.get_element() != 'H']

        return
    

    def set_ligand_charges(self, ligand, standard_protonation_states = 1):
        if standard_protonation_states:
            for atom in ligand.atoms:
                #print('Charge before', atom, atom.charge)
                if atom.name in list(self.sybyl_charges.keys()):
                    atom.charge = self.sybyl_charges[atom.name]
                    #print('Charge', atom, atom.charge)

        else:
            print('Custom protonation state choosen - don\'t know what to do')
        
        return


    def set_charges(self, protein, standard_protonation_states = 1):
        if standard_protonation_states:
            # set side chain charges
            for chain in protein.chains:
                for residue in chain.residues:
                    for atom in residue.atoms:
                        key = '%3s-%s'%(atom.resName, atom.name)
                        if key in list(self.standard_charges.keys()):
                            atom.charge = self.standard_charges[key]
                            #print('Charge', atom, atom.charge)

            # set n-terminal charges
            for chain in protein.chains:
                for residue in chain.residues:
                    if residue.resName.replace(' ','') == 'N+':
                        for atom in residue.atoms:
                            if atom.name == 'N':
                                atom.charge = self.standard_charges['NTERM']
                                #print('Charge', atom, atom.charge)

            # set c-terminal charges
            for chain in protein.chains:
                for residue in chain.residues:
                    if residue.resName.replace(' ','') == 'C-':
                        for atom in residue.atoms:
                            if atom.name in self.my_bond_maker.terminal_oxygen_names: 
                                atom.charge = self.standard_charges['CTERM']
                                #print('Charge', atom, atom.charge)

        else:
            print('Custom protonation state choosen - don\'t know what to do')

        return


    def protonate_atom(self, atom):
        self.set_number_of_protons_to_add(atom)
        self.set_steric_number_and_lone_pairs(atom)
        self.add_protons(atom)

        return


    def set_proton_names(self, heavy_atoms):
        """
        setting an official pdb atom name
        """
        for heavy_atom in heavy_atoms:
            # counting connected hydrogen atoms
            hydrogens = 0; i = 0
            for bonded in heavy_atom.bonded_atoms:
              if bonded.element == 'H':
                hydrogens += 1
            for bonded in heavy_atom.bonded_atoms:
              if bonded.element == 'H':
                name = ""
                if hydrogens > 1:
                  i += 1
                  name += "%d" % (i)
                name += "H"
                name += heavy_atom.name[1:]
                bonded.setProperty(name=name)
                        
        return


    def set_number_of_protons_to_add(self, atom):
        #print('*'*10)
        #print('Setting number of protons to add for',atom)
        atom.number_of_protons_to_add  = 8 
        #print('                  %4d'%8)
        atom.number_of_protons_to_add -= self.valence_electrons[atom.get_element()]
        #print('Valence eletrons: %4d'%-self.valence_electrons[atom.get_element()])
        atom.number_of_protons_to_add -= len(atom.bonded_atoms)
        #print('Number of bonds:  %4d'%- len(atom.bonded_atoms))
        atom.number_of_protons_to_add -= atom.number_of_pi_electrons_in_double_and_triple_bonds
        #print('Pi electrons:     %4d'%-atom.number_of_pi_electrons_in_double_and_triple_bonds)
        atom.number_of_protons_to_add += int(atom.charge)
        #print('Charge:           %4.1f'%atom.charge)

        #print('-'*10)
        #print(atom.number_of_protons_to_add)

        return


    def set_steric_number_and_lone_pairs(self, atom):
        #print('='*10)
        #print('Setting steric number and lone pairs for',atom)

        # costumly set the N backbone atoms up for peptide bond trigonal planer shape
        #if atom.name == 'N' and len(atom.bonded_atoms) == 2:
        #    atom.steric_number = 3
        #    atom.number_of_lone_pairs = 0
        #    print 'Peptide bond: steric number is %d and number of lone pairs is %s'%(atom.steric_number,
         #                                                                             atom.number_of_lone_pairs)
        #    return


        atom.steric_number = 0
        
        #print('%65s: %4d'%('Valence electrons',self.valence_electrons[atom.get_element()]))
        atom.steric_number += self.valence_electrons[atom.get_element()]
        
        #print('%65s: %4d'%('Number of bonds',len(atom.bonded_atoms)))
        atom.steric_number += len(atom.bonded_atoms)

        #print('%65s: %4d'%('Number of hydrogen atoms to add',atom.number_of_protons_to_add))
        atom.steric_number += atom.number_of_protons_to_add

        #print('%65s: %4d'%('Number of pi-electrons in double and triple bonds(-)',atom.number_of_pi_electrons_in_double_and_triple_bonds))
        atom.steric_number -= atom.number_of_pi_electrons_in_double_and_triple_bonds

        #print('%65s: %4d'%('Number of pi-electrons in conjugated double and triple bonds(-)',atom.number_of_pi_electrons_in_conjugate_double_and_triple_bonds))
        atom.steric_number -= atom.number_of_pi_electrons_in_conjugate_double_and_triple_bonds

        #print('%65s: %4d'%('Number of donated co-ordinated bonds',0))
        atom.steric_number += 0

        #print('%65s: %4.1f'%('Charge(-)',atom.charge))
        atom.steric_number -= atom.charge
        
        atom.steric_number = math.floor(atom.steric_number/2.0)

        atom.number_of_lone_pairs = atom.steric_number - len(atom.bonded_atoms) - atom.number_of_protons_to_add

        #print('-'*70)
        #print('%65s: %4d'%('Steric number',atom.steric_number))
        #print('%65s: %4d'%('Number of lone pairs',atom.number_of_lone_pairs))


        return


    def add_protons(self, atom):
        # decide which method to use
        #print('PROTONATING',atom)
        if atom.steric_number in list(self.protonation_methods.keys()):
            self.protonation_methods[atom.steric_number](atom)
        else:
            print('Warning: Do not have a method for protonating',atom,'(steric number: %d)'%atom.steric_number)

        return

    
    def trigonal(self, atom):
        #print('TRIGONAL - %d bonded atoms'%(len(atom.bonded_atoms))) 
        rot_angle = math.radians(120.0)

        c = multi_vector(atom1 = atom)

        # 0 bonds
        if len(atom.bonded_atoms) == 0:
            pass
        
        # 1 bond
        if len(atom.bonded_atoms) == 1 and atom.number_of_protons_to_add > 0:
            # Add another atom with the right angle to the first one
            a = multi_vector(atom1 = atom, atom2 = atom.bonded_atoms[0])
            # use plane of bonded trigonal atom - e.g. arg
            if atom.bonded_atoms[0].steric_number == 3 and len(atom.bonded_atoms[0].bonded_atoms)>1:
                # use other atoms bonded to the neighbour to establish the plane, if possible
                other_atom_indices = []
                for i in range(len(atom.bonded_atoms[0].bonded_atoms)):
                    if atom.bonded_atoms[0].bonded_atoms[i] != atom:
                        other_atom_indices.append(i)

                if len(other_atom_indices)<2:
                    other_atom_indices = [0,1]

                axis = multi_vector(atom1 = atom.bonded_atoms[0], 
                              atom2 = atom.bonded_atoms[0].bonded_atoms[other_atom_indices[0]]
                              )**multi_vector(atom1 = atom.bonded_atoms[0], 
                                              atom2 = atom.bonded_atoms[0].bonded_atoms[other_atom_indices[1]])
            else:
                axis = a.orthogonal()

            a = rotate_multi_vector_around_an_axis(rot_angle, axis, a)
            a = self.set_bond_distance(a, atom.get_element())
            self.add_proton(atom, c+a)

        # 2 bonds
        if len(atom.bonded_atoms) == 2 and atom.number_of_protons_to_add > 0:
            # Add another atom with the right angle to the first two 
            a = multi_vector(atom1 = atom, atom2 = atom.bonded_atoms[1])
            b = multi_vector(atom1 = atom, atom2 = atom.bonded_atoms[0])
            axis = b**a
            new_a = rotate_multi_vector_around_an_axis(rot_angle, axis, a)
            new_a = self.set_bond_distance(new_a, atom.get_element())
            self.add_proton(atom, c+new_a)


        return


    def tetrahedral(self, atom):
        #print('TETRAHEDRAL - %d bonded atoms'%(len(atom.bonded_atoms))) 
        rot_angle = math.radians(109.5)

        # sanity check
        # if atom.number_of_protons_to_add + len(atom.bonded_atoms) != 4:
        # print 'Error: Attempting tetrahedral structure with %d bonds'%(atom.number_of_protons_to_add + 
        #                                                                len(atom.bonded_atoms))
        
        c = multi_vector(atom1 = atom)

        # 0 bonds
        if len(atom.bonded_atoms) == 0:
            pass
        
        # 1 bond
        if len(atom.bonded_atoms) == 1 and atom.number_of_protons_to_add > 0:
            # Add another atom with the right angle to the first one
            a = multi_vector(atom1 = atom, atom2 = atom.bonded_atoms[0])
            axis = a.orthogonal()
            a = rotate_multi_vector_around_an_axis(rot_angle, axis, a)
            a = self.set_bond_distance(a, atom.get_element())
            self.add_proton(atom, c+a)

        # 2 bonds
        if len(atom.bonded_atoms) == 2 and atom.number_of_protons_to_add > 0:
            # Add another atom with the right angle to the first two 
            a = multi_vector(atom1 = atom, atom2 = atom.bonded_atoms[1])
            axis = multi_vector(atom1 = atom.bonded_atoms[0],atom2 = atom)
            new_a = rotate_multi_vector_around_an_axis(math.radians(120), axis, a)
            new_a = self.set_bond_distance(new_a, atom.get_element())
            self.add_proton(atom, c+new_a)

        # 3 bonds
        if len(atom.bonded_atoms) == 3 and atom.number_of_protons_to_add > 0:
            a = multi_vector(atom1 = atom, atom2 = atom.bonded_atoms[2])
            axis = multi_vector(atom1 = atom.bonded_atoms[0],atom2 = atom)
            b = multi_vector(atom1 = atom, atom2 = atom.bonded_atoms[1])
            cross = b**axis
            angle = math.radians(120)
            if angle_degrees(cross.vectors[0],a.vectors[0]) < 90:
                angle = -angle
            new_a = rotate_multi_vector_around_an_axis(angle, axis, a)
            new_a = self.set_bond_distance(new_a, atom.get_element())
            self.add_proton(atom, c+new_a)

      
        return


    def add_proton(self, atom, position):
        residue = atom.residue
        #print(residue)
        # Create the new proton
        new_H = pdb.Atom()
        new_H.setProperty(numb    = None, 
                          name    = 'H', 
                          resName = atom.resName, 
                          chainID = atom.chainID,
                          resNumb = atom.resNumb,
                          x       = None,
                          y       = None,
                          z       = None,
                          occ     = None,
                          beta    = None,
                          element = 'H')

        #print(position)
        # set all the configurations
        for i in range(len(position.keys)):
            #print ('adding',position.keys[i],position.vectors[i])
            new_H.configurations[position.keys[i]] = [position.vectors[i].x, 
                                                      position.vectors[i].y, 
                                                      position.vectors[i].z]
        new_H.setConfiguration(position.keys[0])

        new_H.bonded_atoms = []
        new_H.charge = 0
        new_H.steric_number = 0
        new_H.number_of_lone_pairs = 0
        new_H.number_of_protons_to_add = 0
        new_H.number_of_pi_electrons_in_double_and_triple_bonds = 0
        
        residue.atoms.append(new_H)
        atom.bonded_atoms.append(new_H)
        atom.number_of_protons_to_add -=1
        #print('added',new_H, 'to',atom)

        return


    def set_bond_distance(self, a, element):
        d = 1.0
        if element in list(self.bond_lengths.keys()):
            d = self.bond_lengths[element]
        else:
            print('WARNING: Bond length for %s not found, using the standard value of %f'%(element, d))

        a = a.rescale(d)

        return a


if __name__ == '__main__':
    import protein, pdb, sys,os
    arguments = sys.argv
    if len(arguments) != 2:
        print('Usage: protonate.py <pdb_file>')
        sys.exit(0)

    filename = arguments[1]
    if not os.path.isfile(filename):
        print('Error: Could not find \"%s\"'%filename)
        sys.exit(1)

    
    p = Protonate()
    pdblist = pdb.readPDB(filename)
    my_protein = protein.Protein(pdblist,'test.pdb')
    
    p.remove_all_hydrogen_atoms_from_protein(my_protein)
    my_protein.writePDB('before_protonation.pdb')

    p.protonate_protein(my_protein)

    ## write out protonated file
    my_protein.writePDB('protonated.pdb')
