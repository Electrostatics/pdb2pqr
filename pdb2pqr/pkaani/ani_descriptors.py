# CODE REUSED FROM https://github.com/isayevlab/pKa-ANI
# PLEASE CREDIT THE Isayev Lab if you use pKa-ANI to assign titration states!
# 
# CITATION:
# Gokcan, H.; Isayev, O. Prediction of Protein pKa with Representation Learning. 
# Chemical Science, 2022, 13, 2462–2474. https://doi.org/10.1039/d1sc05610g.

import torch
import torchani
import os
import math
import sys
import numpy as np
import csv
from ase.io import read

device = torch.device('cpu')

###########################################
def pdb_arrays(atoms):
    '''
      READING ARRYAS IN PDB FILE
      USING ASE LIB
    '''
    a=atoms.arrays['residuenames']
    b=atoms.arrays['residuenumbers']
    c=atoms.arrays['atomtypes']             
    d=atoms.get_chemical_symbols()
    e=atoms.get_atomic_numbers()
    e=np.array(e)
    f=atoms.get_positions()
    g=atoms.arrays['chainid']
    h=atoms.arrays['type_atm']
    return a,b,c,d,e,f,g,h
###########################################

def get_titratable(atypes,resnames,resnumbers,chains):

    '''
      FIND RESIDUE NUMBERS OF 
      TITRATABLE RESIDUES WITHIN PDB FILE
    '''

    titratable=('GLU','ASP','HIS','HIE','HID','LYS','TYR')
    restitrate=[]
    chinfo=[]
    for i,a in enumerate(atypes):
        if((str(a)=='CA') and (str(resnames[i]).strip(" ") in titratable)):
            restitrate.append(resnumbers[i])
            chinfo.append(chains[i])

    return restitrate,chinfo

###########################################

'''
    	FUNCTIONS TO FIND 

              * AEV INDICES 
              * NN ACTIVATION INDICES

        TO GET SAME ORDER FOR INDICES
        WE CREATE ARRAY OF THESE INDICES 
        BY THE FOLLOWING ORDER BACKBONE 
           1. BACKBONE N
           2. BACKBONE CA
           3. SIDE CHAIN C ATOMS
           4. IF EXISTS SIDE CHAIN N ATOMS
           5. IF ESISTS SIDE CHAIN O ATOMS
           6. BACKBONE C ATOM
           7. BACKBONE O ATOM


        EXAMPLE: PDB FILE CONTAINING ASP6,GLY7,GLU27
           atom type   NN activation     AEV     Residue    Residue
                          index         index      name        no
             N              0             0        ASP         6
             CA             0             1        ASP         6
             CB             1             2        ASP         6
             CG             2             3        ASP         6
             OD1            0             4        ASP         6
             OD2            1             5        ASP         6
             C              3             6        ASP         6
             O              2             7        ASP         6
             N              1             8        GLY         7
             CA             4             9        GLY         7
             C              5            10        GLY         7
             O              3            11        GLY         7
             N              2            12        GLU        27
             CA             6            13        GLU        27
             CB             7            14        GLU        27
             CG             8            15        GLU        27
             CD             9            16        GLU        27
             OE1            4            17        GLU        27
             OE2            5            18        GLU        27
             C              10           19        GLU        27
             O              6            20        GLU        27
           
           
           ASP: activation indices = 0 0 1 2 0 1 3 2
           GLY: activation indices = 1 4 5 3
           GLU: activation indices = 2 6 7 8 9 4 5 10 6
  
           AEV indices
           ASP : aev indices = 0, 1, 2, 3, 4, 5, 6, 7
           GLU : aev indices = 12, 13, 14, 15, 16, 17, 18, 19, 20


           LETs CHANGE THE ORDER IN PDB FILE 
                FROM   N,CA,SIDE CHAIN,C,O
                TO     N,CA,C,O,SIDE CHAIN (now C index is 2, O index is 3)

                ORDER MUST BE: N,C,CB,CD,OD1,OD2,C,O

                ASP activation indices = 0, 0, 1, 2, 0, 1, 3, 2 
                ASP AEV indices        = 0, 1, 4, 5, 6, 7, 2, 3 
     
                GLU activation indices = 2, 6, 7, 8, 9, 4, 5, 10, 6 
                GLU AEV indices        = 12, 13, 16, 17, 18, 19, 20, 14, 15

                NN ACTIVATIONS DOESN SEEM TO BE EFFECTED, 
                BUT THEY ARE CREATED IN THE SAME FASHIN AS IN AEV
                AND DEPEND ON THE ATOMS CHEMICAL SYMBOL
                  
''' 


def get_SC_index(r,ano,at,atype,k=0):

    #SIDE CHAIN ATOMS
    bb_a=('N','CA','C','O','OXT')
    scai=[]
    scri=[]    
    for i in r:
        if ((ano[i]==at) and (atype[i] not in bb_a)): 
           scai.append(k)
           scri.append(i)
           k=k+1

    return scai,scri,k

def get_BB_index(r,atype,at,k=0):

    #BACKBONE ATOMS
    for i in r:
        if (atype[i]==at): ai=k; ri=i;k=k+1
    
    return ai,ri,k

def get_indices(r,atype,ano,nk,ck,ok):
    #to make sure we have the correct order for every pdb
    #not sure about the side chain
    #it can be modified easily thoug, just like BB atoms
    #but it will requre to write for every
    #titratable residue???
    ai=[]
    ri=[]

    #BB N atom
    nai,nri,nk=get_BB_index(r,atype,'N',nk)
    ai.append(nai)
    ri.append(nri)

    #BB CA atom
    cai,cri,ck=get_BB_index(r,atype,'CA',ck)
    ai.append(cai)
    ri.append(cri)
    #SC C atoms
    cai,cri,ck=get_SC_index(r,ano,6,atype,ck)
    ai.extend(cai)
    ri.extend(cri)
    #SC N atoms
    nai,nri,nk=get_SC_index(r,ano,7,atype,nk)
    if(len(nai)!=0): ai.extend(nai)
    if(len(nri)!=0): ri.extend(nri)
    #SC O atoms
    oai,ori,ok=get_SC_index(r,ano,8,atype,ok)
    ai.extend(oai)
    ri.extend(ori)
    if(len(oai)!=0): ai.extend(nai)
    if(len(ori)!=0): ri.extend(nri)
    #BB C atom
    cai,cri,ck=get_BB_index(r,atype,'C',ck)
    ai.append(cai)
    ri.append(cri)
    #BB O atom
    oai,ori,ok=get_BB_index(r,atype,'O',ok)
    ai.append(oai)
    ri.append(ori)

    return ai,ri,nk,ck,ok


########################################################
def get_desc_arrays(ani,species_coordinates,aev,actindex,aev_indx,elements,atype):
    
    """
       GET DESCRIPTORS (AEV and NN ACTIVATIONS)
    
        AEV shape and index
           i.e. system with 275 atoms
                since we extracted the pdb substructure
                  aev shape --> torch.Size([1, 275, 1008])
                then atom index is on ri
           
           for each atom (ri) in resindex aev is
               np.array(aev[0,ri,:])   
           append each aev to aev_array
    
    
            for each atom (ri) in resindex

                flatten activation_1 (act1) and
                        activation_2 (act2) arrays

                extend  nn_desc1 and nn_desc2 by 
                        activation_1 (act1) and
                        activation_2 (act2) respectively
   
 
        concatenate aev array, nn_desc1, nn_desc2

        EXAMPLE:
                ASP
                AEV  :  [0, 1, 4, 5, 6, 7, 2, 3]
                ELEM :  ['N', 'C', 'C', 'C', 'O', 'O', 'C', 'O']
                NN   :  [0, 0, 1, 2, 0, 1, 3, 2]

                activation  element    nn_act.shape  nn_act2.shape
                index
                    0           N     torch.Size([3, 128]) torch.Size([3, 160])                 
                    0           C     torch.Size([11, 160]) torch.Size([11, 192])
                    1           C     torch.Size([11, 160]) torch.Size([11, 192])
                    2           C     torch.Size([11, 160]) torch.Size([11, 192])
                    0           O     torch.Size([7, 128]) torch.Size([7, 160])
                    1           O     torch.Size([7, 128]) torch.Size([7, 160])
                    3           C     torch.Size([11, 160]) torch.Size([11, 192])
                    2           O     torch.Size([7, 128]) torch.Size([7, 160])

                aev.shape         : (8064,)
                nn_desc1.shape    : (1152,)
                nn_desc2.shape    : (1408,)
                descriptors.shape : (10624,)

                -----------------------------

                GLU
                AEV  :  [12, 13, 16, 17, 18, 19, 20, 14, 15]
                ELEM :  ['N', 'C', 'C', 'C', 'C', 'O', 'O', 'C', 'O']
                NN   :  [0, 0, 1, 2, 3, 0, 1, 4, 2]

                activation  element    nn_act.shape  nn_act2.shape
                index
                    0           N     torch.Size([3, 128]) torch.Size([3, 160])
                    0           C     torch.Size([11, 160]) torch.Size([11, 192])
                    1           C     torch.Size([11, 160]) torch.Size([11, 192])
                    2           C     torch.Size([11, 160]) torch.Size([11, 192])
                    3           C     torch.Size([11, 160]) torch.Size([11, 192])
                    0           O     torch.Size([7, 128]) torch.Size([7, 160])
                    1           O     torch.Size([7, 128]) torch.Size([7, 160])
                    4           C     torch.Size([11, 160]) torch.Size([11, 192])
                    2           O     torch.Size([7, 128]) torch.Size([7, 160])

                aev.shape         : (9072,)
                nn_desc1.shape    : (1312,)
                nn_desc2.shape    : (1600,)
                descriptors.shape : (11984,)

                -----------------------------

    
    """
    # for each species type
    atom_nn=['H','C','N','O','S','F','Cl']


    aev_arr=[]
    nn_desc1=[]
    nn_desc2=[]

    #V0 style
    nn_activations=[]

    aevf=[]
    resatoms=[]
    nn1f=[]
    nn2f=[] 

    #AEV ARRAY
    for ri in aev_indx: #(resindex):
        resatoms.append(atype[ri])
        aa=aev[0,ri,:]
        for j in range(len(aa)):
            aevstr='AEV_'+str(atype[ri])+'_'+str(j)
            aevf.append(aevstr)
  
        aa=torch.flatten(aa).detach().numpy() 
        aev_arr.extend(aa)
       
    aev_arr=np.array(aev_arr).flatten()

    for j,ai in enumerate(actindex): 
        for i, s in enumerate(atom_nn):
            if(s==elements[j]):

               #constuct nn - all layers of first ensemble member, 
               #nn for atoms of specific type
               nn = ani.neural_networks[0][s]
               # apply nn to all atoms of specified type
               nn_act = nn[:-1](aev[species_coordinates[0] == i])
               #if(s=='O'):
               #  print(nn(aev[species_coordinates[0] == i]))
               #  print('---------------------------------------') 

               nn_act2= nn[:-3](aev[species_coordinates[0] == i])
           
               act1=nn_act[ai,:]
               act1=torch.flatten(act1).detach().numpy()

               nn1f_local=[]
               for k in range(len(act1)):
                   nnstr='NN3'+'_'+str(resatoms[j])+'_'+str(k)
                   nn1f_local.append(nnstr)
               
               act2=nn_act2[ai,:]
               act2=torch.flatten(act2).detach().numpy()

               nn2f_local=[]
               for k in range(len(act2)):
                   nnstr='NN2'+'_'+str(resatoms[j])+'_'+str(k)
                   nn2f_local.append(nnstr)
     

               #V0 style 
               #act_s=np.concatenate((act1,act2), axis=0)
               #nn_activations.extend(act_s)
          
               nn_desc1.extend(act1)       
               nn_desc2.extend(act2)
               nn1f.extend(nn1f_local)
               nn2f.extend(nn2f_local)
    
               
    nn_desc1=np.array(nn_desc1).flatten()
    nn_desc2=np.array(nn_desc2).flatten()
    descriptors=np.concatenate((aev_arr,nn_desc1,nn_desc2),axis=0)
    descriptors=np.asarray(descriptors)

    aevf=np.array(aevf)
    nn1f=np.array(nn1f)
    nn2f=np.array(nn2f)
    featurelist=np.concatenate((aevf,nn1f,nn2f),axis=0)
   

    #V0 style
    #nn_activations=np.array(nn_activations)
    #nn_activations=nn_activations.flatten()
    #descriptors=np.concatenate((aev_arr,nn_activations),axis=0)
    #print('AEV,NN1(nn[:-1]), NN2(nn[:-3]),desc,feat',aev_arr.shape,nn_desc1.shape,nn_desc2.shape,descriptors.shape,featurelist.shape)

    return descriptors,featurelist 
#################################################



