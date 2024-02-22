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
import joblib 

from .ani_descriptors import pdb_arrays,get_titratable,get_indices,get_desc_arrays
from .ase_io_proteindatabank_mod import read_proteindatabank



def calculate_pka(pdbfiles,writefile=None):

    
    # device to run the training
    device = torch.device('cpu')

    print("Please cite the original source of pKa-ANI if you use it:")
    print('''# Gokcan, H.; Isayev, O. Prediction of Protein pKa with 
          Representation Learning. Chemical Science, 2022, 13, 2462–2474. 
          https://doi.org/10.1039/d1sc05610g.''')


    print("Loading pKa-ANI Models and ANI-2x...")
    #FEATURES
    tyr_features=joblib.load(os.path.join(os.path.dirname(__file__),'FTYR.joblib'))
    asp_features=joblib.load(os.path.join(os.path.dirname(__file__),'FASP.joblib'))
    glu_features=joblib.load(os.path.join(os.path.dirname(__file__),'FGLU.joblib'))
    lys_features=joblib.load(os.path.join(os.path.dirname(__file__),'FLYS.joblib'))
    his_features=joblib.load(os.path.join(os.path.dirname(__file__),'FHIS.joblib'))

    #MODELS
    asp_model=joblib.load(os.path.join(os.path.dirname(__file__),'ASP_ani2x_FINAL_MODEL_F100.joblib'))
    glu_model=joblib.load(os.path.join(os.path.dirname(__file__),'GLU_ani2x_FINAL_MODEL_F75.joblib'))
    his_model=joblib.load(os.path.join(os.path.dirname(__file__),'HIS_ani2x_FINAL_MODEL_F100.joblib'))
    lys_model=joblib.load(os.path.join(os.path.dirname(__file__),'LYS_ani2x_FINAL_MODEL_F25.joblib'))
    tyr_model=joblib.load(os.path.join(os.path.dirname(__file__),'TYR_ani2x_FINAL_MODEL_F25.joblib'))
    
    #######################################################################        
    #call ani

    ani = torchani.models.ANI2x(periodic_table_index=True)
    print('Finished Loading.')
    
    pkaressize=0
    pkadict={} 
        
    for fpdb in pdbfiles:
        print('Calculating pKa for %s' % fpdb)    

        #pkadict[fpdb]=[]
        pkadict[fpdb]={}
        
        basename=fpdb.rsplit('.', 1)[0]
        infile=str(basename)+".pdb"
        
        if(writefile):
          flog=basename+"_pka.log"
          fo=open(flog,"w") 
          writer = csv.writer(fo,delimiter='\t')
          writer.writerow(['Residue', 'Chain', 'pKa'])
        
        
        #read pdb, and convert to ani type
        atoms=read_proteindatabank(infile, read_arrays=True)
        res,res_no,a_type,a_sym,a_no,pos,chainid,type_atm = pdb_arrays(atoms)

        chainid=np.array(chainid)

        a_no2=np.reshape(a_no, (1, len(a_no)))
        #sptensor=torch.tensor([a_no], device=device) #slow
        sptensor=torch.tensor(a_no2, device=device)
       
        pos=torch.tensor(pos,dtype=torch.float32)
        coords=torch.reshape(pos, (1, len(a_no),3))
        
        species_coordinates = ani.species_converter((sptensor, coords))
        aev = ani.aev_computer(species_coordinates)[1]
            
        #find titratable residues
        pkares,pkach=get_titratable(a_type,res,res_no,chainid)
        pka_res_chain = [(i, j) for i, j in zip(pkares, pkach)]

        #   divide atoms into groups by residue number
        #   this will be used to get activation and aev indices
        #   ca_list: the indices of CA atoms
        #   chlist : the chain IDs
        #   pdball_resi: atom indices per residue
        #   after having thesem we will only chose if the residue
        #   is titratable
        #
        #   i.e.
        #       ca_list  = [6, 7, 8] 
        #       pdball_resi = [array([0, 1, 2, 3, 4, 5, 6, 7]), 
        #                   array([ 8,  9, 10, 11]), 
        #                   array([12, 13, 14, 15, 16, 17, 18, 19, 20])]
        
        atom_list = []
        for i,a in enumerate(a_type):
            if((str(a)=='CA') and type_atm[i].strip()=='ATOM'):
                atom_list.append((res_no[i], chainid[i]))

      
    
        
        pdball_resi=[]
        for i,r in enumerate(atom_list):
            ilist=np.array(np.where((res_no == r[1]) & (chainid == r[1])))
            pdball_resi.append(ilist.flatten())
        
        nk,ck,ok=0,0,0 # counters for activation indices
      
        #GET ACTIVATION AND AEV INDICES
        #FOR ALL ATOMS
        all_acti=[]
        all_aevi=[]
    
        for k,index in enumerate(pdball_resi):
            activation_i,aev_i,nk,ck,ok=get_indices(index,a_type,a_no,nk,ck,ok)
            all_acti.append(activation_i)
            all_aevi.append(aev_i)               
    

        #now we are looping over residues
        #then if the residue is titratable 
        # we get aev, NN activation, and atom indices    
        for i,r in enumerate(atom_list):
            index = pdball_resi[i]
            lres=res[index[0]]
            lchid=str(r[1])

            if(r in pka_res_chain):               
               pkaressize=pkaressize+1
               res_aevi=all_aevi[i]           
               res_acti=all_acti[i]
               
               #lrnum=str(res_no[index[0]])
               lrnum=(chainid[index[0]].strip(), res_no[index[0]].strip())
    
               mychain=lchid
               if not lchid: mychain='A'
      
               a_symbols=[]
               for i in res_aevi: a_symbols.append(a_sym[i])
      
               
               ani_descriptors,features=get_desc_arrays(ani,species_coordinates,aev,res_acti,res_aevi,a_symbols,a_type)
    
    
               checklist=[]
               if(lres=='GLU'):
                  ani_descriptors_model=[] 
                  checklist=glu_features
                  model=glu_model
               if(lres=='ASP'):
                  ani_descriptors_model=[]
                  checklist=asp_features
                  model=asp_model
               if(lres=='LYS'):
                  ani_descriptors_model=[]
                  checklist=lys_features
                  model=lys_model
               if(lres=='HIS' or lres=='HID' or lres=='HIE'):
                  ani_descriptors_model=[]
                  checklist=his_features
                  model=his_model
               if(lres=='TYR'):
                  ani_descriptors_model=[]
                  checklist=tyr_features
                  model=tyr_model
    
    
               for i,fl in enumerate(features):
                   if(str(fl) in checklist):
                       ani_descriptors_model.append(ani_descriptors[i])
               
               ani_descriptors_model=np.array(ani_descriptors_model)
               ani_descriptors=ani_descriptors_model 
               
               X = np.reshape(ani_descriptors,(1, ani_descriptors.size))
      
               estimate_pka=model.predict(X)
               
               if(writefile): 
                 wres=lres+"-"+str(res_no[index[0]])
                 writer.writerow([wres,mychain,'{:2.2f}'.format(estimate_pka[0])])

               #resdict={}
               #resdict[lrnum]=lres,estimate_pka[0]
               #pkadict[fpdb].append(resdict)    
               pkadict[fpdb][lrnum]=lres,estimate_pka[0]
        
        if(writefile):
          fo.close()                  

    return pkadict
       
