#
# Various functions that cluttered pka.py
#
from src.pdb import HETATM, ATOM, ANISOU, SIGUIJ, SIGATM


#
# ----
#

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

#
# ----
#

def dump_protein_no_hydrogens(pdb_list, pdb_out):
    with open(pdb_out,'w') as fd:
        for record in pdb_list:
            if isinstance(record, HETATM):
                continue
            if isinstance(record, (ATOM, ANISOU, SIGUIJ, SIGATM)):
                if record.element == 'H':
                    continue
            fd.write(str(record))
            fd.write('\n')

def remove_hydrogens(pdb_in, pdb_out):
    """Remove hydrogens from the PDB file"""
    with open(pdb_in,'r') as fd:
        l_lines_i = fd.readlines()

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

    with open(pdb_out,'w') as fd:
        fd.writelines(l_lines_o)

    return

#
# ---
#

def is_sameatom(atom1,atom2):
    """Are atom1 and atom2 the same atom?"""
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
# -----
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
