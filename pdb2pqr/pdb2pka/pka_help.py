#
# Various functions that cluttered pka.py
#
import logging
from src.aa import WAT
from src.pdb import HETATM, ATOM, ANISOU, SIGUIJ, SIGATM
import pKaTool.pKa_calc
import pMC_mult

_LOGGER = logging.getLogger(__name__)


def titrate_one_group(name, intpkas, is_charged, acidbase):
    """Titrate a single group and return the pKa value for it"""
    names = [name]
    num_states = len(intpkas)
    state_counter = [num_states]
    linear = [] # The linear matrix
    for _ in range(1):
        for _ in range(num_states):
            for _ in range(num_states):
                linear.append(0.0)
    #
    # Set the MC parameters
    #
    mcsteps = 5000
    phstart = 0.1
    phend = 20.0
    phstep = 0.1
    #
    # Call our little C++ module
    #
    fast = pMC_mult.MC(intpkas, linear, acidbase, state_counter, is_charged)
    fast.set_MCsteps(int(mcsteps))
    _LOGGER.info('Calculating intrinsic pKa value')
    pka_vals = fast.calc_pKas(phstart, phend, phstep)
    count = 0
    intpka = pka_vals[0]
    _LOGGER.info('Simulated intrinsic pKa value: %5.2f', intpka)
    count = 1
    #
    # Get the charges
    #
    charges = {}
    # TODO: 2020/07/06 intendo - ph_start, ph_step, and ph are never used?
    #ph_start = pka_vals[count]
    #ph_step = pka_vals[count + 1]
    num_phs = pka_vals[count + 2]
    count += 2
    phs = []
    charges = []
    #ph = ph_start
    for _ in range(int(num_phs)):
        count += 1
        phs.append(pka_vals[count])
        count += 1
        charges.append(pka_vals[count])
        #ph += ph_step
    if pka_vals[count + 1] == 999.0 and pka_vals[count + 2] == -999.0:
        count += 2
    else:
        _LOGGER.error('Something is wrong')
        _LOGGER.error(pka_vals[count:count+30])
        raise Exception('Incorrect data format from pMC_mult')
    return intpka

#
# ----
#

def dump_protein_no_hydrogens(pdb_list, pdb_out):
    with open(pdb_out, 'w') as wfd:
        for record in pdb_list:
            if isinstance(record, HETATM):
                ##check if the record is not a water in which case we will print a warning
                if not record.res_name in WAT.water_residue_names:
                    _LOGGER.warning("Warning!: HETATM record %s %s that is not a water is being dropped\n  ", record.res_name, record.element)
                    ##raw_input("Press enter to continue...")
                continue
            if isinstance(record, (ATOM, ANISOU, SIGUIJ, SIGATM)):
                if record.element == 'H':
                    continue
            wfd.write(str(record))
            wfd.write('\n')

def remove_hydrogens(pdb_in, pdb_out):
    """Remove hydrogens from the PDB file"""
    with open(pdb_in, 'r') as rfd:
        l_lines_i = rfd.readlines()

    l_lines_o = []
    for s_line in l_lines_i:
        record = s_line[:6].strip()
        if record in ['HETATM']:
            continue
        if record in ['ATOM', 'ANISOU', 'SIGUIJ', 'SIGATM', ]:
            element = s_line[76:78].strip()
            if element == 'H':
                continue
        l_lines_o += [s_line]

    with open(pdb_out, 'w') as wfd:
        wfd.writelines(l_lines_o)

#
# ---
#

def is_sameatom(atom1, atom2):
    """Are atom1 and atom2 the same atom?"""
    #
    # Compare atom1 and atom2
    #
    properties = ['name', 'res_seq', 'chain_id']
    for attr in properties:
        a1_prop = getattr(atom1, attr, None)
        a2_prop = getattr(atom2, attr, None)
        if (attr != 'chain_id' and (not a1_prop or not a2_prop)) or a1_prop != a2_prop:
            return None
    return 1


#
# -----
#

def test_interface():
    """Test the interface with pKaTool"""
    value = pKaTool.pKa_calc.Monte_Carlo_Mult_CPP()
    value.intrinsic_pKa = {':0001:ASP':[0.0, 4.0, 5.0]}
    value.charged_state = {':0001:ASP':[0, 1, 1]}
    value.acid_base = {':0001:ASP':-1}
    value.intene_mult = {':0001:ASP':{':0001:ASP':[[0, 0, 0], [0, 0, 0], [0, 0, 0]]}}
    value._calc_pKas(0.0, 10.0, 0.5)
