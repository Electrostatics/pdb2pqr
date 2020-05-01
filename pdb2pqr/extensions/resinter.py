"""Resinter extension

Print interaction energy between each residue pair in the protein. 

Authors:  Kyle Monson and Emilie Hogan
"""
import logging
from itertools import product, permutations, count
from collections import defaultdict
from ..hydrogens import Optimize
from ..hydrogens import HydrogenRoutines


_LOGGER = logging.getLogger(__name__)


_titrationSets = (('ARG','AR0'),
                  ('ASP', 'ASH'),
                  ('CYS', 'CYX'),
                  ('GLU', 'GLH'),
                  ('HSD', 'HSE', 'HSP'),
                  ('HID', 'HIE', 'HIP'),
                  ('LYN', 'LYS'),
                  ('TYR', 'TYM'),
                  ('NEUTRAL-CTERM', 'CTERM'),
                  ('NEUTRAL-NTERM', 'NTERM'))


_titrationSetsMap = defaultdict(tuple)
for tsSet in _titrationSets:
    for ts in tsSet:
        _titrationSetsMap[ts] = tsSet
#loose ends.
_titrationSetsMap['HIS'] = _titrationSetsMap['HSD']
_titrationSetsMap['CYM'] = _titrationSetsMap['CYS']


_pairEnergyResults = {}
#If the residue pair energy for a specific pair changes less than this ignore it.
PAIR_ENERGY_EPSILON = 1.0e-14


def addExtensionOptions(extensionGroup):
    """
        Add options.
    """
    extensionGroup.add_option('--residue_combinations', 
                              dest='residue_combinations', 
                              action='store_true', 
                              default=False,
                              help=
'''Remap residues to different titration states and rerun resinter appending output.
Consider only the minimum number of whole protein titration combinations needed to
test each possible pairing of residue titration states. Normally used with
--noopt. If a protein titration state combination results in a pair of residue being 
re-tested in the same individual titration states a warning will be generated if the 
re-tested result is different. This warning should not be possible if used with --noopt.''') 
    
    extensionGroup.add_option('--all_residue_combinations', 
                              dest='all_residue_combinations', 
                              action='store_true', 
                              default=False,
                              help=
'''Remap residues to ALL possible titration state combinations and rerun resinter appending output.
Results with --noopt should be the same as --residue_combinations. Runs considerably slower than
--residue_combinations and generates the same type of warnings. 
Use without --noopt to discover how hydrogen optimization affects residue 
interaction energies via the warnings in the output.''') 

def usage():
    """
    Returns usage text for resinter.
    """
    txt = 'Print interaction energy between each residue pair in the protein to {output-path}.resinter.'
    return txt

def _combinations(sublist, remainder):
    """
    combinations function helper
    
    sublist - first list in list of lists for combinations
    remainder - list of remaining lists
    """
    if remainder and sublist:
        for item in sublist:
            for result in _combinations(remainder[0], remainder[1:]):
                yield [item] + result
        return
    elif sublist:
        for item in sublist:
            yield [item]
        return
    elif remainder:
        for result in _combinations(remainder[0], remainder[1:]):
            yield [None] + result
        return
            
    yield [None]

def combinations(initialList):
    """
    Wrapper for main combinations function. 
    
    Iterates over each possible combination of single items 
    from each sub-list. For example:
    
    combinations([[1,2],[3,4]] -> [1,3], [1,4], [2,3], [2,4] in that order.     
    
    
    initialList - list of lists to derive combinations from.
    """
    if not initialList:
        return
    for result in _combinations(initialList[0], initialList[1:]):
        yield result


def pairwiseCombinations(initialList):
    """
    Creates the minimum set of combinations that will make available 
    every possible pair available.
    """
    r = [len(x) for x in initialList]
    m = min(r)
    M = max(r)
    n = len(initialList)
    
    R = set()
    
    for i in range(m):
        t = [initialList[x][i] for x in range(n)]
        R.add(tuple(t))
        
    for i in range(m,M):
        t = [initialList[x][min(r[x]-1,i)] for x in range(n)]
        R.add(tuple(t))
        
    for i in range(m):
        for j in range(n):
            for k in range(i+1,r[j]):
                prejth = [initialList[x][i] for x in range(j)]

                jth = initialList[j][k]
                
                postjth = [initialList[x][i] for x in range(j+1, n)]
                
                t = prejth + [jth] + postjth
                R.add(tuple(t))
    
    for i in range(m,M):
        for j in range(n):
            if r[i] < i:
                continue
            for k in range(i+1,r[j]):
                prejth = [initialList[x][min(r[x]-1,i)] for x in range(j)]

                jth = initialList[j][k]
                
                postjth = [initialList[x][min(r[x]-1,i)] for x in range(j+1, n)]
                
                t = prejth + [jth] + postjth
                R.add(tuple(t))
                
    return R
    
def get_residue_interaction_energy(residue1, residue2):
    """
    Returns to total energy of every atom pair between the two residues.
    
    Uses Optimize.get_pair_energy and it's donor/accepter model 
    to determine energy.
    
    residue1 - "donor" residue
    residue2 - "acceptor" residue
    
    THE RESULTS OF THIS FUNCTION ARE NOT SYMMETRIC. Swapping 
    residue1 and residue2 will not always produce the same result.
    """
    energy = 0.0
    for pair in product(residue1.get_atoms(), residue2.get_atoms()):
        energy += Optimize.get_pair_energy(pair[0], pair[1])
        
    return energy

def save_residue_interaction_energies(residues, output):
    """
    Writes out the residue interaction energy for each possible
    residue pair in the protein.
    """
    residuepairs = permutations(residues, 2)
    
    for pair in residuepairs:
        energy = get_residue_interaction_energy(pair[0], pair[1])
        pairText = str(pair[0]) + ' ' + str(pair[1])
        if pairText in _pairEnergyResults:
            
            oldEnergy = _pairEnergyResults[pairText]
            energyDiff = oldEnergy - energy
            if abs(energyDiff) > PAIR_ENERGY_EPSILON:
                txt = '#%s re-tested' % pairText
                txt += ' with a difference of %s' % repr(energyDiff)
                if (energy != 0):
                    txt += ' and a reference of %s' % repr(energyDiff/energy)
                else:
                    txt += ' and the previous energy was 0'
                txt += '\n'

                output.write(txt)
                 
            continue
        
        _pairEnergyResults[pairText] = energy
                    
        
def get_residue_titration_sets(residues):
    """
    Returns all possible titration states for each residue as a list of lists.
    """
    result = []
    for residue in residues:
        result.append(_titrationSetsMap.get(residue.name, (residue.name,)))
        
    return result

def process_residue_set(residueSet, routines, output, clean = False,
                                                      neutraln = False,
                                                      neutralc = False,
                                                      ligand = None,
                                                      assign_only = False,
                                                      chain = False,
                                                      debump = True,
                                                      opt = True):
    _LOGGER.debug(str(residueSet))
    
    routines.removeHydrogens()
    
    for newResidueName, oldResidue, index in zip(residueSet, routines.protein.get_residues(), count()):
        if newResidueName is None:
            continue
        
        chain = routines.protein.chainmap[oldResidue.chain_id]
        chainIndex = chain.residues.index(oldResidue)
        residueAtoms = oldResidue.atoms
        
        #Create the replacement residue
        newResidue = routines.protein.create_residue(residueAtoms, newResidueName)
        
        #Make sure our names are cleaned up for output.
        newResidue.renameResidue(newResidueName)
        
        #Drop it in
        routines.protein.residues[index] = newResidue
        chain.residues[chainIndex] = newResidue
    
    #Run the meaty bits of PDB2PQR  
    routines.setTermini(neutraln, neutralc)
    routines.updateBonds()
    
    if not clean and not assign_only:
        routines.updateSSbridges()

        if debump:
            routines.debumpProtein()
            
        routines.addHydrogens()

        hydRoutines = HydrogenRoutines(routines)

        if debump:
            routines.debumpProtein()  

        if opt:
            hydRoutines.set_optimizeable_hydrogens()
            hydRoutines.initialize_full_optimization()
            hydRoutines.optimize_hydrogens()
        else:
            hydRoutines.initialize_wat_optimization()
            hydRoutines.optimize_hydrogens()

        # Special for GLH/ASH, since both conformations were added
        hydRoutines.cleanup()
        
    save_residue_interaction_energies(routines.protein.get_residues(), output)        
        

def write_all_residue_interaction_energies_combinations(routines, output, options, all_residue_combinations=False):
    """
    For every titration state combination of residue output the 
    interaction energy for all possible residue pairs. 
    """
    residueNamesList = get_residue_titration_sets(routines.protein.get_residues())
    
    _LOGGER.debug("Testing the following combinations")
    namelist = [r.name for r in routines.protein.get_residues()]
    combinationsData = zip(namelist, residueNamesList)
    for thing in combinationsData:
        _LOGGER.debug(str(thing))
        
    if all_residue_combinations:
        combinationGenerator = combinations(residueNamesList)
    else:
        combinationGenerator = pairwiseCombinations(residueNamesList)
        
    count = 0
    for residueSet in combinationGenerator:
        count += 1
        process_residue_set(residueSet, routines, output, 
                            clean = options.clean,
                            neutraln = options.neutraln,
                            neutralc = options.neutralc,
                            ligand = options.ligand,
                            assign_only = options.assign_only,
                            chain = options.chain,
                            debump = options.debump,
                            opt = options.opt)
        
    for resultKey in sorted(_pairEnergyResults.keys()):
        output.write(resultKey + ' ' + str(_pairEnergyResults[resultKey]) + '\n')
    
    _LOGGER.debug(str(count) + ' residue combinations tried')

def create_resinter_output(routines, outfile, options, 
                           residue_combinations=False,
                           all_residue_combinations=False):
    """
    Output the interaction energy between each possible residue pair.
    """
    _LOGGER.debug("Printing residue interaction energies...")
    
    output = extensions.extOutputHelper(routines, outfile)
    
    if residue_combinations or all_residue_combinations:
        write_all_residue_interaction_energies_combinations(routines, output, options, 
                                                            all_residue_combinations=all_residue_combinations)
    else:
        save_residue_interaction_energies(routines.protein.get_residues(), output)
    

def run_extension(routines, outroot, options):
    outname = outroot + ".resinter"
    with open(outname, "w") as outfile:
        create_resinter_output(routines, outfile, options, 
                               residue_combinations=options.residue_combinations,
                               all_residue_combinations=options.all_residue_combinations)