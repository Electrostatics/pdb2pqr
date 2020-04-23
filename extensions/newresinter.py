"""
    Resinter extension

    Print interaction energy between each residue pair in the protein. 
"""

__date__ = "21 October 2011"
__authors__ = "Kyle Monson and Emile Hogan"

import extensions
from src.hydrogens import Optimize
#itertools FTW!
from itertools import product, permutations, count
from src.hydrogens import hydrogenRoutines

#Here are the Ri -> [Ri0, Ri1] maps:
#ARG -> [{AR0}, {ARG}]
#ASP -> [{ASP}, {ASH}]
#GLU -> [{GLU}, {GLH}]
#CYS -> [{CYM}, {CYS}]
#LYS -> [{LYN}, {LYS}]
#TYR -> [{TYM}, {TYR}]
#CTERM -> [{CTERM}, {NEUTRAL-CTERM}]
#NTERM -> [{NEUTRAL-NTERM}, {NTERM}]
#HIS -> [{HSD, HSE}, {HIS}]
#HIP -> [{HID, HIE}, {HIP}]
#HSP -> [{HSD, HSE}, {HSP}]

_titrationSets = ((('AR0',), 'ARG'),
                  (('ASH',), 'ASP'),
                  (('CYX',), 'CYS'),
                  (('GLU',), 'GLH'),
                  (('HSD', 'HSE'), 'HSP'),
                  (('HID', 'HIE'), 'HIP'),
                  (('LYN',), 'LYS'),
                  (('TYM',), 'TYR'),
                  (('CTERM',), 'NEUTRAL-CTERM'),
                  (('NEUTRAL-NTERM',), 'NTERM'))

_titrationSetsMap = {}

for tsSet in _titrationSets:
    for ts in tsSet[0]:
        _titrationSetsMap[ts] = tsSet
        
    _titrationSetsMap[tsSet[1]] = tsSet
    
#loose ends.
_titrationSetsMap['HIS'] = _titrationSetsMap['HSD']
_titrationSetsMap['CYM'] = _titrationSetsMap['CYS']
        
def usage():
    """
    Returns usage text for newresinter.
    """
    txt = 'Print interaction energy between each residue pair in the protein to {output-path}.newresinter.'
    return txt

def run_extension(routines, outroot, options):
    outname = outroot + ".newresinter"
    with open(outname, "w") as outfile:
        processor = ResInter(routines, outfile, options)
        processor.generate_all()
        processor.write_resinter_output()

class ResInter(object):
    def __init__(self, routines, outfile, options):
        self.pairEnergyResults = {}
        self.combinationCount = 0
        self.totalCombinations = 0
        self.options = options
        self.output = extensions.extOutputHelper(routines, outfile)
        self.routines = routines

    def save_interation_energy(self, first, second):
        energy = get_residue_interaction_energy(first, second)
        pairText = str(first) + ' ' + str(second)
        if pairText in self.pairEnergyResults:            
            txt = '#%s re-tested!!! LOLWAT?\n' % pairText
            self.output.write(txt)             
        else:
            self.pairEnergyResults[pairText] = energy
                    

    def save_all_residue_interaction_energies(self):
        """
        Writes out the residue interaction energy for each possible
        residue pair in the protein.
        """
        residuepairs = permutations(self.routines.protein.getResidues(), 2)
        
        for pair in residuepairs:
            self.save_interation_energy(pair[0], pair[1])
            
    def save_one_with_all_interaction_energies(self, i):
        """
        Writes out the residue interaction energy for each possible
        residue pair in the protein.
        """
        residues = list(self.routines.protein.getResidues())
        target = residues[i]
        del residues[i]        
        
        for residue in residues:
            self.save_interation_energy(target, residue)
            self.save_interation_energy(residue, target)
            
    def save_pair_interaction_energies(self, i, j):
        """
        Writes out the residue interaction energy for each possible
        residue pair in the protein.
        """
        residues = list(self.routines.protein.getResidues())

        self.save_interation_energy(residues[i], residues[j])
        self.save_interation_energy(residues[j], residues[i])
            
    def create_all_protonated(self):
        residueSet = get_residue_titration_set_protonated(self.routines.protein.getResidues())
        self.process_residue_set(residueSet, 
                                clean = self.options.clean,
                                neutraln = self.options.neutraln,
                                neutralc = self.options.neutralc,
                                ligand = self.options.ligand,
                                assign_only = self.options.assign_only,
                                chain = self.options.chain,
                                debump = self.options.debump,
                                opt = self.options.opt)
        
        self.save_all_residue_interaction_energies()
        
    def create_all_single_unprotonated(self):
        combinations = residue_set_single_unprotonated_combinations(self.routines.protein.getResidues())
        for residueSet, i in combinations:
            self.process_residue_set(residueSet, 
                                     clean = self.options.clean,
                                     neutraln = self.options.neutraln,
                                     neutralc = self.options.neutralc,
                                     ligand = self.options.ligand,
                                     assign_only = self.options.assign_only,
                                     chain = self.options.chain,
                                     debump = self.options.debump,
                                     opt = self.options.opt)
        
            self.save_one_with_all_interaction_energies(i)
            
    def create_all_pair_unprotonated(self):
        combinations = residue_set_pair_unprotonated_combinations(self.routines.protein.getResidues())
        for residueSet, i, j in combinations:
            self.process_residue_set(residueSet, 
                                     clean = self.options.clean,
                                     neutraln = self.options.neutraln,
                                     neutralc = self.options.neutralc,
                                     ligand = self.options.ligand,
                                     assign_only = self.options.assign_only,
                                     chain = self.options.chain,
                                     debump = self.options.debump,
                                     opt = self.options.opt)
        
            self.save_pair_interaction_energies(i, j)
            
    def count_combinations(self):
        n = 0 # total iterable residues
        k = 0 # total iterable residues with two possible choices.
        
        allProtonated = get_residue_titration_set_protonated(self.routines.protein.getResidues())
        
        for name in allProtonated:
            if name in _titrationSetsMap:
                n += 1
                
                if len(_titrationSetsMap[name][0]) == 2:
                    k += 1
        
        self.totalCombinations = (((n+k)**2)+(n-k)+2)/2
            
    def generate_all(self):
        """
        For every titration state combination of residue output the 
        interaction energy for all possible residue pairs. 
        """
        self.routines.write("Printing residue interaction energies...\n")
        
        self.count_combinations()
        
        #Phase 1: Everything protonated
        self.create_all_protonated()
        
        #Phase 2: Single unprotonated paired with everything else.
        self.create_all_single_unprotonated()
        
        #Phase 2: Pair unprotonated paired with each other.
        self.create_all_pair_unprotonated()

    def write_resinter_output(self):
        """
        Output the interaction energy between each possible residue pair.
        """
        for resultKey in sorted(self.pairEnergyResults.keys()):
            self.output.write(resultKey + ' ' + str(self.pairEnergyResults[resultKey]) + '\n')
        
        self.routines.write(str(self.combinationCount)+' residue combinations tried\n')
                
    def process_residue_set(self, residueSet,  
                            clean = False,
                            neutraln = False,
                            neutralc = False,
                            ligand = None,
                            assign_only = False,
                            chain = False,
                            debump = True,
                            opt = True):
        
        self.combinationCount += 1
        
        txt = "Running combination {0} of {1}\n".format(self.combinationCount, self.totalCombinations)
        self.routines.write(txt)
        
        self.routines.write(str(residueSet)+'\n')
        
        self.routines.removeHydrogens()
        
        for newResidueName, oldResidue, index in zip(residueSet, self.routines.protein.getResidues(), count()):
            if newResidueName is None:
                continue
            
            chain = self.routines.protein.chainmap[oldResidue.chainID]
            chainIndex = chain.residues.index(oldResidue)
            residueAtoms = oldResidue.atoms
            
            #Create the replacement residue
            newResidue = self.routines.protein.createResidue(residueAtoms, newResidueName)
            
            #Make sure our names are cleaned up for output.
            newResidue.renameResidue(newResidueName)
            
            #Drop it in
            self.routines.protein.residues[index] = newResidue
            chain.residues[chainIndex] = newResidue
        
        #Run the meaty bits of PDB2PQR  
        self.routines.setTermini(neutraln, neutralc)
        self.routines.updateBonds()
        
        if not clean and not assign_only:
            self.routines.updateSSbridges()
    
            if debump:
                self.routines.debumpProtein()
                
            self.routines.addHydrogens()
    
            hydRoutines = hydrogenRoutines(self.routines)
    
            if debump:
                self.routines.debumpProtein()  
    
            if opt:
                hydRoutines.setOptimizeableHydrogens()
                hydRoutines.initializeFullOptimization()
                hydRoutines.optimizeHydrogens()
            else:
                hydRoutines.initializeWaterOptimization()
                hydRoutines.optimizeHydrogens()
    
            # Special for GLH/ASH, since both conformations were added
            hydRoutines.cleanup()
                    
        
def get_residue_titration_set_protonated(residues):
    """
    Returns residue set when everything is protonated.
    """
    result = []
    for residue in residues:
        residueTest = _titrationSetsMap.get(residue.name)
        if residueTest:
            residueTest = residueTest[1]
        else:
            residueTest = residue.name
        result.append(residueTest)
        
    return result

def residue_set_single_unprotonated_combinations(residues):
    """
    Yields pair (residue set, residue index) for 
    every "single unprotonated" combination.
    residue set - set for process_residue_set
    residue index - index of residue that was left unprotonated
    """    
    protonatedNames = get_residue_titration_set_protonated(residues)
    
    for name, i in zip(protonatedNames, count()):
        if not name in _titrationSetsMap:
            continue
        
        tStateSet = _titrationSetsMap[name][0]
        
        for tState in tStateSet:
            result = list(protonatedNames)
            result[i] = tState
            yield result, i
            
def residue_set_pair_unprotonated_combinations(residues):
    """
    Yields pair (residue set, 1rst residue index, 2nd residue index) for 
    every "single unprotonated" combination.
    residue set - set for process_residue_set
    1rst residue index - index of 1rst residue that was left unprotonated
    2nd residue index - index of 2nd residue that was left unprotonated
    """    
    protonatedNames = get_residue_titration_set_protonated(residues)
    
    for i in range(0,len(protonatedNames)):
        firstName = protonatedNames[i]
        if not firstName in _titrationSetsMap:
            continue
        firstStateSet = _titrationSetsMap[firstName][0]
        for j in range(0,i):
            secondName = protonatedNames[j]
            if not secondName in _titrationSetsMap:
                continue            
            
            secondStateSet = _titrationSetsMap[secondName][0]
            
            for firstState in firstStateSet:
                for secondState in secondStateSet:
                    result = list(protonatedNames)
                    result[i] = firstState
                    result[j] = secondState
                    yield result, i, j


def get_residue_interaction_energy(residue1, residue2):
    """
    Returns to total energy of every atom pair between the two residues.
    
    Uses Optimize.getPairEnergy and it's donor/accepter model 
    to determine energy.
    
    residue1 - "donor" residue
    residue2 - "acceptor" residue
    
    THE RESULTS OF THIS FUNCTION ARE NOT SYMMETRIC. Swapping 
    residue1 and residue2 will not always produce the same result.
    """
    energy = 0.0
    for pair in product(residue1.getAtoms(), residue2.getAtoms()):
        energy += Optimize.getPairEnergy(pair[0], pair[1])
        
    return energy
    

