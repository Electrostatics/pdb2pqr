apbs_inputfile=self.apbs_setup.printInput()
self.APBS=runAPBS()
potentials = self.APBS.runAPBS(self.protein, apbs_inputfile, CM)
self.APBS.cleanup()
               

# Here, self.protein is the PDB2PQR instance. you will get this with something like:
dummydef = Definition()
myProtein = Protein(pdblist, dummydef)

# then you do:
myRoutines = Routines(myProtein, verbose) #myDefinition)
myRoutines.updateResidueTypes()
myRoutines.updateSSbridges()
myRoutines.updateBonds()
myRoutines.setTermini()
myRoutines.updateInternalBonds()

myRoutines.applyNameScheme(Forcefield(ff, myDefinition, None))
myRoutines.findMissingHeavy()
myRoutines.addHydrogens()
myRoutines.debumpProtein()

#myRoutines.randomizeWaters()
myProtein.reSerialize()

# to clean up and add hydrogens.

# optimze hydrogens:
from src.hydrogens import hydrogenRoutines
myRoutines.updateInternalBonds()
myRoutines.calculateDihedralAngles()
myhydRoutines = hydrogenRoutines(myRoutines)
#
# Here we should inject the info!!
#
myhydRoutines.setOptimizeableHydrogens()
myhydRoutines.initializeFullOptimization()
myhydRoutines.optimizeHydrogens()
myhydRoutines.cleanup()
myRoutines.setStates()

# Choose the force field:
myForcefield = Forcefield(ff, myDefinition, None)

# To print the charges do:
for chain in protein.getChains():
    for residue in chain.get("residues"):
        for atom in residue.get("atoms"):
            atomname = atom.get("name")
            charge, radius = forcefield.getParams(residue, atomname)
