# These are notes for the current version of PDB2PQR

Please see http://www.poissonboltzmann.org/pdb2pqr/release-history for the complete release history


## NEW FEATURES
* PDB2PQR can no longer be run in Python 2.x.
* PDB2PQR uses PROPKA 3.1 now which add support protein-ligand complexes.
* Added support for reading mmCIF formatted files.
* SCons builder has been updated to 3.1.1

## BUG FIXES
* Miscellaneous bugs were fixed. See closed issues in repo.

## CHANGES
* The networkx library is required for pdb2pka.

## KNOWN BUGS
* If more than one extension is run from the command line and one of the extensions modifies the protein data structure it could affect the output of the other extension. The only included extensions that exhibit this problem are resinter and newresinter.
* PDB2PKA currently leaks memory slowly. Small jobs will use about twice the normally required RAM (ie ~14 titratable residues will use 140MB). Big jobs will use about 5 times the normally required RAM ( 60 titratable residues will use 480MB ). We are working to fix this.
