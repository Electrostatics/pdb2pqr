# These are notes for the current version of PDB2PQR

Please see http://www.poissonboltzmann.org/pdb2pqr/release-history for the complete release history


## NEW FEATURES
* Added alternate method to do visualization using 3dmol.
* Replaced the Monte Carlo method for generating titration curves with graph cut. See http://arxiv.org/abs/1507.07021
* Added compile options to allow for arbitrary flags to be added. Helps work around some platforms where scons does not detect the needed settings correctly.

## BUG FIXES
* Fixed broken links on APBS submission page.
* Added some missing files to querystaus page results.
* Fixed some pages to use the proper CSS file.
* Better error message for --assign-only and HIS residues.
* Fix PROPKA crash for unrecognized residue.
* Debumping routines are now more consistent across platforms. This fixes pdb2pka not giving the same results on different platforms.

## CHANGES
* Added fabric script used to build and test releases.
* The networkx library is now required for pdb2pka.

## KNOWN BUGS
* If more than one extension is run from the command line and one of the extensions modifies the protein data structure it could affect the output of the other extension. The only included extensions that exhibit this problem are resinter and newresinter.
* Running ligands and PDB2PKA at the same time is not currently supported.
* PDB2PKA currently leaks memory slowly. Small jobs will use about twice the normally required RAM (ie ~14 titratable residues will use 140MB). Big jobs will use about 5 times the normally required RAM ( 60 titratable residues will use 480MB ). We are working to fix this.
