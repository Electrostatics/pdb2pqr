# These are notes for the current version of PDB2PQR

Please see http://www.poissonboltzmann.org/pdb2pqr/release-history for the complete release history


## NEW FEATURES
* Option to automatically drop water from pdb file before processing.
* Intergration of PDB2PKA into PDB2PQR as an alternative to PROPKA.
* Support for compiling with VS2008 in Windows.
* Option to build with debug headers.
* PDB2PKA now detects and reports non Henderson-Hasselbalch behavior. 
* PDB2PKA can be instructed whether or not to start from scratch with --pdb2pka-resume
* Can now specify output directory for PDB2PKA.
* Improved error regarding backbone in some cases.
 
## BUG FIXES
* Fixed executable name when creating binaries for Unix based operating systems.
* Fixed potential crash when using --clean with extensions.
* Fixed MAXATOMS display on server home page.
* PDB2PKA now mostly respects the --verbose setting.
* Fixed how hydrogens are added by PDB2PKA for state changes in some cases.
* Fixed psize error check.
* Will now build properly without ligand support if numpy is not installed.
* Removed old automake build files from all tests ported to scons.
 
## CHANGES 
* Command line interface to PROPKA changed to accommodate PDB2PKA. PROPKA is now used with --ph-calc-method=propka. --with-ph now defaults to 7.0 and is only required if a different pH value is required.
* --ph-calc-method to select optional method to calculate pH values used to protonate titratable residues. Possible options are "propka" and "pdb2pka". 
* Dropped support for compilation with mingw. Building on Windows now requires VS 2008 or VS Express 2008 and Windows SDK 6.0A installed in the default location.
* Updated included Scons to 2.3.3
* PDB2PKA can now be run directly (not integrated in PDB2PQR) with pka.py. Arguments are <PDB file> and <Output directory>.
 
## KNOWN BUGS
* If more than one extension is run from the command line and one of the extensions modifies the protein data structure it could affect the output of the other extension. The only included extensions that exhibit this problem are resinter and newresinter.
