# These are notes for the current version of PDB2PQR

Please see http://www.poissonboltzmann.org/pdb2pqr/release-history for the complete release history


## NEW FEATURES
* A ligand file with duplicate atoms will cause pdb2pqr to stop instead of issue a warning.  Trust us, this is a feature, not a bug!
* Improved error reporting.
* Added support for reference command line option for PROPKA.
* Added newresinter plugin to provide alternate methods for calculating interaction energies between residues.
* Mol2 file handling is now case insensitive with atom names. 
* PROPKA with a pH of 7 is now specified by default on the web service.
* Compilation is now done with scons.
* Verbose output now includes information on all patches applied during a run.
* Added stderr and stdout to web error page. 
* Added warning to water optimization when other water is ignored.
* Command line used to generate a pqr is now duplicated in the comments of the output.
* Added support for NUMMDL in parser.
* Added complete commandline feature test. Use complete-test target.
* Added propka support for phosphorous sp3. - Thanks to Dr. Stefan Henrich
* Added a PyInstaller spec file. Standalone pdb2pqr builds are now possible.
 
## BUG FIXES
* Rolled back change that prevented plugins from interfering with each other. Large proteins would cause a stack overflow when trying to do a deep copy.
* Updated INSTALL file to reflect no more need for Fortran.
* Fixed apbs input file to match what web interface produces.
* Fixed user specified mobile ion species not being passed to apbs input file.
* Removed ambiguous A, ADE, C, CYT, G, GUA, T, THY, U, URA as possible residue names.
* Removed eval from pdb parsing routines.
* Updated web links to refer to http://www.poissonboltzmann.org where appropriate. 
* Fixed hbond extension output to include insertion code in residue name.
* Fixed debumping routines not including water in their checks. Fixes bad debump of ASN B 20 in 1gm9 when run with pH 7.0.
* Fixed debumping failing to use best angle for a specific dihedral angle when no tested angles are without conflict.
* Fixed debumping using asymmetrical cutoffs and too large cutoffs in many checks involving hydrogen.
* Fixed debumping accumulating rounding error while checking angles.
* Fixed inconsistencies in pdb parsing. - Thanks to Dr. Stefan Henrich
* Fixed problems with propka handling of aromatic carbon/nitrogen. - Thanks to Dr. Stefan Henrich
* Fixed case where certain apbs compile options would break web visualization. 
* Fixed improper handling of paths with a '.' or filenames with more than one '.' in them.
 
## CHANGES
* Removed numpy from contrib. The user is expected to have numpy installed and available to python at configuration.
* Support for numeric dropped. 
 
## KNOWN BUGS
* If more than one extension is run from the command line and one of the extensions modifies the protein data structure it could affect the output of the other extension. The only included extensions that exhibit this problem are resinter and newresinter.
