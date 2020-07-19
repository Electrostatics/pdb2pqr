===============
Release history
===============

-----------------------
PDB2PQR 2.1.1 (2016-03)
-----------------------

New features
    * Replaced the Monte Carlo method for generating titration curves with Graph Cut. See http://arxiv.org/1507.07021/

Bug fixes
    * Added a check before calculating pKa's for large interaction energies

Known bugs
    * If more than one extension is run from the command line and one of the extensions modifies the protein data structure it could affect the output of the other extension. The only included extensions that exhibit this problem are resinter and newresinter.
    * Running ligands and PDB2PKA at the same time is not currently supported.
    * PDB2PKA currently leaks memory slowly. Small jobs will use about twice the normally required RAM (i.e. ~14 titratable residues will use 140MB). Big jobs will use about 5 times the normally required RAM (60 titratable residues will use 480MB). We are working on this.

-----------------------
PDB2PQR 2.1.0 (2015-12)
-----------------------

New features
    * Added alternate method to do visualization using 3dmol.
    * Replaced the Monte Carlo method for generating titration curves with Graph Cut.  See http://arxiv.org/abs/1507.07021.  If you prefer the Monte Carlo Method, please use http://nbcr-222.ucsd.edu/pdb2pqr_2.0.0/

Bug fixes
    * Added compile options to allow for arbitrary flags to be added.  Helps work around some platforms where scons does not detect the needed settings correctly.
    * Fixed broken links on APBS submission page.
    * Added some missing files to query status page results.
    * Fixed some pages to use the proper CSS file.
    * Better error message for ``--assign-only`` and HIS residues.
    * Fixed PROPKA crash for unrecognized residue.
    * Debumping routines are now more consistent across platforms.  This fixes pdb2pka not giving the same results on different platforms.

Known bugs
    * If more than one extension is run from the command line and one of the extensions modifies the protein data structure it could affect the output of the other extension. The only included extensions that exhibit this problem are resinter and newresinter.
    * Running ligands and PDB2PKA at the same time is not currently supported.
    * PDB2PKA currently leaks memory slowly. Small jobs will use about twice the normally required RAM (i.e. ~14 titratable residues will use 140MB). Big jobs will use about 5 times the normally required RAM (60 titratable residues will use 480MB). We are working on this.

Other comments
    * Added fabric script used to build and test releases.
    * The newtworkx library is now required for pdb2pka.

-----------------------
PDB2PQR 2.0.0 (2014-12)
-----------------------

New features
    * Improved look of web interface.
    * Option to automatically drop water from pdb file before processing.
    * Integration of PDB2PKA  into PDB2PQR as an alternative to PROPKA.
    * Support for compiling with VS2008 in Windows.
    * Option to build with debug headers.
    * PDB2PKA now detects and reports non Henderson-Hasselbalch behavior.
    * PDB2PKA can be instructed whether or not to start from scratch with ``--pdb2pka-resume``.
    * Can now specify output directory for PDB2PKA.
    * Improved error regarding backbone in some cases.
    * Changed time format on query status page.
    * Improved error catching on web interface.

Bug fixes
    * Fixed executable name when creating binaries for Unix based operating systems.
    * Fixed potential crash when using ``--clean`` with extensions.
    * Fixed MAXATOMS display on server home page.
    * PDB2PKA now mostly respects the ``--verbose`` setting.
    * Fixed how hydrogens are added by PDB2PKA for state changes in some cases.
    * Fixed psize error check.
    * Will now build properly without ligand support if numpy is not installed.
    * Removed old automake build files from all test ported to scons.
    * Fixed broken opal backend.

Known bugs
    * If more than one extension is run from the command line and one of the extensions modifies the protein data structure it could affect the output of the other extension.  The only included extensions that exhibit this problem are resinter and newresinter.
    * Running ligands and PDB2PKA at the same time is not currently supported.
    * PDB2PKA currently leaks memory slowly.  Small jobs will use about twice the normally required RAM (i.e. ~14 titratable residues will use 140MB).  Big jobs will use about 5 times the normally required RAM (60 titratable residues will use 480MB).  We are working on this.

Other comments
    * Command line interface to PROPKA changed to accommodate PDB2PKA. PROPKA is now used with ``--ph-calc-method=propka --with-ph`` now defaults to 7.0 and is only required if a different pH value is required.
    * ``--ph-calc-method`` to select optional method to calculate pH values used to protonate titratable residues.  Possible options are "propka" and "pdb2pka".
    * Dropped support for compilation with mingw.  Building on Windows now requires VS 2008 installed in the default location.
    * Updated included Scons to 2.3.3
    * PDB2PKA can now be run directly (not integrated in PDB2PQR) with pka.py. Arguments are PDBfile and Output directory.
    * No longer providing 32-bit binary build.  PDB2PKA support is too memory intensive to make this practical in many cases.

---------------------
PDB2PQR 1.9 (2014-03)
---------------------

New features
    * Binary builds do not require python or numpy be installed to use. Everything needed to run PDB2PQR is included. Just unpack and use.
    * OSX binaries require OSX 10.6 or newer. The OSX binary is 64-bit.
    * Linux binaries require CentOS 6 or newer and have been tested on Ubuntu 12.04 LTS and Linux Mint 13. If you are running 64-bit Linux use the 64-bit libraries. In some cases the needed 32-bit system libraries will not be installed on a 64-bit system.
    * Windows binaries are 32 bit and were built and tested on Windows 7 64-bit but should work on Windows XP, Vista, and 8 both 32 and 64-bit systems.
    * PDB2PQR can now be compiled and run on Windows using MinGW32. See http://mingw.org/ for details.
    * PDB2PQR now uses Scons for compilations. With this comes improved automated testing.
    * A ligand file with duplicate atoms will cause pdb2pqr to stop instead of issue a warning. Trust us, this is a feature, not a bug!
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

Bug fixes
    * Rolled back change that prevented plugins from interfering with each other. Large proteins would cause a stack overflow when trying to do a deep copy
    * Updated INSTALL file to reflect no more need for Fortran.
    * Fixed apbs input file to match what web interface produces.
    * Fixed user specified mobile ion species not being passed to apbs input file.
    * Removed ambiguous A, ADE, C, CYT, G, GUA, T, THY, U, URA as possible residue names.
    * Removed eval from pdb parsing routines.
    * Updated web links where appropriate.
    * Fixed hbond extension output to include insertion code in residue name.
    * Fixed debumping routines not including water in their checks. Fixes bad debump of ASN B 20 in 1gm9 when run with pH 7.0.
    * Fixed debumping failing to use best angle for a specific dihedral angle when no tested angles are without conflict.
    * Fixed debumping using asymmetrical cutoffs and too large cutoffs in many checks involving hydrogen.
    * Fixed debumping accumulating rounding error while checking angles.
    * Fixed inconsistencies in pdb parsing. - Thanks to Dr. Stefan Henrich
    * Fixed problems with propka handling of aromatic carbon/nitrogen. - Thanks to Dr. Stefan Henrich
    * Fixed case where certain apbs compile options would break web visualization.
    * Fixed improper handling of paths with a '.' or filenames with more than one '.' in them.

Known bugs
    * If more than one extension is run from the command line and one of the extensions modifies the protein data structure it could affect the output of the other extension. The only included extensions that exhibit this problem are resinter and newresinter.

Other comments
    * Removed numpy from contrib. The user is expected to have numpy installed and available to python at configuration.
    * Support for numeric dropped.


---------------------
PDB2PQR 1.8 (2012-01)
---------------------

New features
    * Updated PROPKA to version 3.0
    * Added residue interaction energy extension
    * Added protein summary extension
    * Combined hbond and hbondwhatit into one extension (hbond) with new command line parameters
    * Combined rama, phi, psi into one extension (rama) with new command line parameters.
    * Extensions may now add their own command line arguments. Extensions with their own command line arguments will be grouped separately.
    * Improved interface for extensions
    * Added Opal configuration file.

Bug fixes
    * Cleaned up white space in several files and some pydev warnings
    * Creating print output no longer clears the chain id data from atoms in the data. (Affected resinter plugin)
    * Removed possibility of one plug-in affecting the output of another
    * Fixed --protonation=new option for propka30
    * Improved time reporting for apbs jobs
    * Fixed opal runtime reporting
    * Fixed misspelled command line options that prevented the use of PEOEPB and TYL06
    * Fixed error handling when certain data files are missing
    * Fixed LDFLAGS environment variable not being used along with python specific linker flags to link Algorithms.o and _pMC_mult.so
    * Fixed possible Attribute error when applying naming scheme.

---------------------------
PDB2PQR 1.7.1a (2011-09-13)
---------------------------

New features
    * Added force field example.

Bug fixes
    * Fixed ligand command line option.
    * Fixed capitalization of force field in PQR header.
    * Fixed error handling for opal errors.
    * Fixed web logging error when using ligand files, user force fields, and name files.
    * Fixed extension template in documentation.
    * Fixed 1a1p example README to reflect command line changes.

-----------------------
PDB2PQR 1.7.1 (2011-08)
-----------------------

New features
    * Switched Opal service urls from sccne.wustl.edu to NBCR.
    * Added more JMol controls for visualization, JMol code and applets provided by Bob Hanson.
    * Changed default forcefield to PARSE in web interface.

Bug fixes
    * Fixed crash when opal returns an error.
    * Fixed specific combinations of command-line arguments causing :file:`pdb2pqr.py` to crash.
    * Fixed opal job failing when filenames have spaces or dashs.
    * Fixed gap in backbone causing irrationally placed hydrogens.
    * Fixed crash when too many fixes are needed when setting termini.
    * Corrected web and command line error handling in many cases.
    * Fixed ``--username`` command line option.
    * Fixed ambiguous user created forcefield and name handling. Now ``--username`` is required if ``--userff`` is used. 
    * Fixed :file:`querystatus.py` not redirecting to generated error page.

---------------------
PDB2PQR 1.7 (2010-10)
---------------------

* For PDB2PQR web interface users:  the JMol web interface for APBS calculation visualization has been substantially improved, thanks to help from Bob Hanson.  Those performing APBS calculations via the PDB2PQR web interface now have a much wider range of options for visualizing the output online -- as well as downloading for offline analysis.
* For PDB2PQR command-line and custom web interface users:  the Opal service URLs have changed to new NBCR addresses.  Old services hosted at .wustl.edu addresses have been decommissioned.  Please upgrade ASAP to use the new web service.  Thank you as always to the staff at NBCR for their continuing support of APBS/PDB2PQR web servers and services.



---------------------
PDB2PQR 1.6 (2010-04)
---------------------

New features
    * Added Swanson force field based on Swanson et al paper (http://dx.doi.org/10.1021/ct600216k).
    * Modified printAtoms() method. Now "TER" is printed at the end of every chain.
    * Added Google Analytics code to get the statistics on the production server.
    * Modified APBS calculation page layout to hide parameters by default and display PDB ID
    * Added "make test-webserver", which tests a long list of PDBs (246 PDBs) on the production PDB2PQR web server.
    * Removed ``nlev`` from :file:`inputgen.py` and :file:`inputgen_pKa.py` as nlev keyword is now deprecated in APBS.
    * Added PARSE parameters for RNA, data from: Tang C. L., Alexov E, Pyle A. M., Honig B. Calculation of pKas in RNA: On the Structural Origins and Functional Roles of Protonated Nucleotides. Journal of Molecular Biology 366 (5) 1475-1496, 2007.

Bug fixes
    * Fixed a minor bug: when starting :file:`pka.py` from pdb2pka directory using command like "python pka.py [options] inputfile", we need to make sure scriptpath does not end with "/".
    * Fixed a bug which caused "coercing to Unicode: need string or buffer, instance found" when submitting PDB2PQR jobs with user-defined force fields on Opal based web server. 
    * Fixed a bug in :file:`main_cgi.py`, now Opal-based PDB2PQR jobs should also be logged in :file:`usage.txt` file.
    * Updated :file:`src/utilities.py` with a bug fix provided by Greg Cipriano, which prevents infinite loops in analyzing connected atoms in certain cases.
    * Fixed a bug related to neutraln and/or neutralc selections on the web server.
    * Fixed a special case with ``--ffout`` and 1AIK, where the N-terminus is acetylated.
    * Fixed a bug in :file:`psize.py` per Michael Lerner's suggestion. The old version of :file:`psize.py` gives wrong cglen and fglen results in special cases (e.g., all y coordinates are negative values).
    * Fixed a bug in :file:`main_cgi.py`, eliminated input/output file name confusions whether a PDB ID or a pdb file is provided on the web server.
    * Fixed a bug which causes run time error on the web server when user-defined force field and names files are provided.
    * Fixed a bug in :file:`apbs_cgi.py`: pdb file names submitted by users are not always 4 characters long.

---------------------
PDB2PQR 1.5 (2009-10)
---------------------

New features
    * APBS calculations can be executed through the PDB2PQR web interface in the production version of the server
    * APBS-calculated potentials can be visualized via the PDB2PQR web interface thanks to Jmol
    * Disabled Typemap output by default, added --typemap flag to create typemap output if needed.
    * Enabled "Create APBS Input File" by default on the web server, so that APBS calculation and visualization are more obvious to the users.
    * Added warnings to stderr and the REMARK field in the output PQR file regarding multiple occupancy entries in PDB file.
    * Added more informative messages in REMARK field, explaining why PDB2PQR was unable to assign charges to certain atoms.
    * Updated structures.py, now PDB2PQR keeps the insertion codes from PDB files.
    * Added "make test-long", which runs PDB2PQR on a long list (246) of PDBs by default, it is also possible to let it run on specified number of PDBs, e.g.,  export TESTNUM=50; make test-long 
    * Updated NBCR opal service urls from http://ws.nbcr.net/opal/... to http://ws.nbcr.net/opal2/...
    * Compressed APBS OpenDX output files in zip format, so that users can download zip files from the web server.
    * Removed "EXPERIMENTAL" from APBS web solver interface and Jmol visualization interface.
    * Updated all APBS related urls from http://apbs.sourceforge.net/... to http:/apbs.wustl.edu/...
    * Merged PDB2PKA code, PDB2PKA is functional now.
    * Added two new options: --neutraln and --neutralc, so that users can manually make the N-termini or C-termini of their proteins neutral.    
    * Added a local-test, which addresses the issue of Debian-like Linux distros not allowing fetching PDBs from the web.
    * Added deprotonated Arginine form for post-PROPKA routines. This only works for PARSE forcefield as other forcefields lack deprotonated ARG parameters.
    * Updated inputgen.py with --potdx and --istrng options added, original modification code provided by Miguel Ortiz-Lombardía.
    * Changed default Opal service from http://ws.nbcr.net/opal2/services/pdb2pqr_1.4.0 to http://sccne.wustl.edu:8082/opal2/services/pdb2pqr-1.5

Bug fixes
    * Verbosity outputs should be stdouts, not stderrs in web server interface. Corrected this in src/routines.py.
    * Fixed a bug in psize.py: for a pqr file with no ATOM entries but only HETATM entries in it, inputgen.py should still create an APBS input file with reasonable grid lengths. 
    * Added special handling for special mol2 formats (unwanted white spaces or blank lines in ATOM or BOND records).
    * Added template file to doc directory, which fixed a broken link in  programmer guide.

-----------------------
PDB2PQR 1.4.0 (2009-03)
-----------------------

New features
    * Updated html/master-index.html, deleted html/index.php.
    * Updated pydoc by running genpydoc.sh.
    * Added a whitespace option by by putting whitespaces between atom name and residue name, between x and y, and between y and z.
    * Added radius for Chlorine in ligff.py.
    * Added PEOEPB forcefield, data provided by Paul Czodrowski.
    * Updated inputgen.py to write out the electrostatic potential for APBS input file. 
    * Updated CHARMM.DAT with two sets of phosphoserine parameters.
    * Allowed amino acid chains with only one residue, using --assign-only option.
    * Updated server.py.in so that the ligand option is also recorded in usage.txt. 
    * Updated HE21, HE22 coordinates in GLN according to the results from AMBER Leap program.
    * Updated Makefile.am with Manuel Prinz's patch (removed distclean2 and appended its contents to distclean-local).
    * Updated configure.ac, pdb2pqr-opal.py; added AppService_client.py and AppService_types.py with Samir Unni's changes, which fixed earlier problems in invoking Opal services.
    * Applied two patches from Manuel Prinz to pdb2pka/pMC_mult.h and pdb2pka/ligand_topology.py. 
    * Updated PARSE.DAT with the source of parameters. 
    * Created a contrib folder with numpy-1.1.0 package. PDB2PQR will install numpy by default unless any of the following conditions is met:
    * Working version of NumPy dectected by autoconf.
    * User requests no installation with --disable-pdb2pka option.
    * User specifies external NumPy installation.  
    * Merged Samir Unni's branch. Now PDB2PQR Opal and APBS Opal services are available (through --with-opal and/or --with-apbs, --with-apbs-opal options at configure stage).
    * Added error handling for residue name longer than 4 characters.
    * Updated hbond.py with Mike Bradley's definitions for ANGLE_CUTOFF and DIST_CUTOFF by default.
    * Removed PyXML-0.8.4, which is not required for ZSI installation.
    * Updated propka error message for make adv-test -- propka requires a version of Fortran compiler.
    * Updated na.py and PATCHES.xml so that PDB2PQR handles three lettered RNA residue names (ADE, CYT, GUA, THY, and URA) as well.
    * Updated NA.xml with HO2' added as an alternative name for H2'', and H5" added as an alternative name for H5''. 
    * Updated version numbers in html/ and doc/pydoc/ .
    * Updated web server. When selecting user-defined forcefield file from the web server, users should also provide .names file.
    * Removed http://enzyme.ucd.ie/Services/pdb2pqr/ from web server list.
    * Eliminated the need for protein when processing other types (ligands,  nucleic acids).
    * Updated psize.py with Robert Konecny's patch to fix inconsistent assignment of fine grid numbers in some (very) rare cases.
    * Made whitespace option available for both command line and web server versions.
    * Updated inputgen_pKa.py with the latest version.

Bug fixes
    * Fixed a legacy bug with the web server (web server doesn't like ligand files generated on Windows or old Mac OS platforms).
    * Fixed a bug in configure.ac, so that PDB2PQR no longer checks for Numpy.pth at configure stage.
    * Updated pdb2pka/substruct/Makefile.am. 
    * Fixed isBackbone bug in definitions.py.
    * Fixed a bug for Carboxylic residues in hydrogens.py.
    * Fixed a bug in routines.py, which caused hydrogens added in LEU and ILE in eclipsed conformation rather than staggered. 
    * Fixed a bug in configure.ac, now it is OK to configure with double slashes in the prefix path, e.g.,  --prefix=/foo/bar//another/path 
    * Fixed a bug in nucleic acid naming scheme. 
    * Fixed a bug involving MET, GLY as NTERM, CTERM with --ffout option.
    * Fixed a bug for PRO as C-terminus with PARSE forcefield. 
    * Fixed a bug for ND1 in HIS as hacceptor.
    * Fixed the --clean option bug.
    * Fixed a bug in CHARMM naming scheme.
    * Fixed a bug in test.cpp of the simple test (which is related to recent modifications of 1AFS in Protein Data Bank).

-----------------------
PDB2PQR 1.3.0 (2008-01)
-----------------------

New features
    * Added "make test" and "make adv-test"
    * Fixed problems with "make dist"
    * Added integration with Opal for launching jobs as well as querying status
    * The user may use NUMPY to specify the location of NUMPY.
    * Both PDB2PKA and PROPKA are enabled by default.  PDB2PKA is enabled by default since ligand parameterization would fail without this option.
    * For a regular user, "make install" tells the user the exact command the system administrator will use to make the URL viewable.
    * The default value of 7.00 for the pH on the server form is removed due to a problem with browser refershing.
    * Updated warning messages for lines beginning with SITE, TURN, SSBOND and LINK.
    * Switched license from GPL to BSD.
    * Made a new tar ball pdb2pqr-1.3.0-1.tar.gz for Windows users who cannot create pdb2pqr.py through configure process.
    * configure now automatically detects SRCPATH, WEBSITE, and the location of pdb2pqr.cgi.  In version 1.2.1, LOCALPATH(SRCPATH) and WEBSITE were defined in src/server.py and the location of pdb2pqr.cgi was specified in html/server.html (index.html).  Configure now uses variable substitution with new files src/server.py.in and html/server.html.in to create src/server.py and html/server.html (index.html).
    * SRCPATH is automatically set to the current working directory. WEBSITE is automatically set to http://fully_qualified_domain_name/pdb2pqr. Path to CGI is automcailly set to http://fully_qualified_domain_name/pdb2pqr/pdb2pqr.cgi.  
    * In version 1.2.1, there were 3 variables that needed to be changed to set up a server at a location different from agave.wustl.edu.  LOCALPATH, WEBSITE, and the location of the CGI file.  In this version, LOCALPATH has been used to SRCPATH to avoid confusion, since LOCALPATH could be interpreted as the local path for source files or the localpath for the server.
    * Since configure now automatically sets the locations of files/directories based on the machine and configure options, the default  agave.wustl.edu locations are not used anymore.
    * A copy of pdb2pqr.css is included.
    * configure prints out information about parameters such as python flags, srcpath, localpath, website, etc.
    * configure now automatically creates tmp/ with r + w + x permissions.
    * configure now automatically copies pdb2pqr.py to pdb2pqr.cgi.
    * configure now automatically copies html/server.html to index.html after variable substitution.  In src/server.py.in (src/server.py), WEBNAME is changed to index.html. 
    * ${HOME}/pdb2pqr is the default prefix for a regular user
    * /var/www/html is the default prefix for root
    * http://FQDN/pdb2pqr as default website.  
    * "make install" runs "make" first, and the copies the approprite files to --prefix.
    * If root did not specify --prefix and /var/www/html/pdb2pqr already
    * exists, then a warning is issued, and the user may choose to quit or overwrite that directory.  
    * Similary, if a regular user did not specify --prefix and ${HOME}/pdb2pqr already exists, then a warning is issued, and the user may choose to quit or overwrite that directory. 
    * If root does not specify --prefix to be a directory to be inside /var/www/html (for example, --prefix=/share/apps/pdb2pqr), then a symbolic link will be made to /var/www/html/pdb2pqr during "make install".
    * configure option --with-url can be specified either as something like http://sandstone.ucsd.edu/pdb2pqr-test or sandstone.ucsd.edu/pdb2pqr-test.  It also doesn't matter if there's a '/' at the end.
    * If user is root, and the last part of URL and prefix are different, for example, --with-url=athena.nbcr.net/test0 --prefix=/var/www/html/pdb2pqr-test, then a warning will be issued saying the server will be viewable from the URL specified, but not the URL based on pdb2pqr-test.  In other words, the server will be viewable from athena.abcr.net/test0, but not athena.nbcr.net/pdb2pqr-test.  During "make  install", a symbolic link is created to enable users to view the server from --with-url.
    * When making a symbolic link for root, if then link destination already exists as a directory or a symoblic link, then the user may choose to continue with creating the link and overwrite the original directory or quit.
    * If the user changes py_path when running configure for PDB2PQR, then the change also applies to PROPKA.

Bug fixes
    * Fixed the line feed bug. Now PDB2PQR handles different input files (.pdb and .mol2) created or saved on different platforms.
    * Fixed "hbondwhatif" warning at start up.

Known issues
    * The install directory name cannot contain dots.
    * For python 2.2, if PDB2PQR cannot find module sets, then sets needs to be copied from .../python2.2/site-packages/MYSQLdb/sets.py to .../lib/python2.2

-----------------------
PDB2PQR 1.2.1 (2007-04)
-----------------------

New features
    * Updated documentation to include instructions for pdb2pka support, references, more pydoc documents.
    * Added ligand examples to examples/ directory
    * Added native support for the TYL06 forcefield.  For more information on this forcefield please see Tan C, Yang L, Luo R.  How well does Poisson-Boltzmann implicit solvent agree with explicit solvent? A quantitative analysis. Journal of Physical Chemistry B.  110 (37), 18680-7, 2006. 
    * Added a new HTML output page which relays the different atom types between the AMBER and CHARMM forcefields for a generated PQR file (thanks to the anonymous reviewers of the latest PDB2PQR paper).

Bug fixes
    * Fixed bug where a segmentation fault would occur in PropKa if the N atom was not the first atom listed in the residue
    * Fixed error message that occurred when a blank line was found in a parameter file.
    * Better error handling in MOL2 file parsing.
    * Fixed bug where ligands were not supported on PDB files with multiple MODEL fields.

-----------------------
PDB2PQR 1.2.0 (2007-01)
-----------------------

New features
    * Added autoconf support for pdb2pka directory.
    * Added new support for passing in a single ligand residue in MOL2 format via the --ligand command.  Also available from the web server (with link to PRODRG for unsupported ligands).
    * Numerous additions to examples directory (see examples/index.html) and update to User Guide.

Bug fixes
    * Fixed charge assignment error when dealing with LYN in AMBER.
    * Fixed crash when a chain has a single amino acid residue.  The code now reports the offending chain and residue before exiting. 
    * Fixed hydrogen optimization bug where waters with no nearby atoms at certain orientations caused missing hydrogens.

-----------------------
PDB2PQR 1.1.2 (2006-06)
-----------------------

Bug fixes
    * Fixed a bug in the hydrogen bonding routines where PDB2PQR attempted to delete an atom that had already been deleted. (thanks to Rachel Burdge)
    * Fixed a bug in chain detection routines where PDB2PQR was unable to detect multiple chains inside a single unnamed chain (thanks to Rachel Burdge)
    * Fixed a second bug in chain detection routines where HETATM residues with names ending in "3" were improperly chosen for termini (thanks to Reut Abramovich)
    * Fixed a bug where chains were improperly detected when only containing one HETATM residue (thanks to Reut Abramovich)

-----------------------
PDB2PQR 1.1.1 (2006-05)
-----------------------

Bug fixes
    * Fixed a bug which prevented PDB2PQR from recognizing atoms from nucleic acids with "*" in their atom names. (thanks to Jaichen Wang)
    * Fixed a bug in the hydrogen bonding routines where a misnamed object led to a crash for very specific cases. (thanks to Josh Swamidass) 

-----------------------
PDB2PQR 1.1.0 (2006-04)
-----------------------

New features
    * Structural data files have been moved to XML format.  This should make it easier for users and developers to contribute to the project.
    * Added an extensions directory for small scripts.  Scripts in this directory will be automatically loaded into PDB2PQR has command line options for post-processing, and can be easily customized.
    * Code has been greatly cleaned so as to minimize values hard-coded into functions and to allow greater customizability via external XML files.  This includes a more object-oriented hierarchy of structures.
    * Improved detection of the termini of chains.
    * Assign-only now does just that - only assigns parameters to atoms without additions, debumping, or optimizations.
    * Added a --clean command line option which does no additions, optimizations, or forcefield assignment, but simply aligns the PDB columns on output.  Useful for using post-processing scripts like those in the extensions directory without modifying the original input file.
    * The --userff flag has been replaced by opening up the --ff option to user-defined files.
    * Pydoc documentation is now included in html/pydoc.
    * A programmer's guide has been included to explain programming decisions and ease future development.
    * A --ffout flag has been added to allow users to output a PQR file in the naming scheme of the desired forcefield.
    * User guide FAQ updated.
    * The efficiency of the hydrogen bonding detection script (--hbond) has been greatly improved.
    * Increased the number of options available to users via the PDB2PQR web server.

Bug fixes
    * Updated psize.py to use centers and radii when calculating grid sizes (thanks to John Mongan) 
    * Fixed bug where PDB2PQR could not read PropKa results from chains with more than 1000 residues (thanks to Michael Widmann)

-----------------------
PDB2PQR 1.0.2 (2005-12)
-----------------------

New features
   * Added ability for users to add their own forcefield files.  This should be particularly useful for HETATMs.
   * Added sdens keyword to inputgen.py to make PDB2PQR compatibile with APBS 0.4.0. 
   * Added a new examples directory with a basic runthrough on how to use the various features in PDB2PQR.

Bug fixes
   * Fixed a bug that was unable to handle N-Terminal PRO residues with hydrogens already present. 
   * Fixed two instances in the PropKa routines where warnings were improperly handled due to a misspelling.
   * Fixed instance where chain IDs were unable to be assigned to proteins with more than 26 chains.

-----------------------
PDB2PQR 1.0.1 (2005-10)
-----------------------

New features
   * Added citation information to PQR output.

Bug fixes
   * Fixed a bug during hydrogen optimization that left out H2 from water if the oxygen in question had already made 3 hydrogen bonds.

-----------------------
PDB2PQR 1.0.0 (2005-08)
-----------------------

This is the initial version of the PDB2PQR conversion utility.
There are several changes to the various "non-official" versions previously available:

* SourceForge has been chosen as a centralized location for all things related to PDB2PQR, including downloads, mailing lists, and bug reports.
* Several additions to the code have been made, including pKa support via PropKa, a new hydrogen optimization algorithm which should increase both accuracy and speed, and general bug fixes.

