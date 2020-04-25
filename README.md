pdb2pqr
============

This directory is the home for the [PDB2PQR software](http://www.poissonboltzmann.org/docs/structures-ready/).

Some useful information:
* See the [install notes](INSTALL.md) to build pdb2pqr.
* Using pdb2pqr [binaries](BINARY_README.md)
* Latest [changelog](ChangeLog.md)

Please see the [user guide](http://www.poissonboltzmann.org/docs/pdb2pqr-algorithm-description/) for documentation and the [COPYING file](COPYING) for license information.

Please see [programmer's guide](http://www.poissonboltzmann.org/docs/pdb2pqr-programmers/) for information on working with the PDB2PQR code.

# Binary Releases

Binary builds do not require python or numpy be installed to use. Everything needed to run PDB2PQR is included. Just unpack and use.

## OSX
+ Executable is called "pdb2pqr" and should already have the executable flag set
+ OSX binaries require OSX 10.6 or newer. The OSX binary is 64-bit.

## Linux
+ Executable is called "pdb2pqr" and should already have the executable flag set.
+ Linux binaries require CentOS 6 or newer and have been tested on Ubuntu 16.04 LTS and Linux Mint 13. If you are running 64-bit Linux use the 64-bit libraries.
# PDB2PQR Installation Instructions

Installation on most systems is rather straightforward - as the bulk of the PDB2PQR code is written in Python, the PDB2PQR code itself is architecture/compiler independent. PDB2PQR has been tested using Python versions 2.6 through 2.7 - problems will occur with other versions.

If you would like to enable pdb2pka or ligand support, then you must have Numpy installed.
http://numpy.scipy.org/

Networkx is now required for pdb2pka support.
https://networkx.github.io/

# Configuration and Build
PDB2PQR will need to be configured and compiled in order to run:

	$ python scons/scons.py (see python scons/scons.py --help for more options)

This should configure pdb2pqr and compile the pdb2pka wrappers necessary to interface with PDB2PQR and ligand support.

### Configuration File
Compilation and installation can be configured by editing the build_config.py file.
This is the preferred way to configure the program.
Instructions and examples for each setting are included in the file.

### Configuration Command Line Parameters
These will override any setting in build_config.py.

	PREFIX=<DIR>                	Set install directory. Default is ~/pdb2pqr
	APBS=<APBS_BINARY>            Location of APBS binary.
	MAX_ATOMS=<MAX_ATOMS>			Sets the maximum number of atoms in a protein for non-Opal job submission. Only affects web tools. Default is 10000
	BUILD_PDB2PKA=False           Disable pkb2pka compilation. Needed if no C++ compiler is installed or numpy is not installed. pdb2pka is required for ligand support.

## Installation
Installing is only needed if you plan to set up the web service or to make the install available to all users.
The installation location is specified by the <code>PREFIX</code> setting above.

 	$ python scons/scons.py install

## Using pdb2pka
If you are planning to invoke the flag,

	--titration-state-method=pdb2pka

you will need to build the python APBS libraries and set your library path to point to them.  For APBS build instructions see the documention provided [here](../apbs/README.md).

## Using a different python.
pdb2pqr will be configured to use whichever python was used to run the build script.
If you would like to use a different installed python to run pdb2pqr run the scons/scons.py script with that python.
For example:

 	$ /opt/python3/python scons/scons.py

If ligand support is required Numpy must be installed on the python used to build pdb2pqr.

## Windows Support
Compilation of pdb2pka on Windows requires that VS2008 or VS Express 2008 and Windows SDK 10.0A installed in the default location.

## Numpy
If numpy cannot be installed directly on the python used to run pdb2pqr you can either install numpy for your local user account or use virtualenv. Homebrew may be an option on OSX.
http://docs.python.org/2/install/index.html#alternate-installation
http://stackoverflow.com/questions/7465445/how-to-install-python-modules-without-root-access
http://www.virtualenv.org

##Binary Builds
A PyInstaller spec file is included if a stand alone binary build (no python install or compile required to run) is desired.
http://www.pyinstaller.org/

PyInstaller 2.1 or newer installed as python library is required.
To create standalone build, first build the application as normal then run

	pyinstaller pdb2pqr.spec

in the root archive folder. The distributable program will be in the dist/pdb2pqr folder.

Each supported platform has different needs when creating the binary build:

### Linux
The built binaries will dynamically link against the system standard libs.
For this reason the Linux binaries should be built on an older system to increase compatibility if needed.

### OSX
On the Mac you must use homebrew and use the python that it installs with <code>brew python</code>. See http://brew.sh/ and https://github.com/Homebrew/homebrew/wiki/Homebrew-and-Python
Once homebrew python is setup and configured correctly pyinstaller, numpy and networkx can be installed with <code>pip install numpy</code>, <code>pip install pyinstaller</code>, and <code>pip install networkx</code>.
If the global python is used the binary will fail when the <code>--ligand</code> option is used.
Similar to Linux, this binary is dynamically linked against the system libraries. The binary will not work on OSX older than the version built on.
We officially support OSX 10.6 or newer.
