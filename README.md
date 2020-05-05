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

Configuration and build happens from the top-level of the apbs-pdb2pqr repository; e.g.

```
$ mkdir build
$ cd build
$ cmake -DENABLE_PYTHON=ON -DCMAKE_C_FLAGS="-fPIC" ..
$ cmake --build .
```

# Testing

Testing can happen from this directory:

```
$ PYTHONPATH=propka31 python -m pytest
```

# Old documentation is that is probably no longer relevant

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
