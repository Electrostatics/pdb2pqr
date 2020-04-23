#Build Configuration file for pdb2pqr
#While you can specify these on the command line with var=value
# this is the recommended way to setup a build.

#Uncomment the values you would like to change and set new values.


#Installation PREFIX
#Sets the install location of pdb2pqr.
#This defaults to ~/pdb2pqr

#PREFIX="~/pdb2pqr"

#APBS binary
#Change this to specify the location of the APBS binary if installed.
#This is used for the web interface to pdb2pqr. Provide an absolute path. Relative paths and ~ usually will not work correctly.

#APBS=""

#MAX_ATOMS
#Sets the maximum number of atoms in a protein for non-Opal job submission.
#Only affects web tools. Default is 10000.

#MAX_ATOMS=10000


#BUILD_PDB2PKA
#Set to False to skip building ligand and pdb2pka support. Requires numpy.
# Defaults to True

#BUILD_PDB2PKA=False

#DEBUG
#Set to True to build compiled extentions with debug headers.
#Defaults to False

#DEBUG=True

#CXXFLAGS
#Set to add extra CXX flags to the build.
#Defaults to ""

#EXTRA_CXXFLAGS="-fPIC"

#EXTRA_LINKFLAGS
#Set to add extra CXX flags to the build.
#Defaults to ""

#EXTRA_LINKFLAGS=""


#REBUILD_SWIG
#Set to True to rebuild the swig bindings.
# Requires swig on the the user path.
# Defaults to False

#REBUILD_SWIG=True
