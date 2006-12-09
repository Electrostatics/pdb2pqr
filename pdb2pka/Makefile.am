## ###########################################################################
## File:    Makefile.am
##
## Purpose: Automatically generate Makefile.in for autoconf using automake.
##
## Notes:   This is the driver for the pdb2pka directory.
##
## Author:  Todd Dolinsky
## ###########################################################################
AUTOMAKE_OPTIONS = no-dependencies
SUFFIXES         = .H .h .C .cc .c .f .o .cpp .so
AM_CFLAGS        = $(PY_CFLAGS) $(PY_FLAGS)
CCLD             = $(CPP)
CLEANFILES       = _pMCmult.so 
#SUBDIRS = $(PDB2PQR_SUBDIRS)

noinst_PROGRAMS         = _pMCmult.so
_pMCmult_so_SOURCES     = pMC_mult.cpp pMC_mult_wrap.cpp  
_pMCmult_so_LINK        = $(CPP)
_pMCmult_so_LDFLAGS     = -o $@ $(PY_SHARED)
_pMCmult_so_LDADD       = $(PY_LDFLAGS)

distclean-local: clean-local
clean-local: clean-ac
clean-ac:
	-rm -rf autom4te.cache
