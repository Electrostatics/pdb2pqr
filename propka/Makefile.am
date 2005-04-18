## ###########################################################################
## File:    Makefile.am
##
## Purpose: Automatically generate Makefile.in for autoconf using automake.
##
## Author:  Todd Dolinsky
## ###########################################################################

AUTOMAKE_OPTIONS = no-dependencies
SUFFIXES         = .c .C .o .so
INCLUDES         = $(PY_CFLAGS)
CCLD             = $(CC)
CLEANFILES       = propkalibc.so

noinst_PROGRAMS       = propkalibc.so
propkalibc_so_SOURCES   = propkalib.c propka.f propkac.c
propkalibc_so_LINK      = $(F77)
propkalibc_so_LDFLAGS   = -o $@ $(PY_SHARED)
propkalibc_so_LDADD     = $(PY_LDFLAGS)

