## ###########################################################################
## File:    Makefile.am
##
## Purpose: Automatically generate Makefile.in for autoconf using automake.
##
## Notes:   The final makefile produced by autoconf will use the following:
##
##     COMPILE = $(CC) $(DEFS) $(INCLUDES) $(CPPFLAGS) $(CFLAGS)
##     LINK = $(CC) $(CFLAGS) $(LDFLAGS) -o $@
##
## Author:  Nathan Baker
## ###########################################################################
AUTOMAKE_OPTIONS = no-dependencies
SUFFIXES = .H .h .C .cc .c .f .o
SUBDIRS = propka
distclean-local: clean-local
clean-local: clean-ac
clean-ac:
	-rm -f dummy*.o mangle-func* noformain* work.pc*
	-rm -rf autom4te.cache
