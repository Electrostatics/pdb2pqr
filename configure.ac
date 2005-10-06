dnl Some PDB2PQR info

dnl ##########################################################################
dnl File: configure.in
dnl Purpose: configure.in script for PDB2PQR
dnl Note:  For use with autoconf-2.56 and above
dnl Author:  Todd Dolinsky and Nathan Baker
dnl ##########################################################################

dnl ##########################################################################
dnl 0. autoconf initialization
dnl ##########################################################################
AC_INIT([PDB2PQR], [1.0.1], [todd@ccb.wustl.edu])
AC_PREREQ(2.56)
AC_CONFIG_AUX_DIR(config)
AC_PREFIX_DEFAULT(${HOME})
AC_CANONICAL_SYSTEM

dnl ##########################################################################
dnl 1. automake initialization (uncomment only if rebuilding Makefile.in's)
dnl ##########################################################################
dnl # THIS IS THE VERSION INFORMATION FOR THE PDB2PQR LIBRARY.  ACCORDING TO 
dnl # THE LIBTOOL MANUAL:
dnl #   This flag accepts an argument of the form `current[:revision[:age]]'.
dnl #   So, passing `-version-info 3:12:1' sets current to 3, revision to 12,
dnl #   and age to 1.  If either revision or age are omitted, they default to 
dnl #   0. Also note that age must be less than or equal to the current 
dnl #   interface number.  Here are a set of rules to help you update your
dnl #   library version information:
dnl #      1. Start with version information of `0:0:0' for each libtool 
dnl #         library.
dnl #      2. Update the version information only immediately before a public
dnl #         release of your software. More frequent updates are unnecessary,
dnl #         and only guarantee that the current interface number gets larger
dnl #         faster.  
dnl #      3. If the library source code has changed at all since the last
dnl #         update, then increment revision (`c:r:a' becomes `c:r+1:a').
dnl #      4. If any interfaces have been added, removed, or changed since the
dnl #         last update, increment current, and set revision to 0.
dnl #      5. If any interfaces have been added since the last public release,
dnl #         then increment age.
dnl #      6. If any interfaces have been removed since the last public
dnl #         release, then set age to 0.
dnl #
dnl #   Never try to set the interface numbers so that they correspond to the
dnl #   release number of your package. This is an abuse that only fosters
dnl #   misunderstanding of the purpose of library versions.
AM_INIT_AUTOMAKE(PDB2PQR, 1.0.0)
PDB2PQR_VERSION="1:0:0"
AC_SUBST(PDB2PQR_VERSION)
ACLOCAL=":";
AUTOCONF=":";
AUTOMAKE=":";
AUTOHEADER=":";
AC_SUBST(ACLOCAL)
AC_SUBST(AUTOCONF)
AC_SUBST(AUTOMAKE)
AC_SUBST(AUTOHEADER)
dnl ##########################################################################
dnl 2. checks for programs
dnl ##########################################################################

AC_PATH_PROG(py_path, python)
if test -z "${py_path}"; then
   AC_MSG_ERROR([Unable to find python!])
fi

AC_FORTUNES

dnl ##########################################################################
dnl 3. Check for available subdirectories
dnl ##########################################################################

PDB2PQR_SUBDIRS="propka"
AC_MSG_CHECKING([for --disable-propka])
AC_ARG_ENABLE(propka, 
     [AC_HELP_STRING([--disable-propka],[Disable propka compilation])],
     [ with_propka="no" ],
     [ with_propka="" ]
)
if test -z "${with_propka}"; then
     AC_MSG_RESULT([no])
     AC_MSG_RESULT([])
     cd propka
     ./configure 
     cd ..
else
     AC_MSG_RESULT([yes])
     PDB2PQR_SUBDIRS=""
fi
AC_SUBST(PDB2PQR_SUBDIRS,[$PDB2PQR_SUBDIRS])

dnl ##########################################################################
dnl 11. generate the makefiles
dnl ##########################################################################

AC_MSG_RESULT([------------- End Processing PDB2PQR Configure Options -----------])
AC_MSG_RESULT([])
if test -n "${PDB2PQR_SUBDIRS}"; then
   AC_OUTPUT([
        Makefile
   ])
   AC_MSG_RESULT([Please type "make" to complete the installation.])
fi
