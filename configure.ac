AC_INIT([asynch], [1.3.2], [samuel-debionne@uiowa.edu])
AC_CONFIG_MACRO_DIR([m4])

AM_INIT_AUTOMAKE([foreign -Wall -Werror])

# Check for programms
AX_PROG_CC_MPI([], [], [AC_MSG_FAILURE([MPI compiler requested, but couldn't find MPI.])] )
AC_PROG_CC_STDC #C99
AM_PROG_CC_C_O  #automake < 1.14
AM_PROG_AR
AC_PROG_RANLIB

# Checks for libraries
AX_CHECK_ZLIB
AX_LIB_HDF5([serial])
AX_LIB_POSTGRESQL

AM_CONDITIONAL([HAVE_POSTGRESQL], [test "$HAVE_POSTGRESQL" -eq 1])

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile src/Makefile tests/Makefile tools/Makefile])

# Checks for packages
PKG_CHECK_MODULES([CHECK], [check], [HAVE_LIBCHECK=1], [HAVE_LIBCHECK=0])

AC_OUTPUT
