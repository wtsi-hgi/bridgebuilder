#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.61])
AC_INIT([brunel], [1.0], [mp15@sanger.ac.uk])
AM_INIT_AUTOMAKE([-Wall foreign])
AC_CONFIG_SRCDIR([src/main.c])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CC
AC_PROG_CC_C99

# Checks for libraries.
LDFLAGS="-pthread"
AC_CHECK_LIB([hts], [hts_open])
if test "$ac_cl_lib_hts" == no
then
	AC_MSG_WARN([No htslib found])
fi
AC_CHECK_LIB([z], [zlib])

# Checks for header files.
CFLAGS="-DPOSIX_SOURCE -g -I../../htslib -DPOSIX_SOURCE"
LDFLAGS+=" -L../../htslib/htslib -lz -lhts"
LIBS="-lz -lhts"

# Checks for typedefs, structures, and compiler characteristics.

# Checks for library functions.

AC_CONFIG_FILES([Makefile src/Makefile])
AC_OUTPUT
