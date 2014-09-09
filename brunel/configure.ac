# Process this file with autoconf to produce a configure script.
AC_INIT([brunel],
        m4_esyscmd([build-aux/git-version-gen .tarball-version]),
        [mp15@sanger.ac.uk])
AC_PREREQ([2.69])

# Check to make sure that the src dir actually exists (sanity check)
AC_CONFIG_SRCDIR([src/main.c])

# Configure aux dirs
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([m4])
AC_CONFIG_MACRO_DIRS([m4 gl/m4])

# Enable extensions to C or Posix (which gnulib will provide on other systems -- this must be called before compiler is run)
AC_USE_SYSTEM_EXTENSIONS

# Initialize libtool
LT_INIT([])
LT_PREREQ([2.4.2])

# Initialize automake
AC_PROG_MAKE_SET
AM_INIT_AUTOMAKE([foreign 1.13 -Wall -Werror dist-bzip2])

# Check for C compiler
AC_PROG_CC
gl_EARLY
#AC_PROG_CC_C99

# Configure automake to use modulename_CFLAGS, modulename_LDFLAGS, etc.
AM_PROG_CC_C_O

# Test compiler functionality
AC_LANG_C

# Make sure we have install program
AC_PROG_INSTALL

# Checks for stdlib.h, stdarg.h , string.h and float.h, defines STDC_HEADERS on success
AC_HEADER_STDC

# Bring in config headers
AC_CONFIG_HEADERS([config.h])

# Initialise and check gnulib modules (note: gl_EARLY must also be called before this and immediately after AC_PROG_CC)
gl_INIT

# Checks for zlib and adds -lz to LIBS and defined HAVE_LIBZ
AC_ARG_VAR([ZLIB_CFLAGS],[C compiler flags for ZLIB])
AC_ARG_VAR([ZLIB_LDFLAGS],[linker flags for ZLIB])
AC_MSG_CHECKING([for zlib])
AC_CHECK_LIB([z], [zlibVersion], [], [AC_MSG_FAILURE([zlib is required but check for zlibVersion function failed! (is ZLIB_LDFLAGS set correctly?)])], [${ZLIB_LDFLAGS}])

# Check for htslib (which requires zlib)
AC_ARG_VAR([HTSLIB_CFLAGS],[C compiler flags for HTSLIB])
AC_ARG_VAR([HTSLIB_LDFLAGS],[linker flags for HTSLIB])
AC_MSG_CHECKING([for htslib])
AC_CHECK_LIB([hts], [hts_open], [], [AC_MSG_FAILURE([htslib is required but check for hts_open function failed! (is HTSLIB_LDFLAGS set correctly?)])], [${LIBS} ${HTSLIB_LDFLAGS} ${ZLIB_LDFLAGS} ${LTLIBMULTITHREAD}])

# Which files to configure 
AC_CONFIG_FILES([
 Makefile 
 src/Makefile 
 gl/Makefile
])

# Generate all config_files
AC_OUTPUT
