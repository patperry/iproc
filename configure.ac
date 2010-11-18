AC_INIT([libla], [0.1], [patperry@gmail.com])

AC_CONFIG_SRCDIR([libla/vector.h])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([m4])

AM_INIT_AUTOMAKE([-Wall -Werror foreign subdir-objects])
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])

# Checks for programs
AC_PROG_CC
AC_PROG_LIBTOOL
PKG_PROG_PKG_CONFIG

# Checks for libraries
PKG_CHECK_MODULES([GLIB], [glib-2.0 >= 2.14.0])
AC_CHECK_LIBM
AX_CHECK_BLAS
PKG_CHECK_MODULES([CHECK], [check >= 0.9.8])

# Checks for structures
## HACK: not portable
AH_TEMPLATE([f77int], [Define to appropriate type for Fortran integer])
AH_TEMPLATE([F77_INT_MAX], [Define to maximum Fortran integer])
AH_TEMPLATE([F77_INTP],
            [Define as a macro for getting an f77int pointer from an int64_t])
AH_TOP([#include <stdint.h>])
AC_DEFINE([f77int],      [int32_t])
AC_DEFINE([F77_INT_MAX], [INT32_MAX])
AC_DEFINE([F77_INTP(x)], [((int32_t *)(&x))])

# Checks for compiler characteristics
# Checks for library functions
# Checks for system services

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile])
AC_OUTPUT
