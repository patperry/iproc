AC_INIT([iproc], [0.1], [patperry@gmail.com])
AC_SUBST([PACKAGE_DATE], [2011-03-17])
AC_SUBST([RPACKAGE], AC_PACKAGE_NAME[_]AC_PACKAGE_VERSION[.tar.gz])

AC_CONFIG_SRCDIR([src/actors.h])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([m4])

AM_INIT_AUTOMAKE([-Wall -Werror foreign subdir-objects])
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])

# Checks for programs
AC_PROG_CC
AC_PROG_LIBTOOL

# Checks for libraries
AC_CHECK_LIBM
AX_CHECK_BLAS
AX_CHECK_LAPACK

# Checks for structures
AC_HEADER_STDBOOL
AC_TYPE_SSIZE_T
AC_CHECK_SIZEOF([size_t], , [#include <stddef.h> #include <stdio.h>])

## HACK: not portable
AH_TEMPLATE([f77int], [Define to appropriate type for Fortran integer])
AH_TEMPLATE([F77_INT_MAX], [Define to maximum Fortran integer])
AH_TOP([#include <limits.h>])
AC_DEFINE([f77int],      [long int])
AC_DEFINE([F77_INT_MAX], [LONG_MAX])

# Checks for compiler characteristics
# Checks for library functions
# Checks for system services

AC_CONFIG_HEADERS([src/config.h])
AC_CONFIG_FILES([Makefile interfaces/R/DESCRIPTION])
AC_OUTPUT
