AC_INIT([libla], [0.1], [patperry@gmail.com])

AC_CONFIG_SRCDIR([libla/index.h])
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

# Checks for headers
AC_HEADER_ASSERT

# Checks for types
AC_TYPE_SIZE_T

# Checks for structures
# Checks for compiler characteristics
# Checks for library functions
# Checks for system services

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile])
AC_OUTPUT
