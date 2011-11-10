AC_INIT([iproc], [0.1], [patperry@gmail.com])
AC_SUBST([PACKAGE_DATE], [2011-03-17])
AC_SUBST([RPACKAGE], AC_PACKAGE_NAME[_]AC_PACKAGE_VERSION[.tar.gz])

AC_CONFIG_SRCDIR([src/messages.h])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([m4])

AM_INIT_AUTOMAKE([-Wall -Werror foreign subdir-objects])
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])

# Checks for programs
AC_PROG_CC
AC_PROG_RANLIB

# Checks for libraries
AC_CHECK_LIBM
AX_CHECK_BLAS
AX_CHECK_BLAS64
AX_CHECK_LAPACK

# Checks for header files
AC_HEADER_ASSERT

# Checks for typedefs, structures, and compiler characteristics
AC_C_BIGENDIAN
AC_C_INLINE
AC_TYPE_INTPTR_T

# Checks for library functions
AC_CONFIG_LIBOBJ_DIR([lib])

# Checks for system services

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile])
## AC_CONFIG_FILES([interfaces/R/DESCRIPTION])
AC_OUTPUT
