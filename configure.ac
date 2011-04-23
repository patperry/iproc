AC_INIT([iproc], [0.1], [patperry@gmail.com])
AC_SUBST([PACKAGE_DATE], [2011-03-17])
AC_SUBST([RPACKAGE], AC_PACKAGE_NAME[_]AC_PACKAGE_VERSION[.tar.gz])

AC_CONFIG_SRCDIR([src/actors.h])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([m4])

AM_INIT_AUTOMAKE([-Wall -Werror foreign subdir-objects])
m4_ifdef([AM_SILENT_RULES], [AM_SILENT_RULES([yes])])

# Checks for programs
AC_PROG_CC_C99
AC_PROG_LIBTOOL

# Checks for libraries
AC_CHECK_LIBM
AX_CHECK_BLAS
AX_CHECK_BLAS64
AX_CHECK_LAPACK

# Checks for header files
AC_HEADER_ASSERT
AC_HEADER_STDBOOL

# Checks for typedefs, structures, and compiler characteristics
AC_C_BIGENDIAN
AC_C_INLINE
AX_C___ALIGNOF__

AC_TYPE_INT8_T
AC_TYPE_INT16_T
AC_TYPE_INT32_T
AC_TYPE_INT64_T
AC_TYPE_INTMAX_T
AC_TYPE_INTPTR_T
AC_TYPE_SSIZE_T
AC_TYPE_UINT8_T
AC_TYPE_UINT16_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T

AC_CHECK_SIZEOF([int])
AC_CHECK_SIZEOF([intmax_t])
AC_CHECK_SIZEOF([long])
AC_CHECK_SIZEOF([long long])
AC_CHECK_SIZEOF([size_t])
AC_CHECK_SIZEOF([void *])

# Checks for library functions
AC_CONFIG_LIBOBJ_DIR([lib])
AC_REPLACE_FUNCS([imaxabs])

# Checks for system services

AC_CONFIG_HEADERS([config.h])
AC_CONFIG_FILES([Makefile interfaces/R/DESCRIPTION])
AC_OUTPUT
