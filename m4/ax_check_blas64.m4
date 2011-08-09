# SYNOPSIS
#
#   AX_CHECK_BLAS64
#
# DESCRIPTION
#
#   This macro checks whether BLAS supports 64-bit integers.  If so, it defines HAVE_BLAS64.
AC_DEFUN([AX_CHECK_BLAS64], [
AC_REQUIRE([AX_CHECK_BLAS])

if test "x$ax_blas_ok" == xyes; then
    AC_MSG_CHECKING([whether BLAS has 64-bit integers])
    save_LIBS="$LIBS"; LIBS="$BLAS_LIBS $LIBS $FLIBS"
    AC_RUN_IFELSE([AC_LANG_SOURCE([[
    double F77_FUNC(ddot)(int *, double *, int *, double *, int *);

    int main(void)
    {
      double a = 1.0; double b = 1.0; int nn[3] = { -1, 1, -1 }; int inc = 1;
      return (F77_FUNC(ddot)(&nn[1], &a, &inc, &b, &inc) == 1.0) ? 0 : 1;
    }
    ]])],[ac_cv_blas64=no],[ac_cv_blas64=yes],[ac_cv_blas64=no])
    AC_MSG_RESULT($ac_cv_blas64)
    LIBS="$save_LIBS"
    if test "x$ac_cv_blas64" == xyes; then
        AC_DEFINE([HAVE_BLAS64], [1], [Define if BLAS integers are 64-bit.])
    fi
else
    ac_cv_blas64=no
fi
])dnl AX_CHECK_BLAS64
