
#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include "cmockery.h"

#include <assert.h>
#include <float.h> // DBL_DIG, DBL_MANT_DIG
#include <math.h> // fabs
#include <stdlib.h> // rand
#include <string.h> // strlen
#include "ieee754.h" // double_identical
#include "lapack.h"
#include "xalloc.h"
#include "testutil.h"


#define STR(x) #x
#define STRSTR(x) STR(x)
#define LargestRealTypePrintfFormat "%.16g"

void setup_fixture(const char *name)
{
	size_t i, n = strlen(name);
	print_message("%s\n", name);
	for (i = 0; i < n; i++) {
		print_message("-");
	}
	print_message("\n");
}


void teardown_fixture()
{
	print_message("\n\n");
}

double runif(double min, double max)
{
	double u = (double)rand() / RAND_MAX;
	double x =  min + u * (max - min);
	return x;
}


static int reals_approx(const double left, const double right)
{
	const int minprec = DBL_MANT_DIG / 2;
	const int p = double_eqrel(left, right);
	const int rel_approx = p >= minprec;
	const int abs_approx = (fabs(left) <= ROOT4_DBL_EPSILON
				&& fabs(right) <= ROOT4_DBL_EPSILON);
	const int approx = rel_approx || abs_approx;

	return approx;
}


static int reals_approx_display_error(const double left, const double right)
{
	const int approx = reals_approx(left, right);

	if (!approx) {
		const int p = double_eqrel(left, right);
		print_error(LargestRealTypePrintfFormat  " != "
			    LargestRealTypePrintfFormat " (%d bits differ)\n",
			    left, right, DBL_MANT_DIG - p);
	}
	return approx;
}


static int reals_eqrel_display_error(const int precision, const double left, const double right)
{
	const int p = double_eqrel(left, right);
	const int eqrel = p >= precision;
	if (!eqrel) {
		print_error(LargestRealTypePrintfFormat  " != "
			    LargestRealTypePrintfFormat " (%d bits < %d bits)\n",
			    left, right, p, precision);
	}
	return eqrel;
	
}

static int reals_identical_display_error(const double left, const double right)
{
	const int identical = double_identical(left, right);
	if (!identical) {
		print_error(LargestRealTypePrintfFormat  " != "
				LargestRealTypePrintfFormat "\n", left, right);
	}
	return identical;

}


static int vec_approx(const double *expect, const double *actual, size_t n, size_t *idiff)
{
	int approx = 1;
	size_t i;

	for (i = 0; i < n; i++) {
		if (!reals_approx(expect[i], actual[i])) {
			approx = 0;
			*idiff = i;
			break;
		}
	}

	return approx;
}


static int vec_identical(const double *expect, const double *actual, size_t n, size_t *idiff)
{
	int identical = 1;
	size_t i;

	for (i = 0; i < n; i++) {
		if (!double_identical(expect[i], actual[i])) {
			identical = 0;
			*idiff = i;
			break;
		}
	}
	return identical;
}


static int sym_approx(const double *expect, const double *actual,
		      enum blas_uplo uplo, size_t n)
{
	if (n == 0)
		return 1;

	int approx = 1;

	double *expect_copy = xcalloc(n * (n + 1) / 2, sizeof(double));
	double *diff = xcalloc(n * (n + 1) / 2, sizeof(double));
	double *z = xmalloc(n * n * sizeof(double));
	double *w = xmalloc(n * sizeof(double));
	double *err_z = xmalloc(n * sizeof(double));
	const enum blas_uplo fuplo = uplo == BLAS_UPPER ? BLAS_LOWER : BLAS_UPPER;
	size_t lwork, liwork;
	size_t i;

	lwork = lapack_dspevd_lwork(LA_EIG_VEC, n, &liwork);
	double *work = xmalloc(lwork * sizeof(*work));
	ptrdiff_t *iwork = xmalloc(liwork * sizeof(*iwork));

	// diff := expect - actual
	blas_dcopy(n * (n + 1) / 2, expect, 1, diff, 1);
	blas_daxpy(n * (n + 1) / 2, -1.0, actual, 1, diff, 1);

	// compute expect eigendecomp
	blas_dcopy(n * (n + 1) / 2, expect, 1, expect_copy, 1);
	ptrdiff_t info = lapack_dspevd(LA_EIG_VEC, fuplo, n, expect_copy, w, z,
				       n, work, lwork, iwork, liwork);
	assert(info == 0);
	(void)info;

	// check for relative equality of expect and actual
	// on the eigenspaces of expect
	for (i = n; i > 0; i--) {
		const double *zi = z + (i - 1) * n;
		blas_dspmv(fuplo, n, 1.0, diff, zi, 1, 0.0, err_z, 1);
		double z_err_z = blas_ddot(n, zi, 1, err_z, 1);
		if (!(fabs(z_err_z) <=
		      n * ROOT3_DBL_EPSILON * (1 + fabs(w[i - 1])))) {
			approx = 0;
			goto out;
		}
	}

out:
	free(iwork);
	free(work);
	free(err_z);
	free(w);
	free(z);
	free(diff);
	free(expect_copy);

	return approx;
}


static void print_error_vector(const double *x, size_t n)
{

	print_error("{");

	if (n > 0) {
		print_error(LargestRealTypePrintfFormat, x[0]);
	}

	size_t i;
	for (i = 1; i < n; i++) {
		print_error(", "LargestRealTypePrintfFormat, x[i]);
	}

	print_error("}");
}


static int vec_approx_display_error(const double *expect, const double *actual, size_t n)
{
	size_t i;
	const int approx = vec_approx(expect, actual, n, &i);
	if (!approx) {
		print_error_vector(expect, n);
		print_error("\n!=\n");
		print_error_vector(actual, n);
		print_error("\n(index "LargestIntegralTypePrintfFormat": "
			    LargestRealTypePrintfFormat" != "LargestRealTypePrintfFormat")\n",
			    i, expect[i], actual[i]);
	}
	return approx;
}


static int vec_identical_display_error(const double *expect, const double *actual, size_t n)
{
	size_t i;
	const int identical = vec_identical(expect, actual, n, &i);
	if (!identical) {
		print_error_vector(expect, n);
		print_error("\n!=\n");
		print_error_vector(actual, n);
		print_error("\n(index "LargestIntegralTypePrintfFormat": "
			    LargestRealTypePrintfFormat" != "LargestRealTypePrintfFormat")\n",
			    i, expect[i], actual[i]);
	}
	return identical;
}


static int sym_approx_display_error(const double *expect, const double *actual,
				    enum blas_uplo uplo, size_t n)
{
	const int approx = sym_approx(expect, actual, uplo, n);
	if (!approx) {
		print_error_vector(expect, n * (n + 1) / 2);
		print_error("\n!=(sym)\n");
		print_error_vector(actual, n * (n + 1) / 2);
	}
	return approx;
}


void _assert_real_approx(const double a, const double b,
			 const char * const file, const int line)
{
	if (!reals_approx_display_error(a, b)) {
		_fail(file, line);
	}
}

void _assert_real_eqrel(int precision, const double a, const double b,
			const char * const file, const int line)
{
	if (!reals_eqrel_display_error(precision, a, b)) {
		_fail(file, line);
	}
}

void _assert_real_identical(const double a, const double b, const char * const file, const int line)
{
	if (!reals_identical_display_error(a, b)) {
		_fail(file, line);
	}
}

void _assert_vec_approx(const double *a, const double *b, size_t n,
			const char * const file, const int line)
{
	if (!vec_approx_display_error(a, b, n)) {
		_fail(file, line);
	}
}

void _assert_vec_identical(const double *a, const double *b, size_t n,
			   const char * const file, const int line)
{
	if (!vec_identical_display_error(a, b, n)) {
		_fail(file, line);
	}
}

void _assert_sym_approx(const double *a, const double *b, enum blas_uplo uplo,
			size_t n, const char * const file, const int line)
{
	if (!sym_approx_display_error(a, b, uplo, n)) {
		_fail(file, line);
	}
}

