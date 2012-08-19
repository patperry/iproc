#include "port.h"
#include <float.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include <setjmp.h>
#include "blas.h"
#include "cmockery.h"
#include "coreutil.h"
#include "testutil.h"
#include "xalloc.h"

#include "mlogit_glm.h"


struct mlogit_glm MGLM;
double *BETA, *MEAN;
double *X;
size_t N;
size_t P;


static void recompute()
{
	struct mlogit m;
	double *eta = xmalloc(N * sizeof(*eta));
	double *diff = xcalloc(P, sizeof(*diff));
	double p, ptot;
	size_t i;

	mlogit_init(&m, N);
	
	memset(MEAN, 0, P * sizeof(*MEAN));
	
	if (N == 0)
		goto cleanup;

	blas_dgemv(BLAS_NOTRANS, N, P, 1.0, X, N, BETA, 1, 0.0, eta, 1.0);
	mlogit_set_all_eta(&m, eta);

	
	blas_dcopy(P, X, N, MEAN, 1);
	ptot = mlogit_prob(&m, 0);
	
	for (i = 1; i < N; i++) {
		/* diff := x[i,:] - mean */
		blas_dcopy(P, X + i, N, diff, 1);
		blas_daxpy(P, -1.0, MEAN, 1, diff, 1);

		/* ptot += p */
		p = mlogit_prob(&m, i);
		ptot += p;
		
		/* mean := mean + p/ptot * diff */
		blas_daxpy(P, p/ptot, diff, 1, MEAN, 1);
	}
cleanup:
	mlogit_deinit(&m);
	free(diff);
	free(eta);
}


static void test_coefs()
{
	double *beta = mlogit_glm_coefs(&MGLM);
	size_t p = mlogit_glm_dim(&MGLM);
	size_t j;
	
	assert_int_equal(P, p);
	
	for (j = 0; j < p; j++) {
		assert_real_identical(BETA[j], beta[j]);
	}
}


static void test_x()
{
	double *x = mlogit_glm_x(&MGLM);
	size_t n = mlogit_glm_ncat(&MGLM);
	size_t p = mlogit_glm_dim(&MGLM);
	size_t i, j;
	
	assert_int_equal(N, n);
	assert_int_equal(P, p);
	
	for (j = 0; j < P; j++) {
		for (i = 0; i < N; i++) {
			assert_real_identical(MATRIX_ITEM(X, N, i, j), MATRIX_ITEM(x, n, i, j));
		}
	}
}

static void test_mean()
{
	double *mean = mlogit_glm_mean(&MGLM);
	size_t i;
	
	for (i = 0; i < P; i++) {
		assert_real_eqrel(DBL_MANT_DIG / 2, mean[i], MEAN[i]);
	}
}

static void test_inc_x(size_t i, const double *dx, const size_t *jdx, size_t ndx)
{
	size_t j, k;
	
	if (jdx) {
		for (k = 0; k < ndx; k++) {
			MATRIX_ITEM(X, N, i, jdx[k]) += dx[k];
		}
	} else if (ndx) {
		assert(ndx == P);

		for (j = 0; j < P; j++) {
			MATRIX_ITEM(X, N, i, j) += dx[k];
		}
	}
	recompute();
	
	mlogit_glm_inc_x(&MGLM, i, dx, jdx, ndx);
	
	test_x();
	test_mean();
}

static void test_inc_x_rand(double dxmin, double dxmax, size_t ndx)
{
	size_t i = rand() % N;	
	double *dx = xmalloc(ndx * sizeof(*dx));
	size_t *jdx = xmalloc(ndx * sizeof(*jdx));
	size_t k;
	
	for (k = 0; k < ndx; k++) {
		dx[k] = runif(dxmin, dxmax);
		jdx[k] = rand() % P;
	}
	
	test_inc_x(i, dx, jdx, ndx);

	free(jdx);	
	free(dx);
}


static void test_inc_x_small()
{
	test_inc_x_rand(-1.0, 1.0, 1);
}

static void test_inc_x_med()
{
	test_inc_x_rand(-10.0, 10.0, 1);
}

static void test_inc_x_big()
{
	test_inc_x_rand(-100.0, 100.0, 1);
}



static void setup(const double *beta, size_t n, size_t p)
{
	size_t i, j;

	srand(1);
	
	N = n;
	P = p;
	BETA = xmalloc(P * sizeof(*BETA));
	MEAN = xcalloc(P, sizeof(*MEAN));
	
	X = xcalloc(N * P, sizeof(*X));
	
	
	if (beta) {
		memcpy(BETA, beta, P * sizeof(*BETA));
	} else {
		memset(BETA, 0, P * sizeof(*BETA));
	}
	
	for (j = 0; j < P; j++) {
		for (i = 0; i < N; i++) {
			MATRIX_ITEM(X, N, i, j) = runif(-5.0, 5.0);
		}
	}
	
	mlogit_glm_init(&MGLM, N, P);
	_mlogit_glm_check_invariants(&MGLM);
	
	if (beta) {
		mlogit_glm_set_coefs(&MGLM, BETA);
		_mlogit_glm_check_invariants(&MGLM);
	}
	
	mlogit_glm_set_all_x(&MGLM, X);
	recompute();
}

static void teardown()
{
	mlogit_glm_deinit(&MGLM);

	free(X);
	X = NULL;
	
	free(MEAN);
	MEAN = NULL;

	free(BETA);
	BETA = NULL;

	N = 0;
	P = 0;
}


static void empty_setup_fixture()
{
	setup_fixture("Empty");
}

static void empty_setup()
{
	setup(NULL, 50, 0);
}


static void simple_setup_fixture()
{
	setup_fixture("Simple");
}

static void simple_setup()
{
	setup(NULL, 50, 5);
}


static void zeros_setup_fixture()
{
	setup_fixture("Zeros");
}

static void zeros_setup()
{
	setup(NULL, 50, 5);
}


int main()
{
	UnitTest tests[] = {
		unit_test_setup(empty_suite, empty_setup_fixture),
		unit_test_setup_teardown(test_x, empty_setup, teardown),	
		unit_test_setup_teardown(test_coefs, empty_setup, teardown),	
		unit_test_setup_teardown(test_mean, empty_setup, teardown),		
		unit_test_teardown(empty_suite, teardown_fixture),
		
		unit_test_setup(simple_suite, simple_setup_fixture),
		unit_test_setup_teardown(test_x, simple_setup, teardown),
		unit_test_setup_teardown(test_coefs, simple_setup, teardown),
		unit_test_setup_teardown(test_mean, simple_setup, teardown),		
		unit_test_setup_teardown(test_inc_x_small, simple_setup, teardown),
		unit_test_setup_teardown(test_inc_x_med, simple_setup, teardown),
		unit_test_setup_teardown(test_inc_x_big, simple_setup, teardown),		
		unit_test_teardown(simple_suite, teardown_fixture),
		
		unit_test_setup(zeros_suite, zeros_setup_fixture),
		unit_test_setup_teardown(test_x, zeros_setup, teardown),
		unit_test_setup_teardown(test_coefs, zeros_setup, teardown),		
		unit_test_setup_teardown(test_mean, zeros_setup, teardown),		
		unit_test_setup_teardown(test_inc_x_small, zeros_setup, teardown),
		unit_test_setup_teardown(test_inc_x_med, zeros_setup, teardown),
		unit_test_setup_teardown(test_inc_x_big, zeros_setup, teardown),		
		unit_test_teardown(zeros_suite, teardown_fixture),
	};
	
	return run_tests(tests);
}
