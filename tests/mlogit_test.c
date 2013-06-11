#include "port.h"
#include <float.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <setjmp.h>
#include "blas.h"
#include "cmockery.h"
#include "coreutil.h"
#include "testutil.h"
#include "xalloc.h"

#include "mlogit.h"


static struct mlogit MLOGIT;
static double *BETA, *MEAN;
static double *COV;
static enum blas_uplo COV_UPLO;
static double *X;
static size_t N;
static size_t P;


static void recompute()
{
	struct catdist c;
	double *eta = xmalloc(N * sizeof(*eta));
	double *diff = xcalloc(P, sizeof(*diff));
	enum blas_uplo uplo = COV_UPLO == BLAS_LOWER ? BLAS_UPPER : BLAS_LOWER;
	double p, ptot;
	size_t i;

	catdist_init(&c, N);

	memset(MEAN, 0, P * sizeof(*MEAN));
	memset(COV, 0, (P * (P + 1) / 2) * sizeof(*COV));

	if (P == 0)
		goto cleanup;

	blas_dgemv(BLAS_TRANS, P, N, 1.0, X, P, BETA, 1, 0.0, eta, 1.0);
	catdist_set_all_eta(&c, eta);


	blas_dcopy(P, X, 1, MEAN, 1);
	ptot = catdist_prob(&c, 0);

	for (i = 1; i < N; i++) {
		/* diff := x[i,:] - mean */
		blas_dcopy(P, X + i * P, 1, diff, 1);
		blas_daxpy(P, -1.0, MEAN, 1, diff, 1);

		/* ptot += p */
		p = catdist_prob(&c, i);
		ptot += p;

		/* mean := mean + p/ptot * diff */
		if (ptot > 0)
			blas_daxpy(P, p/ptot, diff, 1, MEAN, 1);
	}


	ptot = 0;
	for (i = 0; i < N; i++) {
		/* diff := x[i,:] - mean */
		blas_dcopy(P, X + i * P, 1, diff, 1);
		blas_daxpy(P, -1.0, MEAN, 1, diff, 1);

		/* ptot += p */
		p = catdist_prob(&c, i);
		ptot += p;

		/* COV += sqrt(n) * diff^2 */
		blas_dspr(uplo, P, p, diff, 1, COV);
	}

	/* COV /= ptot */
	blas_dscal(P * (P + 1) / 2, 1.0 / ptot, COV, 1);

cleanup:
	catdist_deinit(&c);
	free(diff);
	free(eta);
}


static void test_coefs()
{
	double *beta = mlogit_coefs(&MLOGIT);
	size_t p = mlogit_dim(&MLOGIT);
	size_t j;

	assert_int_equal(P, p);

	for (j = 0; j < p; j++) {
		assert_real_identical(BETA[j], beta[j]);
	}
}


static void test_x()
{
	double *x = mlogit_x(&MLOGIT);
	size_t n = mlogit_ncat(&MLOGIT);
	size_t p = mlogit_dim(&MLOGIT);
	size_t i, j;

	assert_int_equal(N, n);
	assert_int_equal(P, p);

	for (j = 0; j < P; j++) {
		for (i = 0; i < N; i++) {
			assert_real_identical(X[i * P + j], x[i * p + j]);
		}
	}
}

static void test_mean()
{
	double *mean = mlogit_mean(&MLOGIT);
	size_t i;

	for (i = 0; i < P; i++) {
		assert_real_approx(mean[i], MEAN[i]);
	}
}


static void test_cov()
{
	size_t dim = P * (P + 1) / 2;
	double scale;
	double *cov = mlogit_cov(&MLOGIT, &scale);
	double *cov_copy = xmalloc(dim * sizeof(double));

	memcpy(cov_copy, cov, dim * sizeof(double));
	blas_dscal(dim, 1.0 / scale, cov_copy, 1);

	//assert_false(mlogit_check(&MLOGIT));

	assert_sym_approx(cov_copy, COV, MLOGIT_COV_UPLO, P);

	/*for (i = 0; i < P * (P + 1) / 2; i++) {
		assert_real_approx(cov[i] / scale, COV[i]);
	}*/
}


static void test_inc_x(size_t i, const double *dx, const size_t *jdx, size_t ndx)
{
	size_t j, k;

	if (jdx) {
		for (k = 0; k < ndx; k++) {
			X[i * P + jdx[k]] += dx[k];
		}
	} else if (ndx) {
		assert(ndx == P);

		for (j = 0; j < P; j++) {
			X[i * P + j] += dx[j];
		}
	}
	recompute();

	mlogit_set_x(&MLOGIT, i, 0, P, X + i * P);
	//print_error("\tx...");
	test_x();
	//print_error("ok\n\tmean...");
	test_mean();
	//print_error("ok\n\tcov...");
	test_cov();
	//print_error("ok\n");
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


static void test_many_inc_x_rand(size_t nrep, double dxmin, double dxmax, size_t ndx)
{
	size_t rep;

	for (rep = 0; rep < nrep; rep++) {
		test_inc_x_rand(dxmin, dxmax, ndx);
		//print_message("."); fflush(stdout);

	}
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

static void test_many_inc_x_small()
{
	test_many_inc_x_rand(1000, -1.0, 1.0, 1);
}

static void test_many_inc_x_med()
{
	test_many_inc_x_rand(1000, -10.0, 10.0, 1);
}

static void test_many_inc_x_big()
{
	test_many_inc_x_rand(1000, -100.0, 100.0, 1);
}


static void setup(const double *beta, size_t n, size_t p)
{
	size_t i, j;

	srand(1);

	N = n;
	P = p;
	X = xcalloc(N * P, sizeof(*X));
	BETA = xmalloc(P * sizeof(*BETA));
	MEAN = xcalloc(P, sizeof(*MEAN));
	COV = xcalloc(P * (P + 1) / 2, sizeof(*COV));
	COV_UPLO = MLOGIT_COV_UPLO;

	if (beta) {
		memcpy(BETA, beta, P * sizeof(*BETA));
	} else {
		memset(BETA, 0, P * sizeof(*BETA));
	}

	for (i = 0; i < N; i++) {
		for (j = 0; j < P; j++) {
			X[i * P + j] = runif(-5.0, 5.0);
		}
	}

	mlogit_init(&MLOGIT, N, P, NULL);
	assert_false(mlogit_check(&MLOGIT));

	if (beta) {
		mlogit_set_coefs(&MLOGIT, BETA);
		assert_false(mlogit_check(&MLOGIT));
	}

	mlogit_set_all_x(&MLOGIT, X);
	recompute();
}

static void teardown()
{
	mlogit_deinit(&MLOGIT);

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
	srand(2);
	size_t ncat = 50;
	size_t dim = 5;
	double *beta = xmalloc(dim * sizeof(*beta));
	size_t i;

	for (i = 0; i < dim; i++) {
		beta[i] = runif(-2, 2);
	}

	setup(beta, ncat, dim);

	free(beta);
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
		unit_test_setup_teardown(test_cov, empty_setup, teardown),
		unit_test_teardown(empty_suite, teardown_fixture),

		unit_test_setup(simple_suite, simple_setup_fixture),
		unit_test_setup_teardown(test_x, simple_setup, teardown),
		unit_test_setup_teardown(test_coefs, simple_setup, teardown),
		unit_test_setup_teardown(test_mean, simple_setup, teardown),
		unit_test_setup_teardown(test_cov, simple_setup, teardown),
		unit_test_setup_teardown(test_inc_x_small, simple_setup, teardown),
		unit_test_setup_teardown(test_inc_x_med, simple_setup, teardown),
		unit_test_setup_teardown(test_inc_x_big, simple_setup, teardown),
		unit_test_setup_teardown(test_many_inc_x_small, simple_setup, teardown),
		unit_test_setup_teardown(test_many_inc_x_med, simple_setup, teardown),
		unit_test_setup_teardown(test_many_inc_x_big, simple_setup, teardown),
		unit_test_teardown(simple_suite, teardown_fixture),

		unit_test_setup(zeros_suite, zeros_setup_fixture),
		unit_test_setup_teardown(test_x, zeros_setup, teardown),
		unit_test_setup_teardown(test_coefs, zeros_setup, teardown),
		unit_test_setup_teardown(test_mean, zeros_setup, teardown),
		unit_test_setup_teardown(test_cov, zeros_setup, teardown),
		unit_test_setup_teardown(test_inc_x_small, zeros_setup, teardown),
		unit_test_setup_teardown(test_inc_x_med, zeros_setup, teardown),
		unit_test_setup_teardown(test_inc_x_big, zeros_setup, teardown),
		unit_test_setup_teardown(test_many_inc_x_small, zeros_setup, teardown),
		unit_test_setup_teardown(test_many_inc_x_med, zeros_setup, teardown),
		unit_test_setup_teardown(test_many_inc_x_big, zeros_setup, teardown),
		unit_test_teardown(zeros_suite, teardown_fixture),
	};

	return run_tests(tests);
}
