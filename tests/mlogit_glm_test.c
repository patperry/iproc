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


static void test_mean()
{
	double *mean = mlogit_glm_mean(&MGLM);
	size_t i;
	
	for (i = 0; i < P; i++) {
		assert_real_eqrel(DBL_MANT_DIG / 2, mean[i], MEAN[i]);
	}
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
	
	mlogit_glm_set_x(&MGLM, X);
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
		unit_test_setup_teardown(test_mean, empty_setup, teardown),
		unit_test_teardown(empty_suite, teardown_fixture),
		
		unit_test_setup(simple_suite, simple_setup_fixture),
		unit_test_setup_teardown(test_mean, simple_setup, teardown),
		unit_test_teardown(simple_suite, teardown_fixture),
		
		unit_test_setup(zeros_suite, zeros_setup_fixture),
		unit_test_setup_teardown(test_mean, zeros_setup, teardown),		
		unit_test_teardown(zeros_suite, teardown_fixture),
	};
	
	return run_tests(tests);
}
