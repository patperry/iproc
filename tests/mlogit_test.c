#include "port.h"

#include <stdlib.h> // free
#include <string.h> // memcpy, memset

#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include "cmockery.h"
#include "testutil.h"

#include <math.h>  // INFINITY, isnan
#include <float.h> // DBL_MANT_DIG
#include <stdio.h> // fflush, stdout
#include "ieee754.h" // double_compare
#include "xalloc.h" // xmalloc
#include "mlogit.h"


#define ARRAY_SIZE(x) (sizeof(x) / sizeof((x)[0]))

struct mlogit MLOGIT;
double *ETA;
double PHI;
size_t N;

static double get_phi(const double *eta, size_t n)
{
	if (n == 0)
		return -INFINITY;
	if (!eta)
		return log(n);
	
	double *eta_sort = xmalloc(n * sizeof(*eta_sort));
	memcpy(eta_sort, eta, n * sizeof(*eta_sort));
	qsort(eta_sort, n, sizeof(double), double_compare);

	double etamax = eta_sort[n - 1];
	double sum = 0.0;
	size_t i;
	
	for (i = 0; i < N - 1; i++) {
		sum += exp(eta_sort[i] - etamax);
	}
	
	double phi = etamax + log1p(sum);
	
	assert(isfinite(phi));
	
	free(eta_sort);
	
	return phi;
}

static void setup(const double *eta, size_t n)
{
	srand(1);
	
	N = n;
	ETA = xmalloc(N * sizeof(ETA[0]));

	if (eta) {
		memcpy(ETA, eta, N * sizeof(ETA[0]));
	} else {
		memset(ETA, 0, N * sizeof(ETA[0]));
	}
	
	PHI = get_phi(eta, n);

	mlogit_init(&MLOGIT, N);
	_mlogit_check_invariants(&MLOGIT);

	if (eta) {
		mlogit_set_all_eta(&MLOGIT, ETA);
		_mlogit_check_invariants(&MLOGIT);
	}
}

static void teardown()
{
	mlogit_deinit(&MLOGIT);
	memset(&MLOGIT, 0, sizeof(MLOGIT));
	free(ETA);
	ETA = NULL;
	N = 0;
}

static void zeros_setup_fixture()
{
	setup_fixture("Zeros");
}

static void zeros_setup()
{
	setup(NULL, 100);
}

static void simple_setup_fixture()
{
	setup_fixture("Simple");
}

static void simple_setup()
{
	static double eta[] = { -0.9,  0.8,  0.9,  1.4,  1.2 };
	setup(eta, ARRAY_SIZE(eta));
}

static void empty_setup_fixture()
{
	setup_fixture("Empty");
}

static void empty_setup()
{
	setup(NULL, 0);
}

static void test_ncat()
{
	assert_int_equal(mlogit_ncat(&MLOGIT), N);
}

static void test_eta()
{
	size_t i;

	for (i = 0; i < N; i++) {
		assert_real_identical(mlogit_eta(&MLOGIT, i), ETA[i]);
	}
}

static void test_phi()
{
	assert_real_identical(mlogit_phi(&MLOGIT), PHI);	
}

static void test_lprob()
{
	size_t i;
	
	for (i = 0; i < N; i++) {
		assert_real_identical(mlogit_lprob(&MLOGIT, i), ETA[i] - PHI);
	}
}

static void test_prob()
{
	size_t i;
	
	for (i = 0; i < N; i++) {
		assert_real_identical(mlogit_prob(&MLOGIT, i), exp(ETA[i] - PHI));
	}
}

static void test_set_eta(size_t i, double eta)
{
	assert(i < N);
	assert(!isnan(eta));

	ETA[i] = eta;
	PHI = get_phi(ETA, N);
	
	mlogit_set_eta(&MLOGIT, i, eta);
	_mlogit_check_invariants(&MLOGIT);
	
	test_eta();
	test_phi();
}

static void test_set_eta_small()
{
	size_t i = rand() % N;
	double eta = runif(-1, 1);
	test_set_eta(i, eta);
}

static void test_set_eta_med()
{
	size_t i = rand() % N;
	double eta = runif(-10, 10);
	test_set_eta(i, eta);
}

static void test_set_eta_big()
{
	size_t i = rand() % N;
	double eta = runif(-1000, 1000);
	test_set_eta(i, eta);
}

static void test_many_set_eta(size_t nrep, double min, double max)
{
	size_t rep;
	
	for (rep = 0; rep < nrep; rep++) {
		size_t i = rand() % N;
		double eta = runif(min, max);
		
		ETA[i] = eta;
		PHI = get_phi(ETA, N);
		
		mlogit_set_eta(&MLOGIT, i, eta);
		_mlogit_check_invariants(&MLOGIT);
		//assert_real_eqrel(DBL_MANT_DIG, mlogit_phi(&MLOGIT), PHI);
		print_message(".");
		fflush(stdout);
	}
}

static void test_many_set_eta_small()
{
	test_many_set_eta(10000, -1, 1);
}

static void test_many_set_eta_med()
{
	test_many_set_eta(10000, -10, 10);
}

static void test_many_set_eta_big()
{
	test_many_set_eta(10000, -1000, 1000);
}


int main()
{
	UnitTest tests[] = {
		unit_test_setup(empty_suite, empty_setup_fixture),
		unit_test_setup_teardown(test_ncat, empty_setup, teardown),
		unit_test_setup_teardown(test_eta, empty_setup, teardown),
		unit_test_setup_teardown(test_phi, empty_setup, teardown),
		unit_test_setup_teardown(test_lprob, empty_setup, teardown),
		unit_test_setup_teardown(test_prob, empty_setup, teardown),
		unit_test_teardown(empty_suite, teardown_fixture),

		unit_test_setup(simple_suite, simple_setup_fixture),
		unit_test_setup_teardown(test_ncat, simple_setup, teardown),
		unit_test_setup_teardown(test_eta, simple_setup, teardown),
		unit_test_setup_teardown(test_phi, simple_setup, teardown),
		unit_test_setup_teardown(test_lprob, simple_setup, teardown),
		unit_test_setup_teardown(test_prob, simple_setup, teardown),
		unit_test_setup_teardown(test_set_eta_small, simple_setup, teardown),
		unit_test_setup_teardown(test_many_set_eta_small, simple_setup, teardown),
		unit_test_setup_teardown(test_set_eta_med, simple_setup, teardown),
		unit_test_setup_teardown(test_many_set_eta_med, simple_setup, teardown),
		unit_test_setup_teardown(test_set_eta_big, simple_setup, teardown),
		unit_test_setup_teardown(test_many_set_eta_big, simple_setup, teardown),
		unit_test_teardown(simple_suite, teardown_fixture),

		unit_test_setup(zeros_suite, zeros_setup_fixture),
		unit_test_setup_teardown(test_ncat, zeros_setup, teardown),
		unit_test_setup_teardown(test_eta, zeros_setup, teardown),
		unit_test_setup_teardown(test_phi, zeros_setup, teardown),
		unit_test_setup_teardown(test_lprob, zeros_setup, teardown),
		unit_test_setup_teardown(test_prob, zeros_setup, teardown),
		unit_test_setup_teardown(test_set_eta_small, zeros_setup, teardown),
		unit_test_setup_teardown(test_many_set_eta_small, zeros_setup, teardown),
		unit_test_setup_teardown(test_set_eta_med, zeros_setup, teardown),
		unit_test_setup_teardown(test_many_set_eta_med, zeros_setup, teardown),
		unit_test_setup_teardown(test_set_eta_big, zeros_setup, teardown),
		unit_test_setup_teardown(test_many_set_eta_big, zeros_setup, teardown),
		unit_test_teardown(zeros_suite, teardown_fixture),
	};

	return run_tests(tests);
}
