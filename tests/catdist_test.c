#include "port.h"

#include <stdlib.h>		// free
#include <string.h>		// memcpy, memset

#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include "cmockery.h"
#include "testutil.h"

#include <math.h>		// INFINITY, isnan
#include <float.h>		// DBL_MANT_DIG
#include <stdio.h>		// fflush, stdout
#include "ieee754.h"		// double_compare
#include "xalloc.h"		// xmalloc
#include "catdist.h"

#define ARRAY_SIZE(x) (sizeof(x) / sizeof((x)[0]))

struct catdist CATDIST;
double *ETA;
double PSI;
size_t N;

static double get_psi(const double *eta, size_t n)
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

	for (i = 0; i < n - 1; i++) {
		sum += exp(eta_sort[i] - etamax);
	}

	double psi = etamax + log1p(sum);

	assert(isfinite(psi));

	free(eta_sort);

	return psi;
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

	PSI = get_psi(eta, n);

	catdist_init(&CATDIST, N);
	assert_false(_catdist_check(&CATDIST));

	if (eta) {
		catdist_set_all_eta(&CATDIST, ETA);
		assert_false(_catdist_check(&CATDIST));
	}
}

static void teardown()
{
	catdist_deinit(&CATDIST);
	memset(&CATDIST, 0, sizeof(CATDIST));
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
	static double eta[] = { -0.9, 0.8, 0.9, 1.4, 1.2 };
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
	assert_int_equal(catdist_ncat(&CATDIST), N);
}

static void test_eta()
{
	size_t i;

	for (i = 0; i < N; i++) {
		assert_real_identical(catdist_eta(&CATDIST, i), ETA[i]);
	}
}

static void test_psi()
{
	assert_real_identical(catdist_psi(&CATDIST), PSI);
}

static void test_lprob()
{
	size_t i;

	for (i = 0; i < N; i++) {
		assert_real_identical(catdist_lprob(&CATDIST, i), ETA[i] - PSI);
	}
}

static void test_prob()
{
	size_t i;

	for (i = 0; i < N; i++) {
		assert_real_identical(catdist_prob(&CATDIST, i),
				      exp(ETA[i] - PSI));
	}
}

static void test_set_eta(size_t i, double eta)
{
	assert(i < N);
	assert(!isnan(eta));

	ETA[i] = eta;
	PSI = get_psi(ETA, N);

	catdist_set_eta(&CATDIST, i, eta);
	assert_false(_catdist_check(&CATDIST));

	test_eta();
	test_psi();
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
		PSI = get_psi(ETA, N);

		catdist_set_eta(&CATDIST, i, eta);
		assert_false(_catdist_check(&CATDIST));
		assert_real_eqrel(DBL_MANT_DIG / 2, catdist_psi(&CATDIST), PSI);
		//print_message(".");
		fflush(stdout);
	}
}

static void test_many_set_eta_small()
{
	test_many_set_eta(1000, -1, 1);
}

static void test_many_set_eta_med()
{
	test_many_set_eta(1000, -10, 10);
}

static void test_many_set_eta_big()
{
	test_many_set_eta(1000, -1000, 1000);
}

int main()
{
	UnitTest tests[] = {
		unit_test_setup(empty_suite, empty_setup_fixture),
		unit_test_setup_teardown(test_ncat, empty_setup, teardown),
		unit_test_setup_teardown(test_eta, empty_setup, teardown),
		unit_test_setup_teardown(test_psi, empty_setup, teardown),
		unit_test_setup_teardown(test_lprob, empty_setup, teardown),
		unit_test_setup_teardown(test_prob, empty_setup, teardown),
		unit_test_teardown(empty_suite, teardown_fixture),

		unit_test_setup(simple_suite, simple_setup_fixture),
		unit_test_setup_teardown(test_ncat, simple_setup, teardown),
		unit_test_setup_teardown(test_eta, simple_setup, teardown),
		unit_test_setup_teardown(test_psi, simple_setup, teardown),
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
		unit_test_setup_teardown(test_psi, zeros_setup, teardown),
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
