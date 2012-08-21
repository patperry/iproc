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
#include "mlogit1.h"


#define ARRAY_SIZE(x) (sizeof(x) / sizeof((x)[0]))

static struct mlogit MLOGIT;
static struct mlogit1 MLOGIT1;
static double *ETA1;
static size_t *IND;
static size_t NZ;


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

static void setup(const double *eta0, size_t n)
{
	srand(1);


	mlogit_init(&MLOGIT, n);
	mlogit_set_all_eta(&MLOGIT, eta0);

	mlogit1_init(&MLOGIT1, &MLOGIT);
	assert_false(_mlogit1_check(&MLOGIT1));

	ETA1 = xmalloc(n * sizeof(*ETA1));
	IND = xmalloc(n * sizeof(*IND));
	NZ = 0;
}

static void teardown()
{
	free(IND);
	free(ETA1);
	mlogit1_deinit(&MLOGIT1);
	mlogit_deinit(&MLOGIT);
	memset(&MLOGIT, 0, sizeof(MLOGIT));
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
	assert_int_equal(mlogit1_ncat(&MLOGIT1), mlogit_ncat(&MLOGIT));
}

static void test_eta()
{
	size_t n = mlogit1_ncat(&MLOGIT1);
	double *eta = xmalloc(n * sizeof(*eta));
	size_t i, iz;

	for (i = 0; i < n; i++) {
		eta[i] = mlogit_eta(&MLOGIT, i);
	}

	for (iz = 0; iz < NZ; iz++) {
		eta[IND[iz]] = ETA1[iz];
	}

	for (i = 0; i < n; i++) {
		assert_real_identical(mlogit1_eta(&MLOGIT1, i), eta[i]);
	}

	free(eta);
}

static void test_psi()
{
	size_t n = mlogit1_ncat(&MLOGIT1);
	double *eta = xmalloc(n * sizeof(*eta));
	size_t i;

	for (i = 0; i < n; i++) {
		eta[i] = mlogit1_eta(&MLOGIT1, i);
	}

	double psi = get_psi(eta, n);
	assert_real_approx(mlogit1_psi(&MLOGIT1), psi);

	free(eta);
}

static void test_lprob()
{
	size_t n = mlogit1_ncat(&MLOGIT1);
	double psi = mlogit1_psi(&MLOGIT1);
	size_t i;
	
	for (i = 0; i < n; i++) {
		double eta = mlogit1_eta(&MLOGIT1, i);
		assert_real_approx(mlogit1_lprob(&MLOGIT1, i), eta - psi);
	}
}

static void test_prob()
{
	size_t n = mlogit1_ncat(&MLOGIT1);
	size_t i;
	
	for (i = 0; i < n; i++) {
		double lp = mlogit1_lprob(&MLOGIT1, i);
		assert_real_identical(mlogit1_prob(&MLOGIT1, i), exp(lp));
	}
}

static void test_set_eta(size_t i, double eta)
{
	assert(i < mlogit1_ncat(&MLOGIT1));
	assert(!isnan(eta));

	size_t iz;
	for (iz = 0; iz < NZ; iz++) {
		if (IND[iz] == i)
			break;
	}
	if (iz == NZ)
		IND[NZ++] = iz;

	ETA1[iz] = eta;

	mlogit1_set_eta(&MLOGIT1, i, eta);
	assert_false(_mlogit1_check(&MLOGIT1));

	test_ncat();
	test_eta();
}

static void test_set_eta_small()
{
	size_t n = mlogit1_ncat(&MLOGIT1);
	size_t i = rand() % n;
	double eta = runif(-1, 1);
	test_set_eta(i, eta);
}

static void test_set_eta_med()
{
	size_t n = mlogit1_ncat(&MLOGIT1);
	size_t i = rand() % n;
	double eta = runif(-10, 10);
	test_set_eta(i, eta);
}

static void test_set_eta_big()
{
	size_t n = mlogit1_ncat(&MLOGIT1);
	size_t i = rand() % n;
	double eta = runif(-1000, 1000);
	test_set_eta(i, eta);
}

static void test_many_set_eta(size_t nrep, double min, double max)
{
	size_t n = mlogit1_ncat(&MLOGIT1);
	size_t rep;
	
	for (rep = 0; rep < nrep; rep++) {
		size_t i = rand() % n;
		double eta = runif(min, max);

		test_set_eta(i, eta);
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
