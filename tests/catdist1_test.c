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
#include "catdist.h"
#include "catdist1.h"


#define ARRAY_SIZE(x) (sizeof(x) / sizeof((x)[0]))

static struct catdist CATDIST;
static struct catdist1 CATDIST1;
static double *DETA;
static size_t *IND;
static size_t NZ;


static void clear()
{
	catdist1_set_all_deta(&CATDIST1, NULL, NULL, 0);
	NZ = 0;
}

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


	catdist_init(&CATDIST, n);
	catdist_set_all_eta(&CATDIST, eta0);

	catdist1_init(&CATDIST1, &CATDIST);
	assert_false(catdist1_check(&CATDIST1));

	DETA = xmalloc(n * sizeof(*DETA));
	IND = xmalloc(n * sizeof(*IND));
	NZ = 0;
}

static void teardown()
{
	free(IND);
	free(DETA);
	catdist1_deinit(&CATDIST1);
	catdist_deinit(&CATDIST);
	memset(&CATDIST, 0, sizeof(CATDIST));
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

static void hard_setup_fixture()
{
	setup_fixture("Hard");
}

static void hard_setup()
{
	size_t i, n = 100;
	double eta[n];

	srand(31337);

	for (i = 0; i < 100; i++) {
		eta[i] = runif(-10, 10);
	}

	setup(eta, n);
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
	assert_int_equal(catdist1_ncat(&CATDIST1), catdist_ncat(&CATDIST));
}

static void test_eta()
{
	size_t n = catdist1_ncat(&CATDIST1);
	double *eta = xmalloc(n * sizeof(*eta));
	size_t i, iz;

	for (i = 0; i < n; i++) {
		eta[i] = catdist_eta(&CATDIST, i);
	}

	for (iz = 0; iz < NZ; iz++) {
		eta[IND[iz]] += DETA[iz];
	}

	for (i = 0; i < n; i++) {
		//assert(catdist1_eta(&CATDIST1, i) == eta[i]);
		assert_real_identical(catdist1_eta(&CATDIST1, i), eta[i]);
	}

	free(eta);
}

static void test_psi()
{
	size_t n = catdist1_ncat(&CATDIST1);
	double *eta = xmalloc(n * sizeof(*eta));
	size_t i;

	for (i = 0; i < n; i++) {
		eta[i] = catdist1_eta(&CATDIST1, i);
	}

	double psi = get_psi(eta, n);
	assert_real_approx(catdist1_psi(&CATDIST1), psi);

	free(eta);
}

static void test_lprob()
{
	size_t n = catdist1_ncat(&CATDIST1);
	double psi = catdist1_psi(&CATDIST1);
	size_t i;
	
	for (i = 0; i < n; i++) {
		double eta = catdist1_eta(&CATDIST1, i);
		assert_real_approx(catdist1_lprob(&CATDIST1, i), eta - psi);
	}
}

static void test_prob()
{
	size_t n = catdist1_ncat(&CATDIST1);
	size_t i;
	
	for (i = 0; i < n; i++) {
		double lp = catdist1_lprob(&CATDIST1, i);
		assert_real_identical(catdist1_prob(&CATDIST1, i), exp(lp));
	}
}

static void set_deta(size_t i, double deta)
{
	assert(i < catdist1_ncat(&CATDIST1));
	assert(deta < INFINITY);
	
	size_t iz;
	for (iz = 0; iz < NZ; iz++) {
		if (IND[iz] == i)
			break;
	}
	if (iz == NZ)
		IND[NZ++] = i;
	
	DETA[iz] = deta;
	
	catdist1_set_deta(&CATDIST1, i, deta);
}

static void test_set_deta(size_t i, double deta)
{
	set_deta(i, deta);
	assert_false(catdist1_check(&CATDIST1));
	test_ncat();
	test_eta();
}

static void test_set_deta_rand(double min, double max)
{
	size_t n = catdist1_ncat(&CATDIST1);
	size_t i = rand() % n;
	double deta = runif(min, max);
	test_set_deta(i, deta);
	
}

static void test_set_deta_small()
{
	test_set_deta_rand(-1, 1);
}

static void test_set_deta_med()
{
	test_set_deta_rand(-10, 10);
}

static void test_set_deta_big()
{
	test_set_deta_rand(-1000, 1000);
}



static void test_set_deta_rand_rep(size_t nrep, double min, double max)
{
	size_t rep;
	
	for (rep = 0; rep < nrep; rep++) {
		test_set_deta_rand(min, max);
		print_message(".");
		fflush(stdout);
	}
}

static void test_set_deta_small_rep()
{
	test_set_deta_rand_rep(100, -1, 1);
}

static void test_set_deta_med_rep()
{
	test_set_deta_rand_rep(100, -10, 10);
}

static void test_set_deta_big_rep()
{
	test_set_deta_rand_rep(100, -1000, 1000);
}



static void test_many_set_deta(double min, double max)
{
	size_t n = catdist1_ncat(&CATDIST1);
	size_t iz, nz = 5 * (rand() % n);

	for (iz = 0; iz < nz; iz++) {
		size_t i = rand() % n;
		double deta = runif(min, max);
		set_deta(i, deta);
	}

	assert_false(catdist1_check(&CATDIST1));
	test_ncat();
	test_eta();
}


static void test_many_set_deta_small()
{
	test_many_set_deta(-1, 1);
}

static void test_many_set_deta_med()
{
	test_many_set_deta(-10, 10);
}

static void test_many_set_deta_big()
{
	test_many_set_deta(-1000, 1000);
}


static void test_many_set_deta_rep(size_t nrep, double min, double max)
{
	size_t rep;

	for (rep = 0; rep < nrep; rep++) {
		clear();
		test_many_set_deta(min, max);
		print_message(".");
		fflush(stdout);
	}
}

static void test_many_set_deta_small_rep()
{
	test_many_set_deta_rep(100, -1, 1);
}

static void test_many_set_deta_med_rep()
{
	test_many_set_deta_rep(100, -10, 10);
}

static void test_many_set_deta_big_rep()
{
	test_many_set_deta_rep(100, 0, 500);
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

		unit_test_setup(zeros_suite, zeros_setup_fixture),
		unit_test_setup_teardown(test_ncat, zeros_setup, teardown),
		unit_test_setup_teardown(test_eta, zeros_setup, teardown),
		unit_test_setup_teardown(test_psi, zeros_setup, teardown),
		unit_test_setup_teardown(test_lprob, zeros_setup, teardown),
		unit_test_setup_teardown(test_prob, zeros_setup, teardown),
		unit_test_setup_teardown(test_set_deta_small, zeros_setup, teardown),
		unit_test_setup_teardown(test_set_deta_med, zeros_setup, teardown),		
		unit_test_setup_teardown(test_set_deta_big, zeros_setup, teardown),
		unit_test_setup_teardown(test_set_deta_small_rep, zeros_setup, teardown),
		unit_test_setup_teardown(test_set_deta_med_rep, zeros_setup, teardown),		
		unit_test_setup_teardown(test_set_deta_big_rep, zeros_setup, teardown),		
		unit_test_setup_teardown(test_many_set_deta_small, zeros_setup, teardown),
		unit_test_setup_teardown(test_many_set_deta_med, zeros_setup, teardown),
		unit_test_setup_teardown(test_many_set_deta_big, zeros_setup, teardown),
		unit_test_setup_teardown(test_many_set_deta_small_rep, zeros_setup, teardown),
		unit_test_setup_teardown(test_many_set_deta_med_rep, zeros_setup, teardown),
		unit_test_setup_teardown(test_many_set_deta_big_rep, zeros_setup, teardown),
		unit_test_teardown(zeros_suite, teardown_fixture),

		unit_test_setup(simple_suite, simple_setup_fixture),
		unit_test_setup_teardown(test_ncat, simple_setup, teardown),
		unit_test_setup_teardown(test_eta, simple_setup, teardown),
		unit_test_setup_teardown(test_psi, simple_setup, teardown),
		unit_test_setup_teardown(test_lprob, simple_setup, teardown),
		unit_test_setup_teardown(test_prob, simple_setup, teardown),
		unit_test_setup_teardown(test_set_deta_small, simple_setup, teardown),
		unit_test_setup_teardown(test_set_deta_med, simple_setup, teardown),		
		unit_test_setup_teardown(test_set_deta_big, simple_setup, teardown),
		unit_test_setup_teardown(test_set_deta_small_rep, simple_setup, teardown),
		unit_test_setup_teardown(test_set_deta_med_rep, simple_setup, teardown),		
		unit_test_setup_teardown(test_set_deta_big_rep, simple_setup, teardown),		
		unit_test_setup_teardown(test_many_set_deta_small, simple_setup, teardown),
		unit_test_setup_teardown(test_many_set_deta_med, simple_setup, teardown),
		unit_test_setup_teardown(test_many_set_deta_big, simple_setup, teardown),
		unit_test_setup_teardown(test_many_set_deta_small_rep, simple_setup, teardown),
		unit_test_setup_teardown(test_many_set_deta_med_rep, simple_setup, teardown),
		unit_test_setup_teardown(test_many_set_deta_big_rep, simple_setup, teardown),
		unit_test_teardown(simple_suite, teardown_fixture),

		unit_test_setup(hard_suite, hard_setup_fixture),
		unit_test_setup_teardown(test_ncat, hard_setup, teardown),
		unit_test_setup_teardown(test_eta, hard_setup, teardown),
		unit_test_setup_teardown(test_psi, hard_setup, teardown),
		unit_test_setup_teardown(test_lprob, hard_setup, teardown),
		unit_test_setup_teardown(test_prob, hard_setup, teardown),
		unit_test_setup_teardown(test_set_deta_small, hard_setup, teardown),
		unit_test_setup_teardown(test_set_deta_med, hard_setup, teardown),		
		unit_test_setup_teardown(test_set_deta_big, hard_setup, teardown),
		unit_test_setup_teardown(test_set_deta_small_rep, hard_setup, teardown),
		unit_test_setup_teardown(test_set_deta_med_rep, hard_setup, teardown),		
		unit_test_setup_teardown(test_set_deta_big_rep, hard_setup, teardown),		
		unit_test_setup_teardown(test_many_set_deta_small, hard_setup, teardown),
		unit_test_setup_teardown(test_many_set_deta_med, hard_setup, teardown),
		unit_test_setup_teardown(test_many_set_deta_big, hard_setup, teardown),
		unit_test_setup_teardown(test_many_set_deta_small_rep, hard_setup, teardown),
		unit_test_setup_teardown(test_many_set_deta_med_rep, hard_setup, teardown),
		unit_test_setup_teardown(test_many_set_deta_big_rep, hard_setup, teardown),
		unit_test_teardown(hard_suite, teardown_fixture),
	};

	return run_tests(tests);
}
