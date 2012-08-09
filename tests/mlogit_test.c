#include "port.h"

#include <stdlib.h> // free
#include <string.h> // memcpy, memset

#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include "cmockery.h"
#include "testutil.h"

#include <math.h>  // INFINITY
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
	qsort(eta_sort, n, sizeof(double), double_rcompare);

	double etamax = eta_sort[n - 1];
	double sum = 0.0;
	size_t i;
	
	for (i = 0; i < N - 1; i++) {
		sum += exp(eta_sort[i] - etamax);
	}
	
	double phi = etamax + log1p(sum);
	
	free(eta_sort);
	
	return phi;
}

static void setup(const double *eta, size_t n)
{
	N = n;
	ETA = xmalloc(N * sizeof(ETA[0]));

	if (eta) {
		memcpy(ETA, eta, N * sizeof(ETA[0]));
	} else {
		memset(ETA, 0, N * sizeof(ETA[0]));
	}
	
	PHI = get_phi(eta, n);

	mlogit_init(&MLOGIT, N);

	if (eta) {
		mlogit_set_all_eta(&MLOGIT, ETA);
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
	setup(NULL, 9);
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


int main()
{
	UnitTest tests[] = {
		unit_test_setup(zeros_suite, zeros_setup_fixture),
		unit_test_setup_teardown(test_ncat, zeros_setup, teardown),
		unit_test_setup_teardown(test_eta, zeros_setup, teardown),
		unit_test_setup_teardown(test_phi, zeros_setup, teardown),		
		unit_test_teardown(zeros_suite, teardown_fixture),

		unit_test_setup(simple_suite, simple_setup_fixture),
		unit_test_setup_teardown(test_ncat, simple_setup, teardown),
		unit_test_setup_teardown(test_eta, simple_setup, teardown),
		unit_test_setup_teardown(test_phi, simple_setup, teardown),		
		unit_test_teardown(simple_suite, teardown_fixture),
		
		unit_test_setup(empty_suite, empty_setup_fixture),
		unit_test_setup_teardown(test_ncat, empty_setup, teardown),
		unit_test_setup_teardown(test_eta, empty_setup, teardown),
		unit_test_setup_teardown(test_phi, empty_setup, teardown),		
		unit_test_teardown(empty_suite, teardown_fixture),
		
		
	};

	return run_tests(tests);
}
