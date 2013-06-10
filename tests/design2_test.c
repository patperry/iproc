#include "port.h"
#include <assert.h>
#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <stdlib.h>

#include "../src/history.h"
#include "../src/design.h"
#include "../src/design2.h"

#include "coreutil.h"
#include "testutil.h"
#include "cmockery.h"
#include "ieee754.h"
#include "lapack.h"
#include "matrixutil.h"
#include "sblas.h"
#include "xalloc.h"

#include "fixtures/history.h"
#include "fixtures/design.h"
#include "fixtures/design2.h"


#include "fixtures/history.h"
#include "fixtures/design.h"


struct fixture {
	struct history_fixture h;
	struct design2_fixture d;
};

struct fixture ctx;

#define NRECV	ctx.h.nrecv
#define NSEND	ctx.h.nsend



struct test {
	struct history h;
	struct design2 d;
};

struct test test;

#define DESIGN2	&test.d
#define HISTORY	&test.h



static void fixture_setup_enron()
{
	print_message("Enron employees\n");
	print_message("---------------\n");
	history_fixture_setup_enron(&ctx.h);
	design2_fixture_setup_enron(&ctx.d);
}

static void fixture_teardown()
{
	design2_fixture_teardown(&ctx.d);
	history_fixture_teardown(&ctx.h);
	print_message("\n\n");
}


static void setup()
{
	history_test_setup(&test.h, &ctx.h);
	design2_test_setup(&test.d, &test.h, &ctx.d);
}

static void teardown()
{
	design2_test_teardown(&test.d);
	history_test_teardown(&test.h);
}


static void test_isend()
{
	double window = 1000;
	const struct var2 *v = design2_add_tvar(DESIGN2, "ISend", VAR2_ISEND, window);
	double t0 = 910930020;
	const double *x;
	const size_t *ind;
	size_t nz;

	assert_true(v);

	design2_get_tvar_matrix(DESIGN2, 137, &x, &ind, &nz); /* -inf */
	assert_int_equal(nz, 0);

	history_advance(HISTORY, t0); /* t */
	design2_get_tvar_matrix(DESIGN2, 137, &x, &ind, &nz);
	assert_int_equal(nz, 0);

	history_advance(HISTORY, double_nextup(t0)); /* (t)+ */
	design2_get_tvar_matrix(DESIGN2, 137, &x, &ind, &nz);
	assert_int_equal(nz, 1);
	assert_int_equal(ind[0], 58);
	assert_real_identical(x[0], 1.0);

	design2_get_tvar_matrix(DESIGN2, 58, &x, &ind, &nz);
	assert_int_equal(nz, 0);

	history_advance(HISTORY, t0 + window); /* t + w */
	design2_get_tvar_matrix(DESIGN2, 137, &x, &ind, &nz);
	assert_int_equal(nz, 1);
	assert_int_equal(ind[0], 58);
	assert_real_identical(x[0], 1.0);

	design2_get_tvar_matrix(DESIGN2, 58, &x, &ind, &nz);
	assert_int_equal(nz, 0);

	history_advance(HISTORY, double_nextup(t0 + window)); /* (t + w)+ */
	design2_get_tvar_matrix(DESIGN2, 137, &x, &ind, &nz);
	assert_int_equal(nz, 1);
	assert_int_equal(ind[0], 58);
	assert_real_identical(x[0], 0.0);

	design2_get_tvar_matrix(DESIGN2, 58, &x, &ind, &nz);
	assert_int_equal(nz, 0);
}



int main()
{
	UnitTest tests[] = {
		unit_test_setup(enron_suite, fixture_setup_enron),
		unit_test_setup_teardown(test_isend, setup, teardown),
		unit_test_teardown(enron_suite, fixture_teardown),

	};
	return run_tests(tests);
}


