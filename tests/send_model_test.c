#include "port.h"
#include <assert.h>
#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <stdlib.h>

#include "../src/history.h"
#include "../src/design.h"
#include "../src/send_model.h"


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
#include "fixtures/send_model.h"


struct fixture {
	struct history_fixture h;
	struct design_fixture d;
	struct send_model_fixture m;
};

struct fixture ctx;

#define NRECV	ctx.h.nrecv
#define NSEND	ctx.h.nsend



struct test {
	struct history h;
	struct design d;
	struct send_model m;
};

struct test test;

#define HISTORY	&test.h
#define DESIGN	&test.d
#define SEND_MODEL &test.m



static void fixture_setup_enron()
{
	double width = 6 * 7 * 24 * 60 * 60;
	double intvls[] = { 2 * 60 * 60, 32 * 60 * 60, 3 * 7 * 24 * 60 * 60 };
	size_t nintvl = 3;

	history_fixture_setup_enron(&ctx.h);

	design_fixture_setup_enron(&ctx.d);
	design_fixture_add_traits_enron(&ctx.d);
	design_fixture_add_isendtot(&ctx.d, width);
	design_fixture_add_nrecvtot(&ctx.d, intvls, nintvl);
	design_fixture_add_nsendtot(&ctx.d, intvls, nintvl);

	send_model_fixture_setup(&ctx.m, &ctx.d);
}


static void fixture_teardown()
{
	send_model_fixture_teardown(&ctx.m);
	design_fixture_teardown(&ctx.d);
	history_fixture_teardown(&ctx.h);
}


static void setup()
{
	history_test_setup(&test.h, &ctx.h);
	design_test_setup(&test.d, &test.h, &ctx.d);
	send_model_test_setup(&test.m, &test.d, &ctx.m);
}


static void teardown()
{
	send_model_test_teardown(&test.m);
	design_test_teardown(&test.d);
	history_test_teardown(&test.h);
}



int main()
{
	UnitTest tests[] = {
		unit_test_setup(enron_suite, fixture_setup_enron),
		unit_test_teardown(enron_suite, fixture_teardown),
	};
	return run_tests(tests);
}

