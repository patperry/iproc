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
#define PARAMS	(&ctx.m.params)



struct test {
	struct history h;
	struct design d;
	struct send_model m;
};

struct test test;

#define HISTORY	(&test.h)
#define DESIGN	(&test.d)
#define SEND_MODEL (&test.m)



static void fixture_setup_enron()
{
	double width = 6 * 7 * 24 * 60 * 60;
	double intvls[] = { 2 * 60 * 60, 32 * 60 * 60, 3 * 7 * 24 * 60 * 60 };
	size_t nintvl = 3;

	srand(0);

	history_fixture_setup_enron(&ctx.h);

	design_fixture_setup_enron(&ctx.d);
	design_fixture_add_traits_enron(&ctx.d);
	design_fixture_add_isendtot(&ctx.d, width);
	design_fixture_add_nrecvtot(&ctx.d, intvls, nintvl);
	design_fixture_add_nsendtot(&ctx.d, intvls, nintvl);

	send_model_fixture_setup(&ctx.m, &ctx.d);
	send_model_fixture_set_rand(&ctx.m);
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

static void test_probs()
{
	const struct message *msgs;
	size_t imsg, nmsg;
	size_t i, n = NSEND;
	struct catdist1 *dist;
	double *eta = xmalloc(NSEND * sizeof(double));

	history_get_messages(HISTORY, &msgs, &nmsg);
	for (imsg = 0; imsg < MIN(nmsg, 1000); imsg++) {
		history_advance(HISTORY, msgs[imsg].time);
		design_mul(1.0, DESIGN, &PARAMS->coefs, 0.0, eta);
		dist = send_model_dist(SEND_MODEL);

		for (i = 0; i < n; i++) {
			assert_real_approx(eta[i], catdist1_eta(dist, i));
		}
	}

	free(eta);
}



int main()
{
	UnitTest tests[] = {
		unit_test_setup(enron_suite, fixture_setup_enron),
		unit_test_setup_teardown(test_probs, setup, teardown),
		unit_test_teardown(enron_suite, fixture_teardown),
	};
	return run_tests(tests);
}

