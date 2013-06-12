#include "port.h"
#include <assert.h>
#include <stdarg.h>
#include <stddef.h>
#include <setjmp.h>
#include <stdlib.h>

#include "coreutil.h"
#include "testutil.h"
#include "cmockery.h"
#include "ieee754.h"
#include "lapack.h"
#include "matrixutil.h"
#include "sblas.h"
#include "xalloc.h"

#include "../src/history.h"
#include "../src/design.h"
#include "../src/design2.h"
#include "../src/recv_model.h"

#include "fixtures/history.h"
#include "fixtures/design.h"
#include "fixtures/design2.h"
#include "fixtures/recv_model.h"


struct fixture {
	struct history_fixture h;
	struct design_fixture r;
	struct design2_fixture d;
	struct recv_model_fixture m;
};

struct fixture ctx;

#define NRECV (ctx.h.nrecv)
#define NSEND (ctx.h.nsend)
#define RECV_PARAMS (&ctx.m.params)
#define EXCLUDE_LOOPS (ctx.m.exclude_loops)


struct test {
	struct history h;
	struct design r;
	struct design2 d;
	struct recv_model m;
};

struct test test;

#define HISTORY	(&test.h)
#define RECV_DESIGN (&test.r)
#define DYAD_DESIGN (&test.d)
#define RECV_MODEL (&test.m)


static void fixture_setup_enron()
{
	double width = 6 * 7 * 24 * 60 * 60;
	double intvls[] = { 2 * 60 * 60, 32 * 60 * 60, 3 * 7 * 24 * 60 * 60 };
	size_t nintvl = 3;

	srand(0);

	history_fixture_setup_enron(&ctx.h);

	design_fixture_setup_enron(&ctx.r);
	design_fixture_add_traits_enron(&ctx.r);
	design_fixture_add_isendtot(&ctx.r, width);
	design_fixture_add_nrecvtot(&ctx.r, intvls, nintvl);
	design_fixture_add_nsendtot(&ctx.r, intvls, nintvl);

	design2_fixture_setup_enron(&ctx.d);
	design2_fixture_add_isend(&ctx.d, width);
	design2_fixture_add_nrecv(&ctx.d, intvls, nintvl);
	design2_fixture_add_nsend(&ctx.d, intvls, nintvl);

	recv_model_fixture_setup(&ctx.m, &ctx.r, &ctx.d);
	recv_model_fixture_set_exlude_loops(&ctx.m, 1);
	recv_model_fixture_set_rand(&ctx.m);
}


static void fixture_teardown()
{
	recv_model_fixture_teardown(&ctx.m);
	design2_fixture_teardown(&ctx.d);
	design_fixture_teardown(&ctx.r);
	history_fixture_teardown(&ctx.h);
}


static void setup()
{
	history_test_setup(&test.h, &ctx.h);
	design_test_setup(&test.r, &test.h, &ctx.r);
	design2_test_setup(&test.d, &test.h, &ctx.d);
	recv_model_test_setup(&test.m, &test.r, &test.d, &ctx.m);
}


static void teardown()
{
	recv_model_test_teardown(&test.m);
	design2_test_teardown(&test.d);
	design_test_teardown(&test.r);
	history_test_teardown(&test.h);
}


/*
static void basic_setup()
{
	size_t i;
	double intvls[4] = {
		112.50,  450.00, 1800.00, INFINITY
	};
	int has_loops = 0;
	frame_init(&frame, nsend, nrecv, has_loops, intvls, 4);
	struct design *r = frame_recv_design(&frame);
	design_add_traits(r, trait_names, traits, ntrait);
	
	struct design2 *d = frame_dyad_design(&frame);
	design2_add_tvar(d, "NRecv", VAR2_NRECV);
	recv_coefs_init(&coefs, &frame);
	
	
	for (i = 0; i < coefs.dim; i++) {
		double val = (i % 5 == 0 ? -2.0 :
			      i % 5 == 1 ?  1.0 :
			      i % 5 == 2 ? -1.0 :
			      i % 5 == 3 ?  2.0 : 0.0);
		coefs.all[i] = val;
	}

	recv_model_init(&model, &frame, &coefs);
}

static void teardown()
{
	recv_model_deinit(&model);
	recv_coefs_deinit(&coefs);
	frame_deinit(&frame);
}


static void hard_setup()
{
	double intvls[7] = {
		1800.00,   // 30 min
		7200.00,   //  2 hr
		28800.00,  //  8 hr
		115200.00, // 32 hr
		460800.00, // 5.33 day
		1843200.00, // 21.33 day
		INFINITY
	};

	double beta[] = {
		0.4824403589966597,
		-0.215108203814423,
		-0.4182912193871692,
		-0.04491132567966908,
		0.05698656051879408,
		0.1420793131926562,
		-0.249037515029601,
		-0.02644490380870458,
		0.2987018579760873,
		1.128469426856978,
		1.258299141255399,
		0.3882642769919149,
		0.0984948728445552,
		0.03433916773539594,
		0.01662024071711138,
		-0.007969263858043221,
		0.0006115999551254905,
		3.739814520709765,
		0.9588447963575467,
		0.4521880948047269,
		0.187294928856629,
		0.1757997396356974,
		0.08539282957216279,
		0.05050240510803124,
		0.001111188823465376
	};

	int has_loops = 0;
	frame_init(&frame, nsend, nrecv, has_loops, intvls, 7);
	struct design *r = frame_recv_design(&frame);
	design_add_traits(r, trait_names, traits, ntrait);

	struct design2 *d = frame_dyad_design(&frame);
	design2_add_tvar(d, "IRecv", VAR2_IRECV);
	design2_add_tvar(d, "NRecv", VAR2_NRECV);
	design2_add_tvar(d, "ISend", VAR2_ISEND);
	design2_add_tvar(d, "NSend", VAR2_NSEND);

	recv_coefs_init(&coefs, &frame);
	memcpy(coefs.all, beta, coefs.dim * sizeof(double));

	recv_model_init(&model, &frame, &coefs);
}
*/


static void test_probs()
{
	struct design *r = RECV_DESIGN;
	struct design2 *d = DYAD_DESIGN;
	struct history *h = HISTORY;
	const struct recv_params *p = RECV_PARAMS;
	const struct recv_model *m = RECV_MODEL;
	
	const struct message *msgs;
	size_t imsg, nmsg;

	double *eta = xmalloc(NRECV * sizeof(double));

	history_get_messages(HISTORY, &msgs, &nmsg);
	for (imsg = 0; imsg < MIN(nmsg, 1000); imsg++) {
		const struct message *msg = &msgs[imsg];
		double t = msg->time;
		size_t i = msg->from;

		// fprintf(stderr, "."); fflush(stderr);
		history_advance(h, t);

		design_mul(1.0, r, &p->recv, 0.0, eta);
		design2_mul(1.0, d, i, &p->dyad, 1.0, eta);
			
		if (recv_model_exclude_loops(m))
			eta[i] = -INFINITY;

		struct catdist1 *dist = recv_model_dist(m, i);
		size_t j, n = NRECV;

		for (j = 0; j < n; j++) {
			double eta0 = eta[j];
			double eta1 = catdist1_eta(dist, j);
			assert_real_approx(eta0, eta1);
		}
	}
	
	free(eta);
}


int main()
{
	UnitTest tests[] = {
		unit_test_setup(enron_suite, fixture_setup_enron),
		unit_test_setup_teardown(test_probs, setup, teardown),
		unit_test_teardown(enron_suite, fixture_teardown)
	};
	return run_tests(tests);
}
