#include "port.h"
#include <math.h>
#include <float.h>
#include <assert.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdio.h>
#include <setjmp.h>
#include <stdlib.h>
#include <string.h>
#include "cmockery.h"

#include "coreutil.h"
#include "xalloc.h"
#include "ieee754.h"
#include "util.h"
#include "enron.h"
#include "messages.h"
#include "design.h"
#include "vars.h"
#include "frame.h"
#include "testutil.h"
#include "recv_model.h"


static size_t nsend;
static size_t nrecv;
static size_t ntrait;
static double *traits;
static const char * const *trait_names;
static struct messages messages;
static struct frame frame;
static struct recv_coefs coefs;
static struct recv_model model;


static void enron_setup_fixture()
{
	print_message("Enron\n");
	print_message("-----\n");
	enron_employees_init(&nsend, &traits, &ntrait, &trait_names, ENRON_TERMS_MAX);
	enron_messages_init(&messages, -1);
	nrecv = nsend;
}

static void enron_teardown_fixture()
{
	free(traits);
	messages_deinit(&messages);
	print_message("\n\n");
}

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


static void test_probs()
{
	struct design *r = frame_recv_design(&frame);
	struct design2 *d = frame_dyad_design(&frame);
	struct history *h = frame_history(&frame);

	double *eta, *probs, *logprobs;
	struct messages_iter it;
	const struct message *msg = NULL;
	double t;
	size_t isend, jrecv, nrecv;
	size_t i, n;
	
	nrecv = frame_recv_count(&frame);
	eta = xmalloc(nrecv * sizeof(double));
	probs = xmalloc(nrecv * sizeof(double));
	logprobs = xmalloc(nrecv * sizeof(double));

	MESSAGES_FOREACH(it, &messages) {
		// fprintf(stderr, "."); fflush(stderr);
		t = MESSAGES_TIME(it);
		
		history_advance(h, t);
		
		n = MESSAGES_COUNT(it);
		for (i = 0; i < n; i ++) {
			msg = MESSAGES_VAL(it, i);
			isend = msg->from;
			design_mul(1.0, r, &coefs.recv, 0.0, eta);
			design2_mul(1.0, d, isend, &coefs.dyad, 1.0, eta);
			
			if (!frame_has_loops(&frame))
				eta[isend] = -INFINITY;
			
			blas_dcopy(nrecv, eta, 1, logprobs, 1);
			double max_eta = vector_max(nrecv, eta);
			vector_shift(nrecv, -max_eta, logprobs);
			double log_W = vector_logsumexp(nrecv, logprobs);
			vector_shift(nrecv, -log_W, logprobs);
			
			blas_dcopy(nrecv, logprobs, 1, probs, 1);
			vector_exp(nrecv, probs);

			struct catdist1 *dist = recv_model_dist(&model, isend);
			double psi = catdist1_cached_psi(dist);

			assert(double_eqrel(log_W + max_eta, psi) >= 36);
			assert_in_range(double_eqrel(log_W + max_eta, psi), 36, DBL_MANT_DIG);
			
			for (jrecv = 0; jrecv < nrecv; jrecv++) {
				double lp0 = logprobs[jrecv];
				double lp1 = catdist1_cached_lprob(dist, jrecv);
				assert_real_approx(lp0, lp1);

				double p0 = probs[jrecv];
				double p1 = catdist1_cached_prob(dist, jrecv);
				assert_real_approx(p0, p1);
			}
		}
		
		n = MESSAGES_COUNT(it);
		for (i = 0; i < n; i ++) {
			msg = MESSAGES_VAL(it, i);
			history_add(h, msg);
		}
	}
	
	free(probs);
	free(logprobs);
	free(eta);
}

int main()
{
	UnitTest tests[] = {
		unit_test_setup(enron_suite, enron_setup_fixture),
		unit_test_setup_teardown(test_probs, basic_setup, teardown),
		unit_test_setup_teardown(test_probs, hard_setup, teardown),		
		unit_test_teardown(enron_suite, enron_teardown_fixture),
	};
	return run_tests(tests);
}
