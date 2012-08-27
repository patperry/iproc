#include "port.h"
#include <math.h>
#include <float.h>
#include <assert.h>
#include <stdarg.h>
#include <stddef.h>
#include <stdio.h>
#include <setjmp.h>
#include <stdlib.h>
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
	enron_employees_init(&nsend, &traits, &ntrait, &trait_names);
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
	double intvls[3] = {
		112.50,  450.00, 1800.00,
	};
	int has_loops = 0;
	frame_init(&frame, nsend, nrecv, has_loops, intvls, 3);
	struct design *r = frame_recv_design(&frame);
	design_add_traits(r, trait_names, traits, ntrait);
	
	struct design2 *d = frame_dyad_design(&frame);
	design2_add_tvar(d, "NRecv", DYAD_VAR_NRECV);
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
	size_t i;
	double intvls[3] = {
		112.50,  450.00, 1800.00,
	};
	int has_loops = 0;
	frame_init(&frame, nsend, nrecv, has_loops, intvls, 3);
	struct design *r = frame_recv_design(&frame);
	design_add_traits(r, trait_names, traits, ntrait);
	
	struct design2 *d = frame_dyad_design(&frame);
	design2_add_tvar(d, "NRecv", DYAD_VAR_NRECV);
	
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


static void test_probs()
{
	struct design *r = frame_recv_design(&frame);
	struct design2 *d = frame_dyad_design(&frame);

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

	size_t minprec = DBL_MANT_DIG;
	
	MESSAGES_FOREACH(it, &messages) {
		// fprintf(stderr, "."); fflush(stderr);
		t = MESSAGES_TIME(it);
		
		frame_advance(&frame, t);			
		
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

				if (fabs(lp0) >= 5e-4) {
					//minprec = MIN(minprec, double_eqrel(lp0, lp1));
					assert(double_eqrel(lp0, lp1) >= 36);
					assert_in_range(double_eqrel(lp0, lp1), 36, DBL_MANT_DIG);

				} else {
					assert(fabs(lp0 - lp1) < sqrt(DBL_EPSILON));
					assert_true(fabs(lp0 - lp1) < sqrt(DBL_EPSILON));
				}

				double p0 = probs[jrecv];
				double p1 = catdist1_cached_prob(dist, jrecv);

				if (fabs(p0) >= 5e-4) {
					minprec = MIN(minprec, (size_t)double_eqrel(p0, p1));
					assert(double_eqrel(p0, p1) >= 38);
					assert_in_range(double_eqrel(p0, p1), 38, DBL_MANT_DIG);
				} else {
					assert_true(fabs(p0 - p1) < sqrt(DBL_EPSILON));
				}
			}
		}
		
		n = MESSAGES_COUNT(it);
		for (i = 0; i < n; i ++) {
			msg = MESSAGES_VAL(it, i);
			frame_add(&frame, msg);
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
