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

#include "ieee754.h"
#include "enron.h"
#include "messages.h"
#include "design.h"
#include "vars.h"
#include "frame.h"
#include "recv_model.h"

static struct actors enron_actors;
static struct matrix enron_traits;

static struct actors senders;
static struct actors receivers;
static struct matrix recv_traits;
static struct vector intervals;
static struct messages messages;
static struct design design;
static struct frame frame;
static struct matrix coefs;
static struct recv_model model;


static void enron_setup_fixture(void **state)
{
	print_message("Enron\n");
	print_message("-----\n");
	enron_employees_init(&enron_actors, &enron_traits);
	enron_messages_init(&messages);
	actors_init_copy(&senders, &enron_actors);
	actors_init_copy(&receivers, &enron_actors);	
	matrix_init_copy(&recv_traits, TRANS_NOTRANS, &enron_traits);
}

static void enron_teardown_fixture(void **state)
{
	matrix_deinit(&recv_traits);
	actors_deinit(&receivers);
	actors_deinit(&senders);
	messages_deinit(&messages);
	matrix_deinit(&enron_traits);
	actors_deinit(&enron_actors);	
	print_message("\n\n");
}

static void basic_setup(void **state)
{
	ssize_t i, c;
	double intvls[3] = {
		112.50,  450.00, 1800.00,
	};
	struct vector vintvls = vector_make(intvls, 3);
	bool has_reffects = false;
	bool has_loops = false;
	vector_init(&intervals, 3);
	vector_assign_copy(&intervals, &vintvls);
	design_init(&design, &senders, &receivers, &recv_traits, &intervals);
	design_set_loops(&design, has_loops);
	design_set_recv_effects(&design, has_reffects);
	design_add_recv_var(&design, RECV_VAR_NRECV);
	frame_init(&frame, &design);
	matrix_init(&coefs, design_recv_dim(&design), actors_cohort_count(&senders));
	
	for (c = 0; c < matrix_ncol(&coefs); c++) {
		for (i = 0; i < matrix_nrow(&coefs); i++) {
			double val = (i + 2 * c % 5 == 0 ? -2.0 :
				      i + 3 * c % 5 == 1 ?  1.0 :
				      i + 7 * c % 5 == 2 ? -1.0 :
				      i + 11 * c % 5 == 3 ?  2.0 : 0.0);
			matrix_set_item(&coefs, i, c, val);
		}
	}
	
	recv_model_init(&model, &frame, &senders, &coefs);
}

static void teardown(void **state)
{
	recv_model_deinit(&model);
	matrix_deinit(&coefs);
	frame_deinit(&frame);
	vector_deinit(&intervals);	
	design_deinit(&design);
}


static void hard_setup(void **state)
{
	ssize_t i, c;
	double intvls[3] = {
		112.50,  450.00, 1800.00,
	};
	struct vector vintvls = vector_make(intvls, 3);
	bool has_reffects = false;
	bool has_loops = false;
	vector_init(&intervals, 3);
	vector_assign_copy(&intervals, &vintvls);
	design_init(&design, &senders, &receivers, &recv_traits, &intervals);
	design_set_loops(&design, has_loops);
	design_set_recv_effects(&design, has_reffects);
	design_add_recv_var(&design, RECV_VAR_NRECV);
	frame_init(&frame, &design);
	matrix_init(&coefs, design_recv_dim(&design), actors_cohort_count(&senders));
	
	for (c = 0; c < matrix_ncol(&coefs); c++) {
		for (i = 0; i < matrix_nrow(&coefs); i++) {
			double val = (i + 2 * c % 7 == 0 ?  0.1 :
				      i + 3 * c % 7 == 1 ?  0.3 :
				      i + 5 * c % 7 == 2 ? -0.2 :
				      i + 11 * c % 7 == 4 ? -10 :
				      i % 7 == 6 ? +10 : 0.0);
			matrix_set_item(&coefs, i, c, val);
		}
	}
	
	recv_model_init(&model, &frame, &senders, &coefs);
}


static void test_probs(void **state)
{
	struct vector eta, probs, logprobs, y;
	struct messages_iter it;
	const struct message *msg = NULL;
	double t;
	ssize_t isend, jrecv, nrecv;
	ssize_t i, n;
	
	nrecv = design_recv_count(&design);
	vector_init(&eta, nrecv);	
	vector_init(&probs, nrecv);
	vector_init(&logprobs, nrecv);
	vector_init(&y, nrecv);

	double alpha = 2.0;
	double y0 = 3.14;

	ssize_t minprec = DBL_MANT_DIG;
	
	MESSAGES_FOREACH(it, &messages) {
		// fprintf(stderr, "."); fflush(stderr);
		t = MESSAGES_TIME(it);
		
		frame_advance(&frame, t);			
		
		n = MESSAGES_COUNT(it);
		for (i = 0; i < n; i ++) {
			msg = MESSAGES_VAL(it, i);
			isend = msg->from;
			ssize_t c = actors_items(&senders)[isend].cohort;
			const struct vector col = matrix_col(&coefs, c);
			frame_recv_mul(1.0, TRANS_NOTRANS, &frame, isend, &col, 0.0, &eta);
			
			if (!design_loops(&design))
				vector_set_item(&eta, isend, -INFINITY);
			
			vector_assign_copy(&logprobs, &eta);
			double max_eta = vector_max(&eta);
			vector_shift(&logprobs, -max_eta);
			double log_W = vector_log_sum_exp(&logprobs);
			vector_shift(&logprobs, -log_W);
			
			vector_assign_copy(&probs, &logprobs);
			vector_exp(&probs);
			
			vector_fill(&y, y0);
			recv_model_axpy_probs(alpha, &model, isend, &y);
			
			assert(double_eqrel(log_W + max_eta, recv_model_logsumwt(&model, isend)) >= 36);
			assert_in_range(double_eqrel(log_W + max_eta, recv_model_logsumwt(&model, isend)), 36, DBL_MANT_DIG);
			
			for (jrecv = 0; jrecv < nrecv; jrecv++) {
				double lp0 = vector_item(&logprobs, jrecv);				
				double lp1 = recv_model_logprob(&model, isend, jrecv);

				if (fabs(lp0) >= 5e-4) {
					//minprec = MIN(minprec, double_eqrel(lp0, lp1));
					assert(double_eqrel(lp0, lp1) >= 36);
					assert_in_range(double_eqrel(lp0, lp1), 36, DBL_MANT_DIG);

				} else {
					assert(fabs(lp0 - lp1) < sqrt(DBL_EPSILON));
					assert_true(fabs(lp0 - lp1) < sqrt(DBL_EPSILON));
				}

				double p0 = vector_item(&probs, jrecv);				
				double p1 = recv_model_prob(&model, isend, jrecv);

				if (fabs(p0) >= 5e-4) {
					minprec = MIN(minprec, double_eqrel(p0, p1));
					assert(double_eqrel(p0, p1) >= 38);
					assert_in_range(double_eqrel(p0, p1), 38, DBL_MANT_DIG);
				} else {
					assert_true(fabs(p0 - p1) < sqrt(DBL_EPSILON));
				}
				
				assert_in_range(double_eqrel(alpha * p0 + y0,
							     vector_item(&y, jrecv)),
						45,
						DBL_MANT_DIG);
			}
		}
		
		n = MESSAGES_COUNT(it);
		for (i = 0; i < n; i ++) {
			msg = MESSAGES_VAL(it, i);
			frame_add(&frame, msg);
		}
	}
	
	vector_deinit(&probs);
	vector_deinit(&logprobs);
	vector_deinit(&eta);
	vector_deinit(&y);	
}

int main(int argc, char **argv)
{
	UnitTest tests[] = {
		unit_test_setup(enron_suite, enron_setup_fixture),
		unit_test_setup_teardown(test_probs, basic_setup, teardown),
		unit_test_setup_teardown(test_probs, hard_setup, teardown),		
		unit_test_teardown(enron_suite, enron_teardown_fixture),
	};
	return run_tests(tests);
}
