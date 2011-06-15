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
#include "vnrecv.h"
#include "vrecv.h"
#include "frame.h"
#include "model.h"


static struct actors senders;
static struct actors receivers;
static struct vector intervals;
static struct messages messages;
static struct design design;
static struct frame frame;
static struct vector coefs;
static struct model model;


static void enron_setup_fixture(void **state)
{
	print_message("Enron\n");
	print_message("-----\n");
	enron_employees_init(&senders);
	enron_employees_init(&receivers);
	enron_messages_init(&messages);
}

static void enron_teardown_fixture(void **state)
{
	messages_deinit(&messages);
	actors_deinit(&receivers);
	actors_deinit(&senders);
	print_message("\n\n");
}

static void basic_setup(void **state)
{
	ssize_t i;
	double intvls[3] = {
		112.50,  450.00, 1800.00,
	};
	struct vector vintvls = vector_make(intvls, 3);
	bool has_reffects = false;
	bool has_loops = false;
	vector_init(&intervals, 3);
	vector_assign_copy(&intervals, &vintvls);
	design_init(&design, &senders, &receivers, &intervals);
	design_set_loops(&design, has_loops);
	design_set_reffects(&design, has_reffects);
	design_add_var(&design, VAR_TYPE_RECV);
	frame_init(&frame, &design);
	vector_init(&coefs, design_dim(&design));
	
	for (i = 0; i < vector_dim(&coefs); i++) {
		double val = (i % 5 == 0 ? -2.0 :
			      i % 5 == 1 ?  1.0 :
			      i % 5 == 2 ? -1.0 :
			      i % 5 == 3 ?  2.0 : 0.0);
		vector_set_item(&coefs, i, val);
	}
	
	model_init(&model, &design, &coefs);
}

static void basic_teardown(void **state)
{
	model_deinit(&model);
	vector_deinit(&coefs);
	frame_deinit(&frame);
	vector_deinit(&intervals);	
	design_deinit(&design);
}

static void test_probs(void **state)
{
	struct send_model *sm;
	struct vector eta, probs, logprobs;
	struct messages_iter it;
	const struct message *msg = NULL;
	double t;
	ssize_t isend, jrecv, nrecv;
	ssize_t i, n;
	
	nrecv = design_receiver_count(&design);
	vector_init(&eta, nrecv);	
	vector_init(&probs, nrecv);
	vector_init(&logprobs, nrecv);	

	MESSAGES_FOREACH(it, &messages) {
		t = MESSAGES_TIME(it);
		
		while (frame_next_change(&frame) <= t) {
			model_update(&model, &frame);
			frame_advance(&frame);			
		}
		if (frame_time(&frame) < t) {
			model_update(&model, &frame);
			frame_advance_to(&frame, t);			
		}
		
		n = MESSAGES_COUNT(it);
		for (i = 0; i < n; i ++) {
			msg = MESSAGES_VAL(it, i);
			isend = msg->from;
			sm = model_send_model(&model, &frame, isend);
			
			frame_mul(1.0, TRANS_NOTRANS, &frame, isend, &coefs, 0.0, &eta);
			
			if (!design_loops(&design))
				vector_set_item(&eta, isend, -INFINITY);
			
			vector_assign_copy(&logprobs, &eta);
			double log_W = vector_log_sum_exp(&logprobs);
			vector_shift(&logprobs, -log_W);
			
			vector_assign_copy(&probs, &logprobs);
			vector_exp(&probs);
			
			for (jrecv = 0; jrecv < nrecv; jrecv++) {
				double lp0 = send_model_logprob(sm, jrecv);
				double lp1 = vector_item(&logprobs, jrecv);
				assert_in_range(double_eqrel(lp0, lp1), 49, DBL_MANT_DIG);
				
				double p0 = send_model_prob(sm, jrecv);
				double p1 = vector_item(&probs, jrecv);
				assert_in_range(double_eqrel(p0, p1), 47, DBL_MANT_DIG);
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
}

int main(int argc, char **argv)
{
	UnitTest tests[] = {
		unit_test_setup(enron_suite, enron_setup_fixture),
		unit_test_setup_teardown(test_probs, basic_setup, basic_teardown),
		unit_test_teardown(enron_suite, enron_teardown_fixture),
	};
	return run_tests(tests);
}
