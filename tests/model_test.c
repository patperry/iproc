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
	design_set_recv_effects(&design, has_reffects);
	design_add_recv_var(&design, RECV_VAR_IRECV);
	frame_init(&frame, &design);
	vector_init(&coefs, design_recv_dim(&design));
	
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
	struct recv_model *rm;
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
			rm = model_recv_model(&model, &frame, isend);
			
			frame_recv_mul(1.0, TRANS_NOTRANS, &frame, isend, &coefs, 0.0, &eta);
			
			if (!design_loops(&design))
				vector_set_item(&eta, isend, -INFINITY);
			
			vector_assign_copy(&logprobs, &eta);
			double log_W = vector_log_sum_exp(&logprobs);
			vector_shift(&logprobs, -log_W);
			
			vector_assign_copy(&probs, &logprobs);
			vector_exp(&probs);
			
			vector_fill(&y, y0);
			recv_model_axpy_probs(alpha, rm, &y);
			for (jrecv = 0; jrecv < nrecv; jrecv++) {
				double lp0 = recv_model_logprob(rm, jrecv);
				double lp1 = vector_item(&logprobs, jrecv);
				assert_in_range(double_eqrel(lp0, lp1), 49, DBL_MANT_DIG);
				
				double p0 = recv_model_prob(rm, jrecv);
				double p1 = vector_item(&probs, jrecv);
				assert_in_range(double_eqrel(p0, p1), 47, DBL_MANT_DIG);
				
				assert_true(double_identical(alpha * p0 + y0,
							     vector_item(&y, jrecv)));
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


static void test_mean(void **state)
{
	struct recv_model *rm;
	struct vector probs, mean, mean0;
	struct messages_iter it;
	const struct message *msg = NULL;
	double t;
	ssize_t isend;
	ssize_t i, n;
	ssize_t index, dim = design_recv_dim(&design);
	
	vector_init(&probs, design_recv_count(&design));
	vector_init(&mean, design_recv_dim(&design));
	vector_init(&mean0, design_recv_dim(&design));	
	
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
			rm = model_recv_model(&model, &frame, isend);
			
			vector_fill(&probs, 0.0);
			recv_model_axpy_probs(1.0, rm, &probs);
			
			frame_recv_mul(1.0, TRANS_TRANS, &frame, isend, &probs,
				       0.0, &mean);
			
			vector_fill(&mean0, 0.0);
			recv_model_axpy_mean(1.0, rm, &mean0);
			
			for (index = 0; index < dim; index++) {
				double x0 = vector_item(&mean0, index);
				double x1 = vector_item(&mean, index);				
				assert_in_range(double_eqrel(x0, x1), 44, DBL_MANT_DIG);
			}
		}
		
		n = MESSAGES_COUNT(it);
		for (i = 0; i < n; i ++) {
			msg = MESSAGES_VAL(it, i);
			frame_add(&frame, msg);
		}
	}
	
	vector_deinit(&probs);
	vector_deinit(&mean);
	vector_deinit(&mean0);	
}


static void test_var(void **state)
{
	struct recv_model *rm;
	struct vector mean, y;
	struct svector e_j;
	struct matrix var, var0;
	struct messages_iter it;
	const struct message *msg = NULL;
	double t;
	ssize_t isend, jrecv, nrecv = design_recv_count(&design);
	ssize_t i, n;
	ssize_t index1, index2, dim = design_recv_dim(&design);
	
	vector_init(&mean, dim);
	vector_init(&y, dim);	
	matrix_init(&var, dim, dim);
	matrix_init(&var0, dim, dim);
	svector_init(&e_j, nrecv);
	
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
			rm = model_recv_model(&model, &frame, isend);
			
			vector_fill(&mean, 0.0);
			recv_model_axpy_mean(1.0, rm, &mean);
			matrix_fill(&var, 0.0);
			recv_model_axpy_var(1.0, rm, &frame, &var);

			matrix_fill(&var0, 0.0);
			for (jrecv = 0; jrecv < nrecv; jrecv++) {
				double p = recv_model_prob(rm, jrecv);
				vector_assign_copy(&y, &mean);				
				svector_set_basis(&e_j, jrecv);

				frame_recv_muls(1.0, TRANS_TRANS, &frame, isend, &e_j, -1.0, &y);
				matrix_update1(&var0, p, &y, &y);
			}
			
			for (index2 = 0; index2 < dim; index2++) {
				for (index1 = 0; index1 < dim; index1++) {				
					double v0 = matrix_item(&var0, index1, index2);
					double v1 = matrix_item(&var, index1, index2);
					//printf("v0: %.4f  v1: %.4f (%d)\n", v0, v1, double_eqrel(v0, v1));
					assert(double_eqrel(v0, v1) >= 35
					       || (fabs(v0) < 256 * DBL_EPSILON && fabs(v1) < 256 * DBL_EPSILON));
					if (fabs(v0) >= 256 * DBL_EPSILON) {
						assert_in_range(double_eqrel(v0, v1), 35, DBL_MANT_DIG);
					} else {
						assert_true(fabs(v1) < 256 * DBL_EPSILON);
					}
				}
			}
		}
		
		n = MESSAGES_COUNT(it);
		for (i = 0; i < n; i ++) {
			msg = MESSAGES_VAL(it, i);
			frame_add(&frame, msg);
		}
	}
	
	svector_deinit(&e_j);
	vector_deinit(&mean);
	vector_deinit(&y);
	matrix_deinit(&var);
	matrix_deinit(&var0);	
}



int main(int argc, char **argv)
{
	UnitTest tests[] = {
		unit_test_setup(enron_suite, enron_setup_fixture),
		unit_test_setup_teardown(test_probs, basic_setup, basic_teardown),
		unit_test_setup_teardown(test_mean, basic_setup, basic_teardown),
		unit_test_setup_teardown(test_var, basic_setup, basic_teardown),		
		unit_test_teardown(enron_suite, enron_teardown_fixture),
	};
	return run_tests(tests);
}
