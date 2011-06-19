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
#include "loglik.h"


static struct actors senders;
static struct actors receivers;
static struct vector intervals;
static struct messages messages;
static struct design design;
static struct frame frame;
static struct vector coefs;
static struct model model;
static struct recv_loglik recv_loglik;


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
	design_add_recv_var(&design, RECV_VAR_NRECV);
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
	recv_loglik_init(&recv_loglik, &model);
}

static void basic_teardown(void **state)
{
	recv_loglik_deinit(&recv_loglik);
	model_deinit(&model);
	vector_deinit(&coefs);
	frame_deinit(&frame);
	vector_deinit(&intervals);	
	design_deinit(&design);
}

static void test_dev(void **state)
{
	struct recv_model *rm;
	struct messages_iter it;
	const struct message *msg = NULL;
	double t;
	ssize_t i, itie, ntie, nrecv, n;
	double last_dev0, last_dev1;
	double mean_dev0, mean_dev1, mean_dev_old = 0.0;
	
	nrecv = design_recv_count(&design);

	n = 0;

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
		
		ntie = MESSAGES_COUNT(it);
		for (itie = 0; itie < ntie; itie ++) {
			msg = MESSAGES_VAL(it, itie);
			frame_add(&frame, msg);
			recv_loglik_add(&recv_loglik, &frame, msg);
			n += msg->nto;
			
			rm = model_recv_model(&model, &frame, msg->from);
			
			last_dev0 = 0.0;
			for (i = 0; i < msg->nto; i++) {
				last_dev0 += -2 * recv_model_logprob(rm, msg->to[i]);
			}
			
			last_dev1 = recv_loglik_last_dev(&recv_loglik);

			assert_in_range(double_eqrel(last_dev0, last_dev1), 53, DBL_MANT_DIG);
			
			
			mean_dev0 = mean_dev_old + msg->nto * (((last_dev0 / msg->nto) - mean_dev_old) / n);
			mean_dev1 = recv_loglik_avg_dev(&recv_loglik);
			assert_in_range(double_eqrel(mean_dev0, mean_dev1), 49, DBL_MANT_DIG);
			mean_dev_old = mean_dev1;
		}
	}
}


static void test_score(void **state)
{
	struct recv_model *rm;
	struct vector probs, mean0, mean1, avg_mean0, avg_mean1, diff;
	struct messages_iter it;
	const struct message *msg = NULL;
	double t;
	ssize_t isend;
	ssize_t itie, ntie;
	ssize_t index, dim = design_recv_dim(&design);
	ssize_t n;
	
	vector_init(&probs, design_recv_count(&design));
	vector_init(&mean0, design_recv_dim(&design));
	vector_init(&mean1, design_recv_dim(&design));
	vector_init(&avg_mean0, design_recv_dim(&design));
	vector_init(&avg_mean1, design_recv_dim(&design));
	vector_init(&diff, design_recv_dim(&design));
	
	MESSAGES_FOREACH(it, &messages) {
		printf("."); fflush(stdout);
		t = MESSAGES_TIME(it);
		
		while (frame_next_change(&frame) <= t) {
			model_update(&model, &frame);
			frame_advance(&frame);			
		}
		if (frame_time(&frame) < t) {
			model_update(&model, &frame);
			frame_advance_to(&frame, t);			
		}
		
		ntie = MESSAGES_COUNT(it);
		for (itie = 0; itie < ntie; itie ++) {
			msg = MESSAGES_VAL(it, itie);
			frame_add(&frame, msg);			
			recv_loglik_add(&recv_loglik, &frame, msg);
			n = recv_loglik_count(&recv_loglik);
			isend = msg->from;
			rm = model_recv_model(&model, &frame, isend);
			
			vector_fill(&probs, 0.0);
			recv_model_axpy_probs(1.0, rm, &probs);
			
			frame_recv_mul(msg->nto, TRANS_TRANS, &frame, isend, &probs,
				       0.0, &mean0);
			
			vector_fill(&mean1, 0.0);
			recv_loglik_axpy_last_mean(1.0, &recv_loglik, &mean1);
			
			for (index = 0; index < dim; index++) {
				double x0 = vector_item(&mean0, index);
				double x1 = vector_item(&mean1, index);				
				assert_in_range(double_eqrel(x0, x1), 40, DBL_MANT_DIG);
			}
			
			vector_assign_copy(&diff, &avg_mean0);
			vector_axpy(-1.0/msg->nto, &mean0, &diff);
			vector_axpy(-((double)msg->nto) / n, &diff, &avg_mean0);
			
			vector_fill(&avg_mean1, 0.0);
			recv_loglik_axpy_avg_mean(1.0, &recv_loglik, &avg_mean1);

			for (index = 0; index < dim; index++) {
				double x0 = vector_item(&avg_mean0, index);
				double x1 = vector_item(&avg_mean1, index);
				assert(double_eqrel(x0, x1) >= 40);
				assert_in_range(double_eqrel(x0, x1), 48, DBL_MANT_DIG);
			}

			vector_assign_copy(&avg_mean0, &avg_mean1);
		}
	}
	
	vector_deinit(&diff);
	vector_deinit(&avg_mean1);	
	vector_deinit(&avg_mean0);
	vector_deinit(&mean0);
	vector_deinit(&mean1);	
	vector_deinit(&probs);	
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
		// unit_test_setup_teardown(test_dev, basic_setup, basic_teardown),
		unit_test_setup_teardown(test_score, basic_setup, basic_teardown),
		// unit_test_setup_teardown(test_var, basic_setup, basic_teardown),		
		unit_test_teardown(enron_suite, enron_teardown_fixture),
	};
	return run_tests(tests);
}
