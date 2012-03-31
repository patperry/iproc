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


static size_t nsend;
static size_t nrecv;
static size_t ncohort;
static size_t ntrait;
static size_t *cohorts;
static double *traits;
static const char * const *cohort_names;
static const char * const *trait_names;
static struct vector intervals;
static struct messages messages;
static struct design *design;
static struct frame frame;
static struct matrix coefs;
static struct recv_model model;


static void enron_setup_fixture()
{
	print_message("Enron\n");
	print_message("-----\n");
	enron_employees_init(&nsend, &cohorts, &ncohort, 
			     &cohort_names,
			     &traits, &ntrait, &trait_names);
	enron_messages_init(&messages, -1);
	nrecv = nsend;
}

static void enron_teardown_fixture()
{
	free(cohorts);
	free(traits);
	messages_deinit(&messages);
	print_message("\n\n");
}

static void basic_setup()
{
	size_t i, c;
	double intvls[3] = {
		112.50,  450.00, 1800.00,
	};
	struct vector vintvls = vector_make(intvls, 3);
	int has_effects = 0;
	int has_loops = 0;
	vector_init(&intervals, 3);
	vector_assign_copy(&intervals, &vintvls);
	frame_init(&frame, nsend, nrecv, has_loops, vector_to_ptr(&intervals),
                    vector_dim(&intervals));
	design = frame_recv_design(&frame);
	design_set_has_effects(design, has_effects);
	design_set_traits(design, traits, ntrait, trait_names);
	design_add_dvar(design, RECV_VAR_NRECV, NULL);
	matrix_init(&coefs, design_dim(design), ncohort);
	
	for (c = 0; c < (size_t)matrix_ncol(&coefs); c++) {
		for (i = 0; i < (size_t)matrix_nrow(&coefs); i++) {
			double val = (i + 2 * c % 5 == 0 ? -2.0 :
				      i + 3 * c % 5 == 1 ?  1.0 :
				      i + 7 * c % 5 == 2 ? -1.0 :
				      i + 11 * c % 5 == 3 ?  2.0 : 0.0);
			matrix_set_item(&coefs, i, c, val);
		}
	}
	
	recv_model_init(&model, &frame, ncohort, cohorts, &coefs);
}

static void teardown()
{
	recv_model_deinit(&model);
	matrix_deinit(&coefs);
	frame_deinit(&frame);
	vector_deinit(&intervals);	
}


static void hard_setup()
{
	size_t i, c;
	double intvls[3] = {
		112.50,  450.00, 1800.00,
	};
	struct vector vintvls = vector_make(intvls, 3);
	int has_effects = 0;
	int has_loops = 0;
	vector_init(&intervals, 3);
	vector_assign_copy(&intervals, &vintvls);
	frame_init(&frame, nsend, nrecv, has_loops, vector_to_ptr(&intervals),
                    vector_dim(&intervals));
	design = frame_recv_design(&frame);
	design_set_has_effects(design, has_effects);
	design_set_traits(design, traits, ntrait, trait_names);
	design_add_dvar(design, RECV_VAR_NRECV, NULL);
	matrix_init(&coefs, design_dim(design), ncohort);
	
	for (c = 0; c < (size_t)matrix_ncol(&coefs); c++) {
		for (i = 0; i < (size_t)matrix_nrow(&coefs); i++) {
			double val = (i + 2 * c % 7 == 0 ?  0.1 :
				      i + 3 * c % 7 == 1 ?  0.3 :
				      i + 5 * c % 7 == 2 ? -0.2 :
				      i + 11 * c % 7 == 4 ? -10 :
				      i % 7 == 6 ? +10 : 0.0);
			matrix_set_item(&coefs, i, c, val);
		}
	}
	
	recv_model_init(&model, &frame, ncohort, cohorts, &coefs);
}


static void test_probs()
{
	struct vector eta, probs, logprobs, y;
	struct messages_iter it;
	const struct message *msg = NULL;
	double t;
	size_t isend, jrecv, nrecv;
	size_t i, n;
	
	nrecv = design_count(design);
	vector_init(&eta, nrecv);	
	vector_init(&probs, nrecv);
	vector_init(&logprobs, nrecv);
	vector_init(&y, nrecv);

	double alpha = 2.0;
	double y0 = 3.14;

	size_t minprec = DBL_MANT_DIG;
	
	MESSAGES_FOREACH(it, &messages) {
		// fprintf(stderr, "."); fflush(stderr);
		t = MESSAGES_TIME(it);
		
		frame_advance(&frame, t);			
		
		n = MESSAGES_COUNT(it);
		for (i = 0; i < n; i ++) {
			msg = MESSAGES_VAL(it, i);
			isend = msg->from;
			size_t c = cohorts[isend];
			const struct vector col = vector_make(matrix_col(&coefs, c), matrix_nrow(&coefs));
			frame_recv_mul(1.0, BLAS_NOTRANS, &frame, isend, &col, 0.0, &eta);
			
			if (!frame_has_loops(&frame))
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
					minprec = MIN(minprec, (size_t)double_eqrel(p0, p1));
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
