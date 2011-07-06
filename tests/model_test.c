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
	
	model_init(&model, &frame, &coefs);
}

static void teardown(void **state)
{
	model_deinit(&model);
	vector_deinit(&coefs);
	frame_deinit(&frame);
	vector_deinit(&intervals);	
	design_deinit(&design);
}


static void hard_setup(void **state)
{
	double coefs_data[] = {
	       -2.0078320213241120123852e-16,
	       7.0235805310703240459397e-16,
	       4.0871230686004442771612e-15,
	       3.5061694800573514011041e-15,
	       -9.7506233043687851203037e-15,
	       -1.4542111837460032219838e+00,
	       3.2844263208747004334498e+00,
	       3.4331083526169348107970e-01,
	       7.9362087235055411849061e-01,
	       -1.2066175732167991330179e-01,
	       -1.1378909670106971407932e+00,
	       4.0547542577788908690906e-01,
	       1.1758651538928710511556e+00,
	       -6.7915784985022098485530e-01,
	       6.7407617457656709980540e-01,
	       -2.0047092350113207559481e-01,
	       -5.9427807574077675181745e-03,
	       -3.2814501792059824758496e-01,
	       5.6490694754611092687213e-01,
	       8.4437967320357881773063e-02,
	       -1.3939268854476738468406e+00,
	       1.6815195024824614034031e-01,
	       3.9756417583138609073146e-01,
	       -4.5517302772080470152360e-01,
	       1.5696564637539311970471e+00,
	       4.3464158208588665743832e+00,
	       3.4377415537199795814161e+00,
	       2.3787403030778051515881e+00,
		6.1180077624527857624304e-03 };
	       
	enron_employees_init(&senders);
	enron_employees_init(&receivers);
	enron_messages_init(&messages);
	
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
		vector_set_item(&coefs, i, coefs_data[i]);
	}
	
	model_init(&model, &frame, &coefs);
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

	ssize_t minprec = DBL_MANT_DIG;
	
	MESSAGES_FOREACH(it, &messages) {
		t = MESSAGES_TIME(it);
		
		frame_advance_to(&frame, t);			
		
		n = MESSAGES_COUNT(it);
		for (i = 0; i < n; i ++) {
			msg = MESSAGES_VAL(it, i);
			isend = msg->from;
			rm = model_recv_model(&model, isend);
			
			frame_recv_mul(1.0, TRANS_NOTRANS, &frame, isend, &coefs, 0.0, &eta);
			
			if (!design_loops(&design))
				vector_set_item(&eta, isend, -INFINITY);
			
			vector_assign_copy(&logprobs, &eta);
			double etamax = vector_max(&eta);
			vector_shift(&logprobs, -etamax);
			double log_W = vector_log_sum_exp(&logprobs);
			vector_shift(&logprobs, -log_W);
			
			log_W += etamax;
			
			vector_assign_copy(&probs, &logprobs);
			vector_exp(&probs);
			
			vector_fill(&y, y0);
			recv_model_axpy_probs(alpha, rm, &y);
			
			for (jrecv = 0; jrecv < nrecv; jrecv++) {
				double lp0 = recv_model_logprob(rm, jrecv);
				double lp1 = vector_item(&logprobs, jrecv);
				if (fabs(lp0) >= 5e-4) {
					//minprec = MIN(minprec, double_eqrel(lp0, lp1));
					assert(double_eqrel(lp0, lp1) >= 39);
					assert_in_range(double_eqrel(lp0, lp1), 39, DBL_MANT_DIG);

				} else {
					assert_true(fabs(lp0 - lp1) < sqrt(DBL_EPSILON));
				}
				
				double p0 = recv_model_prob(rm, jrecv);
				double p1 = vector_item(&probs, jrecv);
				if (fabs(p0) >= 5e-4) {
					minprec = MIN(minprec, double_eqrel(p0, p1));
					assert(double_eqrel(p0, p1) >= 47);
					assert_in_range(double_eqrel(p0, p1), 47, DBL_MANT_DIG);
				} else {
					assert_true(fabs(p0 - p1) < sqrt(DBL_EPSILON));
				}
				
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
