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
#include "xalloc.h"

#include "ieee754.h"
#include "enron.h"
#include "messages.h"
#include "design.h"
#include "vars.h"
#include "frame.h"
#include "recv_model.h"
#include "recv_loglik.h"


static size_t enron_nactor;
static size_t enron_ncohort;
static size_t *enron_cohorts;
static struct matrix enron_traits;
static const char * const *enron_cohort_names;
static const char * const *enron_trait_names;


static size_t nsend;
static size_t nrecv;
static size_t ncohort;
static size_t *cohorts;
static struct matrix recv_traits;
static const char * const *recv_trait_names;
static struct vector intervals;
static struct messages messages;
static struct design design;
static struct frame frame;
static struct matrix coefs;
static struct recv_model model;
static struct recv_loglik recv_loglik;


static void enron_setup_fixture()
{
	print_message("Enron\n");
	print_message("-----\n");
	enron_employees_init(&enron_nactor,
			     &enron_ncohort, &enron_cohorts,
			     &enron_cohort_names,
			     &enron_traits,
			     &enron_trait_names);
	nsend = enron_nactor;
	nrecv = enron_nactor;
	ncohort = enron_ncohort;
	cohorts = enron_cohorts;
	matrix_init_copy(&recv_traits, BLAS_NOTRANS, &enron_traits);
	recv_trait_names = enron_trait_names;
	enron_messages_init(&messages, -1);
}

static void enron_teardown_fixture()
{
	free(enron_cohorts);
	messages_deinit(&messages);
	matrix_deinit(&recv_traits);
	matrix_deinit(&enron_traits);
	print_message("\n\n");
}

static void basic_setup()
{
	size_t c, i;
	double intvls[3] = {
		112.50,  450.00, 1800.00,
	};
	struct vector vintvls = vector_make(intvls, 3);
	bool has_reffects = false;
	bool has_loops = false;
	vector_init(&intervals, 3);
	vector_assign_copy(&intervals, &vintvls);
	design_init(&design, nsend, nrecv, &recv_traits, recv_trait_names, &intervals);
	design_set_loops(&design, has_loops);
	design_set_recv_effects(&design, has_reffects);
	design_add_recv_var(&design, RECV_VAR_NRECV, NULL);
	frame_init(&frame, &design);
	matrix_init(&coefs, design_recv_dim(&design), ncohort);
	
	for (c = 0; c < (size_t)matrix_ncol(&coefs); c++) {	
		for (i = 0; i < (size_t)matrix_nrow(&coefs); i++) {
			double val = (i + (c + 1) % 5 == 0 ? -2.0 :
				      i + 2 * (c + 1) % 5 == 1 ?  1.0 :
				      i + 3 * (c + 1) % 5 == 2 ? -1.0 :
				      i + 7 * (c + 1) % 5 == 3 ?  2.0 : 0.0);
			matrix_set_item(&coefs, i, c, val);
		}
	}
	
	recv_model_init(&model, &frame, ncohort, cohorts, &coefs);
	recv_loglik_init(&recv_loglik, &model);
}

static void hard_setup()
{	
	size_t i, c;
	double intvls[3] = {
		112.50,  450.00, 1800.00,
	};
	struct vector vintvls = vector_make(intvls, 3);
	bool has_reffects = false;
	bool has_loops = false;
	vector_init(&intervals, 3);
	vector_assign_copy(&intervals, &vintvls);
	design_init(&design, nsend, nrecv, &recv_traits, recv_trait_names, &intervals);
	design_set_loops(&design, has_loops);
	design_set_recv_effects(&design, has_reffects);
	design_add_recv_var(&design, RECV_VAR_NRECV, NULL);
	frame_init(&frame, &design);
	matrix_init(&coefs, design_recv_dim(&design), ncohort);
	for (c = 0; c < (size_t)matrix_ncol(&coefs); c++) {	
		for (i = 0; i < (size_t)matrix_nrow(&coefs); i++) {
			double val = (i + (c + 1) % 7 == 0 ?  0.1 :
				      i + 2 * (c + 1) % 7 == 1 ?  0.3 :
				      i + 3 * (c + 1) % 7 == 2 ? -0.2 :
				      i + 5 * (c + 1) % 7 == 4 ? -10 :
				      i + 11 * (c + 1) % 7 == 6 ? +10 : 0.0);
			matrix_set_item(&coefs, i, c, val);
		}
	}

	recv_model_init(&model, &frame, ncohort, cohorts, &coefs);
	recv_loglik_init(&recv_loglik, &model);
}

static void teardown()
{
	recv_loglik_deinit(&recv_loglik);
	recv_model_deinit(&model);
	matrix_deinit(&coefs);
	frame_deinit(&frame);
	vector_deinit(&intervals);	
	design_deinit(&design);
}

static void test_dev()
{
	struct messages_iter it;
	const struct message *msg = NULL;
	double t;
	size_t i, itie, ntie, nrecv, nmsg;
	double last_dev0, last_dev1;
	double mean_dev0, mean_dev1, old;
	
	struct vector mean_dev_old;
	vector_init(&mean_dev_old, ncohort);
	
	nrecv = design_recv_count(&design);

	nmsg = 0;

	MESSAGES_FOREACH(it, &messages) {
		printf("."); fflush(stdout);
		t = MESSAGES_TIME(it);
		
		frame_advance(&frame, t);			
		
		ntie = MESSAGES_COUNT(it);
		for (itie = 0; itie < ntie; itie ++) {
			msg = MESSAGES_VAL(it, itie);
			frame_add(&frame, msg);
			recv_loglik_add(&recv_loglik, &frame, msg);
			nmsg += msg->nto;
			
			if (nmsg > 1000)
				goto out;
			
			last_dev0 = 0.0;
			for (i = 0; i < msg->nto; i++) {
				last_dev0 += -2 * recv_model_logprob(&model, msg->from, msg->to[i]);
			}
			
			last_dev1 = recv_loglik_last_dev(&recv_loglik);

			assert_in_range(double_eqrel(last_dev0, last_dev1), 53, DBL_MANT_DIG);
			
			size_t c = cohorts[msg->from];
			size_t n = recv_loglik_count(&recv_loglik, c);
			old = vector_item(&mean_dev_old, c);
			mean_dev0 = old + msg->nto * (((last_dev0 / msg->nto) - old) / n);
			mean_dev1 = recv_loglik_avg_dev(&recv_loglik, c);
			assert_in_range(double_eqrel(mean_dev0, mean_dev1), 48, DBL_MANT_DIG);
			vector_set_item(&mean_dev_old, c, mean_dev1);
		}
	}
out:
	vector_deinit(&mean_dev_old);
	return;
}


static void test_mean()
{
	struct vector probs, mean0, mean1, avg_mean1, diff;
	struct matrix avg_mean0;
	struct messages_iter it;
	const struct message *msg = NULL;
	double t;
	size_t isend;
	size_t itie, ntie;
	size_t index, dim = design_recv_dim(&design);
	size_t nmsg;
	
	vector_init(&probs, design_recv_count(&design));
	vector_init(&mean0, design_recv_dim(&design));
	vector_init(&mean1, design_recv_dim(&design));
	matrix_init(&avg_mean0, design_recv_dim(&design), ncohort);
	vector_init(&avg_mean1, design_recv_dim(&design));
	vector_init(&diff, design_recv_dim(&design));
	
	nmsg = 0;
	
	MESSAGES_FOREACH(it, &messages) {
		t = MESSAGES_TIME(it);

		frame_advance(&frame, t);			
		
		ntie = MESSAGES_COUNT(it);
		for (itie = 0; itie < ntie; itie ++) {
			msg = MESSAGES_VAL(it, itie);
			frame_add(&frame, msg);			
			recv_loglik_add(&recv_loglik, &frame, msg);

			nmsg += msg->nto;
		
			if (nmsg > 1000)
				goto out;

			isend = msg->from;
			size_t c = recv_model_cohort(&model, isend);
			size_t n = recv_loglik_count(&recv_loglik, c);			
			
			vector_fill(&probs, 0.0);
			recv_model_axpy_probs(1.0, &model, isend, &probs);
			
			frame_recv_mul(msg->nto, BLAS_TRANS, &frame, isend, &probs,
				       0.0, &mean0);
			
			vector_fill(&mean1, 0.0);
			recv_loglik_axpy_last_mean(1.0, &recv_loglik, &mean1);
			
			for (index = 0; index < dim; index++) {
				double x0 = vector_item(&mean0, index);
				double x1 = vector_item(&mean1, index);	
				assert(double_eqrel(x0, x1) >= 40);
				assert_in_range(double_eqrel(x0, x1), 40, DBL_MANT_DIG);
			}
			
			struct vector avg_mean0_c = matrix_col(&avg_mean0, c);
			vector_assign_copy(&diff, &avg_mean0_c);
			vector_axpy(-1.0/msg->nto, &mean0, &diff);
			vector_axpy(-((double)msg->nto) / n, &diff, &avg_mean0_c);
			
			vector_fill(&avg_mean1, 0.0);
			recv_loglik_axpy_avg_mean(1.0, &recv_loglik, c, &avg_mean1);

			for (index = 0; index < dim; index++) {
				double x0 = vector_item(&avg_mean0_c, index);
				double x1 = vector_item(&avg_mean1, index);
				
				if (fabs(x0) >= 5e-4) {
					assert(double_eqrel(x0, x1) >= 47);
					assert_in_range(double_eqrel(x0, x1), 47, DBL_MANT_DIG);
					
				} else {
					assert(fabs(x0 - x1) < sqrt(DBL_EPSILON));
					assert_true(fabs(x0 - x1) < sqrt(DBL_EPSILON));
				}
			}

			vector_assign_copy(&avg_mean0_c, &avg_mean1);
		}
	}
out:
	vector_deinit(&diff);
	vector_deinit(&avg_mean1);	
	matrix_deinit(&avg_mean0);
	vector_deinit(&mean0);
	vector_deinit(&mean1);	
	vector_deinit(&probs);	
}


static void test_score()
{
	struct vector score0, score1, avg_score1, diff;
	struct matrix avg_score0;
	struct svector nrecv;
	struct messages_iter it;
	const struct message *msg = NULL;
	double t;
	size_t isend, c, ito;
	size_t itie, ntie;
	size_t index, dim = design_recv_dim(&design);
	size_t n, nmsg;
	
	svector_init(&nrecv, design_recv_count(&design));
	vector_init(&score0, design_recv_dim(&design));
	vector_init(&score1, design_recv_dim(&design));
	matrix_init(&avg_score0, design_recv_dim(&design), recv_model_cohort_count(&model));
	vector_init(&avg_score1, design_recv_dim(&design));
	vector_init(&diff, design_recv_dim(&design));
	
	nmsg = 0;
	
	MESSAGES_FOREACH(it, &messages) {
		t = MESSAGES_TIME(it);
		
		frame_advance(&frame, t);			
		
		ntie = MESSAGES_COUNT(it);
		for (itie = 0; itie < ntie; itie ++) {
			msg = MESSAGES_VAL(it, itie);
			frame_add(&frame, msg);			
			recv_loglik_add(&recv_loglik, &frame, msg);
			nmsg++;

			if (nmsg > 1000)
				goto out;
			
			isend = msg->from;
			c = recv_model_cohort(&model, isend);
			n = recv_loglik_count(&recv_loglik, c);			
			
			svector_clear(&nrecv);
			for (ito = 0; ito < msg->nto; ito++) {
				*svector_item_ptr(&nrecv, msg->to[ito]) += 1.0;
			}

			frame_recv_muls(1.0, BLAS_TRANS, &frame, isend, &nrecv,
					0.0, &score0);
			recv_loglik_axpy_last_mean(-1.0, &recv_loglik, &score0);
			
			vector_fill(&score1, 0.0);
			recv_loglik_axpy_last_score(1.0, &recv_loglik, &score1);
			
			for (index = 0; index < dim; index++) {
				double x0 = vector_item(&score0, index);
				double x1 = vector_item(&score1, index);				
				assert_in_range(double_eqrel(x0, x1), 40, DBL_MANT_DIG);
			}
			
			struct vector avg_score0_c = matrix_col(&avg_score0, c);
			vector_assign_copy(&diff, &avg_score0_c);
			vector_axpy(-1.0/msg->nto, &score0, &diff);
			vector_axpy(-((double)msg->nto) / n, &diff, &avg_score0_c);
			
			vector_fill(&avg_score1, 0.0);
			recv_loglik_axpy_avg_score(1.0, &recv_loglik, c, &avg_score1);
			
			for (index = 0; index < dim; index++) {
				double x0 = vector_item(&avg_score0_c, index);
				double x1 = vector_item(&avg_score1, index);
				if (fabs(x0) >= 5e-4) {
					assert_in_range(double_eqrel(x0, x1), 37, DBL_MANT_DIG);
				} else {
					assert_true(fabs(x0 - x1) < sqrt(DBL_EPSILON));
				}
			}
			
			vector_assign_copy(&avg_score0_c, &avg_score1);
		}
	}
out:
	vector_deinit(&diff);
	vector_deinit(&avg_score1);	
	matrix_deinit(&avg_score0);
	vector_deinit(&score0);
	vector_deinit(&score1);	
	svector_deinit(&nrecv);	
}


static void test_imat()
{
	struct vector mean, y;
	struct svector e_j;
	struct matrix imat0, imat1, diff, *avg_imat0, avg_imat1;
	struct messages_iter it;
	const struct message *msg = NULL;
	double t;
	size_t isend, jrecv, nrecv = design_recv_count(&design);
	size_t itie, ntie, c, n, nmsg;
	size_t index1, index2, dim = design_recv_dim(&design);
	
	vector_init(&mean, dim);
	vector_init(&y, dim);	
	matrix_init(&imat0, dim, dim);	
	matrix_init(&imat1, dim, dim);
	matrix_init(&diff, dim, dim);
	
	avg_imat0 = xcalloc(recv_model_cohort_count(&model), sizeof(*avg_imat0));
	for (c = 0; c < recv_model_cohort_count(&model); c++) {
		matrix_init(&avg_imat0[c], dim, dim);
	}
	matrix_init(&avg_imat1, dim, dim);
	svector_init(&e_j, nrecv);
	
	nmsg = 0;
	
	MESSAGES_FOREACH(it, &messages) {
		t = MESSAGES_TIME(it);
		
		frame_advance(&frame, t);			
		
		ntie = MESSAGES_COUNT(it);
		for (itie = 0; itie < ntie; itie ++) {
			printf("."); fflush(stdout);
			msg = MESSAGES_VAL(it, itie);
			frame_add(&frame, msg);
			recv_loglik_add(&recv_loglik, &frame, msg);

			nmsg++;
			if (nmsg > 1000)
				goto out;

			isend = msg->from;
			c = recv_model_cohort(&model, isend);
			n = recv_loglik_count(&recv_loglik, c);			
			
			vector_fill(&mean, 0.0);
			recv_loglik_axpy_last_mean(1.0 / msg->nto, &recv_loglik, &mean);
			matrix_fill(&imat1, 0.0);
			recv_loglik_axpy_last_imat(1.0, &recv_loglik, &imat1);
			
			matrix_fill(&imat0, 0.0);
			for (jrecv = 0; jrecv < nrecv; jrecv++) {
				double p = recv_model_prob(&model, isend, jrecv);
				vector_assign_copy(&y, &mean);				
				svector_set_basis(&e_j, jrecv);
				
				frame_recv_muls(1.0, BLAS_TRANS, &frame, isend, &e_j, -1.0, &y);
				matrix_update1(&imat0, p, &y, &y);
			}
			matrix_scale(&imat0, msg->nto);
			
			for (index2 = 0; index2 < dim; index2++) {
				for (index1 = 0; index1 < dim; index1++) {				
					double v0 = matrix_item(&imat0, index1, index2);
					double v1 = matrix_item(&imat1, index1, index2);
					//printf("v0: %.4f  v1: %.4f (%d)\n", v0, v1, double_eqrel(v0, v1));
					assert(double_eqrel(v0, v1) >= DBL_MANT_DIG / 2
					       || ((fabs(v0) < 1e-1) && fabs(v0 - v1) < sqrt(DBL_EPSILON)));
					if (abs(v0 - v1) >= sqrt(DBL_EPSILON)) {
						assert_in_range(double_eqrel(v0, v1), 50, DBL_MANT_DIG);
					}
				}
			}
			
			matrix_assign_copy(&diff, BLAS_NOTRANS, &avg_imat0[c]);
			matrix_axpy(-1.0/msg->nto, &imat0, &diff);
			matrix_axpy(-((double)msg->nto) / n, &diff, &avg_imat0[c]);
			
			matrix_fill(&avg_imat1, 0.0);
			recv_loglik_axpy_avg_imat(1.0, &recv_loglik, c, &avg_imat1);
			
			for (index2 = 0; index2 < dim; index2++) {
				for (index1 = 0; index1 < dim; index1++) {				
					double v0 = matrix_item(&avg_imat0[c], index1, index2);
					double v1 = matrix_item(&avg_imat1, index1, index2);
					//assert(double_eqrel(v0, v1) >= 37
					//       || ((fabs(v0) < 1e-1) && fabs(v0 - v1) < sqrt(DBL_EPSILON)));
					if (fabs(v0) >= 5e-4) {
						assert_in_range(double_eqrel(v0, v1), 37, DBL_MANT_DIG);
					} else {
						assert_true(fabs(v0 - v1) < sqrt(DBL_EPSILON));
					}
				}
			}
			
			matrix_assign_copy(&avg_imat0[c], BLAS_NOTRANS, &avg_imat1);

		}
		
	}
	
out:
	svector_deinit(&e_j);
	vector_deinit(&mean);
	vector_deinit(&y);
	matrix_deinit(&imat0);
	matrix_deinit(&imat1);
	matrix_deinit(&diff);
	for (c = 0; c < recv_model_cohort_count(&model); c++) {
		matrix_deinit(&avg_imat0[c]);
	}
	free(avg_imat0);
	matrix_deinit(&avg_imat1);
}


int main()
{
	UnitTest tests[] = {
		unit_test_setup(enron_suite, enron_setup_fixture),
		unit_test_setup_teardown(test_dev, basic_setup, teardown),
		unit_test_setup_teardown(test_dev, hard_setup, teardown),
		unit_test_setup_teardown(test_mean, basic_setup, teardown),
		unit_test_setup_teardown(test_mean, hard_setup, teardown),
		unit_test_setup_teardown(test_score, basic_setup, teardown),
		unit_test_setup_teardown(test_score, hard_setup, teardown),
		unit_test_setup_teardown(test_imat, basic_setup, teardown),
		unit_test_setup_teardown(test_imat, hard_setup, teardown),		
		unit_test_teardown(enron_suite, enron_teardown_fixture),
		
	};
	return run_tests(tests);
}
