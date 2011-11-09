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
#include "blas.h"
#include "matrixutil.h"
#include "lapack.h"
#include "xalloc.h"

#include "ieee754.h"
#include "enron.h"
#include "messages.h"
#include "design.h"
#include "vars.h"
#include "frame.h"
#include "recv_model.h"
#include "recv_loglik.h"




static size_t nsend;
static size_t nrecv;
static size_t ncohort;
static size_t ntrait;
static size_t *cohorts;
static struct dmatrix traits;
static const char * const *trait_names;
static const char * const *cohort_names;
static struct messages messages;
static struct design *design;
static struct frame frame;
static struct dmatrix coefs;
static struct recv_model model;
static struct recv_loglik recv_loglik;


static void enron_setup_fixture()
{
	print_message("Enron\n");
	print_message("-----\n");
	enron_employees_init(&nsend, &cohorts, &ncohort, &cohort_names,
			     &traits.data, &ntrait, &trait_names);
	traits.lda = MAX(1, nsend);
	nrecv = nsend;
	enron_messages_init(&messages, -1);
}

static void enron_teardown_fixture()
{
	free(cohorts);
	free(traits.data);
	messages_deinit(&messages);
	print_message("\n\n");
}

static void basic_setup()
{
	size_t c, i;
	double intvls[3] = {
		112.50,  450.00, 1800.00,
	};
	int has_effects = 0;
	int has_loops = 0;
	frame_init(&frame, nsend, nrecv, has_loops, intvls, 3);
	design = frame_recv_design(&frame);
        design_set_has_effects(design, has_effects);
	design_set_traits(design, ntrait, &traits, trait_names);
        design_add_dvar(design, RECV_VAR_NRECV, NULL);
	coefs.data = xcalloc(design_dim(design) * ncohort, sizeof(double));
	coefs.lda = MAX(1, design_dim(design));
	
	for (c = 0; c < ncohort; c++) {	
		for (i = 0; i < design_dim(design); i++) {
			double val = (i + (c + 1) % 5 == 0 ? -2.0 :
				      i + 2 * (c + 1) % 5 == 1 ?  1.0 :
				      i + 3 * (c + 1) % 5 == 2 ? -1.0 :
				      i + 7 * (c + 1) % 5 == 3 ?  2.0 : 0.0);
			MATRIX_ITEM(&coefs, i, c) =  val;
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
	int has_effects = 0;
	int has_loops = 0;
	frame_init(&frame, nsend, nrecv, has_loops, intvls, 3);
	design = frame_recv_design(&frame);
        design_set_has_effects(design, has_effects);
	design_set_traits(design, ntrait, &traits, trait_names);
        design_add_dvar(design, RECV_VAR_NRECV, NULL);
	coefs.data = xcalloc(design_dim(design) * ncohort, sizeof(double));
	coefs.lda = MAX(1, design_dim(design));

	for (c = 0; c < ncohort; c++) {	
		for (i = 0; i < design_dim(design); i++) {
			double val = (i + (c + 1) % 7 == 0 ?  0.1 :
				      i + 2 * (c + 1) % 7 == 1 ?  0.3 :
				      i + 3 * (c + 1) % 7 == 2 ? -0.2 :
				      i + 5 * (c + 1) % 7 == 4 ? -10 :
				      i + 11 * (c + 1) % 7 == 6 ? +10 : 0.0);
			MATRIX_ITEM(&coefs, i, c) = val;
		}
	}

	recv_model_init(&model, &frame, ncohort, cohorts, &coefs);
	recv_loglik_init(&recv_loglik, &model);
}

static void teardown()
{
	recv_loglik_deinit(&recv_loglik);
	recv_model_deinit(&model);
	free(coefs.data);
	frame_deinit(&frame);
}

static void test_dev()
{
	struct messages_iter it;
	const struct message *msg = NULL;
	double t;
	size_t i, itie, ntie, nrecv, nmsg;
	double last_dev0, last_dev1;
	double mean_dev0, mean_dev1, old;
	
	double *mean_dev_old = xcalloc(ncohort, sizeof(double));
	
	nrecv = design_count(design);

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
			old = mean_dev_old[c];
			mean_dev0 = old + msg->nto * (((last_dev0 / msg->nto) - old) / n);
			mean_dev1 = recv_loglik_avg_dev(&recv_loglik, c);
			assert_in_range(double_eqrel(mean_dev0, mean_dev1), 48, DBL_MANT_DIG);
			mean_dev_old[c] = mean_dev1;
		}
	}
out:
	free(mean_dev_old);
	return;
}


static void test_mean()
{
	double *probs, *mean0, *mean1, *avg_mean1, *diff;
	struct messages_iter it;
	const struct message *msg = NULL;
	double t;
	size_t isend;
	size_t itie, ntie;
	size_t index, dim = design_dim(design);
	size_t nmsg;
	
	probs = xmalloc(design_count(design) * sizeof(double));
	mean0 = xmalloc(dim * sizeof(double));
	mean1 = xmalloc(dim * sizeof(double));
	struct dmatrix avg_mean0 = { xcalloc(design_dim(design) * ncohort, sizeof(double)),
				     MAX(1, design_dim(design)) };
	avg_mean1 = xmalloc(dim * sizeof(double));
	diff = xmalloc(dim * sizeof(double));
	
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
			
			memset(probs, 0, design_count(design) * sizeof(double));
			recv_model_axpy_probs(1.0, &model, isend, probs);
			
			frame_recv_mul(msg->nto, BLAS_TRANS, &frame, isend, probs,
				       0.0, mean0);
			
			memset(mean1, 0, dim * sizeof(double));
			recv_loglik_axpy_last_mean(1.0, &recv_loglik, mean1);
			
			for (index = 0; index < dim; index++) {
				double x0 = mean0[index];
				double x1 = mean1[index];	
				assert(double_eqrel(x0, x1) >= 40);
				assert_in_range(double_eqrel(x0, x1), 40, DBL_MANT_DIG);
			}
			
			double *avg_mean0_c = MATRIX_COL(&avg_mean0, c);
			blas_dcopy(dim, avg_mean0_c, 1, diff, 1);
			blas_daxpy(dim, -1.0/msg->nto, mean0, 1, diff, 1);
			blas_daxpy(dim, -((double)msg->nto) / n, diff, 1, avg_mean0_c, 1);
			
			memset(avg_mean1, 0, dim * sizeof(double));
			recv_loglik_axpy_avg_mean(1.0, &recv_loglik, c, avg_mean1);

			for (index = 0; index < dim; index++) {
				double x0 = avg_mean0_c[index];
				double x1 = avg_mean1[index];
				
				if (fabs(x0) >= 5e-4) {
					assert(double_eqrel(x0, x1) >= 47);
					assert_in_range(double_eqrel(x0, x1), 47, DBL_MANT_DIG);
					
				} else {
					assert(fabs(x0 - x1) < sqrt(DBL_EPSILON));
					assert_true(fabs(x0 - x1) < sqrt(DBL_EPSILON));
				}
			}

			blas_dcopy(dim, avg_mean1, 1, avg_mean0_c, 1);
		}
	}
out:
	free(diff);
	free(avg_mean1);	
	free(avg_mean0.data);
	free(mean0);
	free(mean1);	
	free(probs);	
}


static void test_score()
{
	double *nrecv, *score0, *score1, *avg_score1, *diff;
	struct messages_iter it;
	const struct message *msg = NULL;
	double t;
	size_t isend, c, ito;
	size_t itie, ntie;
	size_t index, dim = design_dim(design);
	size_t n, nmsg;


	nrecv = xcalloc(design_count(design), sizeof(double));
	score0 = xmalloc(dim * sizeof(double));
	score1 = xmalloc(dim * sizeof(double));
	struct dmatrix avg_score0 = { xcalloc(design_dim(design) * recv_model_cohort_count(&model), sizeof(double)),
		MAX(1, design_dim(design)) };
	avg_score1 = xmalloc(dim * sizeof(double));
	diff = xmalloc(dim * sizeof(double));

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

			memset(nrecv, 0, design_count(design) * sizeof(double));
			for (ito = 0; ito < msg->nto; ito++) {
				nrecv[msg->to[ito]] += 1.0;
			}

			frame_recv_mul(1.0, BLAS_TRANS, &frame, isend, nrecv, 0.0, score0);
			recv_loglik_axpy_last_mean(-1.0, &recv_loglik, score0);

			memset(score1, 0, dim * sizeof(double));
			recv_loglik_axpy_last_score(1.0, &recv_loglik, score1);

			for (index = 0; index < dim; index++) {
				double x0 = score0[index];
				double x1 = score1[index];
				assert(double_eqrel(x0, x1) >= 40);
				assert_in_range(double_eqrel(x0, x1), 40, DBL_MANT_DIG);
			}
			
			double *avg_score0_c = MATRIX_COL(&avg_score0, c);
			blas_dcopy(dim, avg_score0_c, 1, diff, 1);
			blas_daxpy(dim, -1.0/msg->nto, score0, 1, diff, 1);
			blas_daxpy(dim, -((double)msg->nto) / n, diff, 1, avg_score0_c, 1);
			
			memset(avg_score1, 0, dim * sizeof(double));
			recv_loglik_axpy_avg_score(1.0, &recv_loglik, c, avg_score1);
			
			for (index = 0; index < dim; index++) {
				double x0 = avg_score0_c[index];
				double x1 = avg_score1[index];
				if (fabs(x0) >= 5e-4) {
					assert_in_range(double_eqrel(x0, x1), 37, DBL_MANT_DIG);
				} else {
					assert_true(fabs(x0 - x1) < sqrt(DBL_EPSILON));
				}
			}
			
			blas_dcopy(dim, avg_score1, 1, avg_score0_c, 1);
		}
	}
out:
	free(diff);
	free(avg_score1);	
	free(avg_score0.data);
	free(score0);
	free(score1);	
	free(nrecv);	
}


static void test_imat()
{
	double *mean, *y;

	struct dmatrix *avg_imat0;
	struct messages_iter it;
	const struct message *msg = NULL;
	double t;
	size_t isend, jrecv, nrecv = design_count(design);
	size_t itie, ntie, c, n, nmsg;
	size_t index1, index2, dim = design_dim(design);
	size_t dim2 = dim * dim;
	double one = 1.0;
	struct vpattern pat_j;
	pat_j.indx = &jrecv;
	pat_j.nz = 1;

	mean = xmalloc(dim * sizeof(double));
	y = xmalloc(dim * sizeof(double));
	struct dmatrix imat0 = { xcalloc(dim2, sizeof(double)), MAX(1, dim) };
	struct dmatrix imat1 = { xcalloc(dim2, sizeof(double)), MAX(1, dim) };
	struct dmatrix diff = { xcalloc(dim2, sizeof(double)), MAX(1, dim) };
	struct dmatrix avg_imat1 = { xcalloc(dim2, sizeof(double)), MAX(1, dim) };	

	avg_imat0 = xcalloc(recv_model_cohort_count(&model), sizeof(*avg_imat0));
	for (c = 0; c < recv_model_cohort_count(&model); c++) {
		avg_imat0[c] = (struct dmatrix) { xcalloc(dim2, sizeof(double)), MAX(1, dim) };
	}
	
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
			
			memset(mean, 0, dim * sizeof(double));
			recv_loglik_axpy_last_mean(1.0 / msg->nto, &recv_loglik, mean);
			matrix_dzero(dim, dim, &imat1);
			recv_loglik_axpy_last_imat(1.0, &recv_loglik, &imat1);
			
			matrix_dzero(dim, dim, &imat0);
			for (jrecv = 0; jrecv < nrecv; jrecv++) {
				double p = recv_model_prob(&model, isend, jrecv);
				blas_dcopy(dim, mean, 1, y, 1);
			
				
				frame_recv_muls(1.0, BLAS_TRANS, &frame, isend, &one, &pat_j, -1.0, y);
				blas_dger(dim, dim, p, y, 1, y, 1, &imat0);
			}
			matrix_dscal(dim, dim, msg->nto, &imat0);
			
			for (index2 = 0; index2 < dim; index2++) {
				for (index1 = 0; index1 < dim; index1++) {				
					double v0 = MATRIX_ITEM(&imat0, index1, index2);
					double v1 = MATRIX_ITEM(&imat1, index1, index2);
					//printf("v0: %.4f  v1: %.4f (%d)\n", v0, v1, double_eqrel(v0, v1));
					assert(double_eqrel(v0, v1) >= DBL_MANT_DIG / 2
					       || ((fabs(v0) < 1e-1) && fabs(v0 - v1) < sqrt(DBL_EPSILON)));
					if (abs(v0 - v1) >= sqrt(DBL_EPSILON)) {
						assert_in_range(double_eqrel(v0, v1), 50, DBL_MANT_DIG);
					}
				}
			}
			
			lapack_dlacpy(LA_COPY_ALL, dim, dim, &avg_imat0[c], &diff);
			matrix_daxpy(dim, dim, -1.0/msg->nto, &imat0, &diff);
			matrix_daxpy(dim, dim, -((double)msg->nto) / n, &diff, &avg_imat0[c]);
			
			matrix_dzero(dim, dim, &avg_imat1);
			recv_loglik_axpy_avg_imat(1.0, &recv_loglik, c, &avg_imat1);
			
			for (index2 = 0; index2 < dim; index2++) {
				for (index1 = 0; index1 < dim; index1++) {				
					double v0 = MATRIX_ITEM(&avg_imat0[c], index1, index2);
					double v1 = MATRIX_ITEM(&avg_imat1, index1, index2);
					//assert(double_eqrel(v0, v1) >= 37
					//       || ((fabs(v0) < 1e-1) && fabs(v0 - v1) < sqrt(DBL_EPSILON)));
					if (fabs(v0) >= 5e-4) {
						assert_in_range(double_eqrel(v0, v1), 37, DBL_MANT_DIG);
					} else {
						assert_true(fabs(v0 - v1) < sqrt(DBL_EPSILON));
					}
				}
			}
			
			lapack_dlacpy(LA_COPY_ALL, dim, dim, &avg_imat1, &avg_imat0[c]);

		}
		
	}
	
out:
	free(mean);
	free(y);
	free(imat0.data);
	free(imat1.data);
	free(diff.data);
	free(avg_imat1.data);	
	for (c = 0; c < recv_model_cohort_count(&model); c++) {
		free(avg_imat0[c].data);
	}
	free(avg_imat0);
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
