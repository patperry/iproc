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
#include "blas.h"
#include "matrixutil.h"
#include "lapack.h"
#include "testutil.h"
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
static size_t ntrait;
static double *traits;
static const char * const *trait_names;
static struct messages messages;
static struct frame frame;
static struct recv_model model;
static struct recv_coefs coefs;
static struct recv_loglik loglik;


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

	size_t i, dim = coefs.dim;
	for (i = 0; i < dim; i++) {
		double val = (i + 1 % 5 == 0 ? -2.0 :
			      i + 2 % 5 == 1 ?  1.0 :
			      i + 3 % 5 == 2 ? -1.0 :
			      i + 7 % 5 == 3 ?  2.0 : 0.0);
		coefs.all[i] =  val;
	}

	recv_model_init(&model, &frame, &coefs);
	recv_loglik_init(&loglik, &model);
}

static void hard_setup()
{	
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

	size_t i, dim = coefs.dim;
	for (i = 0; i < dim; i++) {
		double val = (i + 1 % 7 == 0 ?  0.1 :
			      i + 2 % 7 == 1 ?  0.3 :
			      i + 3 % 7 == 2 ? -0.2 :
			      i + 5 % 7 == 4 ? -10 :
			      i + 11 % 7 == 6 ? +10 : 0.0);
		coefs.all[i] = val;
	}

	recv_model_init(&model, &frame, &coefs);
	recv_loglik_init(&loglik, &model);
}

static void teardown()
{
	recv_loglik_deinit(&loglik);
	recv_model_deinit(&model);
	recv_coefs_deinit(&coefs);
	frame_deinit(&frame);
}


static void test_count()
{
	struct messages_iter it;
	const struct message *msg = NULL;
	double t;
	size_t itie, ntie, nmsg;
	size_t count0, count1;
	size_t last_count0, last_count1;

	count0 = 0;
	nmsg = 0;

	MESSAGES_FOREACH(it, &messages) {
		printf("."); fflush(stdout);
		t = MESSAGES_TIME(it);

		frame_advance(&frame, t);

		ntie = MESSAGES_COUNT(it);
		for (itie = 0; itie < ntie; itie ++) {
			msg = MESSAGES_VAL(it, itie);
			frame_add(&frame, msg);
			recv_loglik_add(&loglik, &frame, msg);
			nmsg += msg->nto;

			if (nmsg > 1000)
				goto out;

			last_count0 = msg->nto;
			last_count1 = recv_loglik_last_count(&loglik);
			assert_int_equal(last_count0, last_count1);

			count0 = count0 + last_count0;
			count1 = recv_loglik_count(&loglik);
			assert_int_equal(count0, count1);
		}
	}
out:
	return;
}


static void test_dev()
{
	struct messages_iter it;
	const struct message *msg = NULL;
	const struct catdist1 *dist = NULL;
	double t;
	size_t i, itie, ntie, nmsg;
	double last_dev0, last_dev1;
	double dev0, dev1, dev_old;
	
	dev_old = 0.0;
	
	nmsg = 0;

	MESSAGES_FOREACH(it, &messages) {
		printf("."); fflush(stdout);
		t = MESSAGES_TIME(it);
		
		frame_advance(&frame, t);			
		
		ntie = MESSAGES_COUNT(it);
		for (itie = 0; itie < ntie; itie ++) {
			msg = MESSAGES_VAL(it, itie);
			frame_add(&frame, msg);
			recv_loglik_add(&loglik, &frame, msg);
			nmsg += msg->nto;
			
			if (nmsg > 1000)
				goto out;
			
			last_dev0 = 0.0;
			dist = recv_model_dist(&model, msg->from);
			for (i = 0; i < msg->nto; i++) {
				last_dev0 += -2 * catdist1_lprob(dist, msg->to[i]);
			}
			
			last_dev1 = recv_loglik_last_dev(&loglik);
			assert_real_identical(last_dev0, last_dev1);
			
			dev0 = dev_old + last_dev0;
			dev1 = recv_loglik_dev(&loglik);
			assert_real_approx(dev0, dev1);

			dev_old = dev1;
		}
	}
out:
	return;
}


static void test_mean()
{
	struct messages_iter it;
	const struct message *msg = NULL;
	const struct catdist1 *dist = NULL;
	double t;
	size_t i, itie, ntie, nmsg;
	double probs[nrecv];
	struct recv_coefs last_mean0, last_mean1;
	struct recv_coefs mean0, mean1;

	recv_coefs_init(&mean0, &frame);
	recv_coefs_init(&mean1, &frame);
	recv_coefs_init(&last_mean0, &frame);
	recv_coefs_init(&last_mean1, &frame);
	size_t dim = mean0.dim;

	const struct design *r = frame_recv_design(&frame);
	const struct design2 *d = frame_dyad_design(&frame);

	nmsg = 0;

	memset(mean0.all, 0, dim * sizeof(*mean0.all));

	MESSAGES_FOREACH(it, &messages) {
		printf("."); fflush(stdout);
		t = MESSAGES_TIME(it);

		frame_advance(&frame, t);

		ntie = MESSAGES_COUNT(it);
		for (itie = 0; itie < ntie; itie ++) {
			msg = MESSAGES_VAL(it, itie);
			frame_add(&frame, msg);
			recv_loglik_add(&loglik, &frame, msg);
			nmsg += msg->nto;

			if (nmsg > 1000)
				goto out;

			dist = recv_model_dist(&model, msg->from);
			for (i = 0; i < nrecv; i++) {
				probs[i] = catdist1_cached_prob(dist, i);
			}
			design_tmul(msg->nto, r, probs, 0.0, &last_mean0.recv);
			design2_tmul(msg->nto, d, msg->from, probs, 0.0, &last_mean0.dyad);

			memset(last_mean1.all, 0, dim * sizeof(*last_mean1.all));
			recv_loglik_axpy_last_mean(1.0, &loglik, &last_mean1);

			for (i = 0; i < dim; i++) {
				assert_real_approx(last_mean0.all[i], last_mean1.all[i]);
			}

			blas_daxpy(dim, 1.0, last_mean0.all, 1, mean0.all, 1);

			memset(mean1.all, 0, dim * sizeof(*mean1.all));
			recv_loglik_axpy_mean(1.0, &loglik, &mean1);

			for (i = 0; i < dim; i++) {
				assert_real_approx(mean0.all[i], mean1.all[i]);
			}

			blas_dcopy(dim, mean1.all, 1, mean0.all, 1);
		}
	}

	recv_coefs_deinit(&last_mean1);
	recv_coefs_deinit(&last_mean0);
	recv_coefs_deinit(&mean1);
	recv_coefs_deinit(&mean0);
out:
	return;
}


static void test_score()
{
	struct messages_iter it;
	const struct message *msg = NULL;
	double t;
	size_t i, itie, ntie, nmsg;
	struct recv_coefs last_score0, last_score1;
	struct recv_coefs score0, score1;

	recv_coefs_init(&score0, &frame);
	recv_coefs_init(&score1, &frame);
	recv_coefs_init(&last_score0, &frame);
	recv_coefs_init(&last_score1, &frame);
	size_t dim = score0.dim;

	const struct design *r = frame_recv_design(&frame);
	const struct design2 *d = frame_dyad_design(&frame);

	nmsg = 0;

	memset(score0.all, 0, dim * sizeof(*score0.all));

	MESSAGES_FOREACH(it, &messages) {
		printf("."); fflush(stdout);
		t = MESSAGES_TIME(it);

		frame_advance(&frame, t);

		ntie = MESSAGES_COUNT(it);
		for (itie = 0; itie < ntie; itie ++) {
			msg = MESSAGES_VAL(it, itie);
			frame_add(&frame, msg);
			recv_loglik_add(&loglik, &frame, msg);
			nmsg += msg->nto;

			if (nmsg > 1000)
				goto out;

			memset(last_score0.all, 0, dim * sizeof(*last_score0.all));
			for (i = 0; i < msg->nto; i++) {
				design_axpy(1.0, r, msg->to[i], &last_score0.recv);
				design2_axpy(1.0, d, msg->from, msg->to[i], &last_score0.dyad);
			}
			recv_loglik_axpy_last_mean(-1.0, &loglik, &last_score0);

			memset(last_score1.all, 0, dim * sizeof(*last_score1.all));
			recv_loglik_axpy_last_score(1.0, &loglik, &last_score1);

			for (i = 0; i < dim; i++) {
				assert_real_approx(last_score0.all[i], last_score1.all[i]);
			}

			blas_daxpy(dim, 1.0, last_score0.all, 1, score0.all, 1);

			memset(score1.all, 0, dim * sizeof(*score1.all));
			recv_loglik_axpy_score(1.0, &loglik, &score1);

			for (i = 0; i < dim; i++) {
				assert_real_approx(score0.all[i], score1.all[i]);
			}

			blas_dcopy(dim, score1.all, 1, score0.all, 1);
		}
	}

	recv_coefs_deinit(&last_score1);
	recv_coefs_deinit(&last_score0);
	recv_coefs_deinit(&score1);
	recv_coefs_deinit(&score0);
out:
	return;
}



/*

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
*/

int main()
{
	UnitTest tests[] = {
		unit_test_setup(enron_suite, enron_setup_fixture),
		unit_test_setup_teardown(test_count, basic_setup, teardown),
		unit_test_setup_teardown(test_count, hard_setup, teardown),
		unit_test_setup_teardown(test_dev, basic_setup, teardown),
		unit_test_setup_teardown(test_dev, hard_setup, teardown),
		unit_test_setup_teardown(test_mean, basic_setup, teardown),
		unit_test_setup_teardown(test_mean, hard_setup, teardown),
		unit_test_setup_teardown(test_score, basic_setup, teardown),
		unit_test_setup_teardown(test_score, hard_setup, teardown),
		//unit_test_setup_teardown(test_imat, basic_setup, teardown),
		//unit_test_setup_teardown(test_imat, hard_setup, teardown),
		unit_test_teardown(enron_suite, enron_teardown_fixture),
		
	};
	return run_tests(tests);
}
