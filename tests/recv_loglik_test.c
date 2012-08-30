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

	size_t dim = recv_model_dim(&model);
	recv_coefs_init(&mean0, &frame);
	recv_coefs_init(&mean1, &frame);
	recv_coefs_init(&last_mean0, &frame);
	recv_coefs_init(&last_mean1, &frame);

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
out:
	recv_coefs_deinit(&last_mean1);
	recv_coefs_deinit(&last_mean0);
	recv_coefs_deinit(&mean1);
	recv_coefs_deinit(&mean0);
}


static void test_score()
{
	struct messages_iter it;
	const struct message *msg = NULL;
	double t;
	size_t i, itie, ntie, nmsg;
	struct recv_coefs last_score0, last_score1;
	struct recv_coefs score0, score1;

	size_t dim = recv_model_dim(&model);
	recv_coefs_init(&score0, &frame);
	recv_coefs_init(&score1, &frame);
	recv_coefs_init(&last_score0, &frame);
	recv_coefs_init(&last_score1, &frame);

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
out:
	recv_coefs_deinit(&last_score1);
	recv_coefs_deinit(&last_score0);
	recv_coefs_deinit(&score1);
	recv_coefs_deinit(&score0);
}


static void test_imat()
{
	struct messages_iter it;
	const struct message *msg = NULL;
	double t;
	size_t i, itie, ntie, nmsg;

	size_t dim = recv_model_dim(&model);
	size_t cov_dim = dim * (dim + 1) / 2;
	struct recv_coefs mean, diff;

	double *cov0 = xmalloc(cov_dim * sizeof(double));
	double *cov1 = xmalloc(cov_dim * sizeof(double));
	double *last_cov0 = xmalloc(cov_dim * sizeof(double));
	double *last_cov1 = xmalloc(cov_dim * sizeof(double));

	recv_coefs_init(&mean, &frame);
	recv_coefs_init(&diff, &frame);

	const struct design *r = frame_recv_design(&frame);
	const struct design2 *d = frame_dyad_design(&frame);
	const struct catdist1 *dist = NULL;

	nmsg = 0;

	memset(cov0, 0, cov_dim * sizeof(*cov0));

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

			/* compute mean */
			memset(mean.all, 0, dim * sizeof(double));
			recv_loglik_axpy_last_mean(1.0, &loglik, &mean);

			/* compute last_cov0 */
			dist = recv_model_dist(&model, msg->from);
			memset(last_cov0, 0, cov_dim * sizeof(double));
			for (i = 0; i < nrecv; i++) {
				memcpy(diff.all, mean.all, dim * sizeof(double));
				design_axpy(-1.0, r, i, &diff.recv);
				design2_axpy(-1.0, d, msg->from, i, &diff.dyad);

				double w = catdist1_cached_prob(dist, i);
				blas_dspr(BLAS_LOWER, dim, w, diff.all, 1, last_cov0);
			}
			blas_dscal(cov_dim, msg->nto, last_cov0, 1);

			/* compute last_cov1 */
			memset(last_cov1, 0, cov_dim * sizeof(double));
			recv_loglik_axpy_last_imat(1.0, &loglik, last_cov1);

			for (i = 0; i < cov_dim; i++) {
				//assert(double_eqrel(last_cov0[i], last_cov1[i]) >= 23);
				assert_real_approx(last_cov0[i], last_cov1[i]);
			}

			/*
			blas_daxpy(cov_dim, 1.0, last_cov0, 1, cov0, 1);

			memset(cov1, 0, cov_dim * sizeof(double));
			recv_loglik_axpy_imat(1.0, &loglik, cov1);

			for (i = 0; i < cov_dim; i++) {
				assert_real_approx(cov0[i], cov1[i]);
			}

			blas_dcopy(dim, cov1, 1, cov0, 1);*/
		}
	}
out:
	recv_coefs_deinit(&diff);
	recv_coefs_deinit(&mean);
	free(last_cov1);
	free(last_cov0);
	free(cov1);
	free(cov0);
}


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
		unit_test_setup_teardown(test_imat, basic_setup, teardown),
		unit_test_setup_teardown(test_imat, hard_setup, teardown),
		unit_test_teardown(enron_suite, enron_teardown_fixture),
		
	};
	return run_tests(tests);
}
