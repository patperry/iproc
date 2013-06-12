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

#include "blas.h"
#include "cmockery.h"
#include "coreutil.h"
#include "ieee754.h"
#include "lapack.h"
#include "matrixutil.h"
#include "testutil.h"
#include "xalloc.h"

#include "../src/history.h"
#include "../src/design.h"
#include "../src/design2.h"
#include "../src/recv_model.h"
#include "../src/recv_loglik.h"

#include "fixtures/history.h"
#include "fixtures/design.h"
#include "fixtures/design2.h"
#include "fixtures/recv_model.h"



struct fixture {
	struct history_fixture h;
	struct design_fixture r;
	struct design2_fixture d;
	struct recv_model_fixture m;
};

struct fixture ctx;

#define NRECV (ctx.h.nrecv)
#define NSEND (ctx.h.nsend)
#define RECV_PARAMS (&ctx.m.params)
#define EXCLUDE_LOOPS (ctx.m.exclude_loops)


struct test {
	struct history h;
	struct design r;
	struct design2 d;
	struct recv_model m;
	struct recv_loglik l;
};

struct test test;

#define HISTORY	(&test.h)
#define RECV_DESIGN (&test.r)
#define DYAD_DESIGN (&test.d)
#define RECV_MODEL (&test.m)
#define RECV_LOGLIK (&test.l)



static void fixture_setup_enron()
{
	double width = 6 * 7 * 24 * 60 * 60;
	double intvls[] = { 2 * 60 * 60, 32 * 60 * 60, 3 * 7 * 24 * 60 * 60 };
	size_t nintvl = 3;

	srand(0);

	history_fixture_setup_enron(&ctx.h);

	design_fixture_setup_enron(&ctx.r);
	design_fixture_add_traits_enron(&ctx.r);
	design_fixture_add_isendtot(&ctx.r, width);
	design_fixture_add_irecvtot(&ctx.r, width);
	design_fixture_add_nrecvtot(&ctx.r, intvls, nintvl);
	design_fixture_add_nsendtot(&ctx.r, intvls, nintvl);

	design2_fixture_setup_enron(&ctx.d);
	design2_fixture_add_irecv(&ctx.d, width);
	design2_fixture_add_isend(&ctx.d, width);
	design2_fixture_add_nrecv(&ctx.d, intvls, nintvl);
	design2_fixture_add_nsend(&ctx.d, intvls, nintvl);

	recv_model_fixture_setup(&ctx.m, &ctx.r, &ctx.d);
	recv_model_fixture_set_exlude_loops(&ctx.m, 1);
	recv_model_fixture_set_rand(&ctx.m);
}


static void fixture_teardown()
{
	recv_model_fixture_teardown(&ctx.m);
	design2_fixture_teardown(&ctx.d);
	design_fixture_teardown(&ctx.r);
	history_fixture_teardown(&ctx.h);
}


static void setup()
{
	history_test_setup(&test.h, &ctx.h);
	design_test_setup(&test.r, &test.h, &ctx.r);
	design2_test_setup(&test.d, &test.h, &ctx.d);
	recv_model_test_setup(&test.m, &test.r, &test.d, &ctx.m);
	recv_loglik_init(&test.l, &test.m);
}


static void teardown()
{
	recv_loglik_deinit(&test.l);
	recv_model_test_teardown(&test.m);
	design2_test_teardown(&test.d);
	design_test_teardown(&test.r);
	history_test_teardown(&test.h);
}


/*
static void basic_setup()
{
	double intvls[4] = {
		112.50,  450.00, 1800.00, INFINITY
	};
	int has_loops = 0;
	frame_init(&frame, nsend, nrecv, has_loops, intvls, 4);
	struct design *r = frame_recv_design(&frame);
	design_add_traits(r, trait_names, traits, ntrait);

	struct design2 *d = frame_dyad_design(&frame);
	design2_add_tvar(d, "NRecv", VAR2_NRECV);
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
	double intvls[7] = {
		1800.00,   // 30 min
		7200.00,   //  2 hr
		28800.00,  //  8 hr
		115200.00, // 32 hr
		460800.00, // 5.33 day
		1843200.00, // 21.33 day
		INFINITY
	};

	double beta[] = {
		0.4824403589966597,
		-0.215108203814423,
		-0.4182912193871692,
		-0.04491132567966908,
		0.05698656051879408,
		0.1420793131926562,
		-0.249037515029601,
		-0.02644490380870458,
		0.2987018579760873,
		1.128469426856978,
		1.258299141255399,
		0.3882642769919149,
		0.0984948728445552,
		0.03433916773539594,
		0.01662024071711138,
		-0.007969263858043221,
		0.0006115999551254905,
		3.739814520709765,
		0.9588447963575467,
		0.4521880948047269,
		0.187294928856629,
		0.1757997396356974,
		0.08539282957216279,
		0.05050240510803124,
		0.001111188823465376
	};

	int has_loops = 0;
	frame_init(&frame, nsend, nrecv, has_loops, intvls, 7);
	struct design *r = frame_recv_design(&frame);
	design_add_traits(r, trait_names, traits, ntrait);

	struct design2 *d = frame_dyad_design(&frame);
	design2_add_tvar(d, "IRecv", VAR2_IRECV);
	design2_add_tvar(d, "NRecv", VAR2_NRECV);
	design2_add_tvar(d, "ISend", VAR2_ISEND);
	design2_add_tvar(d, "NSend", VAR2_NSEND);

	recv_coefs_init(&coefs, &frame);
	memcpy(coefs.all, beta, coefs.dim * sizeof(double));

	recv_model_init(&model, &frame, &coefs);
	recv_loglik_init(&loglik, &model);
}
*/


static void test_count()
{
	const struct message *msgs;
	size_t i, n;

	history_get_messages(HISTORY, &msgs, &n);
	const struct message *msg = NULL;

	size_t count0, count1;
	size_t last_count0, last_count1;

	count0 = 0;

	for (i = 0; i < MIN(1000, n); i++) {
		printf("."); fflush(stdout);
		msg = &msgs[i];

		recv_loglik_add(RECV_LOGLIK, msg);
		assert_real_identical(msg->time, history_time(HISTORY));

		last_count0 = msg->nto;
		last_count1 = recv_loglik_last_count(RECV_LOGLIK);
		assert_int_equal(last_count0, last_count1);

		count0 = count0 + last_count0;
		count1 = recv_loglik_count(RECV_LOGLIK);
		assert_int_equal(count0, count1);
	}
}


static void test_dev()
{
	const struct message *msgs;
	size_t i, n;

	history_get_messages(HISTORY, &msgs, &n);
	const struct message *msg = NULL;

	const struct catdist1 *dist = NULL;
	size_t ito, nmsg;
	double last_dev0, last_dev1;
	double dev0, dev1, dev_old;
	
	dev_old = 0.0;
	
	nmsg = 0;

	for (i = 0; i < MIN(1000, n); i++) {
		printf("."); fflush(stdout);
		msg = &msgs[i];
		
		recv_loglik_add(RECV_LOGLIK, msg);
		assert_real_identical(msg->time, history_time(HISTORY));

		nmsg += msg->nto;
						
		last_dev0 = 0.0;
		dist = recv_model_dist(RECV_MODEL, msg->from);
		for (ito = 0; ito < msg->nto; ito++) {
			last_dev0 += -2 * catdist1_lprob(dist, msg->to[ito]);
		}
			
		last_dev1 = recv_loglik_last_dev(RECV_LOGLIK);
		assert_real_identical(last_dev0, last_dev1);
			
		dev0 = dev_old + last_dev0;
		dev1 = recv_loglik_dev(RECV_LOGLIK);
		assert_real_approx(dev0, dev1);

		dev_old = dev1;
	}
}


static void test_mean()
{
	const struct message *msgs;
	size_t i, j, n;

	history_get_messages(HISTORY, &msgs, &n);
	const struct message *msg = NULL;

	const struct catdist1 *dist = NULL;
	double t;
	double probs[NRECV];
	struct recv_params last_mean0, last_mean1;
	struct recv_params mean0, mean1;

	const struct design *r = RECV_DESIGN;
	const struct design2 *d = DYAD_DESIGN;
	recv_params_init(&mean0, r, d);
	recv_params_init(&mean1, r, d);
	recv_params_init(&last_mean0, r, d);
	recv_params_init(&last_mean1, r, d);


	recv_params_set(&mean0, NULL, r, d);

	for (i = 0; i < MIN(1000, n); i++) {
		printf("."); fflush(stdout);
		msg = &msgs[i];
		t = msg->time;

		history_advance(HISTORY, t);

		recv_loglik_add(RECV_LOGLIK, msg);

		dist = recv_model_dist(RECV_MODEL, msg->from);
		for (j = 0; j < NRECV; j++) {
			probs[j] = catdist1_prob(dist, j);
		}
		design_tmul(msg->nto, r, probs, 0.0, &last_mean0.recv);
		design2_tmul(msg->nto, d, msg->from, probs, 0.0, &last_mean0.dyad);

		recv_params_set(&last_mean1, NULL, r, d);
		recv_loglik_axpy_last_mean(1.0, RECV_LOGLIK, &last_mean1);

		assert_vec_approx(last_mean0.recv.traits, last_mean1.recv.traits, design_trait_dim(r));
		assert_vec_approx(last_mean0.recv.tvars, last_mean1.recv.tvars, design_tvar_dim(r));
		assert_vec_approx(last_mean0.dyad.traits, last_mean1.dyad.traits, design2_trait_dim(d));
		assert_vec_approx(last_mean0.dyad.tvars, last_mean1.dyad.tvars, design2_tvar_dim(d));

		recv_params_axpy(1.0, &last_mean0, &mean0, r, d);

		recv_params_set(&mean1, NULL, r, d);
		recv_loglik_axpy_mean(1.0, RECV_LOGLIK, &mean1);

		assert_vec_approx(mean0.recv.traits, mean1.recv.traits, design_trait_dim(r));
		assert_vec_approx(mean0.recv.tvars, mean1.recv.tvars, design_tvar_dim(r));
		assert_vec_approx(mean0.dyad.traits, mean1.dyad.traits, design2_trait_dim(d));
		assert_vec_approx(mean0.dyad.tvars, mean1.dyad.tvars, design2_tvar_dim(d));

		recv_params_set(&mean1, &mean0, r, d);
	}

	recv_params_deinit(&last_mean1);
	recv_params_deinit(&last_mean0);
	recv_params_deinit(&mean1);
	recv_params_deinit(&mean0);
}


static void test_score()
{
	const struct message *msgs;
	size_t i, ito, n;

	history_get_messages(HISTORY, &msgs, &n);
	const struct message *msg = NULL;

	double t;
	struct recv_params last_score0, last_score1;
	struct recv_params score0, score1;

	const struct design *r = RECV_DESIGN;
	const struct design2 *d = DYAD_DESIGN;
	recv_params_init(&score0, r, d);
	recv_params_init(&score1, r, d);
	recv_params_init(&last_score0, r, d);
	recv_params_init(&last_score1, r, d);

	recv_params_set(&score0, NULL, r, d);

	for (i = 0; i < MIN(1000, n); i++) {
		printf("."); fflush(stdout);
		msg = &msgs[i];
		t = msg->time;

		recv_loglik_add(RECV_LOGLIK, msg);

		recv_params_set(&last_score0, NULL, r, d);
		for (ito = 0; ito < msg->nto; ito++) {
			design_axpy(1.0, r, msg->to[ito], &last_score0.recv);
			design2_axpy(1.0, d, msg->from, msg->to[ito], &last_score0.dyad);
		}
		recv_loglik_axpy_last_mean(-1.0, RECV_LOGLIK, &last_score0);

		recv_params_set(&last_score1, NULL, r, d);
		recv_loglik_axpy_last_score(1.0, RECV_LOGLIK, &last_score1);

		assert_vec_approx(last_score0.recv.traits, last_score1.recv.traits, design_trait_dim(r));
		assert_vec_approx(last_score0.recv.tvars, last_score1.recv.tvars, design_tvar_dim(r));
		assert_vec_approx(last_score0.dyad.traits, last_score1.dyad.traits, design2_trait_dim(d));
		assert_vec_approx(last_score0.dyad.tvars, last_score1.dyad.tvars, design2_tvar_dim(d));

		recv_params_axpy(1.0, &last_score0, &score0, r, d);

		recv_params_set(&score1, NULL, r, d);
		recv_loglik_axpy_score(1.0, RECV_LOGLIK, &score1);

		assert_vec_approx(score0.recv.traits, score1.recv.traits, design_trait_dim(r));
		assert_vec_approx(score0.recv.tvars, score1.recv.tvars, design_tvar_dim(r));
		assert_vec_approx(score0.dyad.traits, score1.dyad.traits, design2_trait_dim(d));
		assert_vec_approx(score0.dyad.tvars, score1.dyad.tvars, design2_tvar_dim(d));

		recv_params_set(&score0, &score1, r, d);
	}

	recv_params_deinit(&last_score1);
	recv_params_deinit(&last_score0);
	recv_params_deinit(&score1);
	recv_params_deinit(&score0);
}


static void test_imat()
{
	const struct message *msgs;
	size_t i, j, n;

	history_get_messages(HISTORY, &msgs, &n);
	const struct message *msg = NULL;

	size_t dim = recv_model_dim(RECV_MODEL);
	size_t cov_dim = dim * (dim + 1) / 2;
	const struct design *r = RECV_DESIGN;
	const struct design2 *d = DYAD_DESIGN;
	struct recv_params mean, diff;

	double *diff_vec = xmalloc(dim * sizeof(double));
	double *cov0 = xmalloc(cov_dim * sizeof(double));
	double *cov1 = xmalloc(cov_dim * sizeof(double));
	double *scaled_cov0 = xmalloc(cov_dim * sizeof(double));
	double *scaled_cov1 = xmalloc(cov_dim * sizeof(double));
	double *last_cov0 = xmalloc(cov_dim * sizeof(double));
	double *last_cov1 = xmalloc(cov_dim * sizeof(double));

	recv_params_init(&mean, r, d);
	recv_params_init(&diff, r, d);

	const struct catdist1 *dist = NULL;
	const enum blas_uplo uplo = MLOGIT_COV_UPLO;
	const enum blas_uplo fuplo = uplo == BLAS_UPPER ? BLAS_LOWER : BLAS_UPPER;

	memset(cov0, 0, cov_dim * sizeof(*cov0));

	for (i = 0; i < MIN(1000, n); i++) {
		printf("."); fflush(stdout);
		msg = &msgs[i];

		recv_loglik_add(RECV_LOGLIK, msg);

		/* compute mean */
		recv_params_set(&mean, NULL, r, d);
		recv_loglik_axpy_last_mean(1.0 / (msg->nto), RECV_LOGLIK, &mean);

		/* compute last_cov0 */
		dist = recv_model_dist(RECV_MODEL, msg->from);
		memset(last_cov0, 0, cov_dim * sizeof(double));
		for (j = 0; j < NRECV; j++) {
			recv_params_set(&diff, &mean, r, d);
			design_axpy(-1.0, r, j, &diff.recv);
			design2_axpy(-1.0, d, msg->from, j, &diff.dyad);

			double w = catdist1_prob(dist, j);
			size_t off = 0;
			memcpy(diff_vec + off, diff.recv.traits, design_trait_dim(r) * sizeof(double)); off += design_trait_dim(r);
			memcpy(diff_vec + off, diff.recv.tvars, design_tvar_dim(r) * sizeof(double)); off += design_tvar_dim(r);
			memcpy(diff_vec + off, diff.dyad.traits, design2_trait_dim(d) * sizeof(double)); off += design2_trait_dim(d);
			memcpy(diff_vec + off, diff.dyad.tvars, design2_tvar_dim(d) * sizeof(double)); off += design2_tvar_dim(d);

			blas_dspr(fuplo, dim, w, diff_vec, 1, last_cov0);
		}
		blas_dscal(cov_dim, msg->nto, last_cov0, 1);

		/* compute last_cov1 */
		memset(last_cov1, 0, cov_dim * sizeof(double));
		recv_loglik_axpy_last_imat(1.0, RECV_LOGLIK, last_cov1);

		assert_sym_approx(last_cov0, last_cov1, uplo, dim);

		/* compute cov0 */
		blas_daxpy(cov_dim, 1.0, last_cov0, 1, cov0, 1);

		/* compute cov1 */
		memset(cov1, 0, cov_dim * sizeof(double));
		recv_loglik_axpy_imat(1.0, RECV_LOGLIK, cov1);

		double n = (double)recv_loglik_count(RECV_LOGLIK);
		memcpy(scaled_cov0, cov0, cov_dim * sizeof(double));
		memcpy(scaled_cov1, cov1, cov_dim * sizeof(double));
		blas_dscal(cov_dim, 1.0/n, scaled_cov0, 1);
		blas_dscal(cov_dim, 1.0/n, scaled_cov1, 1);
		assert_sym_approx(scaled_cov0, scaled_cov1, BLAS_UPPER, dim);

		blas_dcopy(dim, cov1, 1, cov0, 1);
	}

	recv_params_deinit(&diff);
	recv_params_deinit(&mean);
	free(last_cov1);
	free(last_cov0);
	free(scaled_cov1);
	free(scaled_cov0);
	free(cov1);
	free(cov0);
	free(diff_vec);
}



int main()
{
	UnitTest tests[] = {
		unit_test_setup(enron_suite, fixture_setup_enron),
		unit_test_setup_teardown(test_count, setup, teardown),
		unit_test_setup_teardown(test_dev, setup, teardown),
		unit_test_setup_teardown(test_mean, setup, teardown),
		unit_test_setup_teardown(test_score, setup, teardown),
		unit_test_setup_teardown(test_imat, setup, teardown),
		unit_test_teardown(enron_suite, fixture_teardown),
		
	};
	return run_tests(tests);
}
