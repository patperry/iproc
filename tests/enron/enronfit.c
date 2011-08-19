#include "port.h"
#include <assert.h>
#include <stdio.h>
#include <string.h>

#include "ieee754.h"
#include "json.h"
#include "enron.h"
#include "messages.h"
#include "design.h"
#include "vars.h"
#include "frame.h"
#include "recv_model.h"
#include "recv_loglik.h"
#include "recv_fit.h"


static struct actors enron_actors;
static struct matrix enron_traits;

static struct actors senders;
static struct actors receivers;
static struct matrix recv_traits;
static struct vector intervals;
static struct messages messages;
static struct design design;


static void setup(void) {
	enron_employees_init(&enron_actors, &enron_traits);
	enron_messages_init(&messages, 10);	
	
	actors_init_copy(&senders, &enron_actors);
	actors_init_copy(&receivers, &enron_actors);	
	matrix_init_copy(&recv_traits, TRANS_NOTRANS, &enron_traits);
	
	double intvls[] = {
		// 450.00,
		//900.00,    // 15 min
		1800.00,   // 30 min
			   //3600.00,   //  1 hr
		7200.00,   //  2 hr
			   //14400.00,  //  4 hr
		28800.00,  //  8 hr
			   //57600.00,  // 16 hr
		115200.00, // 32 hr
			   // 230400.00, // 2.66 day 
		460800.00, // 5.33 day
			   // 921600.00, // 10.66 day
		1843200.00, // 21.33 day
			    // 3686400.00, // 42.66 day
			    // 7372800.00, // 85.33 day
			    // 14745600.00, // 170.66 day
			    // 29491200.00 // 341.33 day
		// 58982400.00
	};
	ssize_t nintvls = sizeof(intvls) / sizeof(intvls[0]);
	struct vector vintvls = vector_make(intvls, nintvls);
	bool has_reffects = false;
	bool has_loops = false;
	vector_init_copy(&intervals, &vintvls);
	design_init(&design, &senders, &receivers, &recv_traits, &intervals);
	design_set_loops(&design, has_loops);
	design_set_recv_effects(&design, has_reffects);
	design_add_recv_var(&design, RECV_VAR_IRECV, NULL);
	design_add_recv_var(&design, RECV_VAR_NRECV, NULL);
	design_add_recv_var(&design, RECV_VAR_ISEND, NULL);
	design_add_recv_var(&design, RECV_VAR_NSEND, NULL);
	design_add_recv_var(&design, RECV_VAR_IRECV2, NULL);
	design_add_recv_var(&design, RECV_VAR_NRECV2, NULL);
	design_add_recv_var(&design, RECV_VAR_ISEND2, NULL);
	design_add_recv_var(&design, RECV_VAR_NSEND2, NULL);
}

static void teardown(void)
{
	vector_deinit(&intervals);	
	design_deinit(&design);
	matrix_deinit(&recv_traits);
	actors_deinit(&receivers);	
	actors_deinit(&senders);	
	messages_deinit(&messages);
	matrix_deinit(&enron_traits);
	actors_deinit(&enron_actors);
}


#define YG(gen) \
	do { \
		if ((err = gen) != yajl_gen_status_ok) \
			return err; \
	} while (0)
#define YSTR(str) ((const unsigned char *)(str))

#define COEFFICIENTS		"coefficients"
#define COUNT			"count"
#define SCORE			"score"
#define INFORMATION		"information"
#define RANK			"rank"
#define CONSTRAINTS		"constraints"
#define CONSTRAINT_VALUES	"constraint_values"
#define DUALS			"duals"
#define DEVIANCE		"deviance"
#define NULL_DEVIANCE		"null_deviance"
#define DF_RESIDUAL		"df_residual"
#define DF_NULL			"df_null"

yajl_gen_status yaj_gen_recv_fit(yajl_gen hand, const struct recv_fit *fit)
{
	assert(fit);
	yajl_gen_status err = yajl_gen_status_ok;
	const struct recv_loglik *ll = recv_fit_loglik(fit);

	ssize_t ne = recv_fit_constr_count(fit);
	ssize_t dim = recv_model_dim(ll->model);
	ssize_t ic, nc = recv_model_cohort_count(ll->model);

	YG(yajl_gen_map_open(hand));
	{
		/* coefficients */
		YG(yajl_gen_string(hand, YSTR(COEFFICIENTS), strlen(COEFFICIENTS)));
		const struct matrix *coefs = recv_fit_coefs(fit);		
		YG(yajl_gen_matrix(hand, coefs));

		YG(yajl_gen_string(hand, YSTR(COUNT), strlen(COUNT)));
		YG(yajl_gen_array_open(hand));
		for (ic = 0; ic < nc; ic++) {
			ssize_t count = recv_loglik_count(ll, ic);
			YG(yajl_gen_integer(hand, count));
		}
		YG(yajl_gen_array_close(hand));
		
		YG(yajl_gen_string(hand, YSTR(SCORE), strlen(SCORE)));
		YG(yajl_gen_array_open(hand));
		for (ic = 0; ic < nc; ic++) {
			const struct recv_loglik_info *info = recv_loglik_info(ll, ic);
			const struct vector *score = &info->score;
			YG(yajl_gen_vector(hand, score));
		}
		YG(yajl_gen_array_close(hand));
		
		YG(yajl_gen_string(hand, YSTR(INFORMATION), strlen(INFORMATION)));
		YG(yajl_gen_array_open(hand));
		for (ic = 0; ic < nc; ic++) {
			const struct recv_loglik_info *info = recv_loglik_info(ll, ic);
			const struct matrix *imat = &info->imat;
			YG(yajl_gen_matrix(hand, imat));
		}
		YG(yajl_gen_array_close(hand));

		YG(yajl_gen_string(hand, YSTR(CONSTRAINTS), strlen(CONSTRAINTS)));	
		const struct matrix *ce;
		const struct vector *be;
		recv_fit_get_constr(fit, &ce, &be);
		YG(yajl_gen_matrix(hand, ce));
		YG(yajl_gen_string(hand, YSTR(CONSTRAINT_VALUES), strlen(CONSTRAINT_VALUES)));
		YG(yajl_gen_vector(hand, be));
		
		YG(yajl_gen_string(hand, YSTR(DUALS), strlen(DUALS)));
		const struct vector *duals = recv_fit_duals(fit);
		YG(yajl_gen_vector(hand, duals));
		
		/* rank */
		YG(yajl_gen_string(hand, YSTR(RANK), strlen(RANK)));
		ssize_t rank = dim * nc - ne;
		YG(yajl_gen_integer(hand, rank));
		
		/* deviance */
		YG(yajl_gen_string(hand, YSTR(DEVIANCE), strlen(DEVIANCE)));
		double dev = recv_fit_dev(fit);
		YG(yajl_gen_ieee754(hand, dev));
		
		/* df.residual */
		YG(yajl_gen_string(hand, YSTR(DF_RESIDUAL), strlen(DF_RESIDUAL)));
		ssize_t ntot = recv_loglik_count_sum(ll);
		YG(yajl_gen_integer(hand, ntot - rank));
		
		/* df.null */
		YG(yajl_gen_string(hand, YSTR(DF_NULL), strlen(DF_NULL)));
		YG(yajl_gen_integer(hand, ntot));
	}
	YG(yajl_gen_map_close(hand));
		
	/* residuals */
	/* fitted.values */
	/* effects */

	/* qr */	
	/* linear.predictors (eta) */
	/* aic */
		
	/* null.deviance */
	//YG(yajl_gen_string(hand, YSTR(NULL_DEVIANCE), strlen(NULL_DEVIANCE)));
	//YG(yajl_gen_ieee754(hand, info->nrecv * recv_fit_dev0(fit, c)));
		
	/* iter */
	/* weights = wt */
	/* prior.weights = weights */
		
	/* y = y */
	/* converged = conv */
	/* boundary = boundary */
	/* call */
	/* formula */
	/* terms */
	/* data */
	/* offset */
	/* control */
	/* method */
	/* contrasts */
	/* xlevels */
	
	
	return err;
}

static void output(const struct recv_fit *fit)
{
	/* generator config */
	yajl_gen g = yajl_gen_alloc(NULL);
	yajl_gen_config(g, yajl_gen_beautify, 1);
	yaj_gen_recv_fit(g, fit);
	
	
	/* output */
	const unsigned char * buf;
	size_t len;
	yajl_gen_get_buf(g, &buf, &len);
	fwrite(buf, 1, len, stdout);
	yajl_gen_clear(g);
	yajl_gen_free(g);
}


int main()
{
	setup();
	
	ssize_t maxit = 100;
	ssize_t report = 1;
	bool trace = true;
	
	struct recv_fit fit;
	recv_fit_init(&fit, &messages, &design, &senders, NULL);
	ssize_t dim0 = matrix_ncol(&recv_traits);
	ssize_t nc = actors_cohort_count(&senders);	
	ssize_t dim = design_recv_dim(&design);
	struct vector ce;
	vector_init(&ce, dim * nc);
	ssize_t i;
	
	/* constrain dynamic covariates */
	for (i = dim0; i < dim; i++) {
		// 3-way interactions
		vector_fill(&ce, 0.0);
		vector_set_item(&ce, i + dim * ENRON_COHORT_LJF, +1.0);
		vector_set_item(&ce, i + dim * ENRON_COHORT_LJM, -1.0);
		vector_set_item(&ce, i + dim * ENRON_COHORT_LSF, -1.0);
		vector_set_item(&ce, i + dim * ENRON_COHORT_LSM, +1.0);		
		vector_set_item(&ce, i + dim * ENRON_COHORT_OJF, -1.0);
		vector_set_item(&ce, i + dim * ENRON_COHORT_OJM, +1.0);
		vector_set_item(&ce, i + dim * ENRON_COHORT_OSF, +1.0);
		vector_set_item(&ce, i + dim * ENRON_COHORT_OSM, -1.0);
		recv_fit_add_constr(&fit, &ce, 0.0);
		
		vector_fill(&ce, 0.0);
		vector_set_item(&ce, i + dim * ENRON_COHORT_TJF, +1.0);
		vector_set_item(&ce, i + dim * ENRON_COHORT_TJM, -1.0);
		vector_set_item(&ce, i + dim * ENRON_COHORT_TSF, -1.0);
		vector_set_item(&ce, i + dim * ENRON_COHORT_TSM, +1.0);		
		vector_set_item(&ce, i + dim * ENRON_COHORT_OJF, -1.0);
		vector_set_item(&ce, i + dim * ENRON_COHORT_OJM, +1.0);
		vector_set_item(&ce, i + dim * ENRON_COHORT_OSF, +1.0);
		vector_set_item(&ce, i + dim * ENRON_COHORT_OSM, -1.0);
		recv_fit_add_constr(&fit, &ce, 0.0);
				
		// 2-way interactions
		vector_fill(&ce, 0.0);
		vector_set_item(&ce, i + dim * ENRON_COHORT_LJM, +1.0);
		vector_set_item(&ce, i + dim * ENRON_COHORT_LSM, -1.0);
		vector_set_item(&ce, i + dim * ENRON_COHORT_OJM, -1.0);		
		vector_set_item(&ce, i + dim * ENRON_COHORT_OSM, +1.0);
		recv_fit_add_constr(&fit, &ce, 0.0);

		vector_fill(&ce, 0.0);
		vector_set_item(&ce, i + dim * ENRON_COHORT_TJM, +1.0);
		vector_set_item(&ce, i + dim * ENRON_COHORT_TSM, -1.0);
		vector_set_item(&ce, i + dim * ENRON_COHORT_OJM, -1.0);
		vector_set_item(&ce, i + dim * ENRON_COHORT_OSM, +1.0);
		recv_fit_add_constr(&fit, &ce, 0.0);

		vector_fill(&ce, 0.0);
		vector_set_item(&ce, i + dim * ENRON_COHORT_LSF, +1.0);
		vector_set_item(&ce, i + dim * ENRON_COHORT_LSM, -1.0);
		vector_set_item(&ce, i + dim * ENRON_COHORT_OSF, -1.0);
		vector_set_item(&ce, i + dim * ENRON_COHORT_OSM, +1.0);
		recv_fit_add_constr(&fit, &ce, 0.0);

		vector_fill(&ce, 0.0);
		vector_set_item(&ce, i + dim * ENRON_COHORT_TSF, +1.0);
		vector_set_item(&ce, i + dim * ENRON_COHORT_TSM, -1.0);
		vector_set_item(&ce, i + dim * ENRON_COHORT_OSF, -1.0);
		vector_set_item(&ce, i + dim * ENRON_COHORT_OSM, +1.0);
		recv_fit_add_constr(&fit, &ce, 0.0);

		vector_fill(&ce, 0.0);
		vector_set_item(&ce, i + dim * ENRON_COHORT_OJF, +1.0);
		vector_set_item(&ce, i + dim * ENRON_COHORT_OJM, -1.0);
		vector_set_item(&ce, i + dim * ENRON_COHORT_OSF, -1.0);
		vector_set_item(&ce, i + dim * ENRON_COHORT_OSM, +1.0);
		recv_fit_add_constr(&fit, &ce, 0.0);

		// main effects
		vector_fill(&ce, 0.0);
		vector_set_item(&ce, i + dim * ENRON_COHORT_LSM, +1.0);
		vector_set_item(&ce, i + dim * ENRON_COHORT_OSM, -1.0);
		recv_fit_add_constr(&fit, &ce, 0.0);
		
		vector_fill(&ce, 0.0);
		vector_set_item(&ce, i + dim * ENRON_COHORT_TSM, +1.0);
		vector_set_item(&ce, i + dim * ENRON_COHORT_OSM, -1.0);
		recv_fit_add_constr(&fit, &ce, 0.0);
		
		vector_fill(&ce, 0.0);
		vector_set_item(&ce, i + dim * ENRON_COHORT_OJM, +1.0);
		vector_set_item(&ce, i + dim * ENRON_COHORT_OSM, -1.0);
		recv_fit_add_constr(&fit, &ce, 0.0);
		
		vector_fill(&ce, 0.0);
		vector_set_item(&ce, i + dim * ENRON_COHORT_OSF, +1.0);
		vector_set_item(&ce, i + dim * ENRON_COHORT_OSM, -1.0);
		recv_fit_add_constr(&fit, &ce, 0.0);
		
		// intercept
		//vector_fill(&ce, 0.0);
		//vector_set_item(&ce, i + dim * ENRON_COHORT_OSM, +1.0);		
		//recv_fit_add_constr(&fit, &ce, 0.0);
		
	}
	vector_deinit(&ce);
	
	
	/* add constraints to make the model identifiable */
	ssize_t nadd = recv_fit_add_constr_identify(&fit);
	if (nadd > 0)
		fprintf(stderr, "Adding %"SSIZE_FMT" constraints to make parameters identifiable\n", nadd);
	
	enum recv_fit_task task;
	ssize_t it = 0;
	
	//struct matrix coefs0;
	//matrix_init(&coefs0, dim, nc);
	//matrix_fill(&coefs0, 0.0001);
	
	for (it = 0, task = recv_fit_start(&fit, NULL);
	     it < maxit && task == RECV_FIT_STEP;
	     it++, task = recv_fit_advance(&fit)) {
		if (trace && it % report == 0 && task != RECV_FIT_CONV) {
			// const struct recv_loglik *ll = recv_fit_loglik(&fit);
			// const struct recv_loglik_info *info = recv_loglik_info(ll);
			// ssize_t n = info->nrecv;
			// double dev = n * info->dev;
			// const struct vector *score = &info->score;
			// double ngrad = vector_max_abs(score);
			double step = recv_fit_step(&fit);
			double ngrad = recv_fit_grad_norm2(&fit);
			fprintf(stderr, "iter %"SSIZE_FMT"; |grad| = %.16f; step = %.16f\n",
				it, ngrad, step);

			// fprintf(stderr, "iter %"SSIZE_FMT" deviance = %.2f; |grad| = %.16f; step = %.16f\n",
			//	it, dev, ngrad, step);
		}
	}
	
	//matrix_deinit(&coefs0);
	
	int err = 0;
	
	if (task != RECV_FIT_CONV) {
		fprintf(stderr, "ERROR: %s\n", recv_fit_errmsg(&fit));
		err = -1;
	}

	output(&fit);
	recv_fit_deinit(&fit);
	teardown();
	return err;
}
