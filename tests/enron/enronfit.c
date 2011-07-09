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
#include "model.h"
#include "recv_loglik.h"
#include "recv_fit.h"


static struct actors senders;
static struct actors receivers;
static struct vector intervals;
static struct messages messages;
static struct design design;
static struct frame frame;
static struct vector coefs;
static struct model model;
static struct recv_loglik recv_loglik;


static void setup(void) {
	enron_employees_init(&senders);
	enron_employees_init(&receivers);
	enron_messages_init(&messages);
	
	ssize_t i;
	double intvls[] = {
		// 56.25,      
		// 112.50,      
		225.00,      450.00,      900.00,
		1800.00,     3600.00,     7200.00,    14400.00,    28800.00,
		57600.00,    115200.00,   230400.00,   460800.00,   921600.00,
		1843200.00,  3686400.00,  7372800.00, 14745600.00, 29491200.00,
		58982400.00
	};
	ssize_t nintvls = sizeof(intvls) / sizeof(intvls[0]);
	struct vector vintvls = vector_make(intvls, nintvls);
	bool has_reffects = false;
	bool has_loops = false;
	vector_init_copy(&intervals, &vintvls);
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
	recv_loglik_init(&recv_loglik, &model);
}

static void teardown(void)
{
	recv_loglik_deinit(&recv_loglik);
	model_deinit(&model);
	vector_deinit(&coefs);
	frame_deinit(&frame);
	vector_deinit(&intervals);	
	design_deinit(&design);
	messages_deinit(&messages);
	actors_deinit(&receivers);	
	actors_deinit(&senders);	
}


#define YG(gen) \
	do { \
		if ((err = gen) != yajl_gen_status_ok) \
			return err; \
	} while (0)
#define YSTR(str) ((const unsigned char *)(str))

#define COEFFICIENTS		"coefficients"
#define RANK			"rank"
#define CONSTRAINTS		"constraints"
#define CONSTRAINT_VALUES	"constraint_values"
#define INFORMATION		"information"
#define DEVIANCE		"deviance"
#define NULL_DEVIANCE		"null_deviance"
#define DF_RESIDUAL		"df_residual"
#define DF_NULL			"df_null"

yajl_gen_status yaj_gen_recv_fit(yajl_gen hand, const struct recv_fit *fit)
{
	assert(fit);
	yajl_gen_status err = yajl_gen_status_ok;
	const struct recv_loglik *ll = recv_fit_loglik(fit);
	const struct recv_loglik_info *info = recv_loglik_info(ll);
	
	YG(yajl_gen_map_open(hand));
	
	/* coefficients */	
	YG(yajl_gen_string(hand, YSTR(COEFFICIENTS), strlen(COEFFICIENTS)));
	YG(yajl_gen_vector(hand, recv_fit_coefs(fit)));
	
	/* residuals */
	/* fitted.values */
        /* effects */
	
	/* rank */	
	YG(yajl_gen_string(hand, YSTR(RANK), strlen(RANK)));
	YG(yajl_gen_integer(hand, recv_fit_rank(fit)));
	
	/* qr */
	YG(yajl_gen_string(hand, YSTR(CONSTRAINTS), strlen(CONSTRAINTS)));
	YG(yajl_gen_matrix(hand, recv_fit_ce(fit)));

	YG(yajl_gen_string(hand, YSTR(CONSTRAINT_VALUES), strlen(CONSTRAINT_VALUES)));
	YG(yajl_gen_vector(hand, recv_fit_be(fit)));

	YG(yajl_gen_string(hand, YSTR(INFORMATION), strlen(INFORMATION)));
	YG(yajl_gen_matrix(hand, &info->imat));

	/* linear.predictors (eta) */
	
	/* deviance */
	YG(yajl_gen_string(hand, YSTR(DEVIANCE), strlen(DEVIANCE)));
	YG(yajl_gen_ieee754(hand, info->nrecv * info->dev));
	
	/* aic */
	
	/* null.deviance */
	YG(yajl_gen_string(hand, YSTR(NULL_DEVIANCE), strlen(NULL_DEVIANCE)));
	YG(yajl_gen_ieee754(hand, info->nrecv * recv_fit_dev0(fit)));
	
	/* iter */
	/* weights = wt */
	/* prior.weights = weights */
	
	/* df.residual */
	YG(yajl_gen_string(hand, YSTR(DF_RESIDUAL), strlen(DF_RESIDUAL)));
	YG(yajl_gen_integer(hand, info->nrecv - recv_fit_rank(fit)));

	/* df.null */
	YG(yajl_gen_string(hand, YSTR(DF_NULL), strlen(DF_NULL)));
	YG(yajl_gen_integer(hand, info->nrecv));

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
	
	yajl_gen_map_close(hand);
	
	
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


int main(int argc, char **argv)
{
	setup();
	
	/* ssize_t n = messages_recv_count(&messages); */
	//double penalty = 0; // 0.00001; // n / 512.0; // >= 0.00001 works
	ssize_t maxit = 25;
	ssize_t report = 1;
	bool trace = true;
	
	struct recv_fit fit;
	recv_fit_init(&fit, &messages, &design, NULL, NULL);
	enum recv_fit_task task;
	ssize_t it = 0;
	
	do {
		it++;
		task = recv_fit_advance(&fit);

		if (trace && it % report == 0 && task != RECV_FIT_CONV) {
			const struct recv_loglik *ll = recv_fit_loglik(&fit);
			const struct recv_loglik_info *info = recv_loglik_info(ll);
			ssize_t n = info->nrecv;
			double dev = n * info->dev;
			const struct vector *score = &info->score;
			double ngrad = vector_max_abs(score);
			double step = recv_fit_step(&fit);
			fprintf(stderr, "iter %"SSIZE_FMT" deviance = %.2f; |grad| = %.16f; step = %.16f\n",
				it, dev, ngrad, step);
		}
	} while (it < maxit && task == RECV_FIT_STEP);
	
	if (task != RECV_FIT_CONV) {
		printf("ERROR: %s\n", recv_fit_errmsg(&fit));
		return -1;
	} else {
		output(&fit);
	}
	
	teardown();
	return 0;
}
