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
#include "recv_loglik.h"
#include "fit.h"


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


int main(int argc, char **argv)
{
	setup();
	
	double penalty = 1.0;
	ssize_t maxit = 5000;
	ssize_t report = 1;
	bool trace = true;
	
	struct recv_fit fit;
	recv_fit_init(&fit, &messages, &design, NULL, penalty);
	ssize_t it = 0;
	
	do {
		it++;
		recv_fit_step(&fit);
		
		if (trace && it % report == 0) {
			const char *msg = penalty == 0 ? "" : "(penalized) ";
			ssize_t n = recv_loglik_count(&fit.loglik);
			double dev = n * recv_loglik_avg_dev(&fit.loglik);
			double dec = -2 * bfgs_decrement(&fit.opt);
			printf("iter %"SSIZE_FMT" deviance %s%.6f decrement %.6f\n", it,
				msg, dev, dec);
		}
	} while (it < maxit && !recv_fit_converged(&fit));
	
	teardown();
}
