#include "port.h"
#include <assert.h>
#include <stdio.h>

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
		if (trace && it % report == 0) {
			ssize_t n = fit.loglik.info.nrecv;
			double dev = n * recv_loglik_avg_dev(&fit.loglik);
			const struct vector *score = &fit.loglik.info.score;
			double ngrad = vector_max_abs(score);
			printf("iter %"SSIZE_FMT" deviance %.1f |grad| %.22f\n",
			       it, dev, ngrad);
		}
		
		it++;
		task = recv_fit_advance(&fit);
	} while (it < maxit && task == RECV_FIT_STEP);
	
	if (task != RECV_FIT_CONV) {
		printf("ERROR: %s\n", recv_fit_errmsg(&fit));
		return -1;
	}
	
	teardown();
	return 0;
}
