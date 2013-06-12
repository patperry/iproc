#ifndef RECV_FIT_H
#define RECV_FIT_H

#include "constr.h"
#include "design.h"
#include "linesearch.h"
#include "messages.h"

#include "recv_loglik.h"
#include "recv_model.h"



struct recv_fit_params {
	double *all;
	size_t dim;
	struct recv_params params;
	double *duals;
	int owner;
};

struct recv_fit_resid {
	struct recv_fit_params params;
	double norm2;
};

struct recv_fit_eval {
	struct recv_loglik loglik;
	struct recv_fit_params params;
	struct recv_fit_resid resid;
};


struct recv_fit_search {
	struct recv_fit_params params;
};

struct recv_fit_rgrad {
	struct recv_fit_params params;
};

struct recv_fit {
	const struct messages *xmsgs, *ymsgs;

	/* working frame + model */
	struct frame *frame;
	struct recv_model model;

	/* null deviance */
	double dev0;

};

void recv_fit_params_init(struct recv_fit_params *params, const struct frame *f,
		      const struct constr *c);
void recv_fit_params_init_view(struct recv_fit_params *params, const struct frame *f,
			   const struct constr *c, const double *data);
void recv_fit_params_deinit(struct recv_fit_params *params);

size_t constr_add_identify_recv_fit(struct constr *c, struct frame *f,
				    const struct messages *xmsgs,
				    const struct messages *ymsgs);

void recv_fit_init(struct recv_fit *fit,
		   struct frame *f,
		   struct constr *c,
		   const struct messages *xmsgs,
		   const struct messages *ymsgs,
		   const struct recv_fit_ctrl *ctrl);
void recv_fit_deinit(struct recv_fit *fit);



/* fitting */
enum recv_fit_task recv_fit_start(struct recv_fit *fit, const struct recv_fit_params *params0);
enum recv_fit_task recv_fit_advance(struct recv_fit *fit);
const char *recv_fit_errmsg(const struct recv_fit *fit);

/* current values */
const struct recv_loglik *recv_fit_loglik(const struct recv_fit *fit);
double recv_fit_dev(const struct recv_fit *fit);
double recv_fit_dev0(const struct recv_fit *fit);
const struct recv_coefs *recv_fit_coefs(const struct recv_fit *fit);
const double *recv_fit_duals(const struct recv_fit *fit);
double recv_fit_step(const struct recv_fit *fit);
double recv_fit_grad_norm2(const struct recv_fit *fit);


#endif /* RECV_FIT_H */
