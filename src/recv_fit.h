#ifndef RECV_FIT_H
#define RECV_FIT_H

#include "constr.h"
#include "design.h"
#include "newton.h"

#include "recv_loglik.h"
#include "recv_model.h"





struct recv_fit {
	struct newton opt;

	const struct message *msgs;
	size_t nmsg;

	/* working frame + model */
	struct frame *frame;
	struct recv_model model;

	/* null deviance */
	double dev0;

};


/*
size_t constr_add_identify_recv_fit(struct constr *c, struct frame *f,
				    const struct messages *xmsgs,
				    const struct messages *ymsgs);
*/

void recv_fit_init(struct recv_fit *fit, struct design *r, struct design2 *d,
		   const struct message *msgs, size_t nmsg,
		   const struct constr *c, const struct newton_ctrl *ctrl);
void recv_fit_deinit(struct recv_fit *fit);



/* fitting */
enum newton_task recv_fit_start(struct recv_fit *fit,
				const struct recv_params *p,
				const double *duals);
enum newton_task recv_fit_advance(struct recv_fit *fit);


/* current values */
const struct recv_loglik *recv_fit_loglik(const struct recv_fit *fit);
double recv_fit_dev(const struct recv_fit *fit);
double recv_fit_dev0(const struct recv_fit *fit);
const double *recv_fit_duals(const struct recv_fit *fit);


#endif /* RECV_FIT_H */
