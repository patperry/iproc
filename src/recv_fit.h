#ifndef RECV_FIT_H
#define RECV_FIT_H

#include "constr.h"
#include "design.h"
#include "newton.h"

#include "recv_loglik.h"
#include "recv_model.h"

enum recv_fit_task {
	RECV_FIT_CONV = 0,
	RECV_FIT_STEP = 1,
	RECV_FIT_ERR_DOM = -1,
	RECV_FIT_ERR_IMAT = -2,
	RECV_FIT_ERR_LNSRCH = -3,
	RECV_FIT_ERR_XTOL = -4,
};


struct recv_fit {
	struct design *recv;
	struct design2 *dyad;
	const struct message *msgs;
	size_t nmsg;

	struct recv_model model;
	struct recv_loglik loglik;
	double dev;
	double *nscore; /* negative score */
	double *imat;
	enum blas_uplo uplo;

	struct constr constr;
	size_t ncextra;

	struct newton opt;
	enum newton_task task;
};


void recv_fit_init(struct recv_fit *fit, struct design *r, struct design2 *d,
		   int exclude_loops, const struct message *msgs, size_t nmsg,
		   const struct constr *c, const struct newton_ctrl *ctrl);
void recv_fit_deinit(struct recv_fit *fit);


static inline size_t recv_fit_extra_constr_count(const struct recv_fit *fit) { return fit->ncextra; }


/* fitting */
enum recv_fit_task recv_fit_start(struct recv_fit *fit,
				  const struct recv_params *p,
				  const double *duals);
enum recv_fit_task recv_fit_advance(struct recv_fit *fit);


/* current values */
double recv_fit_dev(const struct recv_fit *fit);
double recv_fit_score_norm(const struct recv_fit *fit);
double recv_fit_step_size(const struct recv_fit *fit);

const struct recv_params *recv_fit_params(const struct recv_fit *fit);


double recv_fit_dev0(const struct recv_fit *fit);
const double *recv_fit_duals(const struct recv_fit *fit);


#endif /* RECV_FIT_H */
