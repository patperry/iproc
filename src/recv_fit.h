#ifndef _RECV_FIT_H
#define _RECV_FIT_H

#include "design.h"
#include "linesearch.h"
#include "messages.h"
#include "recv_loglik.h"
#include "recv_model.h"

#define RECV_FIT_GTOL0		(1e-8)
#define RECV_FIT_XTOL0		(1e-8)
#define RECV_FIT_LSMAX0		(10)
#define RECV_FIT_LSCTRL0	LINESEARCH_CTRL0
#define RECV_FIT_VARTOL		(1e-8)
#define RECV_FIT_EIGTOL		(1e-8)
#define RECV_FIT_SVTOL		(1e-8)

#define RECV_FIT_CTRL0 ((struct recv_fit_ctrl) { \
		RECV_FIT_GTOL0, \
		RECV_FIT_XTOL0, \
		RECV_FIT_LSMAX0, \
		RECV_FIT_LSCTRL0, \
		RECV_FIT_VARTOL, \
		RECV_FIT_EIGTOL, \
		RECV_FIT_SVTOL \
	})

struct recv_fit_ctrl {
	double gtol;
	double xtol;
	size_t ls_maxit;
	struct linesearch_ctrl ls;
	double vartol;
	double eigtol;
	double svtol;
};

enum recv_fit_task {
	RECV_FIT_CONV = 0,
	RECV_FIT_STEP = 1,
	RECV_FIT_ERR_LNSRCH = -1,	// linesearch failed to converge
	RECV_FIT_ERR_XTOL = -2,	// step size is smaller than tolerance
};

struct recv_fit_resid {
	double *vector;
	double norm2;
};

struct recv_fit_constr {
	double *wts;
	size_t *wt_inds;
	size_t *wt_cols;
	double *vals;
	char **names;
	size_t n, nmax, wt_nzmax;
};

struct recv_fit_eval {
	double *params;
	struct dmatrix coefs;
	double *duals;
	struct recv_loglik loglik;
	struct recv_fit_resid resid;
};

struct recv_fit_kkt {
	struct dmatrix matrix;
	ptrdiff_t *ldl_ipiv;
	double *ldl_work;
	size_t ldl_lwork;
	enum blas_uplo uplo;
	int factored;
};

struct recv_fit_search {
	double *vector;
};

struct recv_fit_rgrad {
	double *vector;
};

struct recv_fit {
	struct recv_fit_ctrl ctrl;
	const struct messages *xmsgs, *ymsgs;

	/* working frame + model */
	struct frame *frame;
	struct recv_model model;

	/* null deviance */
	double dev0;

	/* optimization constraints */
	struct recv_fit_constr constr;

	/* optimization workspace */
	struct recv_fit_eval eval[2], *cur, *prev;
	struct recv_fit_kkt kkt;
	struct recv_fit_search search;
	struct recv_fit_rgrad rgrad;

	/* additional workspace */
	struct linesearch ls;
	enum recv_fit_task task;
	double step;
};

void recv_fit_init(struct recv_fit *fit,
		   struct frame *f,
		   const struct messages *xmsgs,
		   const struct messages *ymsgs,
		   size_t ncohort,
		   const size_t *cohorts, const struct recv_fit_ctrl *ctrl);
void recv_fit_deinit(struct recv_fit *fit);

/* constraints */
size_t recv_fit_constr_count(const struct recv_fit *fit);
void recv_fit_get_constr(const struct recv_fit *fit, size_t i,
			 const double **weightsp, const size_t **indp,
			 size_t *nzp, double *valp, const char **namep);
void recv_fit_add_constr(struct recv_fit *fit, const double *weights,
			 const size_t *ind, size_t nz, double val,
			 const char *name);
void recv_fit_add_constr_set(struct recv_fit *fit, size_t i, size_t c,
			     double val);
void recv_fit_add_constr_eq(struct recv_fit *fit, size_t i1, size_t c1,
			    size_t i2, size_t c2);
//size_t recv_fit_add_constr_identify(struct recv_fit *fit);

/* fitting */
enum recv_fit_task recv_fit_start(struct recv_fit *fit,
				  const struct dmatrix *coefs0);
enum recv_fit_task recv_fit_advance(struct recv_fit *fit);
const char *recv_fit_errmsg(const struct recv_fit *fit);

/* current values */
const struct recv_loglik *recv_fit_loglik(const struct recv_fit *fit);
double recv_fit_dev(const struct recv_fit *fit);
double recv_fit_dev0(const struct recv_fit *fit);
const struct dmatrix *recv_fit_coefs(const struct recv_fit *fit);
const double *recv_fit_duals(const struct recv_fit *fit);
double recv_fit_step(const struct recv_fit *fit);
double recv_fit_grad_norm2(const struct recv_fit *fit);

/* control parameters */
static inline int recv_fit_ctrl_valid(const struct recv_fit_ctrl *ctrl);

/* inline function definitions */
int recv_fit_ctrl_valid(const struct recv_fit_ctrl *ctrl)
{
	assert(ctrl);

	if (!(ctrl->gtol > 0)) {
		return 0;
	} else if (!(ctrl->ls_maxit > 0)) {
		return 0;
	} else {
		return linesearch_ctrl_valid(&ctrl->ls);
	}
}

#endif /* _RECV_FIT_H */
