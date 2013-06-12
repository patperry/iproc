#ifndef NEWTON_H
#define NEWTON_H

#include <stddef.h>
#include "blas.h"
#include "constr.h"
#include "linesearch.h"


#define NEWTON_HESS_UPLO	BLAS_UPPER
#define NEWTON_HESS_F77UPLO	(NEWTON_HESS_UPLO == BLAS_UPPER ? BLAS_LOWER : BLAS_UPPER)


#define NEWTON_GTOL0	(1e-8)
#define NEWTON_XTOL0	(1e-8)
#define NEWTON_LSMAX0	(10)
#define NEWTON_LSCTRL0	LINESEARCH_CTRL0

#define NEWTON_CTRL0 ((struct newton_ctrl) { \
		NEWTON_GTOL0, \
		NEWTON_XTOL0, \
		NEWTON_LSMAX0, \
		NEWTON_LSCTRL0 \
	})


struct newton_ctrl {
	double gtol;
	double xtol;
	size_t ls_maxit;
	struct linesearch_ctrl ls;
};


enum newton_task {
	NEWTON_CONV = 0,
	NEWTON_STEP = 1,
	NEWTON_HESS = 2,
	NEWTON_ERR_LNSRCH = -1,	/* linesearch failed to converge */
	NEWTON_ERR_XTOL = -2,	/* step size is smaller than tolerance */
	NEWTON_ERR_HESS = -3,   /* Hessian is singular */
	NEWTON_ERR_XDOM = -4    /* starting point is not in domain */
};


struct newton_params {
	double *pri;
	double *dual;
};

struct newton_eval {
	struct newton_params params;
	double val;
	double *grad;

	struct newton_params resid;
	double resid_norm;         /* ||r|| */
	int feasible;
	int indomain;
};

struct newton_kkt {
	double *matrix;
	ptrdiff_t *ldl_ipiv;
	double *ldl_work;
	size_t ldl_lwork;
	int factored;
};



struct newton {
	struct newton_ctrl ctrl;
	size_t dim;
	const struct constr *constr;

	double *hess;

	/* workspace */
	struct newton_eval eval[2], *cur, *next;
	struct newton_kkt kkt;
	struct newton_params search;
	double step;
};


void newton_init(struct newton *opt, size_t dim,
		 const struct constr *constr,
		 const struct newton_ctrl *ctrl);
void newton_deinit(struct newton *opt);

static inline size_t newton_dim(const struct newton *opt) { return opt->dim; }
static inline const struct constr *newton_constr(const struct newton *opt) { return opt->constr; }
static inline const struct newton_ctrl *newton_ctrl(const struct newton *opt) { return &opt->ctrl; }


enum newton_task newton_start(struct newton *opt, const double *x0,
			      double f0, const double *grad0,
			      const double *hess0);

const double *newton_next(const struct newton *opt);

enum newton_task newton_step(struct newton *opt, double f, const double *grad);
enum newton_task newton_set_hess(struct newton *opt, const double *hess);

const char *newton_errmsg(enum newton_task task);

/* current values */
const double *newton_params(const struct newton *opt);
const double *newton_duals(const struct newton *opt);
double newton_val(const struct newton *opt);
const double *newton_grad(const struct newton *opt);
const double *newton_hess(const struct newton *opt);


int newton_ctrl_valid(const struct newton_ctrl *ctrl);



#endif /* NEWTON_H */
