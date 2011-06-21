#ifndef _BFGS_H
#define _BFGS_H

#include "linesearch.h"
#include "matrix.h"
#include "vector.h"


struct bfgs_ctrl {
	struct linesearch_ctrl ls;
};

struct bfgs {
	struct bfgs_ctrl ctrl;
	struct matrix inv_hess;
	struct vector search_dir;
	struct vector *grad0;
	struct vector *grad;
	double step;

	/* workspace */
	struct vector g0;
	struct vector g1;
	struct vector s;
	struct vector y;
	struct vector H_y;
	bool started;
	bool swapped;
};



void bfgs_init(struct bfgs *opt, ssize_t n, const struct bfgs_ctrl *ctrl);
void bfgs_deinit(struct bfgs *opt);
void bfgs_clear(struct bfgs *opt);

const struct vector *bfgs_start(struct bfgs *opt, double f0, const struct vector *grad0);
const struct vector *bfgs_advance(struct bfgs *opt, double f, const struct vector *grad);

static inline bool bfgs_ctrl_valid(const struct bfgs_ctrl *ctrl)
{
	assert(ctrl);
	
	return linesearch_ctrl_valid(&ctrl->ls);
}

#endif /* _BFGS_H */
