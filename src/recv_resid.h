#ifndef _RECV_RESID_H
#define _RECV_RESID_H

#include "messages.h"
#include "design.h"
#include "recv_model.h"
#include "blas.h"

struct recv_resid_count {
	struct dmatrix dyad;
	double *send;
	double *recv;
	size_t tot;

	// internal
	struct dmatrix dyad_trans;
	bool dyad_cached;
};

struct recv_resid {
	struct recv_model model;
	struct recv_resid_count obs;
	struct recv_resid_count exp;
};

void recv_resid_init(struct recv_resid *resid,
		     struct frame *f,
		     const struct messages *msgs,
		     size_t ncohort,
		     const size_t *cohorts, const struct dmatrix *coefs);
void recv_resid_deinit(struct recv_resid *resid);

#endif /* _RECV_RESID_H */
