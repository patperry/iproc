#ifndef RECV_RESID_H
#define RECV_RESID_H

#include "messages.h"
#include "recv_model.h"

struct recv_resid_count {
	double *dyad;
	double *send;
	double *recv;
	size_t tot;
};

struct recv_resid {
	struct recv_model model;
	struct recv_resid_count obs;
	struct recv_resid_count exp;
};

void recv_resid_init(struct recv_resid *resid,
		     struct frame *f,
		     const struct messages *msgs,
		     const struct recv_coefs *coefs);
void recv_resid_deinit(struct recv_resid *resid);

#endif /* RECV_RESID_H */
