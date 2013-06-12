#ifndef RECV_RESID_H
#define RECV_RESID_H

#include "history.h"
#include "recv_model.h"

struct recv_resid_count {
	double *dyad;
	double *send;
	double *recv;
	double tot;
};

struct recv_resid {
	struct recv_resid_count obs;
	struct recv_resid_count fit;
};

void recv_resid_init(struct recv_resid *resid,
		     struct recv_model *model,
		     const struct message *msgs,
		     size_t nmsg);
void recv_resid_deinit(struct recv_resid *resid);

#endif /* RECV_RESID_H */
