#ifndef _RECV_RESID_H
#define _RECV_RESID_H

#include "messages.h"
#include "design.h"
#include "actors.h"
#include "matrix.h"
#include "recv_model.h"


struct recv_resid_count {
	struct matrix dyad;
	struct vector send;
	struct vector recv;
	ssize_t tot;
	
	// internal
	struct matrix dyad_trans;
	bool dyad_cached;
};

struct recv_resid {
	struct frame frame;
	struct recv_model model;
	struct recv_resid_count obs;
	struct recv_resid_count exp;
};

void recv_resid_init(struct recv_resid *resid,
		     const struct messages *msgs,
		     const struct design *design,
		     const struct actors *senders,
		     const struct matrix *coefs);

void recv_resid_deinit(struct recv_resid *resid);


#endif /* _RECV_RESID_H */