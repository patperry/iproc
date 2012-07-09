#ifndef RECV_FMLA_H
#define RECV_FMLA_H

#include "frame.h"

struct recv_fmla {
	const struct frame *frame;
};

struct recv_coefs {
	double *traits;
	double *tvars;
};

struct recv_frame_sender {
	struct vpattern active;
};

struct recv_frame {
	const struct recv_fmla *fmla;
	struct recv_frame_sender *senders;
};

void recv_fmla_init(struct recv_fmla *fmla, const struct frame *f);
void recv_fmla_deinit(struct recv_fmla *fmla);

void recv_fmla_add_kron(struct recv_fmla *fmla, const struct var *s,
			const struct var *r);

void recv_coefs_init(struct recv_coefs *c, const struct recv_fmla *fmla);
void recv_coefs_deinit(struct recv_coefs *c);

void recv_frame_init(struct recv_frame *rf, const struct recv_fmla *fmla);
void recv_frame_deinit(struct recv_frame *rf);

void recv_frame_mul(double alpha, const struct recv_frame *rf, size_t isend,
		    const struct recv_coefs *c, double beta, double *y);
void recv_frame_tmul(double alpha, const struct recv_frame *rf, size_t isend,
		     const double *x, double beta, struct recv_coefs *c);
void recv_frame_axpy(double alpha, const struct recv_frame *rf, size_t isend,
		     size_t jrecv, struct recv_coefs *c);

void recv_frame_mul0(double alpha, const struct recv_frame *rf, size_t isend,
		     size_t jrecv, const double *x, double beta, double *y);
void recv_frame_tmul0(double alpha, const struct recv_frame *rf, size_t isend,
		      size_t jrecv, const double *x, double beta, double *y);
void recv_frame_axpy0(double alpha, const struct recv_frame *rf, size_t isend,
		      size_t jrecv, double *y);


void recv_frame_get_active(const struct recv_frame *rf, size_t isend,
			   const double **dxp, const size_t **activep,
			   size_t *nactivep);



#endif /* RECV_FMLA_H */
