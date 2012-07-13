#ifndef RECV_FMLA_H
#define RECV_FMLA_H

#include "frame.h"

struct recv_fmla {
	struct frame *frame;
	size_t trait_dim;
	size_t tvar_dim;
	size_t ncohort;
	size_t *cohorts;
	size_t *cohort_reps;
};

struct recv_coefs {
	const struct recv_fmla *fmla;
	double *traits;
	double *tvars;
};

//struct recv_frame_sender {
//	struct vpattern active;
//};

struct recv_frame {
	struct recv_fmla *fmla;
	size_t *pat_buf;	
	//	struct recv_frame_sender *senders;
	struct recv_frame_observer *obs;
	size_t nobs, nobs_max;
};

struct recv_frame_callbacks {
	void (*update) (void *udata, struct recv_frame *rf, size_t isend, size_t jrecv,
				const double *delta, const struct vpattern *pat);
	void (*update_all) (void *udata, struct recv_frame *rf, size_t jrecv,
			    const double *delta, const struct vpattern *pat);
	void (*clear) (void *udata, struct recv_frame *rf);
};

struct recv_frame_observer {
	void *udata;
	struct recv_frame_callbacks callbacks;
};


void recv_fmla_init(struct recv_fmla *fmla, const struct frame *f);
void recv_fmla_deinit(struct recv_fmla *fmla);

static inline struct frame *recv_fmla_frame(const struct recv_fmla *fmla);
static inline size_t recv_fmla_trait_dim(const struct recv_fmla *fmla);
static inline size_t recv_fmla_tvar_dim(const struct recv_fmla *fmla);
static inline size_t recv_fmla_cohort(const struct recv_fmla *fmla, size_t isend);
static inline size_t recv_fmla_cohort_rep(const struct recv_fmla *fmla, size_t c);
static inline size_t recv_fmla_cohort_count(const struct recv_fmla *fmla);
static inline void recv_fmla_get_cohorts(const struct recv_fmla *fmla, const size_t **cohortsp, const size_t **repsp, size_t *ncohortp);

//void recv_fmla_add_kron(struct recv_fmla *fmla, const struct var *s,
//			const struct var *r);


void recv_coefs_init(struct recv_coefs *c, const struct recv_fmla *fmla);
void recv_coefs_deinit(struct recv_coefs *c);
void recv_coefs_init_copy(struct recv_coefs *c, const struct recv_coefs *c0);
void recv_coefs_assign_copy(struct recv_coefs *dst, const struct recv_coefs *src);
void recv_coefs_clear(struct recv_coefs *c);


void recv_frame_init(struct recv_frame *rf, const struct recv_fmla *fmla);
void recv_frame_deinit(struct recv_frame *rf);

static inline const struct recv_fmla *recv_frame_fmla(const struct recv_frame *rf);

void recv_frame_mul(double alpha, const struct recv_frame *rf, size_t isend,
		    const struct recv_coefs *c, double beta, double *y);
void recv_frame_tmul(double alpha, const struct recv_frame *rf, size_t isend,
		     const double *x, double beta, struct recv_coefs *c);
void recv_frame_axpy(double alpha, const struct recv_frame *rf, size_t isend,
		     size_t jrecv, struct recv_coefs *c);

void recv_frame_mul0(double alpha, const struct recv_frame *rf, size_t isend,
		     const double *x, double beta, double *y);
void recv_frame_tmul0(double alpha, const struct recv_frame *rf, size_t isend,
		      const double *x, double beta, double *y);
void recv_frame_axpy0(double alpha, const struct recv_frame *rf, size_t isend,
		      size_t jrecv, double *y);

void recv_frame_mul1(double alpha, const struct recv_frame *rf, size_t isend,
		     const double *x, double beta, double *y);
void recv_frame_tmul1(double alpha, const struct recv_frame *rf, size_t isend,
		      const double *x, double beta, double *y);
void recv_frame_axpy1(double alpha, const struct recv_frame *rf, size_t isend,
		      size_t jrecv, double *y);

/* observers */
void recv_frame_add_observer(struct recv_frame *rf, void *udata,
			     const struct recv_frame_callbacks *callbacks);
void recv_frame_remove_observer(struct recv_frame *rf, void *udata);



//void recv_frame_get_active(const struct recv_frame *rf, size_t isend,
//			   const double **dxp, const size_t **activep,
//			   size_t *nactivep);


/* inline function definitions */
struct frame *recv_fmla_frame(const struct recv_fmla *fmla)
{
	return fmla->frame;
}

size_t recv_fmla_trait_dim(const struct recv_fmla *fmla)
{
	return fmla->trait_dim;
}

size_t recv_fmla_tvar_dim(const struct recv_fmla *fmla)
{
	return fmla->tvar_dim;
}


size_t recv_fmla_cohort(const struct recv_fmla *fmla, size_t isend)
{
	assert(isend < frame_recv_count(recv_fmla_frame(fmla)));
	return fmla->cohorts[isend];
}


size_t recv_fmla_cohort_rep(const struct recv_fmla *fmla, size_t c)
{
	assert(c < recv_fmla_cohort_count(fmla));
	return fmla->cohort_reps[c];
}


size_t recv_fmla_cohort_count(const struct recv_fmla *fmla)
{
	return fmla->ncohort;
}


void recv_fmla_get_cohorts(const struct recv_fmla *fmla, const size_t **cohortsp, const size_t **repsp, size_t *ncohortp)
{
	*cohortsp = fmla->cohorts;
	*repsp = fmla->cohort_reps;
	*ncohortp = fmla->ncohort;
}


const struct recv_fmla *recv_frame_fmla(const struct recv_frame *rf)
{
	return rf->fmla;
}



#endif /* RECV_FMLA_H */
