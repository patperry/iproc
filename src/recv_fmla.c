#include "port.h"
#include "xalloc.h"
#include <stdlib.h>
#include <string.h>
#include "coreutil.h"
#include "strata.h"
#include "recv_fmla.h"


void recv_fmla_init(struct recv_fmla *fmla, const struct frame *f)
{
	const struct design *r = frame_recv_design(f);
	const struct design2 *d = frame_dyad_design(f);
	fmla->frame = (struct frame *)f;
	fmla->trait_dim = design_trait_dim(r) + design2_trait_dim(d);
	fmla->tvar_dim = design_tvar_dim(r) + design2_tvar_dim(d);
}


void recv_fmla_deinit(struct recv_fmla *fmla)
{
	(void)fmla;
}


void recv_coefs_init(struct recv_coefs *c, const struct recv_fmla *fmla)
{
	c->fmla = fmla;
	c->traits = xmalloc(recv_fmla_trait_dim(fmla) * sizeof(*c->traits));
	c->tvars = xmalloc(recv_fmla_tvar_dim(fmla) * sizeof(*c->tvars));
}


void recv_coefs_deinit(struct recv_coefs *c)
{
	free(c->tvars);
	free(c->traits);
}


void recv_coefs_init_copy(struct recv_coefs *c, const struct recv_coefs *c0)
{
	recv_coefs_init(c, c0->fmla);
	recv_coefs_assign_copy(c, c0);
}


void recv_coefs_assign_copy(struct recv_coefs *dst, const struct recv_coefs *src)
{
	assert(dst->fmla == src->fmla);
	
	const struct recv_fmla *fmla = dst->fmla;
	memcpy(dst->traits, src->traits, recv_fmla_trait_dim(fmla) * sizeof(*dst->traits));
	memcpy(dst->tvars, src->tvars, recv_fmla_tvar_dim(fmla) * sizeof(*dst->tvars));
}


void recv_coefs_clear(struct recv_coefs *c)
{
	const struct recv_fmla *fmla = c->fmla;
	memset(c->traits, 0, recv_fmla_trait_dim(fmla) * sizeof(*c->traits));
	memset(c->tvars, 0, recv_fmla_tvar_dim(fmla) * sizeof(*c->tvars));
}


static void recv_frame_frame_clear(void *udata, struct frame *f)
{
	struct recv_frame *rf = udata;
	size_t iobs, nobs = rf->nobs;
	for (iobs = 0; iobs < nobs; iobs++) {
		struct recv_frame_observer *obs = &rf->obs[iobs];
		if (obs->callbacks.clear) {
			obs->callbacks.clear(obs->udata, rf);
		}
	}
}


static struct frame_callbacks recv_frame_frame_callbacks = {
	NULL,
	NULL,
	recv_frame_frame_clear
};


static void recv_frame_recv_design_update(void *udata, struct design *d,  const struct var *v, size_t i, const double *delta, const struct vpattern *pat)
{
	struct recv_frame *rf = udata;
	
	if (!rf->nobs)
		return;
	
	const struct recv_fmla *fmla = recv_frame_fmla(rf);
	const struct frame *f = recv_fmla_frame(fmla);
	size_t isend, nsend = frame_send_count(f);
	
	size_t off = v->index;
	size_t iz, nz;
	
	if (pat) {
		nz = pat->nz;
		for (iz = 0; iz < nz; iz++) {
			rf->pat_buf[iz]  = off + pat->indx[iz];
		}
	} else {
		nz = v->dim;
		for (iz = 0; iz < nz; iz++) {
			rf->pat_buf[iz] = off + iz;
		}
	}
	
	struct vpattern my_pat = vpattern_make(rf->pat_buf, nz);
	size_t iobs, nobs = rf->nobs;
	for (iobs = 0; iobs < nobs; iobs++) {
		struct recv_frame_observer *obs = &rf->obs[iobs];
		if (obs->callbacks.update_all) {
			obs->callbacks.update_all(obs->udata, rf, i, delta, &my_pat);
		} else if (obs->callbacks.update) {
			for (isend = 0; isend < nsend; isend++) {
				/* NOTE: this is really inefficient */
				obs->callbacks.update(obs->udata, rf, isend, i, delta, &my_pat);
			}
		}
	}
}


static struct design_callbacks recv_frame_recv_design_callbacks = {
	recv_frame_recv_design_update
};


static void recv_frame_dyad_design_update(void *udata, struct design2 *d,  const struct var2 *v, size_t i, size_t j, const double *delta, const struct vpattern *pat)
{
	struct recv_frame *rf = udata;
	
	if (!rf->nobs)
		return;
	
	const struct recv_fmla *fmla = recv_frame_fmla(rf);
	const struct frame *f = recv_fmla_frame(fmla);
	const struct design *r = frame_recv_design(f);
	
	size_t off = design_tvar_dim(r) + v->index;
	size_t iz, nz;
	
	if (pat) {
		nz = pat->nz;
		for (iz = 0; iz < nz; iz++) {
			rf->pat_buf[iz]  = off + pat->indx[iz];
		}
	} else {
		nz = v->dim;
		for (iz = 0; iz < nz; iz++) {
			rf->pat_buf[iz] = off + iz;
		}
	}
	
	struct vpattern my_pat = vpattern_make(rf->pat_buf, nz);
	size_t iobs, nobs = rf->nobs;
	for (iobs = 0; iobs < nobs; iobs++) {
		struct recv_frame_observer *obs = &rf->obs[iobs];
		if (obs->callbacks.update) {
			obs->callbacks.update(obs->udata, rf, i, j, delta, &my_pat);
		}
	}
}


static struct design2_callbacks recv_frame_dyad_design_callbacks = {
	recv_frame_dyad_design_update
};


void recv_frame_init(struct recv_frame *rf, const struct recv_fmla *fmla)
{
	struct frame *f = recv_fmla_frame(fmla);
	struct design *r = frame_recv_design(f);
	struct design2 *d = frame_dyad_design(f);
	
	rf->fmla = (struct recv_fmla *)fmla;
	rf->pat_buf = xmalloc(fmla->tvar_dim * sizeof(*rf->pat_buf));
	rf->obs = NULL;
	rf->nobs = 0;
	rf->nobs_max = 0;
	
	frame_add_observer(f, rf, &recv_frame_frame_callbacks);
	design_add_observer(r, rf, &recv_frame_recv_design_callbacks);
	design2_add_observer(d, rf, &recv_frame_dyad_design_callbacks);	
}


void recv_frame_deinit(struct recv_frame *rf)
{
	const struct recv_fmla *fmla = recv_frame_fmla(rf);
	struct frame *f = recv_fmla_frame(fmla);
	struct design *r = frame_recv_design(f);
	struct design2 *d = frame_dyad_design(f);
	
	design2_remove_observer(d, rf);	
	design_remove_observer(r, rf);
	frame_remove_observer(f, rf);
	
	free(rf->obs);
	free(rf->pat_buf);
}


void recv_frame_mul(double alpha, const struct recv_frame *rf, size_t isend,
		    const struct recv_coefs *c, double beta, double *y)
{
	recv_frame_mul0(alpha, rf, isend, c->traits, beta, y);
	recv_frame_mul1(alpha, rf, isend, c->tvars, 1.0, y);
}


void recv_frame_tmul(double alpha, const struct recv_frame *rf, size_t isend,
		     const double *x, double beta, struct recv_coefs *c)
{
	recv_frame_tmul0(alpha, rf, isend, x, beta, c->traits);
	recv_frame_tmul1(alpha, rf, isend, x, beta, c->tvars);
}


void recv_frame_axpy(double alpha, const struct recv_frame *rf, size_t isend,
		     size_t jrecv, struct recv_coefs *c)
{
	recv_frame_axpy0(alpha, rf, isend, jrecv, c->traits);
	recv_frame_axpy1(alpha, rf, isend, jrecv, c->tvars);
}


void recv_frame_mul0(double alpha, const struct recv_frame *rf, size_t isend,
		     const double *x, double beta, double *y)
{
	const struct recv_fmla *fmla = recv_frame_fmla(rf);
	const struct frame *f = recv_fmla_frame(fmla);
	const struct design *r = frame_recv_design(f);
	const struct design2 *d = frame_dyad_design(f);
	size_t dimr = design_trait_dim(r);
	const double *xr = x;
	const double *xd = xr + dimr;	
	
	design_traits_mul(alpha, r, xr, beta, y);
	design2_traits_mul(alpha, d, isend, xd, 1.0, y);
}


void recv_frame_tmul0(double alpha, const struct recv_frame *rf, size_t isend,
		      const double *x, double beta, double *y)
{
	const struct recv_fmla *fmla = recv_frame_fmla(rf);
	const struct frame *f = recv_fmla_frame(fmla);
	const struct design *r = frame_recv_design(f);
	const struct design2 *d = frame_dyad_design(f);
	size_t dimr = design_trait_dim(r);
	double *yr = y;
	double *yd = yr + dimr;	
	
	design_traits_tmul(alpha, r, x, beta, yr);
	design2_traits_tmul(alpha, d, isend, x, beta, yd);
}


void recv_frame_axpy0(double alpha, const struct recv_frame *rf, size_t isend, size_t jrecv, double *y)
{
	const struct recv_fmla *fmla = recv_frame_fmla(rf);
	const struct frame *f = recv_fmla_frame(fmla);
	const struct design *r = frame_recv_design(f);
	const struct design2 *d = frame_dyad_design(f);
	size_t dimr = design_trait_dim(r);
	double *yr = y;
	double *yd = yr + dimr;	
	
	design_traits_axpy(alpha, r, isend, yr);
	design2_traits_axpy(alpha, d, isend, jrecv, yd);
}


void recv_frame_mul1(double alpha, const struct recv_frame *rf, size_t isend,
		     const double *x, double beta, double *y)
{
	const struct recv_fmla *fmla = recv_frame_fmla(rf);
	const struct frame *f = recv_fmla_frame(fmla);
	const struct design *r = frame_recv_design(f);
	const struct design2 *d = frame_dyad_design(f);
	size_t dimr = design_tvar_dim(r);
	const double *xr = x;
	const double *xd = xr + dimr;	
	
	design_tvars_mul(alpha, r, xr, beta, y);
	design2_tvars_mul(alpha, d, isend, xd, 1.0, y);
}


void recv_frame_tmul1(double alpha, const struct recv_frame *rf, size_t isend,
		      const double *x, double beta, double *y)
{
	const struct recv_fmla *fmla = recv_frame_fmla(rf);
	const struct frame *f = recv_fmla_frame(fmla);
	const struct design *r = frame_recv_design(f);
	const struct design2 *d = frame_dyad_design(f);
	size_t dimr = design_tvar_dim(r);
	double *yr = y;
	double *yd = yr + dimr;	
	
	design_tvars_tmul(alpha, r, x, beta, yr);
	design2_tvars_tmul(alpha, d, isend, x, beta, yd);
}


void recv_frame_axpy1(double alpha, const struct recv_frame *rf, size_t isend,
		      size_t jrecv, double *y)
{
	const struct recv_fmla *fmla = recv_frame_fmla(rf);
	const struct frame *f = recv_fmla_frame(fmla);
	const struct design *r = frame_recv_design(f);
	const struct design2 *d = frame_dyad_design(f);
	size_t dimr = design_tvar_dim(r);
	double *yr = y;
	double *yd = yr + dimr;	
	
	design_tvars_axpy(alpha, r, isend, yr);
	design2_tvars_axpy(alpha, d, isend, jrecv, yd);
}


static void recv_frame_observers_grow(struct recv_frame *rf, size_t delta)
{
	size_t nmax = array_grow(rf->nobs, rf->nobs_max, delta, SIZE_MAX);
	if (nmax > rf->nobs_max) {
		rf->obs = xrealloc(rf->obs, nmax * sizeof(*rf->obs));
		rf->nobs_max = nmax;
	}
}



void recv_frame_add_observer(struct recv_frame *rf, void *udata,
			     const struct recv_frame_callbacks *callbacks)
{
	recv_frame_observers_grow(rf, 1);
	struct recv_frame_observer *obs = &rf->obs[rf->nobs++];
	obs->udata = udata;
	obs->callbacks = *callbacks;
}


void recv_frame_remove_observer(struct recv_frame *rf, void *udata)
{
	size_t i, n = rf->nobs;
	for (i = n; i > 0; i--) {
		if (rf->obs[i-1].udata == udata) {
			memmove(rf->obs + i - 1, rf->obs + i, (n - i) * sizeof(*rf->obs));
			rf->nobs = n - 1;
			break;
		}
	}
}

