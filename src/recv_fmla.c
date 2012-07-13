#include "port.h"
#include "xalloc.h"
#include <stdlib.h>
#include <string.h>
#include "coreutil.h"
#include "strata.h"
#include "recv_fmla.h"


static void recv_fmla_set_cohorts(struct recv_fmla *fmla)
{
	const struct frame *f = recv_fmla_frame(fmla);
	const struct design *s = frame_send_design(f);
	const struct dmatrix *as = design_traits(s);
	size_t dims = design_trait_dim(s);
	size_t nsend = frame_send_count(f);
	size_t i;
	
	struct strata strat;
	double *buf;
	
	strata_init(&strat, dims);
	buf = xmalloc(dims * sizeof(*buf));	

	for (i = 0; i < nsend; i++) {
		blas_dcopy(dims, MATRIX_PTR(as, i, 0), as->lda, buf, 1);
		fmla->cohorts[i] = strata_add(&strat, buf);
	}
	fmla->ncohort = strata_count(&strat);
	
	free(buf);
	strata_deinit(&strat);
	
	size_t ic, nc = fmla->ncohort;
	
	fmla->cohort_reps = xmalloc(nc * sizeof(*fmla->cohort_reps));	
	
#ifndef NDEBUG
	/* initialize all cohort_reps to invalid values */
	for (ic = 0; ic < nc; ic++) {
		fmla->cohort_reps[ic] = nsend;
	}
#endif

	for (i = nsend; i > 0; i--) {
		size_t isend = i - 1;
		ic = fmla->cohorts[isend];
		fmla->cohort_reps[ic] = isend;
	}
	
#ifndef NDEBUG
	/* ensure all cohorts have reps */
	for (ic = 0; ic < nc; ic++) {
		assert(fmla->cohort_reps[ic] != nsend);
	}
#endif
}



void recv_fmla_init(struct recv_fmla *fmla, const struct frame *f)
{
	const struct design *r = frame_recv_design(f);
	const struct design2 *d = frame_dyad_design(f);
	size_t nsend = frame_send_count(f);
	fmla->frame = (struct frame *)f;
	fmla->trait_dim = design_trait_dim(r) + design2_trait_dim(d);
	fmla->tvar_dim = design_tvar_dim(r) + design2_tvar_dim(d);
	fmla->cohorts = xmalloc(nsend * sizeof(*fmla->cohorts));
	fmla->ncohort = 0;
	
	recv_fmla_set_cohorts(fmla);
}


void recv_fmla_deinit(struct recv_fmla *fmla)
{
	free(fmla->cohort_reps);
	free(fmla->cohorts);
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


static void scale_result(size_t n, double beta, double *y)
{
	if (beta == 1.0) {
		return;
	} else if (beta == 0.0) {
		memset(y, 0, n * sizeof(*y));
	} else {
		blas_dscal(n, beta, y, 1);
	}
}

void recv_frame_mul0(double alpha, const struct recv_frame *rf, size_t isend,
		     const double *x, double beta, double *y)
{
	const struct recv_fmla *fmla = recv_frame_fmla(rf);
	const struct frame *f = recv_fmla_frame(fmla);
	size_t nrecv = frame_recv_count(f);
	
	const struct design *r = frame_recv_design(f);
	const struct dmatrix *ar = design_traits(r);
	size_t dimr = design_trait_dim(r);	
	const double *xr = x;

	blas_dgemv(BLAS_NOTRANS, nrecv, dimr, alpha, ar, xr, 1, beta, y, 1);
	
	const struct design2 *d = frame_dyad_design(f);
	const struct dmatrix *ad = design2_traits(d);
	size_t dimd = design2_trait_dim(d);
	const double *xd = xr + dimr;

	struct dmatrix adsub = *ad;
	adsub.data = MATRIX_PTR(ad, isend * nrecv, 0);
	blas_dgemv(BLAS_NOTRANS, nrecv, dimd, alpha, &adsub, xd, 1, 1.0, y, 1);
}


void recv_frame_tmul0(double alpha, const struct recv_frame *rf, size_t isend,
		      const double *x, double beta, double *y)
{
	const struct recv_fmla *fmla = recv_frame_fmla(rf);
	const struct frame *f = recv_fmla_frame(fmla);
	size_t nrecv = frame_recv_count(f);
	
	const struct design *r = frame_recv_design(f);
	const struct dmatrix *ar = design_traits(r);
	size_t dimr = design_trait_dim(r);	
	double *yr = y;
	
	blas_dgemv(BLAS_TRANS, nrecv, dimr, alpha, ar, x, 1, beta, yr, 1);
	
	const struct design2 *d = frame_dyad_design(f);
	const struct dmatrix *ad = design2_traits(d);
	size_t dimd = design2_trait_dim(d);
	double *yd = yr + dimr;
	
	struct dmatrix adsub = *ad;
	adsub.data = MATRIX_PTR(ad, isend * nrecv, 0);
	blas_dgemv(BLAS_TRANS, nrecv, dimd, alpha, &adsub, x, 1, beta, yd, 1);
}


void recv_frame_axpy0(double alpha, const struct recv_frame *rf, size_t isend, size_t jrecv, double *y)
{
	const struct recv_fmla *fmla = recv_frame_fmla(rf);
	const struct frame *f = recv_fmla_frame(fmla);
	size_t nrecv = frame_recv_count(f);
	
	const struct design *r = frame_recv_design(f);
	const struct dmatrix *ar = design_traits(r);
	size_t dimr = design_trait_dim(r);	
	double *yr = y;
	
	blas_daxpy(dimr, alpha, MATRIX_PTR(ar, jrecv, 0), ar->lda, yr, 1);
	
	const struct design2 *d = frame_dyad_design(f);
	const struct dmatrix *ad = design2_traits(d);
	size_t dimd = design2_trait_dim(d);
	double *yd = yr + dimr;

	blas_daxpy(dimd, alpha, MATRIX_PTR(ad, isend * nrecv + jrecv, 0), ad->lda, yd, 1);
}


void recv_frame_mul1(double alpha, const struct recv_frame *rf, size_t isend,
		     const double *x, double beta, double *y)
{
	const struct recv_fmla *fmla = recv_frame_fmla(rf);
	const struct frame *f = recv_fmla_frame(fmla);
	size_t nrecv = frame_recv_count(f);
	
	const struct design *r = frame_recv_design(f);
	const double *ar;
	const size_t *j;
	size_t nz;
	size_t dimr = design_tvar_dim(r);
	const double *xr = x;
	
	scale_result(nrecv, beta, y);
	
	design_tvars_get(r, &ar, &j, &nz);
	for (; nz != 0; ar += dimr, j++, nz--) {
		y[*j] += alpha * blas_ddot(dimr, ar, 1, xr, 1);
	}

	
	const struct design2 *d = frame_dyad_design(f);
	const double *ad;
	size_t dimd = design2_tvar_dim(d);
	const double *xd = xr + dimr;
	
	design2_tvars_get(d, isend, &ad, &j, &nz);
	
	for (; nz != 0; ad += dimd, j++, nz--) {
		y[*j] += alpha * blas_ddot(dimd, ad, 1, xd, 1);
	}
}


void recv_frame_tmul1(double alpha, const struct recv_frame *rf, size_t isend,
		      const double *x, double beta, double *y)
{
	const struct recv_fmla *fmla = recv_frame_fmla(rf);
	const struct frame *f = recv_fmla_frame(fmla);
	
	const struct design *r = frame_recv_design(f);
	const double *ar;
	const size_t *ix;
	size_t nz;
	size_t dimr = design_tvar_dim(r);
	double *yr = y;
	
	scale_result(dimr, beta, yr);
	
	design_tvars_get(r, &ar, &ix, &nz);
	for (; nz != 0; ar += dimr, ix++, nz--) {
		blas_daxpy(dimr, alpha * x[*ix], ar, 1, yr, 1);
	}
	
	
	const struct design2 *d = frame_dyad_design(f);
	const double *ad;
	const size_t *j;
	size_t dimd = design2_tvar_dim(d);
	double *yd = yr + dimr;

	
	scale_result(dimd, beta, yd);

	design2_tvars_get(d, isend, &ad, &j, &nz);
	for (; nz != 0; ad += dimd, j++, nz--) {
		blas_daxpy(dimd, alpha * x[*j], ad, 1, yd, 1);
	}
}


void recv_frame_axpy1(double alpha, const struct recv_frame *rf, size_t isend,
		      size_t jrecv, double *y)
{
	const struct recv_fmla *fmla = recv_frame_fmla(rf);
	const struct frame *f = recv_fmla_frame(fmla);
	
	const struct design *r = frame_recv_design(f);
	const double *xr = design_tvars(r, jrecv);
	size_t dimr = design_tvar_dim(r);	
	double *yr = y;
	
	if (xr) {
		blas_daxpy(dimr, alpha, xr, 1, yr, 1);
	}
	
	const struct design2 *d = frame_dyad_design(f);
	const double *xd = design2_tvars(d, isend, jrecv);
	size_t dimd = design2_trait_dim(d);
	double *yd = yr + dimr;

	if (xd) {
		blas_daxpy(dimd, alpha, xd, 1, yd, 1);
	}

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

