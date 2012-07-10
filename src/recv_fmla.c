#include "port.h"
#include "xalloc.h"
#include <stdlib.h>
#include "recv_fmla.h"


void recv_fmla_init(struct recv_fmla *fmla, const struct frame *f)
{
	const struct design *r = frame_recv_design(f);
	const struct design *d = frame_dyad_design(f);
	fmla->frame = f;
	fmla->trait_dim = design_trait_dim(r) + design_trait_dim(d);
	fmla->tvar_dim = design_tvar_dim(r) + design_tvar_dim(d);
}


void recv_fmla_deinit(struct recv_fmla *fmla)
{
	(void)fmla; // unused
}


void recv_coefs_init(struct recv_coefs *c, const struct recv_fmla *fmla)
{
	c->traits = xmalloc(recv_fmla_trait_dim(fmla) * sizeof(*c->traits));
	c->tvars = xmalloc(recv_fmla_tvar_dim(fmla) * sizeof(*c->tvars));
}


void recv_coefs_deinit(struct recv_coefs *c)
{
	free(c->tvars);
	free(c->traits);
}


void recv_frame_init(struct recv_frame *rf, const struct recv_fmla *fmla)
{
	rf->fmla = fmla;
}


void recv_frame_deinit(struct recv_frame *rf)
{
	(void)rf; // unused
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
	size_t nrecv = frame_recv_count(f);
	
	const struct design *r = frame_recv_design(f);
	const struct dmatrix *ar = design_traits(r);
	size_t dimr = design_trait_dim(r);	
	const double *xr = x;

	blas_dgemv(BLAS_NOTRANS, nrecv, dimr, alpha, ar, xr, 1, beta, y, 1);
	
	const struct design *d = frame_dyad_design(f);
	const struct dmatrix *ad = design_traits(d);
	size_t dimd = design_trait_dim(d);
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
	
	const struct design *d = frame_dyad_design(f);
	const struct dmatrix *ad = design_traits(d);
	size_t dimd = design_trait_dim(d);
	double *yd = yr + dimr;
	
	struct dmatrix adsub = *ad;
	adsub.data = MATRIX_PTR(ad, isend * nrecv, 0);
	blas_dgemv(BLAS_TRANS, nrecv, dimd, alpha, &adsub, x, 1, 1.0, yd, 1);
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
	
	const struct design *d = frame_dyad_design(f);
	const struct dmatrix *ad = design_traits(d);
	size_t dimd = design_trait_dim(d);
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
	const size_t *ix;
	size_t nz;
	size_t dimr = design_tvar_dim(r);
	const double *xr = x;
	
	design_tvars_get(r, &ar, &ix, &nz);
	for (; nz != 0; ar += dimr, ix++, nz--) {
		y[*ix] += alpha * blas_ddot(dimr, ar, 1, xr, 1);
	}

	
	const struct design *d = frame_dyad_design(f);
	const double *ad0, *ad1;
	const size_t *ix0, *ix1;
	size_t dimd = design_tvar_dim(d);
	const double *xd = xr + dimr;
	size_t off = isend * nrecv;
	
	design_tvar_get_lb(d, off, &ad0, &ix0);
	design_tvar_get_lb(d, off + nrecv, &ad1, &ix1);
	
	for (; ix0 != ix1; ad0 += dimd, ix0++) {
		y[*ix0 - off] += alpha * blas_ddot(dimd, ad0, 1, xd, 1);
	}
}


void recv_frame_tmul1(double alpha, const struct recv_frame *rf, size_t isend,
		      const double *x, double beta, double *y)
{
	const struct recv_fmla *fmla = recv_frame_fmla(rf);
	const struct frame *f = recv_fmla_frame(fmla);
	size_t nrecv = frame_recv_count(f);
	
	const struct design *r = frame_recv_design(f);
	const double *ar;
	const size_t *ix;
	size_t nz;
	size_t dimr = design_tvar_dim(r);
	double *yr = y;
	
	design_tvars_get(r, &ar, &ix, &nz);
	for (; nz != 0; ar += dimr, ix++, nz--) {
		blas_daxpy(dimr, alpha * x[*ix], ar, 1, yr, 1);
	}
	
	
	const struct design *d = frame_dyad_design(f);
	const double *ad0, *ad1;
	const size_t *ix0, *ix1;
	size_t dimd = design_tvar_dim(d);
	double *yd = yr + dimr;
	size_t off = isend * nrecv;
	
	design_tvar_get_lb(d, off, &ad0, &ix0);
	design_tvar_get_lb(d, off + nrecv, &ad1, &ix1);
	
	for (; ix0 != ix1; ad0 += dimd, ix0++) {
		blas_daxpy(dimd, alpha * x[*ix0 - off], ad0, 1, yd, 1);
	}
}


void recv_frame_axpy1(double alpha, const struct recv_frame *rf, size_t isend,
		      size_t jrecv, double *y)
{
	const struct recv_fmla *fmla = recv_frame_fmla(rf);
	const struct frame *f = recv_fmla_frame(fmla);
	size_t nrecv = frame_recv_count(f);
	
	const struct design *r = frame_recv_design(f);
	const double *xr = design_tvars(r, jrecv);
	size_t dimr = design_tvar_dim(r);	
	double *yr = y;
	
	if (xr) {
		blas_daxpy(dimr, alpha, xr, 1, yr, 1);
	}
	
	const struct design *d = frame_dyad_design(f);
	const double *xd = design_tvars(d, isend * nrecv + jrecv);
	size_t dimd = design_trait_dim(d);
	double *yd = yr + dimr;

	if (xd) {
		blas_daxpy(dimd, alpha, xd, 1, yd, 1);
	}

}




#if 0

void frame_recv_mul(double alpha, enum blas_trans trans,
		    const struct frame *f, size_t isend,
		    const double *x, double beta, double *y)
{
	assert(isend < frame_send_count(f));
	
	const struct design *d = frame_recv_design(f);
	size_t off = design_dvars_index(d);
	
	design_mul0(alpha, trans, d, x, beta, y);
	
	if (trans == BLAS_NOTRANS) {
		frame_recv_dmul(alpha, trans, f, isend, x + off, 1.0, y);
	} else {
		frame_recv_dmul(alpha, trans, f, isend, x, 1.0, y + off);
	}
}

void frame_recv_muls(double alpha, enum blas_trans trans,
		     const struct frame *f, size_t isend,
		     const double *x, const struct vpattern *pat,
		     double beta, double *y)
{
	assert(isend < frame_send_count(f));
	
	const struct design *d = frame_recv_design(f);
	size_t off = design_dvars_index(d);
	size_t dim = design_dvars_dim(d);
	
	design_muls0(alpha, trans, d, x, pat, beta, y);
	
	if (trans == BLAS_NOTRANS) {
		ptrdiff_t jx0 = vpattern_find(pat, off);
		ptrdiff_t jx1 = vpattern_find(pat, off + dim);
		size_t jz0 = (jx0 < 0 ? ~jx0 : jx0);
		size_t jz1 = (jx1 < 0 ? ~jx1 : jx1);
		size_t jz;
		
		const struct recv_frame *rf = recv_frames_item(f, isend);
		size_t iz, nz = rf->active.nz;
		
		for (iz = 0; iz < nz; iz++) {
			size_t jrecv = rf->active.indx[iz];
			const double *dx = rf->dx + iz * dim;
			
			double dot = 0;
			for (jz = jz0; jz < jz1; jz++) {
				size_t i = pat->indx[jz] - off;
				dot += x[jz] * dx[i];
			}
			
			y[jrecv] += alpha * dot;
		}
	} else {
		frame_recv_dmuls(alpha, trans, f, isend, x, pat, 1.0, y + off);
	}
}

void frame_recv_dmul(double alpha, enum blas_trans trans,
		     const struct frame *f, size_t isend,
		     const double *x, double beta, double *y)
{
	assert(isend < frame_send_count(f));
	
	const struct design *d = frame_recv_design(f);
	size_t n = design_count(d);
	size_t dim = design_dvars_dim(d);
	size_t ny = (trans == BLAS_NOTRANS ? n : dim);
	
	/* y := beta y */
	if (beta == 0.0) {
		memset(y, 0, ny * sizeof(y[0]));
	} else if (beta != 1.0) {
		blas_dscal(ny, beta, y, 1);
	}
	
	if (dim == 0)
		return;
	
	const struct recv_frame *rf = recv_frames_item(f, isend);
	size_t iz, nz = rf->active.nz;
	size_t jrecv;
	const double *dx;
	
	if (trans == BLAS_NOTRANS) {
		for (iz = 0; iz < nz; iz++) {
			jrecv = rf->active.indx[iz];
			dx = rf->dx + iz * dim;
			
			double dot = blas_ddot(dim, dx, 1, x, 1);
			y[jrecv] += alpha * dot;
		}
	} else {
		for (iz = 0; iz < nz; iz++) {
			jrecv = rf->active.indx[iz];
			dx = rf->dx + iz * dim;
			
			if (x[jrecv] == 0.0)
				continue;
			
			blas_daxpy(dim, alpha * x[jrecv], dx, 1, y, 1);
		}
	}
}

void frame_recv_dmuls(double alpha, enum blas_trans trans,
		      const struct frame *f, size_t isend,
		      const double *x, const struct vpattern *pat,
		      double beta, double *y)
{
	assert(isend < frame_send_count(f));
	
	const struct design *d = frame_recv_design(f);
	size_t n = design_count(d);
	size_t dim = design_dvars_dim(d);
	size_t ny = (trans == BLAS_NOTRANS ? n : dim);
	
	/* y := beta y */
	if (beta == 0.0) {
		memset(y, 0, ny * sizeof(y[0]));
	} else if (beta != 1.0) {
		blas_dscal(ny, beta, y, 1);
	}
	
	if (dim == 0)
		return;
	
	const struct recv_frame *rf = recv_frames_item(f, isend);
	size_t iz, nz = rf->active.nz;
	size_t jrecv;
	const double *dx;
	
	if (trans == BLAS_NOTRANS) {
		
		for (iz = 0; iz < nz; iz++) {
			jrecv = rf->active.indx[iz];
			dx = rf->dx + iz * dim;
			
			double dot = sblas_ddoti(x, pat, dx);
			y[jrecv] += alpha * dot;
		}
	} else {
		size_t jz, mz = pat->nz;
		for (jz = 0; jz < mz; jz++) {
			size_t jrecv = pat->indx[jz];
			ptrdiff_t ix = vpattern_find(&rf->active, jrecv);
			if (ix < 0)
				continue;
			
			dx = rf->dx + ix * dim;
			/* y := y + alpha * x[j] * dx[j] */
			double jscale = alpha * x[jz];
			blas_daxpy(dim, jscale, dx, 1, y, 1);
		}
	}
}

#endif
