#include "port.h"
#include <float.h>		// DBL_MANT_DIG
#include <stdint.h>		// SIZE_MAX
#include <stdio.h>
#include <stdlib.h>		// free
#include <string.h>		// memcpy, memset
#include "blas.h"		// blas_gemv
#include "coreutil.h"
#include "ieee754.h"		// double_eqrel
#include "lapack.h"
#include "sblas.h"
#include "xalloc.h"		// xcalloc
#include "mlogit.h"

#define CHECK(x) \
	do { \
		fail = !(x); \
		assert(!fail); \
		if (fail) \
			goto out; \
	} while(0)

#define CHECK_APPROX(x, y) \
	CHECK(double_eqrel((x), (y)) >= DBL_MANT_DIG / 2 \
	      || (fabs(x) <= SQRT_DBL_EPSILON && fabs(y) <= SQRT_DBL_EPSILON))

#define F77_COV_UPLO (MLOGIT_COV_UPLO == BLAS_LOWER ? BLAS_UPPER : BLAS_LOWER)

#define BLOCK_SIZE	64


static void clear(struct mlogit *m);
static double get_deta(struct mlogit *m, const double *dx, const size_t *jdx,
		       size_t ndx);
static void increment_x(struct mlogit *m, size_t i, const size_t *jdx,
			const double *dx, size_t ndx);
static void recompute_all(struct mlogit *m);
static void recompute_dist(struct mlogit *m);
static void recompute_mean(struct mlogit *m);
static void recompute_cov(struct mlogit *m);
static void update_cov_old(struct mlogit *m, const double *dx, const double *xresid,
			   const double w, const double dw, const double dpsi);
static void update_mean_old(struct mlogit *m, const double *dx,
			    const double *xresid, const double w, const double dw);
static void update(struct mlogit *m);
static void update_x(struct mlogit *m);
static void update_dist(struct mlogit *m);
static void update_mean(struct mlogit *m);
static void update_cov(struct mlogit *m);
static void update_cov1(struct mlogit *m);
static void update_cov2(struct mlogit *m);
static void update_cov3(struct mlogit *m);

static void copy_packed(double *dst, const double *src, size_t dim, enum blas_uplo uplo);
static void grow_ind_array(struct mlogit *m, size_t delta);
static size_t find_ind(const struct mlogit *m, size_t i);
static size_t search_ind(struct mlogit *m, size_t i);
static void raw_inc_x(struct mlogit *m, size_t i, const size_t *jdx,
		      const double *dx, size_t ndx);



void mlogit_work_init(struct mlogit_work *work, size_t ncat, size_t dim)
{
	size_t cov_dim = dim * (dim + 1) / 2;
	size_t cov_dim_full = dim * dim;


	work->cov_diff = xmalloc(cov_dim * sizeof(*work->cov_diff));
	work->cov_diff_full = xmalloc(cov_dim_full * sizeof(*work->cov_diff_full));
	work->cat_buf = xmalloc(ncat * sizeof(*work->cat_buf));
	work->dim_buf1 = xmalloc(dim * sizeof(*work->dim_buf1));
	work->dim_buf2 = xmalloc(dim * sizeof(*work->dim_buf2));

	work->x0 = xmalloc(ncat * dim * sizeof(*work->x0));
	work->one = xmalloc(BLOCK_SIZE * sizeof(*work->one));
	work->w0 = xmalloc(ncat * sizeof(*work->w0));
	work->dw = xmalloc(ncat * sizeof(*work->dw));
	work->mean0 = xmalloc(dim * sizeof(*work->mean0));
	work->dmean = xmalloc(dim * sizeof(*work->dmean));
	work->xbuf1 = xmalloc(BLOCK_SIZE * dim * sizeof(*work->xbuf1));
	work->xbuf2 = xmalloc(BLOCK_SIZE * dim * sizeof(*work->xbuf2));


	size_t k, nk = BLOCK_SIZE;
	for (k = 0; k < nk; k++) {
		work->one[k] = 1.0;
	}
}

void mlogit_work_deinit(struct mlogit_work *work)
{
	free(work->xbuf2);
	free(work->xbuf1);
	free(work->dmean);
	free(work->mean0);
	free(work->dw);
	free(work->w0);
	free(work->one);
	free(work->x0);
	free(work->dim_buf2);
	free(work->dim_buf1);
	free(work->cat_buf);
	free(work->cov_diff_full);
	free(work->cov_diff);
}


void mlogit_init(struct mlogit *m, size_t ncat, size_t dim, struct mlogit_work *work)
{
	assert(ncat == 0 || dim <= SIZE_MAX / sizeof(double) / ncat);
	assert(dim <= SIZE_MAX / sizeof(double) / (dim + 1));

	size_t cov_dim = dim * (dim + 1) / 2;

	catdist_init(&m->dist_, ncat);
	m->x_ = xmalloc(ncat * dim * sizeof(*m->x_));
	m->beta = xmalloc(dim * sizeof(*m->beta));
	m->offset = xmalloc(ncat * sizeof(*m->offset));
	m->mean_ = xmalloc(dim * sizeof(*m->mean_));
	m->cov_ = xmalloc(cov_dim * sizeof(*m->cov_));
	m->log_cov_scale_ = 0;
	m->dim = dim;
	m->moments = 2;

	m->ind = NULL;
	m->dx = NULL;
	m->nz = 0;
	m->nzmax = 0;

	if (work) {
		m->free_work = 0;
	} else {
		work = xmalloc(sizeof(*work));
		mlogit_work_init(work, ncat, dim);
		m->free_work = 1;
	}
	m->work = work;

	clear(m);
}


void mlogit_deinit(struct mlogit *m)
{
	if (m->free_work) {
		mlogit_work_deinit(m->work);
		free(m->work);
	}
	free(m->dx);
	free(m->ind);
	free(m->cov_);
	free(m->mean_);
	free(m->offset);
	free(m->beta);
	free(m->x_);
	catdist_deinit(&m->dist_);
}


void clear(struct mlogit *m)
{
	size_t ncat = mlogit_ncat(m);
	size_t dim = m->dim;
	size_t cov_dim = dim * (dim + 1) / 2;

	catdist_set_all_eta(&m->dist_, NULL);
	memset(m->x_, 0, ncat * dim * sizeof(*m->x_));
	memset(m->beta, 0, dim * sizeof(*m->beta));
	memset(m->offset, 0, ncat * sizeof(*m->offset));
	memset(m->mean_, 0, dim * sizeof(*m->mean_));
	memset(m->cov_, 0, cov_dim * sizeof(*m->cov_));
	m->log_cov_scale_ = 0;
	m->mean_err = 0.0;
	m->cov_err = 0.0;
	m->log_cov_scale_err = 0.0;
	m->nz = 0;
}


void mlogit_set_coefs(struct mlogit *m, const double *beta)
{
	size_t len = mlogit_dim(m) * sizeof(*m->beta);
	if (beta) {
		memcpy(m->beta, beta, len);
	} else {
		memset(m->beta, 0, len);
	}

	recompute_all(m);
}

void mlogit_set_all_offset(struct mlogit *m, const double *offset)
{
	size_t len = mlogit_ncat(m) * sizeof(*m->offset);
	if (offset) {
		memcpy(m->offset, offset, len);
	} else {
		memset(m->offset, 0, len);
	}

	recompute_all(m);
}

void mlogit_set_all_x(struct mlogit *m, const double *x)
{
	size_t len = mlogit_ncat(m) * mlogit_dim(m) * sizeof(*m->x_);
	if (x) {
		memcpy(m->x_, x, len);
	} else {
		memset(m->x_, 0, len);
	}
	m->nz = 0;

	recompute_all(m);
}

void mlogit_set_offset(struct mlogit *m, size_t i, double offset)
{
	assert(i < mlogit_ncat(m));
	assert(offset < INFINITY);

	const size_t dim = mlogit_dim(m);
	const double psi = catdist_psi(&m->dist_);
	const double eta = catdist_eta(&m->dist_, i);
	const double deta = offset - m->offset[i];
	const double eta1 = eta + deta;

	// offset[i] := offset
	m->offset[i] = offset;

	// update eta, psi
	catdist_set_eta(&m->dist_, i, eta1);

	if (m->moments < 1)
		return;

	// compute weight changes
	const double psi1 = catdist_psi(&m->dist_);
	const double dpsi = psi1 - psi;
	const double w1 = exp(eta1 - psi1);
	const double w = exp(eta - psi1);
	const double dw = w1 - w;

	// xresid := x1 - mean
	double *xresid = m->work->dim_buf1;
	blas_dcopy(dim, m->x_ + i * dim, 1, xresid, 1);
	blas_daxpy(dim, -1.0, m->mean_, 1, xresid, 1);

	update_mean_old(m, NULL, xresid, w, dw);
	update_cov_old(m, NULL, xresid, w, dw, dpsi);
}


void mlogit_inc_x(struct mlogit *m, size_t i, const size_t *jdx,
		  const double *dx, size_t ndx)
{
	assert(i < mlogit_ncat(m));
	assert(dx || ndx == 0);

	size_t dim = mlogit_dim(m);
	size_t iz = search_ind(m, i);
	double *dst = m->dx + iz * dim;

	if (jdx) {
#ifndef NDEBUG
		size_t ncat = mlogit_ncat(m);
		size_t k;
		for (k = 0; k < ndx; k++) {
			assert(jdx[k] < ncat);
		}
#endif
		sblas_daxpyi(ndx, 1.0, dx, jdx, dst);
	} else if (ndx) {
		assert(ndx == mlogit_dim(m));
		blas_daxpy(dim, 1.0, dx, 1, dst, 1);
	}
}


void update(struct mlogit *m)
{
	if (m->nz) {
		update_x(m);
		m->nz = 0;
		recompute_all(m);
		//update_dist(m);
		//update_mean(m);
		//update_cov(m);

	}
}


void update_x(struct mlogit *m)
{
	if (!mlogit_dim(m))
		return;

	const size_t dim = mlogit_dim(m);
	double *x0 = m->work->x0;

	// copy x0; compute x
	size_t iz, nz = m->nz;
	for (iz = 0; iz < nz; iz++) {
		size_t i = m->ind[iz];
		double *x_i = m->x_ + i * dim;
		double *x0_i = x0 + iz * dim;
		const double *dx_i = m->dx + iz * dim;

		// copy x0;
		memcpy(x0_i, x_i, dim * sizeof(double));

		// x[i] := x0[i] + dx[i]
		blas_daxpy(dim, 1.0, dx_i, 1, x_i, 1);
	}
}


void update_dist(struct mlogit *m)
{
	if (!mlogit_dim(m))
		return;

	if (m->nz >= 0 * mlogit_ncat(m) / 2) {
		recompute_mean(m);
		return;
	}

	const size_t dim = mlogit_dim(m);
	const double psi = catdist_psi(&m->dist_);
	double *eta = m->work->w0;
	double *eta1 = m->work->dw;

	// copy old values of eta
	size_t iz, nz = m->nz;
	for (iz = 0; iz < nz; iz++) {
		size_t i = m->ind[iz];
		eta[iz] = catdist_eta(&m->dist_, i);
	}

	// compute new values of eta
	blas_dcopy(nz, eta, 1, eta1, 1);
	blas_dgemv(BLAS_TRANS, dim, nz, 1.0, m->dx, dim, m->beta, 1, 1.0, eta1, 1);

	// update dist
	for (iz = 0; iz < nz; iz++) {
		size_t i = m->ind[iz];
		catdist_set_eta(&m->dist_, i, eta1[iz]);
	}

	// compute new value of psi
	const double psi1 = catdist_psi(&m->dist_);
	m->work->dpsi = psi1 - psi;

	// eta (w0) := exp(eta - psi1);
	for (iz = 0; iz < nz; iz++) {
		eta[iz] = exp(eta[iz] - psi1);

	}

	// eta1 := exp(eta1 - psi1)
	for (iz = 0; iz < nz; iz++) {
		eta1[iz] = exp(eta1[iz] - psi1);
	}

	// eta1 (dw) := eta1 - eta (w0)
	blas_daxpy(nz, -1.0, eta, 1, eta1, 1);
}


void update_mean(struct mlogit *m)
{
	if (mlogit_moments(m) < 1 || !mlogit_dim(m))
		return;

	if (m->nz >= 0 * mlogit_ncat(m) / 2) {
		recompute_mean(m);
		return;
	}

	size_t dim = mlogit_dim(m);

	const double *one = m->work->one;
	const double *mean0 = m->mean_;
	const double *dw = m->work->dw;
	const double *w0 = m->work->w0;
	double *dmean = m->work->dmean;
	double *xresid = m->work->xbuf1;

	memset(dmean, 0, dim * sizeof(*dmean));

	size_t iz, nz = m->nz;
	double *xresid_i = xresid;
	size_t nk = 0;

	for (iz = 0; iz < nz; iz++) {
		size_t i = m->ind[iz];
		double dw_i = dw[iz];
		double w0_i = w0[iz];
		const double *x_i = m->x_ + i * dim;
		const double *dx_i = m->dx + iz * dim;

		// xresid[i] := dw[i] * (x1[i] - mu0)
		if (dw_i != 0) {
			blas_dcopy(dim, x_i, 1, xresid_i, 1);
			blas_daxpy(dim, -1.0, mean0, 1, xresid_i, 1);
			blas_dscal(dim, dw_i, xresid_i, 1);
		} else {
			memset(xresid_i, 0, dim * sizeof(*xresid_i));
		}

		// xresid[i] := dw[i] * (x1[i] - mu0) + w0[i] * dx[i]
		if (w0_i != 0) {
			blas_daxpy(dim, w0_i, dx_i, 1, xresid_i, 1);
		}

		m->mean_err += 8 * (fabs(dw_i) + w0_i);

		xresid_i += dim;
		nk++;

		if (nk == BLOCK_SIZE || iz + 1 == nz) {
			blas_dgemv(BLAS_NOTRANS, dim, nk, 1.0, xresid, dim, one, 1, 1.0, dmean, 1);
			xresid_i = xresid;
			nk = 0;
		}
	}

	blas_dcopy(dim, mean0, 1, m->work->mean0, 1);

	m->mean_err += 1.0;
	blas_daxpy(dim, 1.0, dmean, 1, m->mean_, 1);

	const double tol = 1.0 / ROOT4_DBL_EPSILON;
	if (!(m->mean_err < tol)) {
		recompute_mean(m);
		blas_dcopy(dim, m->mean_, 1, dmean, 1);
		blas_daxpy(dim, -1.0, m->work->mean0, 1, dmean, 1);
	}
}


void update_cov(struct mlogit *m)
{
	if (mlogit_moments(m) < 2 || !mlogit_dim(m))
		return;

	if (m->nz >= 0 * mlogit_ncat(m) / 2) {
		if (!(m->mean_err == 0.0))
			recompute_mean(m);
		recompute_cov(m);
		return;
	}

	size_t dim = mlogit_dim(m);
	size_t cov_dim = dim * (dim + 1) / 2;
	const double dpsi = m->work->dpsi;
	const double log_scale = m->log_cov_scale_;
	const double log_scale1 = log_scale + dpsi;
	//const double W = exp(log_scale);
	const double W1 = exp(log_scale1);
	double *dcov_full = m->work->cov_diff_full;
	double *dcov = m->work->cov_diff;

	memset(dcov_full, 0, dim * dim * sizeof(*dcov_full));

	update_cov1(m);
	update_cov2(m);
	update_cov3(m);

	copy_packed(dcov, dcov_full, dim, MLOGIT_COV_UPLO);
	blas_daxpy(cov_dim, W1, dcov, 1, m->cov_, 1);
	m->log_cov_scale_ = log_scale1;

	m->cov_err += 1.0;
	m->log_cov_scale_err += 0.5 * DBL_EPSILON * (fabs(dpsi) + fabs(log_scale1));

	const double tol = 1.0 / ROOT4_DBL_EPSILON;
	double err = (m->cov_err) / exp(log_scale1 - m->log_cov_scale_err);
	
	if (!(err < tol)) {
		blas_dcopy(cov_dim, m->cov_, 1, dcov, 1);

		if (!(m->mean_err == 0.0))
			recompute_mean(m);
		recompute_cov(m);

		blas_dscal(cov_dim, -exp(m->log_cov_scale_ - log_scale1), dcov, 1);
		blas_daxpy(cov_dim, 1.0, m->cov_, 1, dcov, 1);
	}
}


void update_cov1(struct mlogit *m)
{
	size_t dim = mlogit_dim(m);
	const double *x0 = m->work->x0;
	const double *x1 = m->x_;
	const double *w0 = m->work->w0;
	const double *dw = m->work->dw;
	const double *mean0 = m->work->mean0;
	const double *dmean = m->work->dmean;
	const double log_scale = m->log_cov_scale_;
	//const double W0 = exp(log_scale);
	const double dpsi = m->work->dpsi;
	const double W1 = exp(log_scale + dpsi);
	double *dcov = m->work->cov_diff_full;
	double *diff = m->work->xbuf1;

	size_t iz, nz = m->nz;
	size_t nk = 0;
	double *diff_i = diff;

	assert(dim > 0);

	// handle positive dw
	for (iz = 0; iz < nz; iz++) {
		const size_t i = m->ind[iz];
		const double *x0_i = x0 + iz * dim;
		const double *x1_i = x1 + i * dim;
		const double w0_i = w0[iz];
		const double w1_i = catdist_prob(&m->dist_, i);
		const double dw_i = dw[iz];

		if (dw_i > 0) {
			if (w0_i <= w1_i) {
				blas_dcopy(dim, x1_i, 1, diff_i, 1);
			} else {
				blas_dcopy(dim, x0_i, 1, diff_i, 1);
			}
			blas_daxpy(dim, -1.0, mean0, 1, diff_i, 1);
			blas_dscal(dim, sqrt(dw_i), diff_i, 1);
			m->cov_err += 64 * W1 * dw_i;

			diff_i += dim;
			nk++;
		}

		if (nk == BLOCK_SIZE || iz + 1 == nz) {
			blas_dsyrk(F77_COV_UPLO, BLAS_NOTRANS, dim, nk, 1.0, diff, dim, 1.0, dcov, dim);

			diff_i = diff;
			nk = 0;
		}
	}


	assert(nk == 0);
	assert(diff_i == diff);

	// handle dmean and negative dw
	blas_dcopy(dim, dmean, 1, diff_i, 1);
	m->cov_err += W1;

	diff_i += dim;
	nk++;

	for (iz = 0; iz < nz; iz++) {
		const size_t i = m->ind[iz];
		const double *x0_i = x0 + iz * dim;
		const double *x1_i = x1 + i * dim;
		const double w0_i = w0[iz];
		const double w1_i = catdist_prob(&m->dist_, i);
		const double dw_i = dw[iz];

		if (dw_i < 0) {
			if (w0_i <= w1_i) {
				blas_dcopy(dim, x1_i, 1, diff_i, 1);
			} else {
				blas_dcopy(dim, x0_i, 1, diff_i, 1);
			}
			blas_daxpy(dim, -1.0, mean0, 1, diff_i, 1);
			blas_dscal(dim, sqrt(-dw_i), diff_i, 1);
			m->cov_err += 64 * W1 * (-dw_i);

			diff_i += dim;
			nk++;
		}

		if (nk == BLOCK_SIZE || iz + 1 == nz) {
			blas_dsyrk(F77_COV_UPLO, BLAS_NOTRANS, dim, nk, -1.0, diff, dim, 1.0, dcov, dim);

			diff_i = diff;
			nk = 0;
		}
	}
	// note: nz = 0 implies dmean = 0; no need to call dsyrk
}


void update_cov2(struct mlogit *m)
{
	size_t dim = mlogit_dim(m);
	const double *x0 = m->work->x0;
	const double *x1 = m->x_;
	const double *dx = m->dx;
	const double *w0 = m->work->w0;
	const double *mean0 = m->work->mean0;
	const double log_scale = m->log_cov_scale_;
	//const double W0 = exp(log_scale);
	const double dpsi = m->work->dpsi;
	const double W1 = exp(log_scale + dpsi);
	double *dcov = m->work->cov_diff_full;
	double *diff = m->work->xbuf1;

	size_t iz, nz = m->nz;
	size_t nk = 0;
	double *diff_i = diff;

	assert(dim > 0);

	for (iz = 0; iz < nz; iz++) {
		const size_t i = m->ind[iz];
		const double *x0_i = x0 + iz * dim;
		const double *x1_i = x1 + i * dim;
		const double w0_i = w0[iz];
		const double w1_i = catdist_prob(&m->dist_, i);

		if (w0_i <= w1_i) {
			blas_dcopy(dim, x1_i, 1, diff_i, 1);
			blas_daxpy(dim, -1.0, mean0, 1, diff_i, 1);
			blas_dscal(dim, w0_i, diff_i, 1);
			m->cov_err += 2 * 64 * W1 * w0_i;
		} else {
			blas_dcopy(dim, x0_i, 1, diff_i, 1);
			blas_daxpy(dim, -1.0, mean0, 1, diff_i, 1);
			blas_dscal(dim, w1_i, diff_i, 1);
			m->cov_err += 2 * 64 * W1 * w1_i;
		}

		diff_i += dim;
		nk++;

		if (nk == BLOCK_SIZE || iz + 1 == nz) {
			blas_dsyr2k(F77_COV_UPLO, BLAS_NOTRANS, dim, nk, 1.0, diff, dim, dx, dim, 1.0, dcov, dim);
			dx += nk * dim;

			diff_i = diff;
			nk = 0;
		}
	}
}


void update_cov3(struct mlogit *m)
{
	size_t dim = mlogit_dim(m);
	const double *dx = m->dx;
	const double *w0 = m->work->w0;
	const double log_scale = m->log_cov_scale_;
	//const double W0 = exp(log_scale);
	const double dpsi = m->work->dpsi;
	const double W1 = exp(log_scale + dpsi);
	double *dcov = m->work->cov_diff_full;
	double *diff = m->work->xbuf1;

	size_t iz, nz = m->nz;
	size_t nk = 0;
	double *diff_i = diff;

	assert(dim > 0);

	for (iz = 0; iz < nz; iz++) {
		const size_t i = m->ind[iz];
		const double *dx_i = dx + iz * dim;
		const double w0_i = w0[iz];
		const double w1_i = catdist_prob(&m->dist_, i);

		if (w0_i <= w1_i) {
			blas_dcopy(dim, dx_i, 1, diff_i, 1);
			blas_dscal(dim, sqrt(w0_i), diff_i, 1);
			m->cov_err += 64 * W1 * w0_i;

			diff_i += dim;
			nk++;
		}

		if (nk == BLOCK_SIZE || iz + 1 == nz) {
			blas_dsyrk(F77_COV_UPLO, BLAS_NOTRANS, dim, nk, -1.0, diff, dim, 1.0, dcov, dim);

			diff_i = diff;
			nk = 0;
		}
	}

	assert(nk == 0);
	assert(diff_i == diff);

	for (iz = 0; iz < nz; iz++) {
		const size_t i = m->ind[iz];
		const double *dx_i = dx + iz * dim;
		const double w0_i = w0[iz];
		const double w1_i = catdist_prob(&m->dist_, i);

		if (!(w0_i <= w1_i)) {
			blas_dcopy(dim, dx_i, 1, diff_i, 1);
			blas_dscal(dim, sqrt(w1_i), diff_i, 1);
			m->cov_err += 64 * W1 * w1_i;

			diff_i += dim;
			nk++;
		}

		if (nk == BLOCK_SIZE || iz + 1 == nz) {
			blas_dsyrk(F77_COV_UPLO, BLAS_NOTRANS, dim, nk, +1.0, diff, dim, 1.0, dcov, dim);

			diff_i = diff;
			nk = 0;
		}
	}
}


void raw_inc_x(struct mlogit *m, size_t i, const size_t *jdx,
		  const double *dx, size_t ndx)
{
	assert(i < mlogit_ncat(m));
	assert(dx || ndx == 0);
	assert(jdx || ndx == 0 || ndx == mlogit_dim(m));

	size_t dim = mlogit_dim(m);

	if (ndx == 0 || dim == 0)
		return;

	const double psi = catdist_psi(&m->dist_);
	const double eta = catdist_eta(&m->dist_, i);
	const double deta = get_deta(m, dx, jdx, ndx);
	const double eta1 = eta + deta;

	// x := x + dx
	increment_x(m, i, jdx, dx, ndx);

	// update eta, psi
	catdist_set_eta(&m->dist_, i, eta1);

	if (m->moments < 1)
		return;

	// compute weight changes
	const double psi1 = catdist_psi(&m->dist_);
	const double dpsi = psi1 - psi;
	const double w1 = exp(eta1 - psi1);
	const double w = exp(eta - psi1);
	const double dw = w1 - w;

	// xresid := x1 - mean
	double *xresid = m->work->dim_buf1;
	blas_dcopy(dim, m->x_ + i * dim, 1, xresid, 1);
	blas_daxpy(dim, -1.0, m->mean_, 1, xresid, 1);

	// set ddx (dense dx)
	double *ddx = jdx ? m->work->dim_buf2 : (double *)dx;

	if (jdx) {
		memset(ddx, 0, dim * sizeof(*ddx));
		sblas_dsctr(ndx, dx, jdx, ddx);
	}

	update_mean_old(m, ddx, xresid, w, dw);
	update_cov_old(m, ddx, xresid, w, dw, dpsi);
}

static void update_mean_old(struct mlogit *m, const double *dx,
			    const double *xresid, const double w, const double dw)
{
	if (m->moments < 1)
		return;

	size_t dim = mlogit_dim(m);
	double *mean_diff = m->work->dmean;
	const double tol = 1.0 / ROOT4_DBL_EPSILON;

	if (dx) {
		m->mean_err += 1 + 8 * (fabs(dw) + w);	// approximate relative error from update
	} else {
		m->mean_err += 1 + 8 * fabs(dw);
	}

	if (!(m->mean_err < tol)) {
		recompute_mean(m);
		return;
	}
	// mean_diff := dw * (x1 - mean)
	blas_dcopy(dim, xresid, 1, mean_diff, 1);
	blas_dscal(dim, dw, mean_diff, 1);

	if (dx) {
		// mean_diff += w * dx
		blas_daxpy(dim, w, dx, 1, mean_diff, 1);
	}
	// mean += mean_diff
	blas_daxpy(dim, 1.0, mean_diff, 1, m->mean_, 1);
}

static void copy_packed(double *dst, const double *src, size_t dim, enum blas_uplo uplo)
{

	size_t i;

	if (uplo == BLAS_UPPER) {
		for (i = 0; i < dim; i++) {
			size_t rowlen = dim - i;
			memcpy(dst, src + i, rowlen * sizeof(*dst));
			dst += rowlen;
			src += dim;
		}
	} else {
		for (i = 0; i < dim; i++) {
			size_t rowlen = i + 1;
			memcpy(dst, src, rowlen * sizeof(*dst));
			dst += rowlen;
			src += dim;
		}
	}
}

static void update_cov_old(struct mlogit *m, const double *dx, const double *xresid,
			   const double w, const double dw, const double dpsi)
{
	if (m->moments < 2 || mlogit_dim(m) == 0)
		return;

	const size_t dim = mlogit_dim(m);
	const size_t cov_dim = dim * (dim + 1) / 2;
	const size_t cov_dim_full = dim * dim;
	double *cov_diff_full = m->work->cov_diff_full;
	double *cov_diff = m->work->cov_diff;
	const double log_scale = m->log_cov_scale_;
	const double log_scale1 = log_scale + dpsi;
	const double W = exp(log_scale);
	const double W1 = exp(log_scale1);
	const double cov_tol = 1.0 / ROOT5_DBL_EPSILON;

	if (dx) {
		m->cov_err +=
		    1 + 64 * W * (fabs(dw) + w) + 64 * W1 * w * (1 + w);
	} else {
		m->cov_err += 1 + 64 * W * fabs(dw);
	}
	m->log_cov_scale_err +=
	    0.5 * DBL_EPSILON * (fabs(dpsi) + fabs(log_scale1));

	double err_prod = (1 + m->cov_err) * exp(m->log_cov_scale_err) - 1;

	if (!(err_prod < cov_tol)) {
		if (!(m->mean_err == 0.0))
			recompute_mean(m);
		recompute_cov(m);
		return;
	}
	// cov_diff := W * dw * (x1 - mean)^2
	memset(cov_diff_full, 0, cov_dim_full * sizeof(double));
	blas_dsyr(F77_COV_UPLO, dim, W * dw, xresid, 1, cov_diff_full, dim);

	if (dx) {
		// cov_diff += W * w * [ (x1 - mean) * dx + dx * (x1 - mean) ]
		blas_dsyr2(F77_COV_UPLO, dim, W * w, xresid, 1, dx, 1,
			   cov_diff_full, dim);

		// cov_diff += - W1 * w * (1 + w) * (dx)^2
		blas_dsyr(F77_COV_UPLO, dim, -W1 * w * (1 + w), dx, 1,
			  cov_diff_full, dim);
	}
	// cov += cov_diff
	copy_packed(cov_diff, cov_diff_full, dim, MLOGIT_COV_UPLO);
	blas_daxpy(cov_dim, 1.0, cov_diff, 1, m->cov_, 1);
	m->log_cov_scale_ = log_scale1;
}

double get_deta(struct mlogit *m, const double *dx, const size_t *jdx,
		size_t ndx)
{
	double deta;

	if (jdx) {
		deta = sblas_ddoti(ndx, dx, jdx, m->beta);
	} else {
		assert(ndx == m->dim);
		deta = blas_ddot(m->dim, dx, 1, m->beta, 1);
	}

	return deta;
}

void increment_x(struct mlogit *m, size_t i, const size_t *jdx,
		 const double *dx, size_t ndx)
{
	assert(i < mlogit_ncat(m));
	assert(dx || ndx == 0);

	size_t dim = mlogit_dim(m);
	double *x = m->x_ + i * dim;

	if (jdx) {
#ifndef NDEBUG
		size_t ncat = mlogit_ncat(m);
		size_t k;
		for (k = 0; k < ndx; k++) {
			assert(jdx[k] < ncat);
		}
#endif
		sblas_daxpyi(ndx, 1.0, dx, jdx, x);
	} else if (ndx) {
		assert(ndx == mlogit_dim(m));
		blas_daxpy(dim, 1.0, dx, 1, x, 1);
	}
}

void recompute_all(struct mlogit *m)
{
	recompute_dist(m);
	recompute_mean(m);
	recompute_cov(m);
}

void recompute_dist(struct mlogit *m)
{
	size_t ncat = mlogit_ncat(m);
	size_t dim = mlogit_dim(m);

	const double *beta = m->beta;
	double *eta = m->work->cat_buf;

	blas_dcopy(ncat, m->offset, 1, eta, 1);

	// eta := offset + x * beta
	if (dim > 0)
		blas_dgemv(BLAS_TRANS, dim, ncat, 1.0, m->x_, dim, beta, 1, 1.0, eta, 1);

	catdist_set_all_eta(&m->dist_, eta);
}

void recompute_mean(struct mlogit *m)
{
	size_t ncat = mlogit_ncat(m);
	size_t dim = mlogit_dim(m);

	if (dim == 0)
		return;

	double *mean = m->mean_;
	double *prob = m->work->cat_buf;
	size_t i;

	for (i = 0; i < ncat; i++) {
		prob[i] = catdist_prob(&m->dist_, i);
	}

	// mean := t(X) * prob
	blas_dgemv(BLAS_NOTRANS, dim, ncat, 1.0, m->x_, dim, prob, 1, 0, mean,
		   1);
	m->mean_err = 0.0;
}

void recompute_cov(struct mlogit *m)
{
	size_t ncat = mlogit_ncat(m);
	size_t dim = mlogit_dim(m);

	if (dim == 0)
		return;

	const double *x = m->x_;
	const double *mean = m->mean_;
	double *diff = m->work->dim_buf1;
	double *cov = m->cov_;
	double *cov_full = m->work->cov_diff_full;
	enum blas_uplo uplo = F77_COV_UPLO;

	size_t i;
	double p, ptot;

	// cov := 0; ptot := 0
	memset(cov_full, 0, dim * dim * sizeof(*cov_full));
	ptot = 0;

	for (i = 0; i < ncat; i++) {
		/* diff := mean - x[i,:] */
		memcpy(diff, x + i * dim, dim * sizeof(*diff));
		blas_daxpy(dim, -1.0, mean, 1, diff, 1);

		/* ptot += p[i] */
		p = catdist_prob(&m->dist_, i);
		ptot += p;

		/* cov += p[i] * diff^2 */
		blas_dsyr(uplo, dim, p, diff, 1, cov_full, dim);
	}
	copy_packed(cov, cov_full, dim, MLOGIT_COV_UPLO);

	m->log_cov_scale_ = log(ptot);
	m->cov_err = 0.0;
	m->log_cov_scale_err = 0.0;
}


void mlogit_set_moments(struct mlogit *m, int k)
{
	assert(k >= 0);

	int k0 = m->moments;
	m->moments = k;

	if (k0 < k) {
		if (k0 <= 0)
			recompute_mean(m);
		if (k0 <= 1)
			recompute_cov(m);
	}
}


double *mlogit_x(const struct mlogit *m)
{
	update((struct mlogit *)m);
	return m->x_;
}


struct catdist *mlogit_dist(const struct mlogit *m)
{
	update((struct mlogit *)m);
	return &((struct mlogit *)m)->dist_;
}


double *mlogit_mean(const struct mlogit *m)
{
	assert(m->moments >= 1);
	update((struct mlogit *)m);
	return m->mean_;
}


double *mlogit_cov(const struct mlogit *m, double *cov_scale)
{
	assert(m->moments >= 2);
	update((struct mlogit *)m);
	*cov_scale = exp(m->log_cov_scale_);
	return m->cov_;
}


int mlogit_check(const struct mlogit *m)
{
	int fail = 0;
	size_t n = mlogit_ncat(m);
	size_t p = mlogit_dim(m);
	const double *x = mlogit_x(m);
	const double *beta = mlogit_coefs(m);
	const double *mean = mlogit_mean(m);
	double cov_scale;
	const double *cov = mlogit_cov(m, &cov_scale);
	size_t i, j;

	double *eta0 = xcalloc(n, sizeof(*eta0));
	double *prob0 = xcalloc(n, sizeof(*prob0));
	double *mean0 = xcalloc(p, sizeof(*mean0));
	double *cov0 = xcalloc(p * (p + 1) / 2, sizeof(*cov0));
	double *cov0_copy = xcalloc(p * (p + 1) / 2, sizeof(*cov0_copy));
	double *cov_err = xcalloc(p * (p + 1) / 2, sizeof(*cov_err));
	double *z = xmalloc(p * p * sizeof(*z));
	double *w = xmalloc(p * sizeof(*w));
	double *diff = xcalloc(p, sizeof(*diff));
	double *err_z = xmalloc(p * sizeof(*err_z));
	const enum blas_uplo uplo = F77_COV_UPLO;
	double probtot0 = 0;
	size_t lwork, liwork;

	lwork = lapack_dspevd_lwork(LA_EIG_VEC, p, &liwork);
	double *work = xmalloc(lwork * sizeof(*work));
	ptrdiff_t *iwork = xmalloc(liwork * sizeof(*iwork));

	CHECK(!catdist_check(&m->dist_));

	blas_dcopy(n, m->offset, 1, eta0, 1);

	if (p > 0)
		blas_dgemv(BLAS_TRANS, p, n, 1.0, x, p, beta, 1, 1.0, eta0, 1);

	for (i = 0; i < n; i++) {
		CHECK_APPROX(eta0[i], catdist_eta(&m->dist_, i));
	}

	for (i = 0; i < n; i++) {
		prob0[i] = catdist_prob(&m->dist_, i);
		probtot0 += prob0[i];
	}

	if (m->moments < 1)
		goto out;

	if (p > 0)
		blas_dgemv(BLAS_NOTRANS, p, n, 1.0, x, p, prob0, 1, 0.0, mean0,
			   1);

	for (j = 0; j < p; j++) {
		CHECK_APPROX(mean0[j], mean[j]);
	}

	if (m->moments < 2)
		goto out;

	for (i = 0; i < n; i++) {
		blas_dcopy(p, x + i * p, 1, diff, 1);
		blas_daxpy(p, -1.0, mean0, 1, diff, 1);
		blas_dspr(uplo, p, prob0[i], diff, 1, cov0);
	}
	blas_dscal(p * (p + 1) / 2, 1.0 / probtot0, cov0, 1);

	// cov_err := cov0 - (1/cov_scale) * cov
	blas_dcopy(p * (p + 1) / 2, cov0, 1, cov_err, 1);
	blas_daxpy(p * (p + 1) / 2, -1.0 / cov_scale, cov, 1, cov_err, 1);

	// compute cov0 eigendecomp
	blas_dcopy(p * (p + 1) / 2, cov0, 1, cov0_copy, 1);
	ptrdiff_t info = lapack_dspevd(LA_EIG_VEC, uplo, p, cov0_copy, w, z,
				       MAX(1, p), work, lwork, iwork, liwork);
	assert(info == 0);
	(void)info;

	// check for relative equality of cov0 and (1/cov_scale) * cov
	// on the eigenspaces of cov0
	for (i = p; i > 0; i--) {
		const double *zi = z + (i - 1) * p;
		blas_dspmv(uplo, p, 1.0, cov_err, zi, 1, 0.0, err_z, 1);
		double z_err_z = blas_ddot(p, zi, 1, err_z, 1);
		CHECK(fabs(z_err_z) <=
		      2 * p * SQRT_DBL_EPSILON * (1 + fabs(w[i - 1])));
	}

out:
	free(iwork);
	free(work);
	free(err_z);
	free(diff);
	free(w);
	free(z);
	free(cov_err);
	free(cov0_copy);
	free(cov0);
	free(mean0);
	free(prob0);
	free(eta0);

	return fail;
}

void grow_ind_array(struct mlogit *m, size_t delta)
{
	size_t nz = m->nz;
	size_t nz1 = nz + delta;
	size_t nzmax = m->nzmax;

	if (nz1 <= nzmax)
		return;

	size_t nzmax1 = array_grow(nz, nzmax, delta, SIZE_MAX);
	size_t dim = m->dim;
	assert(nzmax1 >= nz1);

	m->ind = xrealloc(m->ind, nzmax1 * sizeof(*m->ind));
	m->dx = xrealloc(m->dx, nzmax1 * sizeof(*m->dx) * dim);
	m->nzmax = nzmax1;
}


size_t search_ind(struct mlogit *m, size_t i)
{
	const size_t *base = m->ind, *ptr;
	size_t nz;

	for (nz = m->nz; nz != 0; nz >>= 1) {
		ptr = base + (nz >> 1);
		if (i == *ptr) {
			return ptr - m->ind;
		}
		if (i > *ptr) {
			base = ptr + 1;
			nz--;
		}
	}

	/* not found */
	size_t iz = base - m->ind;
	size_t ntail = m->nz - iz;
	size_t dim = m->dim;

	grow_ind_array(m, 1);
	memmove(m->ind + iz + 1, m->ind + iz, ntail * sizeof(*m->ind));
	memmove(m->dx + (iz + 1) * dim, m->dx + iz * dim, ntail * sizeof(*m->dx) * dim);

	m->ind[iz] = i;
	m->offset[iz] = 0.0;
	memset(m->dx + iz * dim, 0, sizeof(*m->dx) * dim);
	m->nz++;

	return iz;
}


size_t find_ind(const struct mlogit *m, size_t i)
{
	const size_t *base = m->ind, *ptr;
	size_t nz;

	for (nz = m->nz; nz != 0; nz >>= 1) {
		ptr = base + (nz >> 1);
		if (i == *ptr) {
			return ptr - m->ind;
		}
		if (i > *ptr) {
			base = ptr + 1;
			nz--;
		}
	}

	/* not found */
	return m->nz;
}
