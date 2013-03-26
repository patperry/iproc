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
#include "matrixutil.h"
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
static void recompute_all(struct mlogit *m);
static void recompute_dist(struct mlogit *m);
static void recompute_mean(struct mlogit *m);
static void recompute_cov(struct mlogit *m);
static void update(struct mlogit *m);
static void update_version(struct mlogit *m);
static void update_dx(struct mlogit *m);
static void update_dist(struct mlogit *m);
static void update_mean(struct mlogit *m);
static void update_cov(struct mlogit *m);
static void compute_cov_diff(struct mlogit *m);

static void grow_ind_array(struct mlogit *m, size_t delta);
static size_t search_ind(struct mlogit *m, size_t i);


void mlogit_work_init(struct mlogit_work *work, size_t ncat, size_t dim)
{
	size_t cov_dim = dim * (dim + 1) / 2;
	size_t cov_dim_full = dim * dim;


	work->cov_diff = xmalloc(cov_dim * sizeof(*work->cov_diff));
	work->cov_diff_full = xmalloc(cov_dim_full * sizeof(*work->cov_diff_full));
	work->cat_buf = xmalloc(ncat * sizeof(*work->cat_buf));
	work->dim_buf1 = xmalloc(dim * sizeof(*work->dim_buf1));
	work->dim_buf2 = xmalloc(dim * sizeof(*work->dim_buf2));

	work->dx = xmalloc(ncat * dim * sizeof(*work->dx));
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
	free(work->dx);
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
	m->x = xmalloc(ncat * dim * sizeof(*m->x));
	m->beta = xmalloc(dim * sizeof(*m->beta));
	m->offset = xmalloc(ncat * sizeof(*m->offset));
	m->mean_ = xmalloc(dim * sizeof(*m->mean_));
	m->cov_ = xmalloc(cov_dim * sizeof(*m->cov_));
	m->log_cov_scale_ = 0;
	m->dim = dim;
	m->moments = 2;

	m->ind = NULL;
	m->x0 = NULL;
	m->nz = 0;
	m->nzmax = 0;

	m->version = 0;
	m->cp = NULL;
	m->ncp = 0;
	m->ncpmax = 0;

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
	free(m->cp);
	free(m->x0);
	free(m->ind);
	free(m->cov_);
	free(m->mean_);
	free(m->offset);
	free(m->beta);
	free(m->x);
	catdist_deinit(&m->dist_);
}

void update_version(struct mlogit *m)
{
	if (m->version == SIZE_MAX) {
		size_t i, n = m->ncp;
		for (i = 0; i < n; i++) {
			m->cp[i]->version = 0;
		}
		m->version = 0;
	}

	m->version++;
}




void mlogit_add_checkpoint(struct mlogit *m, struct mlogit_checkpoint *cp)
{
	if (needs_grow(m->ncp + 1, &m->ncpmax)) {
		m->cp = xrealloc(m->cp, m->ncpmax * sizeof(m->cp[0]));
	}
	m->cp[m->ncp++] = cp;
	cp->version = m->version;
}

void mlogit_remove_checkpoint(struct mlogit *m, struct mlogit_checkpoint *cp)
{
	size_t i, n = m->ncp;
	for (i = n; i > 0; i--) {
		if (m->cp[i - 1] == cp) {
			memmove(m->cp + i - 1, m->cp + i, (n - i) * sizeof(m->cp[0]));
			break;
		}
	}
}

int mlogit_checkpoint_passed(const struct mlogit_checkpoint *cp, const struct mlogit *m)
{
	return cp->version < m->version;
}

void mlogit_checkpoint_set(struct mlogit_checkpoint *cp, const struct mlogit *m)
{
	cp->version = m->version;
}



void clear(struct mlogit *m)
{
	size_t ncat = mlogit_ncat(m);
	size_t dim = m->dim;
	size_t cov_dim = dim * (dim + 1) / 2;

	catdist_set_all_eta(&m->dist_, NULL);
	memset(m->x, 0, ncat * dim * sizeof(*m->x));
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
	update_version(m);
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
	update_version(m);
}

void mlogit_set_all_x(struct mlogit *m, const double *x)
{
	size_t len = mlogit_ncat(m) * mlogit_dim(m) * sizeof(*m->x);
	if (x) {
		memcpy(m->x, x, len);
	} else {
		memset(m->x, 0, len);
	}
	m->nz = 0;

	recompute_all(m);
	update_version(m);
}


void mlogit_set_offset(struct mlogit *m, size_t i, double offset)
{
	assert(i < mlogit_ncat(m));
	assert(offset < INFINITY);

	// offset[i] := offset
	m->offset[i] = offset;

	if (m->nz) {
		/* can't handle updates to both offset and x at same time;
		 * flush updates to x right now */
		m->nz = 0;
	}

	recompute_all(m);
	update_version(m);
}


void mlogit_set_x(struct mlogit *m, size_t i, const size_t *jx,
		  const double *x, size_t nx)
{
	assert(i < mlogit_ncat(m));
	assert(x || nx == 0);

	size_t dim = mlogit_dim(m);
	search_ind(m, i); /* store x0[i] if not already present */
	double *dst = m->x + i * dim;

	if (jx) {
#ifndef NDEBUG
		size_t k;
		for (k = 0; k < nx; k++) {
			assert(jx[k] < dim);
		}
#endif
		sblas_dsctr(nx, x, jx, dst);
	} else if (nx) {
		assert(nx == mlogit_dim(m));
		memcpy(dst, x, dim * sizeof(dst[0]));
	}

	update_version(m);

}



void mlogit_inc_x(struct mlogit *m, size_t i, const size_t *jdx,
		  const double *dx, size_t ndx)
{
	assert(i < mlogit_ncat(m));
	assert(dx || ndx == 0);

	size_t dim = mlogit_dim(m);
	double *x = m->x + i * dim;
	double *x1 = m->work->dim_buf1;

	if (jdx) {
#ifndef NDEBUG
		size_t k;
		for (k = 0; k < ndx; k++) {
			assert(jdx[k] < dim);
		}
#endif
		sblas_dgthr(ndx, x, x1, jdx);
		blas_daxpy(ndx, 1.0, dx, 1, x1, 1);
	} else if (ndx) {
		assert(ndx == mlogit_dim(m));

		memcpy(x1, x, dim * sizeof(double));
		blas_daxpy(dim, 1.0, dx, 1, x1, 1);
	}

	mlogit_set_x(m, i, jdx, x1, ndx);
	update_version(m);
}


void update(struct mlogit *m)
{
	if (m->nz) {
		update_dx(m);
		update_dist(m);
		update_mean(m);
		update_cov(m);
		m->nz = 0;
	}
}


void update_dx(struct mlogit *m)
{
	if (!mlogit_dim(m))
		return;

	const size_t dim = mlogit_dim(m);
	double *dx = m->work->dx;

	// compute dx
	size_t iz, nz = m->nz;
	for (iz = 0; iz < nz; iz++) {
		size_t i = m->ind[iz];
		const double *x0_i = m->x0 + iz * dim;
		const double *x_i = m->x + i * dim;

		double *dx_i = dx + iz * dim;

		/* copy x */
		memcpy(dx_i, x_i, dim * sizeof(double));

		/* dx[i] := x[i] - x0[i] */
		blas_daxpy(dim, -1.0, x0_i, 1, dx_i, 1);
	}
}


void update_dist(struct mlogit *m)
{
	if (!mlogit_dim(m))
		return;

	if (m->nz >= 0 * mlogit_ncat(m) / 2) {
		recompute_dist(m);
		return;
	}

	const size_t dim = mlogit_dim(m);
	const double psi = catdist_psi(&m->dist_);
	const double *x = m->x;
	double *eta = m->work->w0;
	double *eta1 = m->work->dw;

	// copy old values of eta and compute new eta values
	size_t iz, nz = m->nz;
	for (iz = 0; iz < nz; iz++) {
		size_t i = m->ind[iz];
		eta[iz] = catdist_eta(&m->dist_, i);
		eta1[iz] = blas_ddot(dim, m->beta, 1, x + i * dim, 1);
	}

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
	const double *dx = m->work->dx;
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
		const double *x_i = m->x + i * dim;
		const double *dx_i = dx + iz * dim;

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
	double *dcov = m->work->cov_diff;

	compute_cov_diff(m);

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
		printf("!");

		blas_dscal(cov_dim, -exp(m->log_cov_scale_ - log_scale1), dcov, 1);
		blas_daxpy(cov_dim, 1.0, m->cov_, 1, dcov, 1);
	}
}


void compute_cov_diff(struct mlogit *m)
{
	size_t dim = mlogit_dim(m);
	const double *x0 = m->x0;
	const double *x1 = m->x;
	const double *w0 = m->work->w0;
	const double *mean0 = m->work->mean0;
	const double *dmean = m->work->dmean;
	const double log_scale = m->log_cov_scale_;
	const double dpsi = m->work->dpsi;
	const double W1 = exp(log_scale + dpsi);
	double *dcov = m->work->cov_diff_full;
	double *diff = m->work->xbuf1;

	size_t iz, nz = m->nz;
	size_t nk = 0;
	double *diff_i = diff;

	assert(dim > 0);

	memset(dcov, 0, dim * dim * sizeof(*dcov));

	// handle positive terms
	for (iz = 0; iz < nz; iz++) {
		const size_t i = m->ind[iz];
		const double *x1_i = x1 + i * dim;
		const double w1_i = catdist_prob(&m->dist_, i);

		if (w1_i > 0) {
			blas_dcopy(dim, x1_i, 1, diff_i, 1);
			blas_daxpy(dim, -1.0, mean0, 1, diff_i, 1);
			blas_dscal(dim, sqrt(w1_i), diff_i, 1);
			m->cov_err += 64 * W1 * w1_i;

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

	// handle dmean and negative terms
	blas_dcopy(dim, dmean, 1, diff_i, 1);
	m->cov_err += W1;

	diff_i += dim;
	nk++;

	for (iz = 0; iz < nz; iz++) {
		const double *x0_i = x0 + iz * dim;
		const double w0_i = w0[iz];

		if (w0_i > 0) {
			blas_dcopy(dim, x0_i, 1, diff_i, 1);
			blas_daxpy(dim, -1.0, mean0, 1, diff_i, 1);
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
	// note: nz = 0 implies dmean = 0; no need to call dsyrk

	packed_dgthr(F77_COV_UPLO, dim, dcov, dim, m->work->cov_diff);
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
		blas_dgemv(BLAS_TRANS, dim, ncat, 1.0, m->x, dim, beta, 1, 1.0, eta, 1);

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
	blas_dgemv(BLAS_NOTRANS, dim, ncat, 1.0, m->x, dim, prob, 1, 0, mean, 1);
	m->mean_err = 0.0;
}

void recompute_cov(struct mlogit *m)
{
	size_t i, ncat = mlogit_ncat(m);
	size_t dim = mlogit_dim(m);

	if (dim == 0)
		return;

	const double *x = m->x;
	const double *mean = m->mean_;
	double *cov = m->cov_;
	double *cov_full = m->work->cov_diff_full;
	double *diff = m->work->xbuf1;
	double ptot;

	size_t nk = 0;
	double *diff_i = diff;

	// cov := 0; ptot := 0
	memset(cov_full, 0, dim * dim * sizeof(*cov_full));
	ptot = 0;

	for (i = 0; i < ncat; i++) {
		double p = catdist_prob(&m->dist_, i);

		if (p) {
			/* diff := mean - x[i,:] */
			memcpy(diff_i, x + i * dim, dim * sizeof(*diff_i));
			blas_daxpy(dim, -1.0, mean, 1, diff_i, 1);
			blas_dscal(dim, sqrt(p), diff_i, 1);
			ptot += p;

			diff_i += dim;
			nk++;
		}

		if (nk == BLOCK_SIZE || i + 1 == ncat) {
			blas_dsyrk(F77_COV_UPLO, BLAS_NOTRANS, dim, nk, 1.0, diff, dim, 1.0, cov_full, dim);

			diff_i = diff;
			nk = 0;
		}
	}

	assert(nk == 0);
	assert(diff_i == diff);

	packed_dgthr(F77_COV_UPLO, dim, cov_full, dim, cov);

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

	update_version(m);
}


double *mlogit_x(const struct mlogit *m)
{
	return m->x;
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
	size_t dim = m->dim;
	size_t nz = m->nz;
	size_t nz1 = nz + delta;

	if (needs_grow(nz1, &m->nzmax)) {
		m->ind = xrealloc(m->ind, m->nzmax * sizeof(m->ind[0]));
		m->x0 = xrealloc(m->x0, m->nzmax * dim * sizeof(m->x0[0]));
	}
}


size_t search_ind(struct mlogit *m, size_t i)
{
	ptrdiff_t siz = find_index(i, m->ind, m->nz);
	size_t iz;

	if (siz < 0) {
		iz = ~siz;
		
		size_t ntail = m->nz - iz;
		size_t dim = m->dim;

		grow_ind_array(m, 1);
		memmove(m->ind + iz + 1, m->ind + iz, ntail * sizeof(*m->ind));
		memmove(m->x0 + (iz + 1) * dim, m->x0 + iz * dim, ntail * sizeof(*m->x0) * dim);


		m->ind[iz] = i;
		m->offset[iz] = 0.0;

		/* store the initial value of x in x0 */
		memcpy(m->x0 + iz * dim, m->x + i * dim, dim * sizeof(*m->x0));
		m->nz++;
	} else {
		iz = siz;
	}

	return iz;
}
