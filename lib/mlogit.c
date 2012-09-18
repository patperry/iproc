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

static void clear(struct mlogit *m);
static double get_deta(struct mlogit *m, const double *dx, const size_t *jdx,
		       size_t ndx);
static void increment_x(struct mlogit *m, size_t i, const size_t *jdx,
			const double *dx, size_t ndx);
static void recompute_all(struct mlogit *m);
static void recompute_dist(struct mlogit *m);
static void recompute_mean(struct mlogit *m);
static void recompute_cov(struct mlogit *m);
static void update_cov(struct mlogit *m, const double *dx, const double *xresid,
		       const double w, const double dw, const double dpsi);
static void update_mean(struct mlogit *m, const double *dx,
			const double *xresid, const double w, const double dw);


void mlogit_work_init(struct mlogit_work *work, size_t ncat, size_t dim)
{
	size_t cov_dim = dim * (dim + 1) / 2;
	work->mean_diff = xmalloc(dim * sizeof(*work->mean_diff));
	work->cov_diff = xmalloc(cov_dim * sizeof(*work->cov_diff));
	work->cat_buf = xmalloc(ncat * sizeof(*work->cat_buf));
	work->dim_buf1 = xmalloc(dim * sizeof(*work->dim_buf1));
	work->dim_buf2 = xmalloc(dim * sizeof(*work->dim_buf2));

}

void mlogit_work_deinit(struct mlogit_work *work)
{
	free(work->dim_buf2);
	free(work->dim_buf1);
	free(work->cat_buf);
	free(work->cov_diff);
	free(work->mean_diff);
}


void mlogit_init(struct mlogit *m, size_t ncat, size_t dim, struct mlogit_work *work)
{
	assert(ncat == 0 || dim <= SIZE_MAX / sizeof(double) / ncat);
	assert(dim <= SIZE_MAX / sizeof(double) / (dim + 1));

	size_t cov_dim = dim * (dim + 1) / 2;

	catdist_init(&m->dist, ncat);
	m->x = xmalloc(ncat * dim * sizeof(*m->x));
	m->beta = xmalloc(dim * sizeof(*m->beta));
	m->offset = xmalloc(ncat * sizeof(*m->offset));
	m->mean = xmalloc(dim * sizeof(*m->mean));
	m->cov = xmalloc(cov_dim * sizeof(*m->cov));
	m->log_cov_scale = 0;
	m->dim = dim;
	m->moments = 2;

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
	free(m->cov);
	free(m->mean);
	free(m->offset);
	free(m->beta);
	free(m->x);
	catdist_deinit(&m->dist);
}


void clear(struct mlogit *m)
{
	size_t ncat = mlogit_ncat(m);
	size_t dim = m->dim;
	size_t cov_dim = dim * (dim + 1) / 2;

	catdist_set_all_eta(&m->dist, NULL);
	memset(m->x, 0, ncat * dim * sizeof(*m->x));
	memset(m->beta, 0, dim * sizeof(*m->beta));
	memset(m->offset, 0, ncat * sizeof(*m->offset));
	memset(m->mean, 0, dim * sizeof(*m->mean));
	memset(m->cov, 0, cov_dim * sizeof(*m->cov));
	m->log_cov_scale = 0;
	m->mean_err = 0.0;
	m->cov_err = 0.0;
	m->log_cov_scale_err = 0.0;
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
	size_t len = mlogit_ncat(m) * mlogit_dim(m) * sizeof(*m->x);
	if (x) {
		memcpy(m->x, x, len);
	} else {
		memset(m->x, 0, len);
	}

	recompute_all(m);
}

void mlogit_set_offset(struct mlogit *m, size_t i, double offset)
{
	assert(i < mlogit_ncat(m));
	assert(offset < INFINITY);

	const size_t dim = mlogit_dim(m);
	const double psi = catdist_psi(&m->dist);
	const double eta = catdist_eta(&m->dist, i);
	const double deta = offset - m->offset[i];
	const double eta1 = eta + deta;

	// offset[i] := offset
	m->offset[i] = offset;

	// update eta, psi
	catdist_set_eta(&m->dist, i, eta1);

	if (m->moments < 1)
		return;

	// compute weight changes
	const double psi1 = catdist_psi(&m->dist);
	const double dpsi = psi1 - psi;
	const double w1 = exp(eta1 - psi1);
	const double w = exp(eta - psi1);
	const double dw = w1 - w;

	// xresid := x1 - mean
	double *xresid = m->work->dim_buf1;
	blas_dcopy(dim, m->x + i * dim, 1, xresid, 1);
	blas_daxpy(dim, -1.0, m->mean, 1, xresid, 1);

	update_mean(m, NULL, xresid, w, dw);
	update_cov(m, NULL, xresid, w, dw, dpsi);
}

void mlogit_inc_x(struct mlogit *m, size_t i, const size_t *jdx,
		  const double *dx, size_t ndx)
{
	assert(i < mlogit_ncat(m));
	assert(dx || ndx == 0);
	assert(jdx || ndx == 0 || ndx == mlogit_dim(m));

	size_t dim = mlogit_dim(m);

	if (ndx == 0 || dim == 0)
		return;

	const double psi = catdist_psi(&m->dist);
	const double eta = catdist_eta(&m->dist, i);
	const double deta = get_deta(m, dx, jdx, ndx);
	const double eta1 = eta + deta;

	// x := x + dx
	increment_x(m, i, jdx, dx, ndx);

	// update eta, psi
	catdist_set_eta(&m->dist, i, eta1);

	if (m->moments < 1)
		return;

	// compute weight changes
	const double psi1 = catdist_psi(&m->dist);
	const double dpsi = psi1 - psi;
	const double w1 = exp(eta1 - psi1);
	const double w = exp(eta - psi1);
	const double dw = w1 - w;

	// xresid := x1 - mean
	double *xresid = m->work->dim_buf1;
	blas_dcopy(dim, m->x + i * dim, 1, xresid, 1);
	blas_daxpy(dim, -1.0, m->mean, 1, xresid, 1);

	// set ddx (dense dx)
	double *ddx = jdx ? m->work->dim_buf2 : (double *)dx;

	if (jdx) {
		memset(ddx, 0, dim * sizeof(*ddx));
		sblas_dsctr(ndx, dx, jdx, ddx);
	}

	update_mean(m, ddx, xresid, w, dw);
	update_cov(m, ddx, xresid, w, dw, dpsi);
}

static void update_mean(struct mlogit *m, const double *dx,
			const double *xresid, const double w, const double dw)
{
	if (m->moments < 1)
		return;

	size_t dim = mlogit_dim(m);
	double *mean_diff = m->work->mean_diff;
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
	blas_daxpy(dim, 1.0, mean_diff, 1, m->mean, 1);
}

static void update_cov(struct mlogit *m, const double *dx, const double *xresid,
		       const double w, const double dw, const double dpsi)
{
	if (m->moments < 2)
		return;

	const size_t dim = mlogit_dim(m);
	const size_t cov_dim = dim * (dim + 1) / 2;
	double *cov_diff = m->work->cov_diff;
	const double log_scale = m->log_cov_scale;
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
	memset(cov_diff, 0, cov_dim * sizeof(double));
	blas_dspr(F77_COV_UPLO, dim, W * dw, xresid, 1, cov_diff);

	if (dx) {
		// cov_diff += W * w * [ (x1 - mean) * dx + dx * (x1 - mean) ]
		blas_dspr2(F77_COV_UPLO, dim, W * w, xresid, 1, dx, 1,
			   cov_diff);

		// cov_diff += - W1 * w * (1 + w) * (dx)^2
		blas_dspr(F77_COV_UPLO, dim, -W1 * w * (1 + w), dx, 1,
			  cov_diff);
	}
	// cov += cov_diff
	blas_daxpy(cov_dim, 1.0, cov_diff, 1, m->cov, 1);
	m->log_cov_scale = log_scale1;
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
	double *x = m->x + i * dim;

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
		blas_dgemv(BLAS_TRANS, dim, ncat, 1.0, m->x, dim, beta, 1, 1.0, eta, 1);

	catdist_set_all_eta(&m->dist, eta);
}

void recompute_mean(struct mlogit *m)
{
	size_t ncat = mlogit_ncat(m);
	size_t dim = mlogit_dim(m);

	if (dim == 0)
		return;

	double *mean = m->mean;
	double *prob = m->work->cat_buf;
	size_t i;

	for (i = 0; i < ncat; i++) {
		prob[i] = catdist_prob(&m->dist, i);
	}

	// mean := t(X) * prob
	blas_dgemv(BLAS_NOTRANS, dim, ncat, 1.0, m->x, dim, prob, 1, 0, mean,
		   1);
	m->mean_err = 0.0;
}

void recompute_cov(struct mlogit *m)
{
	size_t ncat = mlogit_ncat(m);
	size_t dim = mlogit_dim(m);

	if (dim == 0)
		return;

	const double *x = m->x;
	const double *mean = m->mean;
	double *diff = m->work->dim_buf1;
	double *cov = m->cov;
	enum blas_uplo uplo = F77_COV_UPLO;

	size_t i;
	double p, ptot;

	// cov := 0; ptot := 0
	memset(cov, 0, dim * (dim + 1) / 2 * sizeof(*cov));
	ptot = 0;

	for (i = 0; i < ncat; i++) {
		/* diff := mean - x[i,:] */
		memcpy(diff, x + i * dim, dim * sizeof(*diff));
		blas_daxpy(dim, -1.0, mean, 1, diff, 1);

		/* ptot += p[i] */
		p = catdist_prob(&m->dist, i);
		ptot += p;

		/* cov += p[i] * diff^2 */
		blas_dspr(uplo, dim, p, diff, 1, cov);
	}

	m->log_cov_scale = log(ptot);
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

	CHECK(!catdist_check(&m->dist));

	blas_dcopy(n, m->offset, 1, eta0, 1);

	if (p > 0)
		blas_dgemv(BLAS_TRANS, p, n, 1.0, x, p, beta, 1, 1.0, eta0, 1);

	for (i = 0; i < n; i++) {
		CHECK_APPROX(eta0[i], catdist_eta(&m->dist, i));
	}

	for (i = 0; i < n; i++) {
		prob0[i] = catdist_prob(&m->dist, i);
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
