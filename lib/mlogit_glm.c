#include "port.h"
#include <float.h>  // DBL_MANT_DIG
#include <stdint.h> // SIZE_MAX
#include <stdio.h>
#include <stdlib.h> // free
#include <string.h> // memcpy, memset
#include "blas.h"   // blas_gemv
#include "coreutil.h"
#include "ieee754.h" // double_eqrel
#include "sblas.h"
#include "xalloc.h" // xcalloc
#include "mlogit_glm.h"


static double get_deta(struct mlogit_glm *m, size_t i, const double *dx, const size_t *jdx, size_t ndx);
static void increment_x(struct mlogit_glm *m, size_t i, const double *dx, const size_t *jdx, size_t ndx);
static void recompute_all(struct mlogit_glm *m);
static void recompute_values(struct mlogit_glm *m);
static void recompute_mean(struct mlogit_glm *m);
static void recompute_cov(struct mlogit_glm *m);


#define assert_approx(x, y) \
	assert(double_eqrel((x), (y)) >= DBL_MANT_DIG / 2 \
	       || (fabs(x) <= SQRT_DBL_EPSILON && fabs(y) <= SQRT_DBL_EPSILON))


void mlogit_glm_init(struct mlogit_glm *m, size_t ncat, size_t dim)
{
	assert(ncat == 0 || dim <= SIZE_MAX / sizeof(double) / ncat);
	assert(dim <= SIZE_MAX / sizeof(double) / (dim + 1));
	
	size_t cov_dim = dim * (dim + 1) / 2;
	
	mlogit_init(&m->values, ncat);
	m->x = xcalloc(ncat * dim, sizeof(*m->x));
	m->beta = xcalloc(dim, sizeof(*m->beta));
	m->mean = xcalloc(dim, sizeof(*m->mean));
	m->mean_diff = xcalloc(dim, sizeof(*m->mean));
	
	m->cov = xcalloc(cov_dim, sizeof(*m->cov));
	m->log_cov_scale = 0;
	m->cov_diff = xcalloc(cov_dim, sizeof(*m->cov_diff));
	
	m->cat_buf = xcalloc(ncat, sizeof(*m->cat_buf));
	m->dim_buf = xcalloc(dim, sizeof(*m->dim_buf));	
	m->dim = dim;
	m->mean_err = 0.0;
}

void mlogit_glm_deinit(struct mlogit_glm *m)
{
	free(m->dim_buf);	
	free(m->cat_buf);
	free(m->cov_diff);
	free(m->cov);
	free(m->mean_diff);
	free(m->mean);	
	free(m->beta);
	free(m->x);
	mlogit_deinit(&m->values);
}


void mlogit_glm_set_coefs(struct mlogit_glm *m, const double *beta)
{
	size_t len = mlogit_glm_dim(m) * sizeof(*m->beta);
	if (beta) {
		memcpy(m->beta, beta, len);
	} else {
		memset(m->beta, 0, len);
	}

	recompute_all(m);
}

void mlogit_glm_set_all_x(struct mlogit_glm *m, const double *x)
{
	size_t len = mlogit_glm_ncat(m) * mlogit_glm_dim(m) * sizeof(*m->x);
	if (x) {
		memcpy(m->x, x, len);
	} else {
		memset(m->x, 0, len);
	}
	
	recompute_all(m);
}


void mlogit_glm_inc_x(struct mlogit_glm *m, size_t i, const double *dx, const size_t *jdx, size_t ndx)
{
	size_t dim = mlogit_glm_dim(m);
	
	if (ndx == 0 || dim == 0)
		return;

	const double psi = mlogit_psi(&m->values);
	const double eta = mlogit_eta(&m->values, i);
	const double deta = get_deta(m, i, dx, jdx, ndx);
	const double eta1 = eta + deta;
	double *mean_diff = m->mean_diff;
	double *diff = m->dim_buf;
	double *mean = m->mean;
	double *x = m->x + i * dim;

	// x := x + dx
	increment_x(m, i, dx, jdx, ndx);
	
	// set ddx
	const double *ddx;
	
	if (jdx) {
		memset(m->dim_buf, 0, dim * sizeof(*ddx));
		sblas_dsctr(ndx, dx, jdx, m->dim_buf);
		ddx = m->dim_buf;
	} else {
		assert(ndx == dim);
		ddx = dx;
	}
	
	// update eta, psi
	mlogit_set_eta(&m->values, i, eta1);

	const double psi1 = mlogit_psi(&m->values);
	const double w1 = exp(eta1 - psi1);
	const double w = exp(eta - psi1);
	const double dw = w1 - w;
	const double tol = 1.0 / ROOT4_DBL_EPSILON;

	m->mean_err += 1 + 8 * (fabs(dw) + w); // approximate relative error from update

	// diff := x1 - mean
	blas_dcopy(dim, x, 1, mean_diff, 1);
	blas_daxpy(dim, -1.0, mean, 1, mean_diff, 1);

	
	if (m->mean_err <= tol) {
		// mean_diff := dw * (x1 - mean)
		blas_dcopy(dim, diff, 1, mean_diff, 1);
		blas_dscal(dim, dw, mean_diff, 1);

		// mean_diff += w * dx
		blas_daxpy(ndx, w, ddx, 1, mean_diff, 1);

		// mean += mean_diff
		blas_daxpy(dim, 1.0, mean_diff, 1, mean, 1);
	} else {
		printf("!");
		fflush(stdout);
		recompute_mean(m);
	}

	size_t cov_dim = dim * (dim + 1) / 2;	
	enum blas_uplo uplo = MLOGIT_GLM_COV_UPLO == BLAS_LOWER ? BLAS_UPPER : BLAS_LOWER;
	double *cov_diff = m->cov_diff;
	const double log_cov_scale = m->log_cov_scale;	
	const double W = exp(psi - log_cov_scale);
	const double W1 = exp(psi1 - log_cov_scale);	
	double *cov = m->cov;

	// cov_diff := W * dw * (x1 - mean)^2
	memset(cov_diff, 0, cov_dim * sizeof(double));
	blas_dspr(uplo, dim, W * dw, diff, 1, cov_diff);
	
	// cov_diff += W * w * [ (x1 - mean) * dx + dx * (x1 - mean) ]
	blas_dspr2(uplo, dim, W * w, diff, 1, ddx, 1, cov_diff);
	
	// cov_diff += - W1 * w * (1 + w) * (dx)^2
	blas_dspr(uplo, dim, -W1 * w * (1 + w), ddx, 1, cov_diff);
	
	// cov += cov_diff
	blas_daxpy(cov_dim, 1.0, cov_diff, 1, cov, 1);
	
	/* update cov */
	//recompute_cov(m);
	
	_mlogit_glm_check_invariants(m); // DEBUG
					 //free(x0); // DEBUG
					 //free(mean0); // DEBUG
}


static double get_deta(struct mlogit_glm *m, size_t i, const double *dx, const size_t *jdx, size_t ndx)
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


void increment_x(struct mlogit_glm *m, size_t i, const double *dx, const size_t *jdx, size_t ndx)
{
	assert(i < mlogit_glm_ncat(m));
	assert(dx || ndx == 0);
	
	size_t ncat = mlogit_glm_ncat(m);
	size_t dim = mlogit_glm_dim(m);
	double *x = m->x + i * dim;
	
	if (jdx) {
#ifndef NDEBUG
		size_t k;
		for (k = 0; k < ndx; k++) {
			assert(jdx[k] < ncat);
		}
#endif
		sblas_daxpyi(ndx, 1.0, dx, jdx, x);
	} else if (ndx) {
		assert(ndx == mlogit_glm_dim(m));
		blas_daxpy(dim, 1.0, dx, 1, x, 1);
	}
}


void recompute_all(struct mlogit_glm *m)
{
	recompute_values(m);
	recompute_mean(m);
	recompute_cov(m);
}


void recompute_values(struct mlogit_glm *m)
{
	size_t ncat = mlogit_glm_ncat(m);
	size_t dim = mlogit_glm_dim(m);
	
	if (dim == 0)
		return;
	
	const double *beta = m->beta;
	double *eta = m->cat_buf;	
	
	// tmp := x * beta
	blas_dgemv(BLAS_TRANS, dim, ncat, 1.0, m->x, dim, beta, 1, 0.0, eta, 1);
	mlogit_set_all_eta(&m->values, eta);
}


void recompute_mean(struct mlogit_glm *m)
{
	size_t ncat = mlogit_glm_ncat(m);
	size_t dim = mlogit_glm_dim(m);
	
	if (dim == 0)
		return;
	
	double *mean = m->mean;
	double *prob = m->cat_buf;	
	size_t i;
	
	for (i = 0; i < ncat; i++) {
		prob[i] = mlogit_prob(&m->values, i);
	}
	
	// mean := t(X) * prob
	blas_dgemv(BLAS_NOTRANS, dim, ncat, 1.0, m->x, dim, prob, 1, 0, mean, 1);
	memset(m->mean_diff, 0, dim * sizeof(*m->mean_diff));
	m->mean_err = 0.0;	
}


void recompute_cov(struct mlogit_glm *m)
{
	size_t ncat = mlogit_glm_ncat(m);
	size_t dim = mlogit_glm_dim(m);
	size_t cov_dim = dim * (dim + 1) / 2;
	
	if (dim == 0)
		return;

	const double *x = m->x;
	const double *mean = m->mean;
	double *diff = m->dim_buf;
	double *cov = m->cov;
	enum blas_uplo uplo = MLOGIT_GLM_COV_UPLO == BLAS_LOWER ? BLAS_UPPER : BLAS_LOWER;
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
		p = mlogit_prob(&m->values, i);
		ptot += p;
		
		/* cov += p[i] * diff^2 */
		blas_dspr(uplo, dim, p, diff, 1, cov);
	}
	
	m->log_cov_scale = log(ptot);
	memset(m->cov_diff, 0, cov_dim * sizeof(*m->mean_diff));
}


void _mlogit_glm_check_invariants(const struct mlogit_glm *m)
{
	size_t n = mlogit_glm_ncat(m);
	size_t p = mlogit_glm_dim(m);
	const double *x = mlogit_glm_x(m);	
	const double *beta = mlogit_glm_coefs(m);
	const double *mean = mlogit_glm_mean(m);
	double cov_scale;
	const double *cov = mlogit_glm_cov(m, &cov_scale);
	size_t i, j, k;

	double *eta0 = xcalloc(n, sizeof(*eta0));
	double *prob0 = xcalloc(n, sizeof(*prob0));
	double *mean0 = xcalloc(p, sizeof(*mean0));
	double *cov0 = xcalloc(p * (p + 1) / 2, sizeof(*cov0));
	double *diff = xcalloc(p, sizeof(*diff));
	const enum blas_uplo uplo = (MLOGIT_GLM_COV_UPLO == BLAS_LOWER ? BLAS_UPPER : BLAS_LOWER);	
	double probtot0 = 0;
	
	if (p > 0)
		blas_dgemv(BLAS_TRANS, p, n, 1.0, x, p, beta, 1, 0.0, eta0, 1);
	
	for (i = 0; i < n; i++) {
		assert_approx(eta0[i], mlogit_eta(&m->values, i));
	}
	
	for (i = 0; i < n; i++) {
		prob0[i] = mlogit_prob(&m->values, i);
		probtot0 += prob0[i];
	}
	
	if (p > 0)
		blas_dgemv(BLAS_NOTRANS, p, n, 1.0, x, p, prob0, 1, 0.0, mean0, 1);
	
	for (j = 0; j < p; j++) {
		assert_approx(mean0[j], mean[j]);
	}
	
	for (i = 0; i < n; i++) {
		blas_dcopy(p, x + i * p, 1, diff, 1);
		blas_daxpy(p, -1.0, mean0, 1, diff, 1);
		blas_dspr(uplo, p, prob0[i], diff, 1, cov0);
	}
	blas_dscal(p * (p + 1) / 2, 1.0 / probtot0, cov0, 1);
	
	for (k = 0; k < p * (p + 1) / 2; k++) {
		assert_approx(cov0[k], cov[k] / cov_scale);
	}
	
	free(diff);
	free(cov0);
	free(mean0);
	free(prob0);
	free(eta0);
}
