#include "port.h"
#include <stdint.h> // SIZE_MAX
#include <stdlib.h> // free
#include <string.h> // memcpy, memset
#include "blas.h"   // blas_gemv
#include "xalloc.h" // xcalloc
#include "mlogit_glm.h"

static void increment_x(struct mlogit_glm *m, size_t i, const double *dx, const size_t *jdx, size_t ndx);
static void recompute(struct mlogit_glm *m);
static void recompute_mean(struct mlogit_glm *m);
static void recompute_cov(struct mlogit_glm *m);



void mlogit_glm_init(struct mlogit_glm *m, size_t ncat, size_t dim)
{
	assert(ncat == 0 || dim <= SIZE_MAX / sizeof(double) / ncat);
	assert(dim <= SIZE_MAX / sizeof(double) / (dim + 1));
	
	mlogit_init(&m->values, ncat);
	m->x = xcalloc(ncat * dim, sizeof(*m->x));
	m->beta = xcalloc(dim, sizeof(*m->beta));
	m->mean = xcalloc(dim, sizeof(*m->mean));
	m->cov = xcalloc(dim * (dim + 1) / 2, sizeof(*m->cov));
	m->cat_buf = xcalloc(ncat, sizeof(*m->cat_buf));
	m->dim_buf = xcalloc(dim, sizeof(*m->dim_buf));	
	m->dim = dim;
	
}

void mlogit_glm_deinit(struct mlogit_glm *m)
{
	free(m->dim_buf);	
	free(m->cat_buf);
	free(m->cov);
	free(m->mean);	
	free(m->beta);
	free(m->x);
	mlogit_deinit(&m->values);
}


void mlogit_glm_set_coefs(struct mlogit_glm *m, const double *beta)
{
	size_t len = mlogit_glm_dim(m);
	if (beta) {
		memcpy(m->beta, beta, len);
	} else {
		memset(m->beta, 0, len);
	}

	recompute(m);
}

void mlogit_glm_set_all_x(struct mlogit_glm *m, const double *x)
{
	size_t len = mlogit_glm_ncat(m) * mlogit_glm_dim(m) * sizeof(*m->x);
	if (x) {
		memcpy(m->x, x, len);
	} else {
		memset(m->x, 0, len);
	}
	
	recompute(m);
}


void mlogit_glm_inc_x(struct mlogit_glm *m, size_t i, const double *dx, const size_t *jdx, size_t ndx)
{
	increment_x(m, i, dx, jdx, ndx);
	recompute(m);
}


void increment_x(struct mlogit_glm *m, size_t i, const double *dx, const size_t *jdx, size_t ndx)
{
	size_t j, k;
	
	assert(i < mlogit_glm_ncat(m));
	assert(dx || ndx == 0);
	
	size_t ncat = mlogit_glm_ncat(m);
	size_t dim = mlogit_glm_dim(m);
	double *x = m->x;
	
	if (jdx) {
		for (k = 0; k < ndx; k++) {
			assert(jdx[k] < ncat);
			x[i * dim + jdx[k]] += dx[k];
		}
	} else if (ndx) {
		assert(ndx == mlogit_glm_dim(m));
		for (j = 0; j < ndx; j++) {
			x[i * dim + j] += dx[j];
		}
	}
}


void recompute(struct mlogit_glm *m)
{
	recompute_mean(m);
	recompute_cov(m);
}


void recompute_mean(struct mlogit_glm *m)
{
	size_t ncat = mlogit_glm_ncat(m);
	size_t dim = mlogit_glm_dim(m);
	
	if (dim == 0)
		return;
	
	const double *beta = m->beta;
	double *mean = m->mean;
	double *tmp = m->cat_buf;	
	size_t i;

	// tmp := x * beta
	blas_dgemv(BLAS_TRANS, dim, ncat, 1.0, m->x, dim, beta, 1, 0.0, tmp, 1);
	mlogit_set_all_eta(&m->values, tmp);
	
	// tmp := p
	for (i = 0; i < ncat; i++) {
		tmp[i] = mlogit_prob(&m->values, i);
	}
	
	// mean := t(X) * p
	blas_dgemv(BLAS_NOTRANS, dim, ncat, 1.0, m->x, dim, tmp, 1, 0, mean, 1);
}


void recompute_cov(struct mlogit_glm *m)
{
	size_t ncat = mlogit_glm_ncat(m);
	size_t dim = mlogit_glm_dim(m);
	
	if (dim == 0)
		return;

	const double *x = m->x;
	const double *mean = m->mean;
	double *tmp = m->dim_buf;	
	double *cov = m->cov;
	enum blas_uplo uplo = MLOGIT_GLM_COV_UPLO == BLAS_LOWER ? BLAS_UPPER : BLAS_LOWER;
	size_t i;
	double ptot;
	double p;

	/* cov := 0; ptot := 0 */
	memset(cov, 0, dim * (dim + 1) / 2 * sizeof(*cov));
	ptot = 0;
	
	for (i = 0; i < ncat; i++) {
		/* tmp := mean - x[i,:] */
		memcpy(tmp, x + i * dim, dim * sizeof(*tmp));
		blas_daxpy(dim, -1.0, mean, 1, tmp, 1);
		
		/* ptot += p[i] */
		p = mlogit_prob(&m->values, i);
		ptot += p;
		
		/* cov += p[i] * tmp^2 */
		blas_dspr(uplo, dim, p, tmp, 1, cov);
	}
	
	/* cov /= ptot */
	blas_dscal(dim * (dim + 1) / 2, 1.0 / ptot, cov, 1);
}


void _mlogit_glm_check_invariants(const struct mlogit_glm *m)
{
	(void)m;
}


#if 0



void mlogit_mean_update(struct mlogit_mean *m, const struct mlogit *mlogit,
			const double *x1, const double *dx,
			const struct vpattern *ix)
{
	size_t n = m->dim;
	double *buf = m->xbuf;
	double eta0 = mlogit->eta0;
	double eta_max = mlogit->eta_max;
	double phi = mlogit->phi;
	double expm1_deta = mlogit->expm1_deta;
	
	// buf := expm1(deta) * (x1 - mean0)
	memcpy(buf, x1, n * sizeof(*buf));
	blas_daxpy(n, -1.0, m->mean, 1, buf, 1); 
	blas_dscal(n, expm1_deta, buf, 1);
	
	// buf += dx
	if (ix) {
		sblas_daxpyi(1.0, dx, ix, buf);
	} else {
		blas_daxpy(n, 1.0, dx, 1, buf, 1);
	}
	
	blas_daxpy(n, exp(eta0 - eta_max - phi), buf, 1, m->mean, 1);
}

#endif


