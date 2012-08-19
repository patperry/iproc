#include "port.h"
#include <stdint.h> // SIZE_MAX
#include <stdlib.h> // free
#include <string.h> // memcpy, memset
#include "blas.h"   // struct dmatrix, blas_gemv
#include "xalloc.h" // xcalloc
#include "mlogit_glm.h"

static void recompute(struct mlogit_glm *m);


void mlogit_glm_init(struct mlogit_glm *m, size_t ncat, size_t dim)
{
	assert(ncat == 0 || dim <= SIZE_MAX / ncat);
	
	mlogit_init(&m->values, ncat);
	m->x = xcalloc(ncat * dim, sizeof(*m->x));
	m->beta = xcalloc(dim, sizeof(*m->beta));
	m->cat_buf = xcalloc(ncat, sizeof(*m->cat_buf));
	m->mean = xcalloc(dim, sizeof(*m->mean));
	m->dim = dim;
	
}

void mlogit_glm_deinit(struct mlogit_glm *m)
{
	free(m->mean);
	free(m->cat_buf);
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
}

void mlogit_glm_set_x(struct mlogit_glm *m, const double *x)
{
	size_t len = mlogit_glm_ncat(m) * mlogit_glm_dim(m) * sizeof(*m->x);
	if (x) {
		memcpy(m->x, x, len);
	} else {
		memset(m->x, 0, len);
	}
	
	recompute(m);
}

void mlogit_glm_inc_x(struct mlogit_glm *m, size_t i, const double *dx, const size_t *idx, size_t ndx)
{
	size_t j;
	
	assert(i < mlogit_glm_ncat(m));
	assert(dx || ndx == 0);
	assert(ndx < mlogit_glm_dim(m));
	
	size_t ncat = mlogit_glm_ncat(m);
	size_t dim = mlogit_glm_dim(m);
	
	if (idx) {
		for (j = 0; j < ndx; j++) {
			assert(idx[j] < ncat);
			MATRIX_ITEM(m->x, ncat, i, idx[j]) += dx[j];
		}
	} else {
		for (j = 0; j < dim; j++) {
			MATRIX_ITEM(m->x, ncat, i, j) = 0;
		}
	}
	
	recompute(m);
}


void recompute(struct mlogit_glm *m)
{
	size_t ncat = mlogit_glm_ncat(m);
	size_t dim = mlogit_glm_dim(m);
	
	if (ncat == 0 || dim == 0)
		return;
	
	const double *beta = m->beta;
	double *mean = m->mean;
	double *tmp = m->cat_buf;
	size_t i;

	// tmp := x * beta
	blas_dgemv(BLAS_NOTRANS, ncat, dim, 1.0, m->x, ncat, beta, 1, 0.0, tmp, 1);
	mlogit_set_all_eta(&m->values, tmp);
	
	// tmp := p
	for (i = 0; i < ncat; i++) {
		tmp[i] = mlogit_prob(&m->values, i);
	}
	
	// mean := t(X) * p
	blas_dgemv(BLAS_TRANS, ncat, dim, 1.0, m->x, ncat, tmp, 1, 0, mean, 1);
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


