#include "port.h"
#include <assert.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "ieee754.h"
#include "timsort.h"
#include "xalloc.h"
#include "mlogit.h"

static int pdouble_rcompar(const void *x, const void *y);
static void sort_eta(struct mlogit *m);
static void set_eta_tail(struct mlogit *m);


void mlogit_init(struct mlogit *m, size_t ncat)
{
	m->ncat = ncat;
	
	m->eta = xmalloc(ncat * sizeof(*m->eta));
	m->eta_order = xmalloc(ncat * sizeof(*m->eta_order));
	m->eta_rank = xmalloc(ncat * sizeof(*m->eta_rank));
	
	mlogit_clear(m);
}

void mlogit_clear(struct mlogit *m)
{
	size_t i, n = m->ncat;
	memset(m->eta, 0, n * sizeof(*m->eta));
	
	for (i = 0; i < n; i++) {
		m->eta_order[i] = i;
		m->eta_rank[i] = i;
	}
	
	m->eta_max = 0;
	m->eta_tail = n == 0 ? NAN : n - 1;
	m->phi_shift = log(n);
}

void mlogit_deinit(struct mlogit *m)
{
	free(m->eta_rank);
	free(m->eta_order);
	free(m->eta);
}

void mlogit_set_all_eta(struct mlogit *m, const double *eta)
{
	size_t n = mlogit_ncat(m);

#ifndef NDEBUG
	size_t i;
	for (i = 0; i < n; i++) {
		assert(!isnan(eta[i]));
	}
#endif

	memcpy(m->eta, eta, n * sizeof(*m->eta));
	
	
	sort_eta(m);
	set_eta_tail(m);
	m->phi_shift = log1p(m->eta_tail);
}




int pdouble_rcompar(const void *x, const void *y)
{
	return double_rcompare(*(const double **)x, *(const double **)y);
}


void sort_eta(struct mlogit *m)
{
	if (m->ncat == 0)
		return;

	size_t i, n = m->ncat;
	const double *eta = m->eta;
	const double **peta = xmalloc(n * sizeof(*peta));
	size_t *restrict order = m->eta_order;
	size_t *restrict rank = m->eta_rank;
	
	for (i = 0; i < n; i++) {
		peta[i] = eta + i;
	}
	
	if (timsort(peta, n, sizeof(*peta), pdouble_rcompar)) {
		// fallback to qsort if memory allocation fails
		qsort(peta, n, sizeof(*peta), pdouble_rcompar);
	}
	
	for (i = 0; i < n; i++) {
		order[i] = peta[i] - eta;
		rank[order[i]] = i;
	}
	
	m->eta_max = m->eta[order[0]];
	
	free(peta);
}


void set_eta_tail(struct mlogit *m)
{
	if (m->ncat == 0)
		return;

	size_t *order = m->eta_order;
	double eta_max = m->eta_max;
	double eta_tail1 = 0;
	
	size_t i, n = m->ncat;
	
	assert(n > 0);
	for (i = n - 1; i > 0; i--) {
		eta_tail1 += exp(m->eta[order[i]] - eta_max);
	}
	
	m->eta_tail = eta_tail1;
}

#if 0


#define IS_ROOT(i) (rank[i] == 0)

#define HAS_PARENT(i) (rank[i] != 0)
#define PARENT(i) (order[rank[i] / 2])

#define HAS_LEFT(i) (rank[i] < (n - 1) / 2)
#define LEFT(i) (order[2 * rank[i] + 1])

#define HAS_RIGHT(i) (rank[i] < n / 2 - 1)
#define RIGHT(i) (order[2 * rank[i] + 2])

#define SWAP(i,j) \
	tmp = order[rank[i]]; \
	order[rank[i]] = order[rank[j]]; \
	order[rank[j]] = tmp; \
	tmp = rank[i]; \
	rank[i] = rank[j]; \
	rank[j] = tmp;


static void update_eta(struct mlogit *m, size_t i, double deta)
{
	double eta1 = m->eta[i] + deta;
	size_t j, left, right, tmp;
	size_t n = m->n;
	size_t *restrict order = m->eta_order;
	size_t *restrict rank = m->eta_rank;
	
	m->eta[i] = eta1;
	
	if (deta >= 0) {
		// swap with parent until heap property is restored
		while (HAS_PARENT(i)) {
			j = PARENT(i);
			if (m->eta[i] <= m->eta[j])
				break;
			SWAP(i,j);
		}
	} else {
		// swap with greatest child until heap property is restored
		while (HAS_LEFT(i)) {
			left = LEFT(i);
			
			if (HAS_RIGHT(i)) {
				right = RIGHT(i);
				// rightmost child has greatest key
				if (m->eta[left] <= m->eta[right]) {
					if (m->eta[i] >= m->eta[right])
						break;
					SWAP(i, right);
				}
			}
			// leftmost child has greatest val
			if (m->eta[i] >= m->eta[left])
				break;
			SWAP(i, left);
		}
	}
}

void mlogit_update(struct mlogit *m, size_t i, double deta)
{
	size_t *order = m->eta_order;
	size_t *rank = m->eta_rank;
	double eta = m->eta[i];
	double eta1 = eta + deta;
	double eta_max = m->eta_max;
	double eta_tail = m->eta_tail;
	double expm1_deta = expm1(deta);
	
	if (IS_ROOT(i)) {
		update_eta(m, i, deta);
		if (IS_ROOT(i)) {
			m->eta_max = eta1;
			m->eta_tail = eta_tail * exp(-deta);
		} else {
			size_t root = order[0];
			m->eta_max = m->eta[root];
			compute_eta_tail(m);
		}
	} else {
		update_eta(m, i, deta);
		if (IS_ROOT(i)) {
			m->eta_max = eta1;
			compute_eta_tail(m);
		} else {
			m->eta_tail = eta_tail + exp(eta - eta_max) * expm1_deta;
		}
	}
	m->phi = log1p(m->eta_tail);
	
	// record parameters of update
	m->eta0 = eta;
	m->deta = deta;
	m->expm1_deta = expm1_deta;
}

#undef IS_ROOT
#undef HAS_PARENT
#undef PARENT
#undef HAS_LEFT
#undef LEFT
#undef HAS_RIGHT
#undef RIGHT
#undef SWAP


void mlogit_mean_init(struct mlogit_mean *m, size_t dim, const double *mean0)
{
	m->dim = dim;
	m->mean = xmalloc(dim * sizeof(*m->mean));
	m->xbuf = xmalloc(dim * sizeof(*m->xbuf));
	
	memcpy(m->mean, mean0, dim * sizeof(*m->mean));
}


void mlogit_mean_deinit(struct mlogit_mean *m)
{
	free(m->xbuf);
	free(m->mean);
}


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


