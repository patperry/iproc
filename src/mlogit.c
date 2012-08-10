#include "port.h"
#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "coreutil.h" // MAX
#include "ieee754.h"
#include "timsort.h"
#include "xalloc.h"
#include "mlogit.h"

static int pdouble_rcompar(const void *x, const void *y);
static void replace_eta(struct mlogit *m, size_t i, double eta1);
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
	m->eta_tail_err = 0;
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


#define IS_ROOT(i) (rank[i] == 0)

#define HAS_PARENT(i) (rank[i] != 0)
#define PARENT(i) (order[(rank[i] - 1) / 2])

#define HAS_LEFT(i) (rank[i] < n / 2)
#define LEFT(i) (order[2 * rank[i] + 1])

#define HAS_RIGHT(i) (rank[i] < (n - 1) / 2)
#define RIGHT(i) (order[2 * rank[i] + 2])

#define SWAP(i,j) \
	tmp = order[rank[i]]; \
	order[rank[i]] = order[rank[j]]; \
	order[rank[j]] = tmp; \
	tmp = rank[i]; \
	rank[i] = rank[j]; \
	rank[j] = tmp;


void mlogit_set_eta(struct mlogit *m, size_t i, double eta1)
{
	assert(i < mlogit_ncat(m));
	assert(!isnan(eta1));

	size_t *restrict order = m->eta_order;
	size_t *restrict rank = m->eta_rank;
	double eta = m->eta[i];
	double deta = eta1 - eta;
	double eta_max = m->eta_max;
	double eta_tail = m->eta_tail;
	// double expm1_deta = expm1(deta);
	
	if (IS_ROOT(i)) {
		replace_eta(m, i, eta1);
		if (IS_ROOT(i)) {
			m->eta_max = eta1;
			m->eta_tail = eta_tail * exp(-deta);
			if (isnan(m->eta_tail) || m->eta_max == INFINITY) // overflow
				set_eta_tail(m);
			assert(m->eta_tail >= 0);
		} else {
			size_t root = order[0];
			m->eta_max = m->eta[root];
			set_eta_tail(m);
		}
	} else {
		replace_eta(m, i, eta1);
		if (IS_ROOT(i)) {
			m->eta_max = eta1;
			set_eta_tail(m);
		} else {
			// m->eta_tail = eta_tail + exp(eta - eta_max) * expm1_deta;
			m->eta_tail = (eta_tail - exp(eta - eta_max)) + exp(eta1 - eta_max);
			if (!(m->eta_tail >= 0)) {
				m->eta_tail = 0;
			}
			assert(m->eta_tail >= 0);
		}
	}
	m->phi_shift = log1p(m->eta_tail);
}



int pdouble_rcompar(const void *x, const void *y)
{
	return double_rcompare(*(const double **)x, *(const double **)y);
}



static void replace_eta(struct mlogit *m, size_t i, double eta1)
{
	double eta0 = m->eta[i];
	size_t j, left, right, tmp;
	size_t n = m->ncat;
	size_t *restrict order = m->eta_order;
	size_t *restrict rank = m->eta_rank;
	const double *eta = m->eta;
	
	m->eta[i] = eta1;
	
	if (eta1 >= eta0) {
		// swap with parent until heap property is restored
		while (HAS_PARENT(i)) {
			j = PARENT(i);
			if (eta[i] <= eta[j])
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
				if (eta[left] <= eta[right]) {
					if (eta[i] >= eta[right])
						break;
					SWAP(i, right);
					continue;
				}
			}
			// leftmost child has greatest val
			if (eta[i] >= eta[left])
				break;
			SWAP(i, left);
		}
	}
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


static double twosum(double a, double b, double *err)
{
	double a1, b1, da, db;
	double res = a + b;
	a1 = res - b;
	b1 = res - a1;
	da = a - a1;
	db = b - b1;
	*err = da + db;
	return res;
}

static double fast_twosum(double a, double b, double *err)
{
	double res = a + b;
	double z = res - a;
	*err = b - z;
	return res;
}

static double exp_e(double x, double *err)
{
	double res = exp(x);
	*err = 2.0 * DBL_EPSILON * fabs(res);
	return res;
}

static double exp_err_e(double x, double dx, double *err)
{
	double adx = fabs(dx);
	double ex  = exp(x);
	double edx = exp(adx);
	double res  = ex;
	*err  = ex * MAX(DBL_EPSILON, edx - 1.0/edx);
	*err += 2.0 * DBL_EPSILON * fabs(res);
	return res;
}

struct sum_t {
	size_t n;
	double val;
	double comp;
	double acomp;
	double err;
};

void sum_init(struct sum_t *s)
{
	s->n = 0;
	s->val = 0;
	s->comp = 0;
	s->acomp = 0;
	s->err = 0;
}

void sum_add(struct sum_t *s, double x)
{
	double err;
	s->n++;
	s->val = twosum(s->val, x, &err);
	s->comp += err;
	s->acomp += fabs(err);
}

/* Ogita, Rump, & Oishi (2005). "Accurate sum and dot product." SIAM Journal
 * on Scientific Computing (SISC), 26(6):1955-1988.
 *
 * Cor. 4.7.
 */
double sum_get(const struct sum_t *s, double *err)
{
	assert((s->n) < 1.0 / DBL_EPSILON);

	double eps = DBL_EPSILON / 2;
	double eta = DBL_MIN * DBL_EPSILON;
	double res = s->val + s->comp;
	double ares = fabs(res);

	double eps2 = eps * eps;
	double two_n_eps = 2 * s->n * eps;
	double beta = two_n_eps / (1 - two_n_eps) * s->acomp;
	*err = eps * ares + (beta + (2 * eps2 * ares + 3 * eta));
	return res;
}

void set_eta_tail(struct mlogit *m)
{
	if (m->ncat == 0)
		return;

	size_t *order = m->eta_order;
	double eta_max = m->eta_max;
	double diff, diff_err;
	double ediff, ediff_err;
	struct sum_t eta_tail, eta_tail_err;
	double extra;
	
	sum_init(&eta_tail);
	sum_init(&eta_tail_err);
	
	size_t i, n = m->ncat;
	
	assert(n > 0);
	for (i = n - 1; i > 0; i--) {
		// eta_tail1 += exp(m->eta[order[i]] - eta_max);
		diff = twosum(m->eta[order[i]], -eta_max, &diff_err);
		ediff = exp_err_e(diff, diff_err, &ediff_err);
		sum_add(&eta_tail, ediff);
		sum_add(&eta_tail_err, ediff_err);
	}
	
	m->eta_tail = sum_get(&eta_tail, &m->eta_tail_err);
	m->eta_tail_err += sum_get(&eta_tail_err, &extra);
	m->eta_tail_err += extra;
	
	printf("\neta_tail = %.10e +/- %.10e ", m->eta_tail, m->eta_tail_err);
}



void _mlogit_check_invariants(const struct mlogit *m)
{
	size_t i, n = m->ncat;
	const double *eta = m->eta;
	const size_t *order = m->eta_order;
	const size_t *rank = m->eta_rank;
	
	for (i = 0; i < n; i++) {
		assert(order[rank[i]] == i);
		assert(rank[order[i]] == i);
	}
	
	for (i = 0; i < n; i++) {
		if (HAS_LEFT(i))
			assert(eta[i] >= eta[LEFT(i)]);
		if (HAS_RIGHT(i))
			assert(eta[i] >= eta[RIGHT(i)]);
	}
	
	for (i = 0; i < n; i++) {
		assert(m->eta_max >= eta[i]);
	}
	
	if (m->ncat > 0) {
		assert(!isnan(m->eta_tail));
		assert(m->eta_tail != INFINITY);
		// assert(m->phi_shift == log1p(m->eta_tail));
	}
}


	

#if 0


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


