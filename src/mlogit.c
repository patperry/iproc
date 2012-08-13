#include "port.h"
#include <assert.h>
#include <float.h>
#include <math.h>
//#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "coreutil.h" // MAX
#include "ieee754.h"
#include "timsort.h"
#include "xalloc.h"
#include "mlogit.h"

#define ETA (DBL_MIN * DBL_EPSILON)
#define EPS (DBL_EPSILON / 2)

static int pdouble_rcompar(const void *x, const void *y);
static void replace_eta(struct mlogit *m, size_t i, double eta1);
static void sort_eta(struct mlogit *m);
static void set_eta_tail(struct mlogit *m);

struct valerr_t {
	double val;
	double err;
};

static struct valerr_t twosum(double a, double b);
static struct valerr_t twosum_err(double a, struct valerr_t b);
static struct valerr_t fast_twosum(double a, double b);
static struct valerr_t twomult(double a, double b);
static struct valerr_t twomult_err(double a, struct valerr_t b);
static struct valerr_t exp_e(double x);
static struct valerr_t exp_err_e(struct valerr_t x);
static struct valerr_t expm1_e(double x);

struct sum_t {
	size_t n;
	double val;
	double comp;
	double acomp;
	double err;
};

static void sum_init(struct sum_t *s);
static void sum_add(struct sum_t *s, double x);
static struct valerr_t sum_get(const struct sum_t *s);



void mlogit_init(struct mlogit *m, size_t ncat)
{
	m->ncat = ncat;
	
	m->eta = xmalloc(ncat * sizeof(*m->eta));
	m->eta_order = xmalloc(ncat * sizeof(*m->eta_order));
	m->eta_rank = xmalloc(ncat * sizeof(*m->eta_rank));
	m->tol = sqrt(EPS) * sqrt(sqrt(EPS));
	
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
	const double eta = m->eta[i];
	const double eta_max = m->eta_max;
	const double eta_tail = m->eta_tail;
	const double eta_tail_err = m->eta_tail_err;
	
	if (IS_ROOT(i)) {
		replace_eta(m, i, eta1);
		if (IS_ROOT(i)) {
			m->eta_max = eta1;
			// m->eta_tail = m->eta_tail * exp(-(eta1 - eta));
			
			struct valerr_t ndeta = twosum(eta, -eta1);
			struct valerr_t endeta = exp_err_e(ndeta);
			struct valerr_t eta_tail1 = twomult_err(eta_tail, endeta);

			m->eta_tail = eta_tail1.val;
			m->eta_tail_err = (eta_tail1.err + (m->eta_tail_err * (endeta.val + endeta.err) + ETA)) / (1 - 4 * EPS);
			
		} else {
			size_t root = order[0];
			double eta_max1 = m->eta[root];
			// m->eta_tail = eta_tail * exp(eta_max - eta_max1) + exp(eta1 - eta_max1) - 1;
			
			struct valerr_t ndeta_max = twosum(eta_max, -eta_max1);
			struct valerr_t endeta_max = exp_err_e(ndeta_max);
			struct valerr_t eta_tail_scale = twomult_err(eta_tail, endeta_max);
			eta_tail_scale.err += eta_tail_err * (endeta_max.val + endeta_max.err) + ETA;
			eta_tail_scale.err /= (1 - 4 * EPS);
			
			struct valerr_t diff = twosum(eta1, -eta_max1);
			struct valerr_t ediff = exp_err_e(diff);
			
			struct valerr_t eta_tail1;
			eta_tail1.val = eta_tail_scale.val + ediff.val - 1;
			eta_tail1.err = 3 * EPS / (1 - 3 * EPS) * eta_tail1.val + ETA;
			eta_tail1.err += (eta_tail_scale.err + ediff.err) / (1 - 2 * EPS);
			
			m->eta_max = eta_max1;
			m->eta_tail = eta_tail1.val;
			m->eta_tail_err = eta_tail1.err;
		}
	} else {
		replace_eta(m, i, eta1);
		if (IS_ROOT(i)) {
			double eta_max1 = eta1;
			
			// m->eta_tail = (eta_tail - exp(eta - eta_max) + 1) * exp(eta_max - eta_max1) 
			
			struct valerr_t diff = twosum(eta, -eta_max);
			struct valerr_t ediff = exp_err_e(diff);

			struct valerr_t tail_adj;
			tail_adj.val = (eta_tail - ediff.val) + 1;
			tail_adj.err = 3 * EPS / (1 - 3 * EPS) * tail_adj.val + ETA;
			tail_adj.err += (eta_tail_err + ediff.err) / (1 + 3 * EPS);
			
			struct valerr_t ndeta_max = twosum(eta_max, -eta_max1);
			struct valerr_t endeta_max = exp_err_e(ndeta_max);
			struct valerr_t eta_tail1 = twomult_err(tail_adj.val, endeta_max);

			m->eta_max = eta_max1;
			m->eta_tail = eta_tail1.val;
			m->eta_tail_err = (m->eta_tail_err * (endeta_max.val + endeta_max.err) + ETA
					   + eta_tail1.err) / (1 - 4 * EPS);
		} else {
			//m->eta_tail = (eta_tail - exp(eta - eta_max)) + exp(eta1 - eta_max);
			struct valerr_t diff = twosum(eta, -eta_max);
			struct valerr_t ediff = exp_err_e(diff);
			struct valerr_t diff1 = twosum(eta1, -eta_max);
			struct valerr_t ediff1 = exp_err_e(diff1);
			
			struct sum_t eta_tail_sum;
			struct sum_t eta_tail_err_sum;
			
			sum_init(&eta_tail_sum);
			sum_add(&eta_tail_sum, m->eta_tail);
			
			sum_init(&eta_tail_err_sum);
			sum_add(&eta_tail_err_sum, m->eta_tail_err);

			sum_add(&eta_tail_sum, -(ediff.val));
			sum_add(&eta_tail_err_sum, ediff.err);
			sum_add(&eta_tail_sum, ediff1.val);
			sum_add(&eta_tail_err_sum, ediff1.err);
			
			struct valerr_t eta_tail = sum_get(&eta_tail_sum);
			struct valerr_t eta_tail_err = sum_get(&eta_tail_err_sum);
			
			m->eta_tail = eta_tail.val;
			m->eta_tail_err = (eta_tail.err + eta_tail_err.err + eta_tail_err.val) / (1 - 3 * EPS);
		}
	}
	
	
	
	if (!(m->eta_tail_err < m->tol * (1 + m->eta_tail))) {
		//assert(isnan(m->eta_tail_err));
		set_eta_tail(m);
	} else if (!(m->eta_tail >= 0)) {
		assert(m->eta_tail >= - m->eta_tail_err);
		m->eta_tail_err = (m->eta_tail_err + m->eta_tail) / (1 - DBL_EPSILON);
		m->eta_tail = 0;
	}
	
	//printf("\neta_tail = %.10e +/- %.10e ", m->eta_tail, m->eta_tail_err);
	assert(m->eta_tail >= 0);	
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


struct valerr_t twosum(double a, double b)
{
	struct valerr_t res;
	double a1, b1, da, db;
	res.val = a + b;
	a1 = res.val - b;
	b1 = res.val - a1;
	da = a - a1;
	db = b - b1;
	res.err = da + db;
	return res;
}

struct valerr_t fast_twosum(double a, double b)
{
	struct valerr_t res;
	res.val = a + b;
	double z = res.val - a;
	res.err = b - z;
	return res;
}

static struct valerr_t twosum_err(double a, struct valerr_t b)
{
	struct valerr_t res = twosum(a, b.val);
	res.err = (res.err + b.err) / (1 - 2 * EPS);
	return res;
}

struct valerr_t twomult(double a, double b)
{
	struct valerr_t res;
	res.val = a * b;
	res.err = fma(a, b, -res.val);
	return res;
}

static struct valerr_t twomult_err(double a, struct valerr_t b)
{
	struct valerr_t res = twomult(a, b.val);
	res.err = fma(a, b.err, res.err) / (1 - 2 * EPS);
	return res;
}

struct valerr_t exp_e(double x)
{
	struct valerr_t res;
	res.val = exp(x);
	res.err = 2.0 * EPS * res.val + ETA;
	return res;
}


struct valerr_t exp_err_e(struct valerr_t x)
{
	struct valerr_t res;
	double adx = fabs(x.err);
	double ex  = exp(x.val);
	double edx = exp(adx);
	res.val = ex;
	res.err  = 2.0 * ((ex + ETA) * (edx + EPS) + ETA);
	return res;
}


static struct valerr_t expm1_e(double x)
{
	struct valerr_t res;
	res.val = expm1(x);
	res.err = 2.0 * EPS * fabs(res.val) + ETA;
	return res;
}


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
	struct valerr_t t;
	s->n++;
	t = twosum(s->val, x);
	s->val = t.val;
	s->comp += t.err;
	s->acomp += fabs(t.err);
}

/* Ogita, Rump, & Oishi (2005). "Accurate sum and dot product." SIAM Journal
 * on Scientific Computing (SISC), 26(6):1955-1988.
 *
 * Cor. 4.7.
 */
struct valerr_t sum_get(const struct sum_t *s)
{
	assert((s->n) < 1.0 / DBL_EPSILON);

	struct valerr_t res;
	double eps2 = EPS * EPS;
	double two_n_eps = 2 * s->n * EPS;
	
	res.val = s->val + s->comp;
	double ares = fabs(res.val);
	double beta = two_n_eps / (1 - two_n_eps) * s->acomp;
	res.err = EPS * ares + (beta + (2 * eps2 * ares + 3 * ETA));
	return res;
}

void set_eta_tail(struct mlogit *m)
{
	if (m->ncat == 0)
		return;

	size_t *order = m->eta_order;
	double eta_max = m->eta_max;
	struct sum_t eta_tail_sum, eta_tail_err_sum;
	
	sum_init(&eta_tail_sum);
	sum_init(&eta_tail_err_sum);
	
	size_t i, n = m->ncat;
	
	assert(n > 0);
	for (i = n - 1; i > 0; i--) {
		// m->eta_tail += exp(m->eta[order[i]] - eta_max);
		struct valerr_t diff = twosum(m->eta[order[i]], -eta_max);
		struct valerr_t ediff = exp_err_e(diff);
		sum_add(&eta_tail_sum, ediff.val);
		sum_add(&eta_tail_err_sum, ediff.err);
	}
	
	struct valerr_t eta_tail = sum_get(&eta_tail_sum);
	struct valerr_t eta_tail_err = sum_get(&eta_tail_err_sum);
	m->eta_tail = eta_tail.val;
	m->eta_tail_err = (eta_tail.err + eta_tail_err.val + eta_tail_err.err) / (1 - 3 * EPS);
	
	//printf("\neta_tail = %.10e +/- %.10e ", m->eta_tail, m->eta_tail_err);
	//printf("!"); fflush(stdout);
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
		
		struct sum_t eta_tail_sum;
		struct sum_t eta_tail_err_sum;
		
		sum_init(&eta_tail_sum);
		sum_init(&eta_tail_err_sum);
		
		for (i = 0; i < n; i++) {
			if (IS_ROOT(i))
				continue;
			    
			struct valerr_t diff = twosum(eta[i], -(m->eta_max));
			struct valerr_t ediff = exp_err_e(diff);
				
			sum_add(&eta_tail_sum, ediff.val);
			sum_add(&eta_tail_err_sum, ediff.err);
		}
		
		double eps = DBL_EPSILON / 2;
		struct valerr_t eta_tail = sum_get(&eta_tail_sum);
		struct valerr_t eta_tail_err = sum_get(&eta_tail_err_sum);
		double val = eta_tail.val;
		double err = (eta_tail_err.err + eta_tail.err + eta_tail_err.val) / (1 - 3 * eps);
		
		if (m->eta_tail < val) {
			// eta_tail + eta_tail_err >= val - err
			assert(val <= (m->eta_tail + err + m->eta_tail_err) / (1 - 3 * eps));
		} else {
			// eta_tail - eta_tail_err <= val + err
			assert(m->eta_tail <= (val + err + m->eta_tail_err) / (1 - 3 * eps));
		}

		
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


