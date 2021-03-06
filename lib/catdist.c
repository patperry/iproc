#include "port.h"
#include <assert.h>
#include <float.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "coreutil.h"
#include "ieee754.h"		// double_rcompare
#include "timsort.h"
#include "xalloc.h"
#include "catdist.h"

#define ETA (DBL_MIN * DBL_EPSILON)
#define EPS (DBL_EPSILON / 2)

#define CHECK(x) \
	do { \
		fail = !(x); \
		assert(!fail); \
		if (fail) \
			goto out; \
	} while(0)

struct valerr_t {
	double val;
	double err;
};

static struct valerr_t twosum(double a, double b);
static struct valerr_t mult_err_e(struct valerr_t x, struct valerr_t y);
static struct valerr_t exp_err_e(struct valerr_t x);
static struct valerr_t exp_sum_e(double x, double y);

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

static int pdouble_rcompar(const void *x, const void *y);

static void clear(struct catdist *c);
static struct valerr_t compute_eta_tail(const struct catdist *c);
static void replace_eta(struct catdist *c, size_t i, double eta1);
static void sort_eta(struct catdist *c);
static void set_eta_tail(struct catdist *c);
static void update_version(struct catdist *c);

void catdist_init(struct catdist *c, size_t ncat)
{
	c->ncat = ncat;

	c->eta = xmalloc(ncat * sizeof(*c->eta));
	c->eta_order = xmalloc(ncat * sizeof(*c->eta_order));
	c->eta_rank = xmalloc(ncat * sizeof(*c->eta_rank));
	c->tol = sqrt(EPS) * sqrt(sqrt(EPS));
	version_init(&c->version);

	clear(c);
}

void clear(struct catdist *c)
{
	size_t i, n = c->ncat;
	memset(c->eta, 0, n * sizeof(*c->eta));

	for (i = 0; i < n; i++) {
		c->eta_order[i] = i;
		c->eta_rank[i] = i;
	}

	c->eta_max = 0;
	c->eta_tail = n == 0 ? NAN : n - 1;
	c->eta_tail_err = 0;
	c->psi_shift = log(n);
}

void catdist_deinit(struct catdist *c)
{
	version_deinit(&c->version);
	free(c->eta_rank);
	free(c->eta_order);
	free(c->eta);
}


void update_version(struct catdist *c)
{
	version_update(&c->version);
}


void catdist_set_all_eta(struct catdist *c, const double *eta)
{
	if (!eta) {
		clear(c);
		return;
	}

	size_t n = catdist_ncat(c);

#ifndef NDEBUG
	size_t i;
	for (i = 0; i < n; i++) {
		assert(!isnan(eta[i]));
	}
#endif

	memcpy(c->eta, eta, n * sizeof(*c->eta));

	sort_eta(c);
	set_eta_tail(c);
	c->psi_shift = log1p(c->eta_tail);
	update_version(c);
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

void catdist_set_eta(struct catdist *c, size_t i, double eta1)
{
	assert(i < catdist_ncat(c));
	assert(!isnan(eta1));

	size_t *restrict order = c->eta_order;
	size_t *restrict rank = c->eta_rank;
	const double eta = c->eta[i];
	const double eta_max = c->eta_max;
	struct valerr_t eta_tail = { c->eta_tail, c->eta_tail_err };
	double eta_max1;
	struct valerr_t eta_tail1;

	if (IS_ROOT(i)) {
		replace_eta(c, i, eta1);
		if (IS_ROOT(i)) {
			eta_max1 = eta1;
			// eta_tail1 = c->eta_tail * exp(-(eta1 - eta));

			struct valerr_t endeta = exp_sum_e(eta, -eta1);
			eta_tail1 = mult_err_e(eta_tail, endeta);
		} else {
			eta_max1 = c->eta[order[0]];
			// eta_tail1 = eta_tail * exp(eta_max - eta_max1) + exp(eta1 - eta_max1) - 1;

			struct valerr_t endeta_max =
			    exp_sum_e(eta_max, -eta_max1);
			struct valerr_t eta_tail_scale =
			    mult_err_e(eta_tail, endeta_max);
			struct valerr_t ediff = exp_sum_e(eta1, -eta_max1);

			eta_tail1.val = eta_tail_scale.val + ediff.val - 1;
			eta_tail1.err =
			    3 * EPS / (1 - 3 * EPS) * eta_tail1.val + ETA;
			eta_tail1.err +=
			    (eta_tail_scale.err + ediff.err) / (1 - 2 * EPS);

			c->eta_tail = eta_tail1.val;
			c->eta_tail_err = eta_tail1.err;
		}
	} else {
		replace_eta(c, i, eta1);
		if (IS_ROOT(i)) {
			eta_max1 = eta1;
			// eta_tail1 = (eta_tail - exp(eta - eta_max) + 1) * exp(eta_max - eta_max1) 

			struct valerr_t ediff = exp_sum_e(eta, -eta_max);
			struct valerr_t endeta_max =
			    exp_sum_e(eta_max, -eta_max1);
			struct valerr_t tail_adj;

			tail_adj.val = (eta_tail.val - ediff.val);
			tail_adj.err = abs(tail_adj.val);

			tail_adj.val += 1.0;
			tail_adj.err += abs(tail_adj.val);

			tail_adj.err =
			    EPS * tail_adj.err + (eta_tail.err + ediff.err);
			tail_adj.err = tail_adj.err / (1 - 5 * EPS) + 2 * ETA;

			eta_tail1 = mult_err_e(tail_adj, endeta_max);
		} else {
			eta_max1 = eta_max;
			//eta_tail1 = (eta_tail - exp(eta - eta_max)) + exp(eta1 - eta_max);

			struct valerr_t ediff = exp_sum_e(eta, -eta_max);
			struct valerr_t ediff1 = exp_sum_e(eta1, -eta_max);

			eta_tail1.val = eta_tail.val - ediff.val;
			eta_tail1.err = fabs(eta_tail1.val);

			eta_tail1.val += ediff1.val;
			eta_tail1.err += abs(eta_tail1.val);

			eta_tail1.err =
			    EPS * eta_tail1.err + (eta_tail.err + ediff.err +
						   ediff1.err);
			eta_tail1.err = eta_tail1.err / (1 - 5 * EPS) + 2 * ETA;
		}
	}

	c->eta_max = eta_max1;
	c->eta_tail = eta_tail1.val;
	c->eta_tail_err = eta_tail1.err;

	if (!(c->eta_tail_err < c->tol * (1 + c->eta_tail))) {
		//assert(isnan(c->eta_tail_err));
		set_eta_tail(c);
	} else if (!(c->eta_tail >= 0)) {
		assert(c->eta_tail >= -c->eta_tail_err);
		c->eta_tail_err =
		    (c->eta_tail_err + c->eta_tail) / (1 - DBL_EPSILON);
		c->eta_tail = 0;
	}
	//printf("\neta_tail = %.10e +/- %.10e ", c->eta_tail, c->eta_tail_err);
	assert(c->eta_tail >= 0);
	c->psi_shift = log1p(c->eta_tail);
	update_version(c);
}

int pdouble_rcompar(const void *x, const void *y)
{
	return double_rcompare(*(const double **)x, *(const double **)y);
}

static void replace_eta(struct catdist *c, size_t i, double eta1)
{
	double eta0 = c->eta[i];
	size_t j, left, right, tmp;
	size_t n = c->ncat;
	size_t *restrict order = c->eta_order;
	size_t *restrict rank = c->eta_rank;
	const double *eta = c->eta;

	c->eta[i] = eta1;

	if (eta1 >= eta0) {
		// swap with parent until heap property is restored
		while (HAS_PARENT(i)) {
			j = PARENT(i);
			if (eta[i] <= eta[j])
				break;
			SWAP(i, j);
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

void sort_eta(struct catdist *c)
{
	if (c->ncat == 0)
		return;

	size_t i, n = c->ncat;
	const double *eta = c->eta;
	const double **peta = xmalloc(n * sizeof(*peta));
	size_t *restrict order = c->eta_order;
	size_t *restrict rank = c->eta_rank;

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

	c->eta_max = c->eta[order[0]];

	free(peta);
}

struct valerr_t twosum(double a, double b)
{
	assert(FLT_ROUNDS == 1);	/* require round-to-nearest */

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

struct valerr_t mult_err_e(struct valerr_t x, struct valerr_t y)
{
	struct valerr_t res;
	double adx = fabs(x.err);
	double ady = fabs(y.err);

	res.val = x.val * y.val;
	res.err = EPS * fabs(res.val);
	res.err += adx * fabs(y.val) + ady * fabs(x.val) + adx * ady;
	res.err = res.err / (1 - 6 * EPS) + 4 * ETA;

	return res;
}

struct valerr_t exp_err_e(struct valerr_t x)
{
	struct valerr_t res;
	double adx = fabs(x.err);
	double ex = exp(x.val);
	double edxm1 = expm1(adx);
	res.val = ex;
	res.err =
	    ((2 * EPS * ex + 3 * ETA) + (ex + ETA) * (edxm1 + ETA)) / (1 -
								       7 * EPS);
	return res;
}

struct valerr_t exp_sum_e(double x, double y)
{
	struct valerr_t s = twosum(x, y);
	return exp_err_e(s);
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

void set_eta_tail(struct catdist *c)
{
	struct valerr_t eta_tail = compute_eta_tail(c);

	c->eta_tail = eta_tail.val;
	c->eta_tail_err = eta_tail.err;

	//printf("\neta_tail = %.10e +/- %.10e ", c->eta_tail, c->eta_tail_err);
	//printf("!"); fflush(stdout); 
}

struct valerr_t compute_eta_tail(const struct catdist *c)
{
	struct valerr_t eta_tail = { NAN, 0 };

	if (c->ncat == 0)
		return eta_tail;

	size_t *order = c->eta_order;
	double eta_max = c->eta_max;
	struct sum_t eta_tail_sum;
	double eta_tail_err = 0.0;

	sum_init(&eta_tail_sum);

	size_t i, n = c->ncat;

	assert(n > 0);
	for (i = n - 1; i > 0; i--) {
		// c->eta_tail += exp(c->eta[order[i]] - eta_max);
		struct valerr_t ediff = exp_sum_e(c->eta[order[i]], -eta_max);
		sum_add(&eta_tail_sum, ediff.val);
		eta_tail_err += ediff.err;
	}

	eta_tail = sum_get(&eta_tail_sum);
	eta_tail.err += eta_tail_err;
	eta_tail.err /= (1 - n * EPS);

	return eta_tail;
}

int catdist_check(const struct catdist *c)
{
	int fail = 0;
	size_t i, n = c->ncat;
	const double *eta = c->eta;
	const size_t *order = c->eta_order;
	const size_t *rank = c->eta_rank;

	for (i = 0; i < n; i++) {
		CHECK(order[rank[i]] == i);
		CHECK(rank[order[i]] == i);
	}

	for (i = 0; i < n; i++) {
		if (HAS_LEFT(i))
			CHECK(eta[i] >= eta[LEFT(i)]);
		if (HAS_RIGHT(i))
			CHECK(eta[i] >= eta[RIGHT(i)]);
	}

	for (i = 0; i < n; i++) {
		CHECK(c->eta_max >= eta[i]);
	}

	if (c->ncat > 0) {
		// assert(c->psi_shift == log1p(c->eta_tail));

		struct valerr_t eta_tail = compute_eta_tail(c);

		if (c->eta_tail < eta_tail.val) {
			// eta_tail + eta_tail_err >= val - err
			CHECK(eta_tail.val <=
			      (c->eta_tail + eta_tail.err +
			       c->eta_tail_err) / (1 - 3 * EPS));
		} else {
			// eta_tail - eta_tail_err <= val + err
			CHECK(c->eta_tail <=
			      (eta_tail.val + eta_tail.err +
			       c->eta_tail_err) / (1 - 3 * EPS));
		}
	}
out:
	return fail;
}
