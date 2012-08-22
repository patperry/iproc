#include "port.h"
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "coreutil.h"
#include "ieee754.h"
#include "mlogit1.h"
#include "xalloc.h"

#define ETA (DBL_MIN * DBL_EPSILON)
#define EPS (DBL_EPSILON / 2)


#define CHECK(x) \
	do { \
		fail = !(x); \
		assert(!fail); \
		if (fail) \
			goto out; \
	} while(0)

#define CHECK_APPROX(x, y) \
	CHECK(double_eqrel((x), (y)) >= DBL_MANT_DIG / 2 \
	      || (fabs(x) <= ROOT3_DBL_EPSILON && fabs(y) <= ROOT3_DBL_EPSILON))


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



static double get_deta(const struct mlogit1 *m1, size_t i);
static double get_dpsi(const struct mlogit1 *m1);
static double get_dpsi_safe(const struct mlogit1 *m1);
static double get_dpsi_safer(const struct mlogit1 *m1);
static void grow_deta_array(struct mlogit1 *m1, size_t delta);
static size_t find_ind(const struct mlogit1 *m1, size_t i);
static size_t search_ind(struct mlogit1 *m1, size_t i);


void mlogit1_init(struct mlogit1 *m1, const struct mlogit *parent)
{
	assert(parent);

	m1->parent = parent;
	m1->deta = NULL;
	m1->ind = NULL;
	m1->nz = 0;
	m1->nzmax = 0;

	mlogit1_clear(m1);
}


void mlogit1_deinit(struct mlogit1 *m1)
{
	free(m1->ind);
	free(m1->deta);
}

void mlogit1_clear(struct mlogit1 *m1)
{
	m1->nz = 0;

}

double get_deta(const struct mlogit1 *m1, size_t i)
{
	assert(i < mlogit1_ncat(m1));

	double deta = 0.0;
	size_t iz = find_ind(m1, i);

	if (iz != m1->nz) {
		deta = m1->deta[iz];
	}

	return deta;
}

double mlogit1_eta(const struct mlogit1 *m1, size_t i)
{
	assert(i < mlogit1_ncat(m1));

	double eta = mlogit_eta(m1->parent, i);
	double deta = get_deta(m1, i);
	return eta + deta;
}

double mlogit1_prob(const struct mlogit1 *m1, size_t i)
{
	assert(i < mlogit1_ncat(m1));
	double lp = mlogit1_lprob(m1, i);
	return exp(lp);
}


double mlogit1_lprob(const struct mlogit1 *m1, size_t i)
{
	assert(i < mlogit1_ncat(m1));
	double lp0 = mlogit_lprob(m1->parent, i);
	double deta = get_deta(m1, i);
	double dpsi = get_dpsi(m1);
	double lp = lp0 + (deta - dpsi);
	return lp;
}


double get_dpsi(const struct mlogit1 *m1)
{
	printf("_"); fflush(stdout);
	
	size_t iz, nz = m1->nz;
	double sum = 0.0;

	for (iz = 0; iz < nz; iz++) {
		double eta = mlogit_lprob(m1->parent, m1->ind[iz]);
		double deta = m1->deta[iz];
		double eta1 = eta + deta;
		double dw = exp(eta1) - exp(eta);

		sum += dw;
	}

	double dpsi = log1p(sum);
	
	if (!isfinite(dpsi))
		dpsi = get_dpsi_safe(m1);
	
	return dpsi;
}


double get_dpsi_safe(const struct mlogit1 *m1)
{
	printf("s"); fflush(stdout);
	
	size_t nz = m1->nz;
	double dpsi = NAN;
	
	if (!nz)
		return 0.0;

	size_t iz;
	double shift = -INFINITY;
	double shift1 = -INFINITY;
	double sum1 = 0.0;
	double sum = 0.0;
	
	for (iz = 0; iz < nz; iz++) {
		double eta = mlogit_lprob(m1->parent, m1->ind[iz]);
		double deta = m1->deta[iz];
		double eta1 = eta + deta;

		shift = MAX(shift, eta);		
		shift1 = MAX(shift1, eta1);
	}

	if (shift1 < shift) {
		printf("1"); fflush(stdout);
		//double sum1 = 0.0;	
		//double sum = 0.0;
		
		for (iz = 0; iz < nz; iz++) {
			double eta = mlogit_lprob(m1->parent, m1->ind[iz]);
			double deta = m1->deta[iz];
			double eta1 = eta + deta;
			
			sum1 += exp(eta1 - shift1);
			sum += exp(eta);
		}
	
		dpsi = shift1 + log(exp(-shift1) * (1 - sum) + sum1);
	} else if (shift1 - shift < log(DBL_MAX)) {
		printf("2"); fflush(stdout);
		//double sum = 0.0;
		
		for (iz = 0; iz < nz; iz++) {
			double eta = mlogit_lprob(m1->parent, m1->ind[iz]);
			double deta = m1->deta[iz];
			double eta1 = eta + deta;
			double dw = exp(eta1 - shift) - exp(eta);
			
			sum += dw;
		}
		
		dpsi = shift + log(exp(-shift) + sum);
	} else {
		dpsi = get_dpsi_safer(m1);
	}
	
	assert(isfinite(dpsi));
	
	return dpsi;
}


double get_dpsi_safer(const struct mlogit1 *m1)
{
	printf("S"); fflush(stdout);
	
	size_t iz, nz = m1->nz;
	size_t i, n = mlogit1_ncat(m1);	
	
	if (!nz || !n)
		return 0.0;
	
	double etamax = -INFINITY;
	size_t imax = 0;
	
	for (i = 0, iz = 0; i < n; i++) {
		double eta = mlogit_eta(m1->parent, i);
		
		if (iz < nz && m1->ind[iz] == i) {
			eta += m1->deta[iz];
			iz++;
		}
		
		if (eta > etamax) {
			etamax = eta;
			imax = i;
		}
	}
	
	double sum = 0.0;

	for (i = 0, iz = 0; i < n; i++) {
		double eta = mlogit_eta(m1->parent, i);
		
		if (iz < nz && m1->ind[iz] == i) {
			eta += m1->deta[iz];
			iz++;
		}
		
		if (i == imax)
			continue;

		sum += exp(eta - etamax);
	}

	double psi = mlogit_psi(m1->parent);
	double psi1 = etamax + log1p(sum);
	double dpsi = psi1 - psi;
	
	assert(isfinite(dpsi));
	return dpsi;
}



double mlogit1_psi(const struct mlogit1 *m1)
{
	double psi = mlogit_psi(m1->parent);
	double dpsi = get_dpsi(m1);
	double psi1 = psi + dpsi;
	return psi1;
}


void mlogit1_set_deta(struct mlogit1 *m1, size_t i, double deta)
{
	assert(i < mlogit1_ncat(m1));
	assert(deta < INFINITY);

	size_t iz = search_ind(m1, i);
	m1->deta[iz] = deta;
}

int _mlogit1_check(const struct mlogit1 *m1)
{
	int fail = 0;
	size_t i, n = mlogit1_ncat(m1);
	size_t iz, nz = m1->nz;
	const size_t *ind = m1->ind;
	double *eta = xmalloc(n * sizeof(*eta));

	for (iz = 1; iz < nz; iz++) {
		CHECK(ind[iz - 1] < ind[iz]);
	}

	for (i = 0; i < n; i++) {
		eta[i] = mlogit_eta(m1->parent, i);
	}

	for (iz = 0; iz < nz; iz++) {
		eta[ind[iz]] += m1->deta[iz];
	}

	for (i = 0; i < n; i++) {
		CHECK(mlogit1_eta(m1, i) == eta[i]);
	}

	double etamax = -INFINITY;
	size_t imax = 0;

	for (i = 0; i < n; i++) {
		if (eta[i] > etamax) {
			etamax = eta[i];
			imax = i;
		}
	}

	double sum = 0.0;
	for (i = 0; i < n; i++) {
		if (i != imax)
			sum += exp(eta[i] - etamax);
	}
	double psi = etamax + log1p(sum);
	CHECK_APPROX(mlogit1_psi(m1), psi);

	for (i = 0; i < n; i++) {
		CHECK_APPROX(mlogit1_lprob(m1, i), eta[i] - psi);
	}

	for (i = 0; i < n; i++) {
		CHECK_APPROX(mlogit1_prob(m1, i), exp(eta[i] - psi));
	}

out:
	free(eta);
	return fail;
}

void grow_deta_array(struct mlogit1 *m1, size_t delta)
{
	size_t nz = m1->nz;
	size_t nz1 = nz + delta;
	size_t nzmax = m1->nzmax;

	if (nz1 <= nzmax)
		return;

	size_t nzmax1 = array_grow(nz, nzmax, delta, SIZE_MAX);
	assert(nzmax1 >= nz1);

	m1->deta = xrealloc(m1->deta, nzmax1 * sizeof(*m1->deta));
	m1->ind = xrealloc(m1->ind, nzmax1 * sizeof(*m1->ind));
	m1->nzmax = nzmax1;
}

size_t search_ind(struct mlogit1 *m1, size_t i)
{
	const size_t *base = m1->ind, *ptr;
	size_t nz;

	for (nz = m1->nz; nz != 0; nz >>= 1) {
		ptr = base + (nz >> 1);
		if (i == *ptr) {
			return ptr - m1->ind;
		}
		if (i > *ptr) {
			base = ptr + 1;
			nz--;
		}
	}

	/* not found */
	size_t iz = base - m1->ind;
	size_t ntail = m1->nz - iz;

	grow_deta_array(m1, 1);
	memmove(m1->deta + iz + 1, m1->deta + iz, ntail * sizeof(*m1->deta));
	memmove(m1->ind + iz + 1, m1->ind + iz, ntail * sizeof(*m1->ind));
	m1->nz++;
	m1->ind[iz] = i;

	return iz;
}

size_t find_ind(const struct mlogit1 *m1, size_t i)
{
	const size_t *base = m1->ind, *ptr;
	size_t nz;

	for (nz = m1->nz; nz != 0; nz >>= 1) {
		ptr = base + (nz >> 1);
		if (i == *ptr) {
			return ptr - m1->ind;
		}
		if (i > *ptr) {
			base = ptr + 1;
			nz--;
		}
	}

	/* not found */
	return m1->nz;
}




struct valerr_t twosum(double a, double b)
{
	assert(FLT_ROUNDS == 1); /* require round-to-nearest */
	
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
	double ex  = exp(x.val);
	double edxm1 = expm1(adx);
	res.val = ex;
	res.err  = ((2 * EPS * ex + 3 * ETA) + (ex + ETA) * (edxm1 + ETA)) / (1 - 7 * EPS);
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

