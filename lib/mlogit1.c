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

double mlogit_deta(const struct mlogit1 *m1, size_t i)
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
	double deta = mlogit_deta(m1, i);
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
	double deta = mlogit_deta(m1, i);
	double dpsi = mlogit_dpsi(m1);
	double lp = lp0 + (deta - dpsi);
	return lp;
}

double mlogit_dpsi(const struct mlogit1 *m1)
{
	size_t iz, nz = m1->nz;
	double sum = 0.0;
	double sum1 = 0.0;

	for (iz = 0; iz < nz; iz++) {
		double eta = mlogit_lprob(m1->parent, m1->ind[iz]);
		double deta = m1->deta[iz];
		double eta1 = eta + deta;
		double w1 = exp(eta1);
		double w = exp(eta);
		sum += w;
		sum1 += w1;
	}

	double dpsi = log((1 - sum) + sum1);

	if (!isfinite(dpsi))
		dpsi = get_dpsi_safe(m1);

	return dpsi;
}

double get_dpsi_safe(const struct mlogit1 *m1)
{
	size_t nz = m1->nz;
	double dpsi = NAN;

	if (!nz)
		return 0.0;

	size_t iz;
	double max = -INFINITY;
	double max1 = -INFINITY;
	double sum1 = 0.0;
	double sum = 0.0;

	for (iz = 0; iz < nz; iz++) {
		double eta = mlogit_lprob(m1->parent, m1->ind[iz]);
		double deta = m1->deta[iz];
		double eta1 = eta + deta;

		max = MAX(max, eta);
		max1 = MAX(max1, eta1);
	}

	double shift = MIN(max, max1);

	for (iz = 0; iz < nz; iz++) {
		double eta = mlogit_lprob(m1->parent, m1->ind[iz]);
		double deta = m1->deta[iz];
		double eta1 = eta + deta;

		sum1 += exp(eta1 - shift);
		sum += exp(eta);
	}

	dpsi = shift + log(exp(-shift) * (1 - sum) + sum1);

	if (!isfinite(dpsi))
		dpsi = get_dpsi_safer(m1);

	return dpsi;
}

double get_dpsi_safer(const struct mlogit1 *m1)
{
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
	double dpsi = mlogit_dpsi(m1);
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
