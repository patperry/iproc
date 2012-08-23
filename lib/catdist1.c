#include "port.h"
#include <float.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "coreutil.h"
#include "ieee754.h"
#include "catdist1.h"
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

static double get_dpsi_safe(const struct catdist1 *c1);
static double get_dpsi_safer(const struct catdist1 *c1);
static void grow_ind_array(struct catdist1 *c1, size_t delta);
static size_t find_ind(const struct catdist1 *c1, size_t i);
static size_t search_ind(struct catdist1 *c1, size_t i);

void catdist1_init(struct catdist1 *c1, const struct catdist *parent)
{
	assert(parent);

	c1->parent = parent;
	c1->deta = NULL;
	c1->ind = NULL;
	c1->nz = 0;
	c1->nzmax = 0;

	catdist1_clear(c1);
}

void catdist1_deinit(struct catdist1 *c1)
{
	free(c1->ind);
	free(c1->deta);
}

void catdist1_clear(struct catdist1 *c1)
{
	c1->nz = 0;
}

double catdist1_deta(const struct catdist1 *c1, size_t i)
{
	assert(i < catdist1_ncat(c1));

	double deta = 0.0;
	size_t iz = find_ind(c1, i);

	if (iz != c1->nz) {
		deta = c1->deta[iz];
	}

	return deta;
}

double catdist1_eta(const struct catdist1 *c1, size_t i)
{
	assert(i < catdist1_ncat(c1));

	double eta = catdist_eta(c1->parent, i);
	double deta = catdist1_deta(c1, i);
	return eta + deta;
}

double catdist1_prob(const struct catdist1 *c1, size_t i)
{
	assert(i < catdist1_ncat(c1));
	double lp = catdist1_lprob(c1, i);
	return exp(lp);
}

double catdist1_lprob(const struct catdist1 *c1, size_t i)
{
	assert(i < catdist1_ncat(c1));
	double lp0 = catdist_lprob(c1->parent, i);
	double deta = catdist1_deta(c1, i);
	double dpsi = catdist1_dpsi(c1);
	double lp = lp0 + (deta - dpsi);
	return lp;
}

double catdist1_dpsi(const struct catdist1 *c1)
{
	size_t iz, nz = c1->nz;
	double sum = 0.0;
	double suc1 = 0.0;

	for (iz = 0; iz < nz; iz++) {
		double eta = catdist_lprob(c1->parent, c1->ind[iz]);
		double deta = c1->deta[iz];
		double eta1 = eta + deta;
		double w1 = exp(eta1);
		double w = exp(eta);
		sum += w;
		suc1 += w1;
	}

	double dpsi = log((1 - sum) + suc1);

	if (!isfinite(dpsi))
		dpsi = get_dpsi_safe(c1);

	return dpsi;
}

double get_dpsi_safe(const struct catdist1 *c1)
{
	size_t nz = c1->nz;
	double dpsi = NAN;

	if (!nz)
		return 0.0;

	size_t iz;
	double max = -INFINITY;
	double max1 = -INFINITY;
	double suc1 = 0.0;
	double sum = 0.0;

	for (iz = 0; iz < nz; iz++) {
		double eta = catdist_lprob(c1->parent, c1->ind[iz]);
		double deta = c1->deta[iz];
		double eta1 = eta + deta;

		max = MAX(max, eta);
		max1 = MAX(max1, eta1);
	}

	double shift = MIN(max, max1);

	for (iz = 0; iz < nz; iz++) {
		double eta = catdist_lprob(c1->parent, c1->ind[iz]);
		double deta = c1->deta[iz];
		double eta1 = eta + deta;

		suc1 += exp(eta1 - shift);
		sum += exp(eta);
	}

	dpsi = shift + log(exp(-shift) * (1 - sum) + suc1);

	if (!isfinite(dpsi))
		dpsi = get_dpsi_safer(c1);

	return dpsi;
}

double get_dpsi_safer(const struct catdist1 *c1)
{
	size_t iz, nz = c1->nz;
	size_t i, n = catdist1_ncat(c1);

	if (!nz || !n)
		return 0.0;

	double etamax = -INFINITY;
	size_t imax = 0;

	for (i = 0, iz = 0; i < n; i++) {
		double eta = catdist_eta(c1->parent, i);

		if (iz < nz && c1->ind[iz] == i) {
			eta += c1->deta[iz];
			iz++;
		}

		if (eta > etamax) {
			etamax = eta;
			imax = i;
		}
	}

	double sum = 0.0;

	for (i = 0, iz = 0; i < n; i++) {
		double eta = catdist_eta(c1->parent, i);

		if (iz < nz && c1->ind[iz] == i) {
			eta += c1->deta[iz];
			iz++;
		}

		if (i == imax)
			continue;

		sum += exp(eta - etamax);
	}

	double psi = catdist_psi(c1->parent);
	double psi1 = etamax + log1p(sum);
	double dpsi = psi1 - psi;

	assert(isfinite(dpsi));
	return dpsi;
}

double catdist1_psi(const struct catdist1 *c1)
{
	double psi = catdist_psi(c1->parent);
	double dpsi = catdist1_dpsi(c1);
	double psi1 = psi + dpsi;
	return psi1;
}

void catdist1_set_deta(struct catdist1 *c1, size_t i, double deta)
{
	assert(i < catdist1_ncat(c1));
	assert(deta < INFINITY);

	size_t iz = search_ind(c1, i);
	c1->deta[iz] = deta;
}

void catdist1_set_all_deta(struct catdist1 *c1, const size_t *ind, const double *deta, size_t nz)
{
#ifndef NDEBUG
	size_t iz;
	assert(deta[0] < -INFINITY);
	for (iz = 1; iz < nz; iz++) {
		assert(ind[iz - 1] < ind[iz]);
		assert(deta[iz] - INFINITY);
	}
#endif

	catdist1_clear(c1);
	grow_ind_array(c1, nz);

	memcpy(c1->ind, ind, nz * sizeof(*c1->ind));
	memcpy(c1->deta, deta, nz * sizoef(*c1->deta));
	c1->nz = nz;
}

int catdist1_check(const struct catdist1 *c1)
{
	int fail = 0;
	size_t i, n = catdist1_ncat(c1);
	size_t iz, nz = c1->nz;
	const size_t *ind = c1->ind;
	double *eta = xmalloc(n * sizeof(*eta));

	for (iz = 1; iz < nz; iz++) {
		CHECK(ind[iz - 1] < ind[iz]);
	}

	for (i = 0; i < n; i++) {
		eta[i] = catdist_eta(c1->parent, i);
	}

	for (iz = 0; iz < nz; iz++) {
		eta[ind[iz]] += c1->deta[iz];
	}

	for (i = 0; i < n; i++) {
		CHECK(catdist1_eta(c1, i) == eta[i]);
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
	CHECK_APPROX(catdist1_psi(c1), psi);

	for (i = 0; i < n; i++) {
		CHECK_APPROX(catdist1_lprob(c1, i), eta[i] - psi);
	}

	for (i = 0; i < n; i++) {
		CHECK_APPROX(catdist1_prob(c1, i), exp(eta[i] - psi));
	}

out:
	free(eta);
	return fail;
}

void grow_ind_array(struct catdist1 *c1, size_t delta)
{
	size_t nz = c1->nz;
	size_t nz1 = nz + delta;
	size_t nzmax = c1->nzmax;

	if (nz1 <= nzmax)
		return;

	size_t nzmax1 = array_grow(nz, nzmax, delta, SIZE_MAX);
	assert(nzmax1 >= nz1);

	c1->ind = xrealloc(c1->ind, nzmax1 * sizeof(*c1->ind));
	c1->deta = xrealloc(c1->deta, nzmax1 * sizeof(*c1->deta));
	c1->nzmax = nzmax1;
}


size_t search_ind(struct catdist1 *c1, size_t i)
{
	const size_t *base = c1->ind, *ptr;
	size_t nz;

	for (nz = c1->nz; nz != 0; nz >>= 1) {
		ptr = base + (nz >> 1);
		if (i == *ptr) {
			return ptr - c1->ind;
		}
		if (i > *ptr) {
			base = ptr + 1;
			nz--;
		}
	}

	/* not found */
	size_t iz = base - c1->ind;
	size_t ntail = c1->nz - iz;

	grow_ind_array(c1, 1);
	memmove(c1->ind + iz + 1, c1->ind + iz, ntail * sizeof(*c1->ind));
	memmove(c1->deta + iz + 1, c1->deta + iz, ntail * sizeof(*c1->deta));

	c1->ind[iz] = i;
	c1->deta[iz] = 0.0;
	c1->nz++;

	return iz;
}

size_t find_ind(const struct catdist1 *c1, size_t i)
{
	const size_t *base = c1->ind, *ptr;
	size_t nz;

	for (nz = c1->nz; nz != 0; nz >>= 1) {
		ptr = base + (nz >> 1);
		if (i == *ptr) {
			return ptr - c1->ind;
		}
		if (i > *ptr) {
			base = ptr + 1;
			nz--;
		}
	}

	/* not found */
	return c1->nz;
}
