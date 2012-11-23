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

static void diff_init(struct catdist1_diff *diff);
static void diff_deinit(struct catdist1_diff *diff);
static void diff_clear(struct catdist1_diff *diff);
static void diff_set(struct catdist1_diff *diff, size_t i, double deta);
static void diff_set_all(struct catdist1_diff *diff, const size_t *ind, const double *deta, size_t nz);
static ptrdiff_t diff_find(const struct catdist1_diff *diff, size_t i);
static size_t diff_search(struct catdist1_diff *diff, size_t i);
static void diff_grow(struct catdist1_diff *diff, size_t delta);


static void catdist1_update_(struct catdist1 *c1);
static void update_diff(struct catdist1 *c1);
static void update_dpsi(struct catdist1 *c1);
static double get_dpsi_safe(const struct catdist1 *c1);
static double get_dpsi_safer(const struct catdist1 *c1);




void catdist1_init(struct catdist1 *c1, const struct catdist *parent)
{
	assert(parent);
	c1->parent = parent;
	diff_init(&c1->diff);
	diff_init(&c1->pending);
	c1->cached_dpsi = 0;
	c1->cleared = 0;
}


void catdist1_deinit(struct catdist1 *c1)
{
	diff_deinit(&c1->pending);
	diff_deinit(&c1->diff);
}


void catdist1_update_(struct catdist1 *c1)
{
	if (c1->cleared || c1->pending.nz) {
		update_diff(c1);
		update_dpsi(c1);
	}
}


void update_diff(struct catdist1 *c1)
{
	if (c1->cleared) {
		diff_set_all(&c1->diff, c1->pending.ind, c1->pending.deta, c1->pending.nz);
		c1->cleared = 0;
	} else {
		size_t iz, nz = c1->pending.nz;
		for (iz = 0; iz < nz; iz++) {
			size_t i = c1->pending.ind[iz];
			double deta = c1->pending.deta[iz];
			size_t iz1 = diff_search(&c1->diff, i);
			c1->diff.deta[iz1] = deta;
		}
	}
	diff_clear(&c1->pending);
}


double catdist1_deta(const struct catdist1 *c1, size_t i)
{
	assert(i < catdist1_ncat(c1));

	catdist1_update_((struct catdist1 *)c1);

	double deta = 0.0;
	ptrdiff_t iz = diff_find(&c1->diff, i);

	if (iz >= 0) {
		deta = c1->diff.deta[iz];
	}

	return deta;
}


double catdist1_dpsi(const struct catdist1 *c1)
{
	catdist1_update_((struct catdist1 *)c1);
	return c1->cached_dpsi;
}


double catdist1_eta(const struct catdist1 *c1, size_t i)
{
	assert(i < catdist1_ncat(c1));

	double eta = catdist_eta(c1->parent, i);
	double deta = catdist1_deta(c1, i);
	return eta + deta;
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


double catdist1_prob(const struct catdist1 *c1, size_t i)
{
	assert(i < catdist1_ncat(c1));
	double lp = catdist1_lprob(c1, i);
	return exp(lp);
}


double catdist1_psi(const struct catdist1 *c1)
{
	double psi = catdist_psi(c1->parent);
	double dpsi = catdist1_dpsi(c1);
	double psi1 = psi + dpsi;
	return psi1;
}


void catdist1_get_deta(const struct catdist1 *c1,
		       const size_t **ind,
		       const double **deta,
		       size_t *nz)
{
	catdist1_update_((struct catdist1 *)c1);
	*ind = c1->diff.ind;
	*deta = c1->diff.deta;
	*nz = c1->diff.nz;
}


void update_dpsi(struct catdist1 *c1)
{
	assert(!c1->pending.nz);

	size_t iz, nz = c1->diff.nz;
	double sum = 0.0;
	double sum1 = 0.0;

	for (iz = 0; iz < nz; iz++) {
		double eta = catdist_lprob(c1->parent, c1->diff.ind[iz]);
		double deta = c1->diff.deta[iz];
		double eta1 = eta + deta;
		double w1 = exp(eta1);
		double w = exp(eta);
		sum += w;
		sum1 += w1;
	}

	double dpsi = log((1 - sum) + sum1);

	if (!isfinite(dpsi))
		dpsi = get_dpsi_safe(c1);

	c1->cached_dpsi = dpsi;
}


double get_dpsi_safe(const struct catdist1 *c1)
{
	assert(!c1->pending.nz);

	size_t nz = c1->diff.nz;
	double dpsi = NAN;

	if (!nz)
		return 0.0;

	size_t iz;
	double max = -INFINITY;
	double max1 = -INFINITY;
	double suc1 = 0.0;
	double sum = 0.0;

	for (iz = 0; iz < nz; iz++) {
		double eta = catdist_lprob(c1->parent, c1->diff.ind[iz]);
		double deta = c1->diff.deta[iz];
		double eta1 = eta + deta;

		max = MAX(max, eta);
		max1 = MAX(max1, eta1);
	}

	double shift = MIN(max, max1);

	for (iz = 0; iz < nz; iz++) {
		double eta = catdist_lprob(c1->parent, c1->diff.ind[iz]);
		double deta = c1->diff.deta[iz];
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
	assert(!c1->pending.nz);

	size_t iz, nz = c1->diff.nz;
	size_t i, n = catdist1_ncat(c1);

	if (!nz || !n)
		return 0.0;

	double etamax = -INFINITY;
	size_t imax = 0;

	for (i = 0, iz = 0; i < n; i++) {
		double eta = catdist_eta(c1->parent, i);

		if (iz < nz && c1->diff.ind[iz] == i) {
			eta += c1->diff.deta[iz];
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

		if (iz < nz && c1->diff.ind[iz] == i) {
			eta += c1->diff.deta[iz];
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


void catdist1_set_deta(struct catdist1 *c1, size_t i, double deta)
{
	assert(i < catdist1_ncat(c1));
	assert(deta < INFINITY);
	diff_set(&c1->pending, i, deta);

}

void catdist1_set_all_deta(struct catdist1 *c1, const size_t *ind, const double *deta, size_t nz)
{
#ifndef NDEBUG
	size_t iz;
	assert(nz == 0 || ind[nz - 1] < catdist1_ncat(c1));

	for (iz = 0; iz < nz; iz++) {
		assert(deta[iz] < INFINITY);
	}
#endif
	diff_set_all(&c1->pending, ind, deta, nz);
	c1->cleared = 1;
}


int catdist1_check(const struct catdist1 *c1)
{
	int fail = 0;
	size_t i, n = catdist1_ncat(c1);
	size_t iz, nz = c1->diff.nz;
	const size_t *ind = c1->diff.ind;
	double *eta = xmalloc(n * sizeof(*eta));

	for (iz = 1; iz < nz; iz++) {
		CHECK(ind[iz - 1] < ind[iz]);
	}

	for (i = 0; i < n; i++) {
		eta[i] = catdist_eta(c1->parent, i);
	}

	for (iz = 0; iz < nz; iz++) {
		eta[ind[iz]] += c1->diff.deta[iz];
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
	double psi1 = catdist_psi(c1->parent) + c1->cached_dpsi;
	CHECK_APPROX(psi1, psi);
out:
	free(eta);
	return fail;
}



void diff_init(struct catdist1_diff *diff)
{
	diff->deta = NULL;
	diff->ind = NULL;
	diff->nzmax = 0;
	diff_clear(diff);
}


void diff_deinit(struct catdist1_diff *diff)
{
	free(diff->deta);
	free(diff->ind);
}


void diff_clear(struct catdist1_diff *diff)
{
	diff->nz = 0;
}


void diff_set(struct catdist1_diff *diff, size_t i, double deta)
{
	size_t iz = diff_search(diff, i);
	diff->deta[iz] = deta;
}


void diff_set_all(struct catdist1_diff *diff, const size_t *ind, const double *deta, size_t nz)
{
#ifndef NDEBUG
	size_t iz;
	for (iz = 1; iz < nz; iz++) {
		assert(ind[iz - 1] < ind[iz]);
		assert(deta[iz] < INFINITY);
	}
#endif
	diff->nz = 0;
	diff_grow(diff, nz);

	memcpy(diff->ind, ind, nz * sizeof(diff->ind[0]));
	memcpy(diff->deta, deta, nz * sizeof(diff->deta[0]));
	diff->nz = nz;
}


size_t diff_search(struct catdist1_diff *diff, size_t i)
{
	ptrdiff_t siz = find_index(i, diff->ind, diff->nz);
	size_t iz;

	/* not found */
	if (siz < 0) {
		iz = ~siz;
		size_t ntail = diff->nz - iz;

		diff_grow(diff, 1);
		memmove(diff->ind + iz + 1, diff->ind + iz, ntail * sizeof(*diff->ind));
		memmove(diff->deta + iz + 1, diff->deta + iz, ntail * sizeof(*diff->deta));

		diff->ind[iz] = i;
		diff->deta[iz] = 0.0;
		diff->nz++;
	} else {
		iz = siz;
	}

	return iz;
}


ptrdiff_t diff_find(const struct catdist1_diff *diff, size_t i)
{
	return find_index(i, diff->ind, diff->nz);
}


void diff_grow(struct catdist1_diff *diff, size_t delta)
{
	size_t nz = diff->nz;
	size_t nz1 = nz + delta;

	if (needs_grow(nz1, &diff->nzmax)) {
		diff->ind = xrealloc(diff->ind, diff->nzmax * sizeof(diff->ind[0]));
		diff->deta = xrealloc(diff->deta, diff->nzmax * sizeof(diff->deta[0]));
	}
}


