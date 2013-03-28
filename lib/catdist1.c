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


static int needs_update(const struct catdist1 *c1);
static void update_cache(struct catdist1 *c1);

static void update_nzmax(struct catdist1 *c1);

static double get_dpsi(const struct catdist1 *c1);
static double get_dpsi_safe(const struct catdist1 *c1);
static double get_dpsi_safer(const struct catdist1 *c1);


void catdist1_init(struct catdist1 *c1, struct catdist *parent)
{
	assert(parent);

	c1->parent = parent;
	uintset_init(&c1->ind);
	c1->deta = NULL;
	c1->nzmax = 0;
	c1->dpsi = 0.0;
	c1->cached = 1;
	catdist_add_checkpoint(c1->parent, &c1->parent_cp);
}

void catdist1_deinit(struct catdist1 *c1)
{
	catdist_remove_checkpoint(c1->parent, &c1->parent_cp);
	free(c1->deta);
	uintset_deinit(&c1->ind);
}

double catdist1_deta(const struct catdist1 *c1, size_t i)
{
	assert(i < catdist1_ncat(c1));

	double deta = 0.0;
	size_t iz;

	if (uintset_find(&c1->ind, i, &iz)) {
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


double get_dpsi(const struct catdist1 *c1)
{
	const size_t *ind;
	size_t iz, nz;
	double sum = 0.0;
	double sum1 = 0.0;

	uintset_get_vals(&c1->ind, &ind, &nz);

	for (iz = 0; iz < nz; iz++) {
		double eta = catdist_lprob(c1->parent, ind[iz]);
		double deta = c1->deta[iz];
		double eta1 = eta + deta;
		double w1 = exp(eta1);
		double w = exp(eta);
		sum += w;
		sum1 += w1;
	}

	double dpsi = log((1 - sum) + sum1);

	if (!isfinite(dpsi))
		dpsi = get_dpsi_safe(c1);

	return dpsi;
}


int needs_update(const struct catdist1 *c1)
{
	return !c1->cached || catdist_checkpoint_passed(&c1->parent_cp, c1->parent);
}


void update_cache(struct catdist1 *c1)
{
	double dpsi = get_dpsi(c1);
	c1->dpsi = dpsi;
	c1->cached = 1;
	catdist_checkpoint_set(&c1->parent_cp, c1->parent);
}


double get_dpsi_safe(const struct catdist1 *c1)
{
	const size_t *ind;
	size_t iz, nz;
	double dpsi = NAN;

	uintset_get_vals(&c1->ind, &ind, &nz);

	if (!nz)
		return 0.0;

	double max = -INFINITY;
	double max1 = -INFINITY;
	double suc1 = 0.0;
	double sum = 0.0;

	for (iz = 0; iz < nz; iz++) {
		double eta = catdist_lprob(c1->parent, ind[iz]);
		double deta = c1->deta[iz];
		double eta1 = eta + deta;

		max = MAX(max, eta);
		max1 = MAX(max1, eta1);
	}

	double shift = MIN(max, max1);

	for (iz = 0; iz < nz; iz++) {
		double eta = catdist_lprob(c1->parent, ind[iz]);
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
	const size_t *ind;
	size_t iz, nz;
	size_t i, n = catdist1_ncat(c1);

	uintset_get_vals(&c1->ind, &ind, &nz);

	if (!nz || !n)
		return 0.0;

	double etamax = -INFINITY;
	size_t imax = 0;

	for (i = 0, iz = 0; i < n; i++) {
		double eta = catdist_eta(c1->parent, i);

		if (iz < nz && ind[iz] == i) {
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

		if (iz < nz && ind[iz] == i) {
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


double catdist1_dpsi(const struct catdist1 *c1)
{
	if (needs_update(c1))
		update_cache((struct catdist1 *)c1);

	return c1->dpsi;
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

	struct uintset *ind = &c1->ind;
	size_t iz;

	if (!uintset_find(ind, i, &iz)) {
		size_t ntail = uintset_insert(ind, iz, i);

		update_nzmax(c1);
		memmove(c1->deta + iz + 1, c1->deta + iz, ntail * sizeof(double));
	}

	c1->deta[iz] = deta;
	c1->cached = 0;
}

void catdist1_set_all_deta(struct catdist1 *c1, const size_t *ind, const double *deta, size_t nz)
{
#ifndef NDEBUG
	size_t iz;
	assert(nz == 0 || deta[0] < INFINITY);
	for (iz = 1; iz < nz; iz++) {
		assert(ind[iz - 1] < ind[iz]);
		assert(deta[iz] < INFINITY);
	}
#endif
	uintset_assign_array(&c1->ind, ind, nz, 1);
	update_nzmax(c1);
	memcpy(c1->deta, deta, nz * sizeof(double));

	c1->cached = 0;
}

int catdist1_check(const struct catdist1 *c1)
{
	int fail = 0;
	size_t i, n = catdist1_ncat(c1);
	size_t iz, nz;
	const size_t *ind;
	double *eta = xmalloc(n * sizeof(*eta));

	uintset_get_vals(&c1->ind, &ind, &nz);

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


void update_nzmax(struct catdist1 *c1)
{
	size_t nzmax = uintset_capacity(&c1->ind);
	if (c1->nzmax != nzmax) {
		c1->nzmax = nzmax;
		c1->deta = xrealloc(c1->deta, nzmax * sizeof(double));
	}
}
