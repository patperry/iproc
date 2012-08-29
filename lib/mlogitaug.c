#include "port.h"
#include <stdint.h>
#include <string.h>
#include "coreutil.h"
#include "sblas.h"
#include "xalloc.h"
#include "mlogitaug.h"


#define CHECK(x) \
	do { \
		fail = !(x); \
		assert(!fail); \
		if (fail) \
			goto out; \
	} while(0)


static void clear(struct mlogitaug *m1);
static void clear_x(struct mlogitaug *m1);
static void clear_offset(struct mlogitaug *m1);
static void recompute_dist(struct mlogitaug *m1);
static void recompute_mean(struct mlogitaug *m1);
static void set_x_iz(struct mlogitaug *m1, size_t iz, const double *x);
static void grow_ind_array(struct mlogitaug *m1, size_t delta);
static size_t find_ind(const struct mlogitaug *m1, size_t i);
static size_t search_ind(struct mlogitaug *m1, size_t i);


void mlogitaug_init(struct mlogitaug *m1, const struct mlogit *base,
		size_t dim)
{
	m1->base = base;
	catdist1_init(&m1->dist, mlogit_dist(base));

	m1->beta = xmalloc(dim * sizeof(*m1->beta));
	m1->dim = dim;

	m1->ind = NULL;
	m1->offset = NULL;
	m1->x = NULL;
	m1->deta = NULL;
	m1->nz = 0;
	m1->nzmax = 0;
	m1->mean = xmalloc(dim * sizeof(*m1->mean));
	m1->base_mean = xmalloc(mlogit_dim(m1->base) * sizeof(*m1->base_mean));
	m1->xbuf = xmalloc(dim * sizeof(*m1->xbuf));

	clear(m1);
}


void clear(struct mlogitaug *m1)
{
	size_t dim = m1->dim;
	catdist1_set_all_deta(&m1->dist, NULL, NULL, 0);
	memset(m1->beta, 0, dim * sizeof(*m1->beta));
	m1->nz = 0;
}


void mlogitaug_deinit(struct mlogitaug *m1)
{
	free(m1->xbuf);
	free(m1->base_mean);
	free(m1->mean);
	free(m1->deta);
	free(m1->x);
	free(m1->offset);
	free(m1->ind);
	free(m1->beta);
	catdist1_deinit(&m1->dist);
}


double *mlogitaug_coefs(const struct mlogitaug *m1)
{
	return m1->beta;
}


double mlogitaug_offset(const struct mlogitaug *m1, size_t i)
{
	size_t iz = find_ind(m1, i);
	if (iz == m1->nz)
		return 0.0;

	return m1->offset[iz];
}


double *mlogitaug_x(const struct mlogitaug *m1, size_t i)
{
	size_t iz = find_ind(m1, i);
	if (iz == m1->nz)
		return NULL;

	return m1->x + iz * m1->dim;
}


void mlogitaug_set_coefs(struct mlogitaug *m1, const double *beta)
{
	size_t len = m1->dim * sizeof(*m1->beta);
	
	if (!beta) {
		memset(m1->beta, 0, len);
	} else {
		memcpy(m1->beta, beta, len);
	}
}


void mlogitaug_set_offset(struct mlogitaug *m1, size_t i, double offset)
{
	size_t iz = search_ind(m1, i);
	m1->offset[iz] = offset;
}


void mlogitaug_set_all_offset(struct mlogitaug *m1, const size_t *i, const double *offset, size_t nz)
{
	size_t iz;
	clear_offset(m1);
	for (iz = 0; iz < nz; iz++) {
		mlogitaug_set_offset(m1, i[iz], offset[iz]);
	}
}


void set_x_iz(struct mlogitaug *m1, size_t iz, const double *x)
{
	size_t dim = m1->dim;
	double *dst = m1->x + iz * dim;
	size_t len = dim * sizeof(*dst);

	if (x) {
		memcpy(dst, x, len);
	} else {
		memset(dst, 0, len);
	}

}


void mlogitaug_set_x(struct mlogitaug *m1, size_t i, const double *x)
{
	size_t iz = search_ind(m1, i);
	set_x_iz(m1, iz, x);
}


void mlogitaug_inc_x(struct mlogitaug *m1, size_t i, const size_t *jdx,
		const double *dx, size_t ndx)
{
	assert(i < mlogitaug_ncat(m1));

	size_t dim = m1->dim;
	size_t iz = search_ind(m1, i);
	double *x = m1->x + iz * dim;

	if (jdx) {
		sblas_daxpyi(ndx, 1.0, dx, jdx, x);
	} else {
		assert(ndx == dim);
		blas_daxpy(dim, 1.0, dx, 1, x, 1);
	}
}


/* Note: this implementation leaves the ind array unchanged; we could
 * potentially reset the ind array, but we would have to leave the indices
 * of the nonzero offsets.
 */
void clear_x(struct mlogitaug *m1)
{
	size_t dim = m1->dim;
	size_t nz = m1->nz;

	memset(m1->x, 0, dim * nz * sizeof(*m1->x));
}

/* Note from clear_x applies here */
void clear_offset(struct mlogitaug *m1)
{
	size_t nz = m1->nz;
	memset(m1->offset, 0, nz * sizeof(*m1->offset));
}


void mlogitaug_set_all_x(struct mlogitaug *m1, const size_t *i, const double *x, size_t nz)
{
	size_t dim = m1->dim;
	size_t iz;

	clear_x(m1);

	for (iz = 0; iz < nz; iz++) {
		mlogitaug_set_x(m1, i[iz], x + iz * dim);
	}
}


void recompute_dist(struct mlogitaug *m1)
{
	size_t nz = m1->nz;
	size_t dim = m1->dim;
	const double *beta = m1->beta;
	const double *offset = m1->offset;
	double *deta = m1->deta;
	const double *x = m1->x;

	blas_dcopy(nz, offset, 1, deta, 1);
	if (dim)
		blas_dgemv(BLAS_TRANS, dim, nz, 1.0, x, dim, beta, 1, 1.0, deta, 1);

	catdist1_set_all_deta(&m1->dist, m1->ind, m1->deta, m1->nz);
}


struct catdist1 *mlogitaug_dist(const struct mlogitaug *m1)
{
	recompute_dist((struct mlogitaug *)m1);
	return &((struct mlogitaug *)m1)->dist;
}


static void recompute_mean(struct mlogitaug *m1)
{
	recompute_dist(m1);

	const struct catdist1 *dist;
	size_t dim = m1->dim;
	double *mean = m1->mean;
	double *diff = m1->xbuf;
	size_t iz, nz = m1->nz;
	double w, wtot = 0.0;

	if (nz == 0) {
		memset(mean, 0, dim * sizeof(*mean));
		return;
	}

	blas_dcopy(dim, m1->x, 1, mean, 1);
	wtot += catdist1_cached_prob(dist, m1->ind[0]);

	for (iz = 1; iz < nz; iz++) {
		blas_dcopy(dim, m1->x + iz * dim, 1, diff, 1);
		blas_daxpy(dim, -1.0, mean, 1, diff, 1);

		w = catdist1_cached_prob(dist, m1->ind[iz]);
		wtot += w;

		blas_daxpy(dim, w/wtot, diff, 1, mean, 1);
	}

	blas_dscal(dim, wtot, mean, 1);
}


double *mlogitaug_mean(const struct mlogitaug *m1)
{
	recompute_mean((struct mlogitaug *)m1);
	return m1->mean;
}

double *mlogitaug_base_mean(const struct mlogitaug *m1)
{
	return m1->base_mean;
}



int mlogitaug_check(const struct mlogitaug *m1)
{
	int fail = 0;
	CHECK(!catdist1_check(&m1->dist));
out:
	return fail;
}


void grow_ind_array(struct mlogitaug *m1, size_t delta)
{
	size_t nz = m1->nz;
	size_t nz1 = nz + delta;
	size_t nzmax = m1->nzmax;

	if (nz1 <= nzmax)
		return;

	size_t nzmax1 = array_grow(nz, nzmax, delta, SIZE_MAX);
	size_t dim = m1->dim;
	assert(nzmax1 >= nz1);

	m1->ind = xrealloc(m1->ind, nzmax1 * sizeof(*m1->ind));
	m1->offset = xrealloc(m1->offset, nzmax1 * sizeof(*m1->offset));
	m1->x = xrealloc(m1->x, nzmax1 * sizeof(*m1->x) * dim);
	m1->deta = xrealloc(m1->deta, nzmax1 * sizeof(*m1->deta));
	m1->nzmax = nzmax1;
}


size_t search_ind(struct mlogitaug *m1, size_t i)
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
	size_t dim = m1->dim;

	grow_ind_array(m1, 1);
	memmove(m1->ind + iz + 1, m1->ind + iz, ntail * sizeof(*m1->ind));
	memmove(m1->offset + iz + 1, m1->offset + iz, ntail * sizeof(*m1->offset));
	memmove(m1->x + (iz + 1) * dim, m1->x + iz * dim, ntail * sizeof(*m1->x) * dim);
	memmove(m1->deta + iz + 1, m1->deta + iz, ntail * sizeof(*m1->deta));

	m1->ind[iz] = i;
	m1->offset[iz] = 0.0;
	memset(m1->x + iz * dim, 0, sizeof(*m1->x) * dim);
	m1->deta[iz] = 0.0;
	m1->nz++;

	return iz;
}


size_t find_ind(const struct mlogitaug *m1, size_t i)
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
