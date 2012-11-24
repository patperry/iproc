#include "port.h"
#include <stdint.h>
#include <string.h>
#include "coreutil.h"
#include "matrixutil.h"
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

#define F77_COV_UPLO (MLOGIT_COV_UPLO == BLAS_LOWER ? BLAS_UPPER : BLAS_LOWER)


#define BLOCK_SIZE	64

static int mlogitaug_moments(const struct mlogitaug *m1);
static void mlogitaug_update_cache(struct mlogitaug *m1);

static int needs_update(const struct mlogitaug *m1);
static void clear(struct mlogitaug *m1);
static void clear_x(struct mlogitaug *m1);
static void clear_offset(struct mlogitaug *m1);
static void recompute_dist(struct mlogitaug *m1);
static void recompute_mean(struct mlogitaug *m1);
static void recompute_base_mean(struct mlogitaug *m1);
static void recompute_cov(struct mlogitaug *m1);
static void recompute_base_cov(struct mlogitaug *m1);
static void recompute_cross_cov(struct mlogitaug *m1);
static void set_x_iz(struct mlogitaug *m1, size_t iz, const double *x);
static void grow_ind_array(struct mlogitaug *m1, size_t delta);
static ptrdiff_t find_ind(const struct mlogitaug *m1, size_t i);
static size_t search_ind(struct mlogitaug *m1, size_t i);


void mlogitaug_work_init(struct mlogitaug_work *work, size_t base_dim,
			 size_t aug_dim)
{
	size_t dim = MAX(base_dim, aug_dim);
	size_t dim2 = MAX(dim * dim, BLOCK_SIZE);
	work->cov_full = xmalloc(dim2 * sizeof(double));
	work->xbuf = xmalloc(dim * BLOCK_SIZE * sizeof(double));
}

void mlogitaug_work_deinit(struct mlogitaug_work *work)
{
	free(work->cov_full);
	free(work->xbuf);
}


void mlogitaug_init(struct mlogitaug *m1, struct mlogit *base,
		size_t dim, struct mlogitaug_work *work)
{
	size_t base_dim = mlogit_dim(base);
	size_t cov_dim = dim * (dim + 1) / 2;
	size_t base_cov_dim = base_dim * (base_dim + 1) / 2;

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
	m1->cov = xmalloc(cov_dim * sizeof(*m1->cov));

	m1->base_mean = xmalloc(base_dim * sizeof(*m1->base_mean));
	m1->base_cov = xmalloc(base_cov_dim * sizeof(*m1->base_cov));
	m1->base_dmean = xmalloc(base_dim * sizeof(*m1->base_dmean));

	m1->cross_cov = xmalloc(base_dim * dim * sizeof(*m1->cross_cov));

	if (work) {
		m1->free_work = 0;
	} else {
		work = xmalloc(sizeof(*work));
		mlogitaug_work_init(work, base_dim, dim);
		m1->free_work = 1;
	}
	m1->work = work;

	clear(m1);
	m1->cached = 0;

	mlogit_add_checkpoint(m1->base, &m1->base_cp);
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
	mlogit_remove_checkpoint(m1->base, &m1->base_cp);
	
	if (m1->free_work) {
		mlogitaug_work_deinit(m1->work);
		free(m1->work);
	}
	free(m1->cross_cov);
	free(m1->base_dmean);
	free(m1->base_cov);
	free(m1->base_mean);
	free(m1->cov);
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
	ptrdiff_t iz = find_ind(m1, i);
	if (iz < 0)
		return 0.0;

	return m1->offset[iz];
}


double *mlogitaug_x(const struct mlogitaug *m1, size_t i)
{
	ptrdiff_t iz = find_ind(m1, i);
	if (iz < 0)
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
	m1->cached = 0;
}


void mlogitaug_set_offset(struct mlogitaug *m1, size_t i, double offset)
{
	size_t iz = search_ind(m1, i);
	m1->offset[iz] = offset;
	m1->cached = 0;
}


void mlogitaug_set_all_offset(struct mlogitaug *m1, const size_t *i, const double *offset, size_t nz)
{
	size_t iz;
	clear_offset(m1);
	for (iz = 0; iz < nz; iz++) {
		mlogitaug_set_offset(m1, i[iz], offset[iz]);
	}
	m1->cached = 0;
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
	m1->cached = 0;
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

	m1->cached = 0;
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

	m1->cached = 0;
}


void recompute_dist(struct mlogitaug *m1)
{
	size_t nz = m1->nz;
	size_t dim = m1->dim;
	const double *beta = m1->beta;
	const double *offset = m1->offset;
	double *deta = m1->deta;
	const double *x = m1->x;

	(void)mlogit_dist(m1->base); // update pareent dist;

	blas_dcopy(nz, offset, 1, deta, 1);
	if (dim)
		blas_dgemv(BLAS_TRANS, dim, nz, 1.0, x, dim, beta, 1, 1.0, deta, 1);

	catdist1_set_all_deta(&m1->dist, m1->ind, m1->deta, m1->nz);
}


int needs_update(const struct mlogitaug *m1)
{
	return !m1->cached || mlogit_checkpoint_passed(&m1->base_cp, m1->base);
}

struct catdist1 *mlogitaug_dist(const struct mlogitaug *m1)
{
	if (needs_update(m1))
		mlogitaug_update_cache((struct mlogitaug *)m1);
	return &((struct mlogitaug *)m1)->dist;
}


void recompute_mean(struct mlogitaug *m1)
{
	if (mlogitaug_moments(m1) < 1)
		return;

	const struct catdist1 *dist = &m1->dist;
	size_t dim = m1->dim;
	double *mean = m1->mean;
	double *diff = m1->work->xbuf;
	double *prob = m1->work->cov_full;

	if (dim == 0)
		return;

	memset(mean, 0, dim * sizeof(*mean));

	size_t iz, nz = m1->nz;
	double *diff_i = diff;
	size_t nk = 0;

	for (iz = 0; iz < nz; iz++) {
		double w = catdist1_prob(dist, m1->ind[iz]);

		if (w > 0) {
			prob[nk] = w;
			blas_dcopy(dim, m1->x + iz * dim, 1, diff_i, 1);
			diff_i += dim;
			nk++;
		}

		if (nk == BLOCK_SIZE || iz + 1 == nz) {
			blas_dgemv(BLAS_NOTRANS, dim, nk, 1.0, diff, dim, prob, 1, 1.0, mean, 1);
			diff_i = diff;
			nk = 0;
		}
	}
}


double *mlogitaug_mean(const struct mlogitaug *m1)
{
	assert(mlogitaug_moments(m1) >= 1);
	if (needs_update(m1))
		mlogitaug_update_cache((struct mlogitaug *)m1);
	return m1->mean;
}


void recompute_base_mean(struct mlogitaug *m1)
{
	if (mlogitaug_moments(m1) < 1)
		return;

	size_t dim = mlogit_dim(m1->base);

	if (dim == 0)
		return;

	const struct catdist1 *dist = &m1->dist;
	const struct catdist *dist0 = catdist1_parent(dist);

	double *mean = m1->base_mean;
	double *dmean = m1->base_dmean;
	double *diff = m1->work->xbuf;
	double *dw = m1->work->cov_full;

	double psi = catdist1_psi(dist);
	const double *x = mlogit_x(m1->base);
	const double *mean0 = mlogit_mean(m1->base);


	memset(dmean, 0, dim * sizeof(*dmean));

	size_t iz, nz = m1->nz;
	double *diff_i = diff;
	size_t nk = 0;

	for (iz = 0; iz < nz; iz++) {
		size_t i = m1->ind[iz];
		double eta0 = catdist_eta(dist0, i);
		double w0 = exp(eta0 - psi);
		double w = catdist1_prob(dist, i);
		double dw_i = w - w0;

		if (dw_i != 0) {
			dw[nk] = dw_i;
			blas_dcopy(dim, x + i * dim, 1, diff_i, 1);
			blas_daxpy(dim, -1.0, mean0, 1, diff_i, 1);
			diff_i += dim;
			nk++;
		}

		if (nk == BLOCK_SIZE || iz + 1 == nz) {
			blas_dgemv(BLAS_NOTRANS, dim, nk, 1.0, diff, dim, dw, 1, 1.0, dmean, 1);
			diff_i = diff;
			nk = 0;
		}
	}

	blas_dcopy(dim, mean0, 1, mean, 1);
	blas_daxpy(dim, 1.0, dmean, 1, mean, 1);
}


int mlogitaug_moments(const struct mlogitaug *m1)
{
	const struct mlogit *m = mlogitaug_base(m1);
	return mlogit_moments(m);
}



double *mlogitaug_base_mean(const struct mlogitaug *m1)
{
	assert(mlogitaug_moments(m1) >= 1);
	if (needs_update(m1))
		mlogitaug_update_cache((struct mlogitaug *)m1);
	return m1->base_mean;
}


void recompute_cov(struct mlogitaug *m1)
{
	if (mlogitaug_moments(m1) < 2)
		return;

	const struct catdist1 *dist = &m1->dist;
	size_t dim = m1->dim;
	const double *mean = m1->mean;
	double *diff = m1->work->xbuf;
	double *cov = m1->work->cov_full;
	double w, wtot = 0.0;

	if (dim == 0)
		return;

	memset(cov, 0, dim * dim * sizeof(*cov));

	size_t iz, nz = m1->nz;
	size_t nk = 0;
	double *diff_i = diff;

	for (iz = 0; iz < nz; iz++) {
		w = catdist1_prob(dist, m1->ind[iz]);

		if (w != 0) {
			blas_dcopy(dim, m1->x + iz * dim, 1, diff_i, 1);
			blas_daxpy(dim, -1.0, mean, 1, diff_i, 1);
			blas_dscal(dim, sqrt(w), diff_i, 1);
			wtot += w;

			diff_i += dim;
			nk++;
		}

		if (nk == BLOCK_SIZE) {
			blas_dsyrk(F77_COV_UPLO, BLAS_NOTRANS, dim, nk, 1.0, diff, dim, 1.0, cov, dim);

			nk = 0;
			diff_i = diff;
		}
	}

	assert(nk < BLOCK_SIZE);
	if (wtot < 1) { // note: wtot might be > 1 due to numerical error
		blas_dcopy(dim, mean, 1, diff_i, 1);
		blas_dscal(dim, sqrt(1.0 - wtot), diff_i, 1);
		nk++;
	}
	blas_dsyrk(F77_COV_UPLO, BLAS_NOTRANS, dim, nk, 1.0, diff, dim, 1.0, cov, dim);

	packed_dgthr(F77_COV_UPLO, dim, cov, dim, m1->cov);
}


double *mlogitaug_cov(const struct mlogitaug *m1)
{
	assert(mlogitaug_moments(m1) >= 2);
	if (needs_update(m1))
		mlogitaug_update_cache((struct mlogitaug *)m1);
	return m1->cov;
}


static void recompute_base_cov(struct mlogitaug *m1)
{
	if (mlogitaug_moments(m1) < 2)
		return;

	const struct catdist1 *dist = &m1->dist;
	const struct catdist *dist0 = catdist1_parent(dist);
	size_t dim = mlogit_dim(m1->base);
	size_t cov_dim = dim * (dim + 1) / 2;

	double *dmean = m1->base_dmean;
	double *diff = m1->work->xbuf;
	double *cov = m1->work->cov_full;

	double psi0 = catdist_psi(dist0);
	double psi = catdist1_psi(dist);
	double dpsi = psi - psi0;
	const double *x = mlogit_x(m1->base);
	const double *mean0 = mlogit_mean(m1->base);
	//const double *cov0 = m1->base->cov;
	//double log_scale0 = m1->base->log_cov_scale;
	double scale0;
	const double *cov0 = mlogit_cov(m1->base, &scale0);
	const double log_scale0 = log(scale0);
	const double log_scale = log_scale0 + dpsi;
	const double W = exp(log_scale);
	const double scale = isfinite(W) ? W : 1.0;

	if (dim == 0)
		return;

	memset(cov, 0, dim * dim * sizeof(*cov));

	size_t iz, nz = m1->nz;
	size_t nk = 0;
	double *diff_i = diff;

	// handle positive dw
	for (iz = 0; iz < nz; iz++) {
		size_t i = m1->ind[iz];
		double eta0 = catdist_eta(dist0, i);
		double w0 = exp(eta0 - psi);
		double w = catdist1_prob(dist, i);
		double dw = w - w0;

		if (dw > 0) {
			blas_dcopy(dim, x + i * dim, 1, diff_i, 1);
			blas_daxpy(dim, -1.0, mean0, 1, diff_i, 1);
			blas_dscal(dim, sqrt(dw), diff_i, 1);

			diff_i += dim;
			nk++;
		}

		if (nk == BLOCK_SIZE || iz + 1 == nz) {
			blas_dsyrk(F77_COV_UPLO, BLAS_NOTRANS, dim, nk, scale, diff, dim, 1.0, cov, dim);

			diff_i = diff;
			nk = 0;
		}
	}

	nk = 0;
	diff_i = diff;

	// handle dmean and negative dw
	blas_dcopy(dim, dmean, 1, diff_i, 1);
	diff_i += dim;
	nk++;

	for (iz = 0; iz < nz; iz++) {
		size_t i = m1->ind[iz];
		double eta0 = catdist_eta(dist0, i);
		double w0 = exp(eta0 - psi);
		double w = catdist1_prob(dist, i);
		double dw = w - w0;

		if (dw < 0) {
			blas_dcopy(dim, x + i * dim, 1, diff_i, 1);
			blas_daxpy(dim, -1.0, mean0, 1, diff_i, 1);
			blas_dscal(dim, sqrt(-dw), diff_i, 1);

			diff_i += dim;
			nk++;
		}

		if (nk == BLOCK_SIZE || iz + 1 == nz) {
			blas_dsyrk(F77_COV_UPLO, BLAS_NOTRANS, dim, nk, -scale, diff, dim, 1.0, cov, dim);

			diff_i = diff;
			nk = 0;
		}
	}

	// note: nz = 0 implies dmean = 0; no need to call dsyrk

	packed_dgthr(F77_COV_UPLO, dim, cov, dim, m1->base_cov);

	if (isfinite(W))
		blas_daxpy(cov_dim, 1.0, cov0, 1, m1->base_cov, 1);
	blas_dscal(cov_dim, 1.0/scale, m1->base_cov, 1);
}



double *mlogitaug_base_cov(const struct mlogitaug *m1)
{
	assert(mlogitaug_moments(m1) >= 2);
	if (needs_update(m1))
		mlogitaug_update_cache((struct mlogitaug *)m1);
	return m1->base_cov;
}


void recompute_cross_cov(struct mlogitaug *m1)
{
	if (mlogitaug_moments(m1) < 2)
		return;

	size_t dim = m1->dim;
	size_t base_dim = mlogit_dim(m1->base);
	const struct catdist1 *dist = &m1->dist;
	const double *base_mean = m1->base_mean;
	const double *base_x = mlogit_x(m1->base);
	const double *x = m1->x;
	double *cov = m1->cross_cov;
	double *diff = m1->work->xbuf;

	if (dim == 0 || base_dim == 0)
		return;

	memset(cov, 0, dim * base_dim * sizeof(*cov));

	size_t iz, nz = m1->nz;
	size_t nk = 0;
	double *diff_i = diff;

	for (iz = 0; iz < nz; iz++) {
		size_t i = m1->ind[iz];
		double w = catdist1_prob(dist, i);

		if (w != 0) {
			blas_dcopy(base_dim, base_x + i * base_dim, 1, diff_i, 1);
			blas_daxpy(base_dim, -1.0, base_mean, 1, diff_i, 1);
			blas_dscal(base_dim, w, diff_i, 1);
		} else {
			memset(diff_i, 0, base_dim * sizeof(double));
		}

		diff_i += base_dim;
		nk++;

		if (nk == BLOCK_SIZE || iz + 1 == nz) {
			if (MLOGIT_COV_UPLO == BLAS_LOWER) {
				blas_dgemm(BLAS_NOTRANS, BLAS_TRANS, base_dim, dim, nk,
					   1.0, diff, base_dim, x, dim, 1.0, cov, base_dim);
			} else {
				blas_dgemm(BLAS_NOTRANS, BLAS_TRANS, dim, base_dim, nk,
					   1.0, x, dim, diff, base_dim, 1.0, cov, dim);
			}
			x += nk * dim;
			
			nk = 0;
			diff_i = diff;
		}
	}
}


double *mlogitaug_cross_cov(const struct mlogitaug *m1)
{
	assert(mlogitaug_moments(m1) >= 2);
	if (needs_update(m1))
		mlogitaug_update_cache((struct mlogitaug *)m1);
	return m1->cross_cov;
}


void mlogitaug_update_cache(struct mlogitaug *m1)
{
	recompute_dist((struct mlogitaug *)m1);
	recompute_mean((struct mlogitaug *)m1);
	recompute_base_mean((struct mlogitaug *)m1);
	recompute_cov((struct mlogitaug *)m1);
	recompute_base_cov((struct mlogitaug *)m1);
	recompute_cross_cov((struct mlogitaug *)m1);
	m1->cached = 1;
	mlogit_checkpoint_set(&m1->base_cp, m1->base);
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

	if (needs_grow(nz1, &m1->nzmax)) {
		size_t nzmax = m1->nzmax;
		size_t dim = m1->dim;

		m1->ind = xrealloc(m1->ind, nzmax * sizeof(*m1->ind));
		m1->offset = xrealloc(m1->offset, nzmax * sizeof(*m1->offset));
		m1->x = xrealloc(m1->x, nzmax * sizeof(*m1->x) * dim);
		m1->deta = xrealloc(m1->deta, nzmax * sizeof(*m1->deta));
		m1->nzmax = nzmax;
	}
}


size_t search_ind(struct mlogitaug *m1, size_t i)
{
	ptrdiff_t siz = find_ind(m1, i);

	if (siz >= 0)
		return (size_t)(siz);

	/* not found */
	size_t iz = ~siz;
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


ptrdiff_t find_ind(const struct mlogitaug *m1, size_t i)
{
	return find_index(i, m1->ind, m1->nz);
}
