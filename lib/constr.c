#include "port.h"
#include "coreutil.h"
#include "lapack.h"
#include "matrixutil.h"
#include "xalloc.h"
#include <inttypes.h>
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "constr.h"


#define SVTOL	MATRIX_DRANK_SVTOL0
#define VARTOL	(1e-8)
#define EIGTOL	(1e-8)


static int constr_add_check(struct constr *c);
static void constr_grow(struct constr *c, size_t delta);
static size_t constr_rank(const struct constr *c);


void constr_init(struct constr *c, size_t dim)
{
	c->dim = dim;
	c->wts = NULL;
	c->vals = NULL;
	c->n = 0;
	c->nmax = 0;
}

void constr_deinit(struct constr *c)
{
	free(c->vals);
	free(c->wts);
}


int constr_add(struct constr *c, const double *wts, double val)
{
	assert(wts || !constr_dim(c));
	
	size_t dim = c->dim;
	size_t n = c->n;
	
	constr_grow(c, 1);

	memcpy(c->wts + n * dim, wts, dim * sizeof(c->wts[0]));
	c->vals[n] = val;
	c->n = n + 1;

	return constr_add_check(c);
}


int constr_add_set(struct constr *c, size_t i, double val)
{
	assert(i < constr_dim(c));
	assert(isfinite(val));

	size_t dim = c->dim;
	size_t n = c->n;

	constr_grow(c, 1);

	memset(c->wts + n * dim, 0, dim * sizeof(c->wts[0]));
	c->wts[i + n * dim] = 1.0;
	c->vals[n] = val;
	c->n = n + 1;

	return constr_add_check(c);
}


int constr_add_eq(struct constr *c, size_t i1, size_t i2)
{
	assert(i1 != i2);
	assert(i1 < constr_dim(c));
	assert(i2 < constr_dim(c));

	size_t dim = c->dim;
	size_t n = c->n;

	constr_grow(c, 1);

	memset(c->wts + n * dim, 0, dim * sizeof(c->wts[0]));
	c->wts[i1 + n * dim] = +1.0;
	c->wts[i2 + n * dim] = -1.0;
	c->vals[n] = 0.0;
	c->n = n + 1;

	return constr_add_check(c);
}


size_t constr_add_identify(struct constr *c, const double *imatp,
			   enum blas_uplo uplo)
{
	enum blas_uplo f77uplo = uplo == BLAS_LOWER ? BLAS_UPPER : BLAS_LOWER;
	size_t iadd, nadd = 0;
	size_t i, dim = c->dim;
	size_t n = c->n;

	if (dim == 0)
		goto out;

	/* copy the information matrix */
	double *imat = xmalloc(dim * dim * sizeof(imat[0]));
	packed_dsctr(f77uplo, dim, imatp, imat, dim);

	/* compute the scale */
	double *scale = xmalloc(dim * sizeof(scale[0]));
	for (i = 0; i < dim; i++) {
		double s2 = imat[i * dim + i];
		scale[i] = 1.0 / (s2 < VARTOL ? 1.0 : sqrt(s2));
	}

	/* scale the information matrix */
	matrix_dscalr(dim, dim, scale, 1, imat, dim);
	matrix_dscalc(dim, dim, scale, 1, imat, dim);

	/* compute the eigendecomposition of the information matrix */
	enum lapack_eigjob jobz = LA_EIG_VEC;
	size_t lwork, liwork;
	lwork = lapack_dsyevd_lwork(jobz, dim, &liwork);
	double *work = xmalloc(lwork * sizeof(work[0]));
	ptrdiff_t *iwork = xmalloc(liwork * sizeof(iwork[0]));
	double *evals = xmalloc(dim * sizeof(evals[0]));
	ptrdiff_t info = lapack_dsyevd(jobz, f77uplo, dim, imat, dim, evals,
				       work, lwork, iwork, liwork);
	assert(info == 0);

	/* change back to the original scale */
	matrix_dscalr(dim, dim, scale, 1, imat, dim);

	/* compute the dimension of the nullspace */
	size_t nulldim;
	for (nulldim = 0; nulldim < dim; nulldim++) {
		if (evals[nulldim] > EIGTOL)
			break;
	}
	/* at this point, the first nulldim cols of imat store its nullspace */

	if (nulldim == 0)
		goto cleanup_eig_imat;

	if (n == 0) {
		/* add the nullspace of the information matrix as constraints */
		nadd = nulldim;
		for (iadd = 0; iadd < nadd; iadd++) {
			int ok = constr_add(c, imat + iadd * dim, 0.0);
			assert(ok);
			(void)ok;
		}
	} else {
		/* project the constraints onto the nullspace */
		double *proj = xmalloc(nulldim * n * sizeof(proj[0]));
		blas_dgemm(BLAS_TRANS, BLAS_NOTRANS, nulldim, n, dim, 1.0,
			   imat, dim, c->wts, dim, 0.0, proj, nulldim);

		/* compute the nullspace of the projection */
		double *s = xmalloc(MIN(n, nulldim) * sizeof(s[0]));
		double *u = xmalloc(nulldim * nulldim * sizeof(u[0]));
		double *vt = xmalloc(n * n * sizeof(vt[0]));

		enum lapack_svdjob svd_jobz = LA_SVD_ALL;
		lwork = lapack_dgesdd_lwork(svd_jobz, nulldim, n, &liwork);
		work = xrealloc(work, lwork * sizeof(work[0]));
		iwork = xrealloc(iwork, liwork * sizeof(iwork[0]));
		ptrdiff_t svd_info = lapack_dgesdd(svd_jobz, nulldim, n, proj,
						   nulldim, s, u, nulldim, vt,
						   n, work, lwork, iwork);
		assert(svd_info == 0);

		/* compute the rank */
		size_t rank, maxrank = MIN(n, nulldim);
		for (rank = maxrank; rank > 0; rank--) {
			if (s[rank - 1] > SVTOL * fabs(s[maxrank - 1]))
				break;
		}

		/* at this point the last nulldim - rank columns of u hold the
		   coefficients of the orthogonal complement to the constraint
		   matrix in the nullspace */

		/* compute a basis for the orthogonal complement;
		   these are the constraints */
		nadd = nulldim - rank;
		double *compl = xmalloc(dim * nadd * sizeof(compl[0]));
		blas_dgemm(BLAS_NOTRANS, BLAS_NOTRANS, dim, nadd, nulldim, 1.0,
			   imat, dim, u + rank, nulldim, 0.0, compl, dim);

		for (iadd = 0; iadd < nadd; iadd++) {
			int ok = constr_add(c, compl + iadd * dim, 0.0);
			assert(ok);
			(void)ok;
		}

		free(compl);
		free(vt);
		free(u);
		free(s);
		free(proj);
	}

cleanup_eig_imat:
	free(evals);
	free(iwork);
	free(work);
	free(scale);
	free(imat);
out:
	return nadd;

}


int constr_add_check(struct constr *c)
{
	assert(c->n);
	size_t n = c->n;
	size_t rank = constr_rank(c);

	if (rank < n) {
		c->n = n - 1;
		return 0;
	} else {
		return 1;
	}
}

void constr_grow(struct constr *c, size_t delta)
{
	if (needs_grow(c->n + delta, &c->nmax)) {
		size_t dim = c->dim;
		c->wts = xrealloc(c->wts, c->nmax * dim * sizeof(*c->wts));
		c->vals = xrealloc(c->vals, c->nmax * sizeof(*c->vals));
	}
}

static size_t constr_rank(const struct constr *c)
{
	if (!c->dim)
		return 0;

	size_t dim = c->dim;
	size_t n = c->n;
	size_t liwork = 0;
	size_t lwork = matrix_drank_lwork(dim, n, &liwork);

	double *work = xmalloc(lwork * sizeof(double));
	ptrdiff_t *iwork = xmalloc(liwork * sizeof(ptrdiff_t));
	double *wts = xmalloc(n * dim * sizeof(double));

	memcpy(wts, c->wts, n * dim * sizeof(double));
	size_t rank = matrix_drank(dim, n, wts, dim, SVTOL, work, lwork, iwork);

	free(wts);
	free(iwork);
	free(work);

	return rank;
}


