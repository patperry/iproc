#include "port.h"
#include <assert.h>
#include <stdlib.h>
#include "lapack.h"
#include "xalloc.h"
#include "linalg.h"

ssize_t chol_solve(enum blas_uplo uplo, struct matrix *a, struct matrix *b)
{
	assert(matrix_nrow(a) == matrix_ncol(a));
	assert(matrix_nrow(b) == matrix_nrow(a));

	size_t n = matrix_nrow(b);
	size_t nrhs = matrix_ncol(b);
	double *pa = matrix_to_ptr(a);
	size_t lda = matrix_lda(a);
	double *pb = matrix_to_ptr(b);
	size_t ldb = matrix_lda(b);
	ptrdiff_t info;

	info = lapack_dposv(uplo, n, nrhs, pa, lda, pb, ldb);
	assert(info >= 0);

	return (ssize_t)info;
}

static void ldlfac_get_worksize(ssize_t n, size_t *lwork)
{
	assert(lwork);
	assert(0 <= n);

	*lwork = lapack_dsysv_lwork(n);
}

void ldlfac_init(struct ldlfac *fac, ssize_t n)
{
	assert(fac);
	assert(0 <= n);

	fac->n = 0;
	fac->work = NULL;
	fac->lwork = 0;
	fac->ipiv = NULL;
	fac->info = -1;

	if (n > 0)
		ldlfac_reinit(fac, n);
}

void ldlfac_reinit(struct ldlfac *fac, ssize_t n)
{
	assert(fac);
	assert(0 <= n);

	fac->n = n;
	ldlfac_get_worksize(n, &fac->lwork);
	fac->work = xrealloc(fac->work, fac->lwork * sizeof(fac->work[0]));
	fac->ipiv = xrealloc(fac->ipiv, n * sizeof(fac->ipiv[0]));
	fac->info = -1;
}

void ldlfac_deinit(struct ldlfac *fac)
{
	assert(fac);
	free(fac->ipiv);
	free(fac->work);
}

ptrdiff_t ldlfac_solve(struct ldlfac *fac, enum blas_uplo uplo,
		       struct matrix *a, struct matrix *b)
{
	assert(fac);
	assert(a);
	assert(b);
	assert(matrix_nrow(a) == ldlfac_dim(fac));
	assert(matrix_ncol(a) == ldlfac_dim(fac));
	assert(matrix_nrow(b) == ldlfac_dim(fac));

	size_t n = matrix_nrow(b);
	size_t nrhs = matrix_ncol(b);
	double *pa = matrix_to_ptr(a);
	size_t lda = matrix_lda(a);
	ptrdiff_t *ipiv = fac->ipiv;
	double *pb = matrix_to_ptr(b);
	size_t ldb = matrix_lda(b);
	double *work = fac->work;
	size_t lwork = fac->lwork;

	fac->info = lapack_dsysv(uplo, n, nrhs, pa, lda, ipiv, pb, ldb, work,
			  	 lwork);

	return fac->info;
}

static void symeig_get_worksize(ssize_t n, enum eig_job job,
				size_t *lwork, size_t *liwork)
{
	assert(0 <= n);
	assert(lwork);
	assert(liwork);

	*lwork = lapack_dsyevd_lwork(job, n, liwork);
}

void symeig_init(struct symeig *eig, ssize_t n, enum eig_job job)
{
	assert(eig);
	assert(0 <= n);

	eig->n = 0;
	eig->job = EIG_NOVEC;
	eig->work = NULL;
	eig->lwork = 0;
	eig->iwork = NULL;
	eig->liwork = 0;
	eig->info = 0;
	symeig_reinit(eig, n, job);
}

void symeig_reinit(struct symeig *eig, ssize_t n, enum eig_job job)
{
	assert(eig);
	assert(0 <= n);

	size_t lwork, liwork;

	symeig_get_worksize(n, job, &lwork, &liwork);

	eig->n = n;
	eig->job = job;
	eig->work = xrealloc(eig->work, lwork * sizeof(eig->work[0]));
	eig->lwork = lwork;
	eig->iwork = xrealloc(eig->iwork, liwork * sizeof(eig->iwork[0]));
	eig->liwork = liwork;
	eig->info = 0;
}

void symeig_deinit(struct symeig *eig)
{
	assert(eig);

	free(eig->work);
	free(eig->iwork);
}

bool symeig_factor(struct symeig *eig, enum blas_uplo uplo,
		   struct matrix *a, struct vector *w)
{
	assert(eig);
	assert(matrix_nrow(a) == symeig_dim(eig));
	assert(matrix_ncol(a) == symeig_dim(eig));
	assert(vector_dim(w) == symeig_dim(eig));

	size_t n = symeig_dim(eig);
	double *pa = matrix_to_ptr(a);
	size_t lda = (size_t)matrix_lda(a);
	double *pw = vector_to_ptr(w);
	double *work = eig->work;
	size_t lwork = eig->lwork;
	ptrdiff_t *iwork = eig->iwork;
	size_t liwork = eig->liwork;

	if (n <= 1) {
		assert(lwork >= 1);
		assert(liwork >= 1);
	} else if (eig->job == EIG_NOVEC) {
		assert(lwork >= 2 * n + 1);
		assert(liwork >= 1);
	} else {
		assert(lwork >= 1 + 6 * n + 2 * n * n);
		assert(liwork >= 3 + 5 * n);
	}
	
	eig->info = lapack_dsyevd(eig->job, uplo, n, pa, lda, pw, work, lwork,
				  iwork, liwork);
	assert(eig->info >= 0);
	return (eig->info == 0);
}

static void svdfac_get_worksize(ssize_t m, ssize_t n, enum svd_job job,
				size_t *lwork, size_t *liwork)
{
	assert(0 <= m);
	assert(0 <= n);
	assert(lwork);
	assert(liwork);
	
	*lwork = lapack_dgesdd_lwork(job, m, n, liwork);
}

void svdfac_init(struct svdfac *svd, ssize_t m, ssize_t n, enum svd_job job)
{
	assert(svd);
	assert(0 <= m && m <= F77INT_MAX);	
	assert(0 <= n && n <= F77INT_MAX);
	
	svd->work = NULL;
	svd->iwork = NULL;
	svdfac_reinit(svd, m, n, job);
}

void svdfac_reinit(struct svdfac *svd, ssize_t m, ssize_t n, enum svd_job job)
{
	assert(svd);
	assert(0 <= m && m <= F77INT_MAX);	
	assert(0 <= n && n <= F77INT_MAX);
	
	size_t lwork, liwork;
	
	svdfac_get_worksize(m, n, job, &lwork, &liwork);
	
	svd->m = m;
	svd->n = n;
	svd->job = job;
	svd->work = xrealloc(svd->work, lwork * sizeof(svd->work[0]));
	svd->lwork = lwork;
	svd->iwork = xrealloc(svd->iwork, liwork * sizeof(svd->iwork[0]));
	svd->info = -1;
}

void svdfac_deinit(struct svdfac *svd)
{
	assert(svd);
	
	free(svd->work);
	free(svd->iwork);
}

bool svdfac_factor(struct svdfac *svd, struct matrix *a, struct vector *s, struct matrix *u, struct matrix *vt)
{
	enum svd_job job = svdfac_job(svd);	
	size_t m = svdfac_row_dim(svd);
	size_t n = svdfac_col_dim(svd);
#ifndef NDEBUG
	size_t mn = MIN(m, n);
#endif
	
	assert(svd);
	assert(a);
	assert((size_t)matrix_nrow(a) == m);
	assert((size_t)matrix_ncol(a) == n);
	assert(s);
	assert((size_t)vector_dim(s) == mn);
	switch(job) {
	case SVD_ALL:
		assert(u);
		assert((size_t)matrix_nrow(u) == m);
		assert((size_t)matrix_ncol(u) == m);
		assert(vt);
		assert((size_t)matrix_nrow(vt) == n);
		assert((size_t)matrix_ncol(vt) == n);
		break;
	case SVD_SEPARATE:
		assert(u);
		assert((size_t)matrix_nrow(u) == m);
		assert((size_t)matrix_ncol(u) == mn);
		assert(vt);
		assert((size_t)matrix_nrow(vt) == mn);
		assert((size_t)matrix_ncol(vt) == n);
		break;
	case SVD_OVERWRITE:
		if (m < n) {
			assert(u);
			assert((size_t)matrix_nrow(u) == m);
			assert((size_t)matrix_ncol(u) == m);
			assert(!vt);
		} else {
			assert(!u);
			assert(vt);
			assert((size_t)matrix_nrow(vt) == n);
			assert((size_t)matrix_ncol(vt) == n);
		}
		break;
	case SVD_NOVEC:
		assert(!u);
		assert(!vt);
		break;
	}

	double *pa = matrix_to_ptr(a);
	size_t lda = (size_t)matrix_lda(a);
	double *ps = vector_to_ptr(s);
	double *pu = u ? matrix_to_ptr(u) : NULL;
	size_t ldu = u ? (size_t)matrix_lda(u) : 1;
	double *pvt = vt ? matrix_to_ptr(vt) : NULL;
	size_t ldvt = vt ? (size_t)matrix_lda(vt) : 1;
	double *work = svd->work;
	size_t lwork = svd->lwork;
	ptrdiff_t *iwork = svd->iwork;
	
	svd->info = lapack_dgesdd(job, m, n, pa, lda, ps, pu, ldu, pvt, ldvt,
				  work, lwork, iwork);
	return (svd->info == 0);
}


