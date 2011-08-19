#include "port.h"
#include <assert.h>
#include "lapack-private.h"
#include "util.h"
#include "linalg.h"

ssize_t chol_solve(enum matrix_uplo uplo, struct matrix *a, struct matrix *b)
{
	assert(matrix_nrow(a) == matrix_ncol(a));
	assert(matrix_nrow(b) == matrix_nrow(a));

	char *cuplo = uplo == UPLO_LOWER ? "L" : "U";
	f77int n = (f77int)matrix_nrow(b);
	f77int nrhs = (f77int)matrix_ncol(b);
	double *pa = matrix_to_ptr(a);
	f77int lda = (f77int)matrix_lda(a);
	double *pb = matrix_to_ptr(b);
	f77int ldb = (f77int)matrix_lda(b);
	f77int info = 0;

	F77_FUNC(dposv) (cuplo, &n, &nrhs, pa, &lda, pb, &ldb, &info);
	assert(info >= 0);

	return (ssize_t)info;
}

static void ldlfac_get_worksize(ssize_t n, f77int *lwork)
{
	assert(lwork);
	assert(0 <= n && n <= F77INT_MAX);

	const char *uplo = "U";
	f77int n1 = (f77int)n;
	f77int nrhs = 1;
	double *a = NULL;
	f77int lda = MAX(1, n1);
	f77int *ipiv = NULL;
	double *b = NULL;
	f77int ldb = lda;
	double work = 0;
	f77int info = 0;

	*lwork = -1;
	F77_FUNC(dsysv) (uplo, &n1, &nrhs, a, &lda, ipiv, b, &ldb, &work,
			 lwork, &info);
	*lwork = (f77int)work;
	assert(info == 0);
}

void ldlfac_init(struct ldlfac *fac, ssize_t n)
{
	assert(fac);
	assert(0 <= n && n <= F77INT_MAX);

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
	assert(0 <= n && n <= F77INT_MAX);

	fac->n = n;
	ldlfac_get_worksize(n, &fac->lwork);
	fac->work = xrealloc(fac->work, fac->lwork * sizeof(fac->work[0]));
	fac->ipiv = xrealloc(fac->ipiv, n * sizeof(fac->ipiv[0]));
	fac->info = -1;
}

void ldlfac_deinit(struct ldlfac *fac)
{
	assert(fac);
	xfree(fac->ipiv);
	xfree(fac->work);
}

ssize_t ldlfac_solve(struct ldlfac *fac, enum matrix_uplo uplo,
		     struct matrix *a, struct matrix *b)
{
	assert(fac);
	assert(a);
	assert(b);
	assert(matrix_nrow(a) == ldlfac_dim(fac));
	assert(matrix_ncol(a) == ldlfac_dim(fac));
	assert(matrix_nrow(b) == ldlfac_dim(fac));

	const char *cuplo = uplo == UPLO_LOWER ? "L" : "U";
	f77int n = (f77int)matrix_nrow(b);
	f77int nrhs = (f77int)matrix_ncol(b);
	double *pa = matrix_to_ptr(a);
	f77int lda = (f77int)matrix_lda(a);
	f77int *ipiv = fac->ipiv;
	double *pb = matrix_to_ptr(b);
	f77int ldb = (f77int)matrix_lda(b);
	double *work = fac->work;
	f77int lwork = fac->lwork;

	F77_FUNC(dsysv) (cuplo, &n, &nrhs, pa, &lda, ipiv, pb, &ldb, work,
			 &lwork, &fac->info);

	return (ssize_t)fac->info;
}

static void symeig_get_worksize(ssize_t n, enum eig_job job,
				f77int *lwork, f77int *liwork)
{
	assert(lwork);
	assert(liwork);
	assert(0 <= n && n <= F77INT_MAX);

	f77int n1 = (f77int)n;
	const char *jobz = job == EIG_NOVEC ? "N" : "V";
	const char *uplo = "U";
	double *a = NULL;
	f77int lda = MAX(1, n1);
	double *w = NULL;
	double work = 0;
	f77int inone = -1;
	f77int info = 0;

	F77_FUNC(dsyevd) (jobz, uplo, &n1, a, &lda, w, &work, &inone, liwork,
			  &inone, &info);
	assert(info == 0);
	*lwork = (f77int)work;

	if (n <= 1) {
		assert(*lwork >= 1);
		assert(*liwork >= 1);
	} else if (job == EIG_NOVEC) {
		assert(*lwork >= 2 * n + 1);
		assert(*liwork >= 1);
	} else {
		assert(*lwork >= 1 + 6 * n + 2 * n * n);
		assert(*liwork >= 3 + 5 * n);
	}
}

void symeig_init(struct symeig *eig, ssize_t n, enum eig_job job)
{
	assert(eig);
	assert(0 <= n && n <= F77INT_MAX);

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
	assert(0 <= n && n <= F77INT_MAX);

	f77int lwork, liwork;

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

	xfree(eig->work);
	xfree(eig->iwork);
}

bool symeig_factor(struct symeig *eig, enum matrix_uplo uplo,
		   struct matrix *a, struct vector *w)
{
	assert(eig);
	assert(matrix_nrow(a) == symeig_dim(eig));
	assert(matrix_ncol(a) == symeig_dim(eig));
	assert(vector_dim(w) == symeig_dim(eig));
	assert(matrix_lda(a) <= F77INT_MAX);

	f77int n = (f77int)symeig_dim(eig);
	const char *jobz = eig->job == EIG_NOVEC ? "N" : "V";
	const char *suplo = uplo == UPLO_UPPER ? "U" : "L";
	double *pa = matrix_to_ptr(a);
	f77int lda = (f77int)matrix_lda(a);
	double *pw = vector_to_ptr(w);
	double *work = eig->work;
	f77int lwork = eig->lwork;
	f77int *iwork = eig->iwork;
	f77int liwork = eig->liwork;

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
	
	F77_FUNC(dsyevd) (jobz, suplo, &n, pa, &lda, pw, work, &lwork, iwork,
			  &liwork, &eig->info);
	assert(eig->info >= 0);

	return (eig->info == 0);
}

static void svdfac_get_worksize(ssize_t m, ssize_t n, enum svd_job job,
				f77int *lwork)
{
	assert(lwork);
	assert(0 <= m && m <= F77INT_MAX);	
	assert(0 <= n && n <= F77INT_MAX);
	
	f77int m1 = (f77int)m;
	f77int n1 = (f77int)n;
	f77int mn1 = MIN(m1, n1);
	const char *jobz = (job == SVD_ALL ? "A"
			    : job == SVD_SEPARATE ? "S"
			    : job == SVD_OVERWRITE ? "O"
			    : "N");
	double *a = NULL;
	f77int lda = MAX(1, m1);
	double *s = NULL;
	double *u = NULL;
	f77int ldu = MAX(1, m1);
	double *vt = NULL;
	f77int ldvt = (job == SVD_ALL ? MAX(1, n1)
		       : job == SVD_OVERWRITE && m >= n ? MAX(1, n1)
		       : job == SVD_SEPARATE ? MAX(1, mn1)
		       : 1);
	double work = 0;
	f77int inone = -1;
	f77int *iwork = NULL;
	f77int info = 0;
	
	F77_FUNC(dgesdd) (jobz, &m1, &n1, a, &lda, s, u, &ldu, vt, &ldvt, &work, &inone, iwork, &info);
	assert(info == 0);
	*lwork = (f77int)work;
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
	
	f77int lwork, liwork;
	
	svdfac_get_worksize(m, n, job, &lwork);
	liwork = 8 * (f77int)MIN(m, n);
	
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
	
	xfree(svd->work);
	xfree(svd->iwork);
}

bool svdfac_factor(struct svdfac *svd, struct matrix *a, struct vector *s, struct matrix *u, struct matrix *vt)
{
	f77int m = (f77int)svdfac_row_dim(svd);
	f77int n = (f77int)svdfac_col_dim(svd);
	f77int mn = MIN(m, n);
	enum svd_job job = svdfac_job(svd);
	
	assert(svd);
	assert(a);
	assert(matrix_nrow(a) == m);
	assert(matrix_ncol(a) == n);
	assert(s);
	assert(vector_dim(s) == mn);
	switch(job) {
	case SVD_ALL:
		assert(u);
		assert(matrix_nrow(u) == m);
		assert(matrix_ncol(u) == m);
		assert(vt);
		assert(matrix_nrow(vt) == n);
		assert(matrix_ncol(vt) == n);
		break;
	case SVD_SEPARATE:
		assert(u);
		assert(matrix_nrow(u) == m);
		assert(matrix_ncol(u) == mn);
		assert(vt);
		assert(matrix_nrow(vt) == mn);
		assert(matrix_ncol(vt) == n);
		break;
	case SVD_OVERWRITE:
		if (m < n) {
			assert(u);
			assert(matrix_nrow(u) == m);
			assert(matrix_ncol(u) == m);
			assert(!vt);
		} else {
			assert(!u);
			assert(vt);
			assert(matrix_nrow(vt) == n);
			assert(matrix_ncol(vt) == n);
		}
		break;
	case SVD_NOVEC:
		assert(!u);
		assert(!vt);
		break;
	}

	const char *jobz = (job == SVD_ALL ? "A"
			    : job == SVD_SEPARATE ? "S"
			    : job == SVD_OVERWRITE ? "O"
			    : "N");
	double *pa = matrix_to_ptr(a);
	f77int lda = (f77int)matrix_lda(a);
	double *ps = vector_to_ptr(s);
	double *pu = u ? matrix_to_ptr(u) : NULL;
	f77int ldu = u ? (f77int)matrix_lda(u) : 1;
	double *pvt = vt ? matrix_to_ptr(vt) : NULL;
	f77int ldvt = vt ? (f77int)matrix_lda(vt) : 1;
	double *work = svd->work;
	f77int lwork = svd->lwork;
	f77int *iwork = svd->iwork;
	f77int *info = &svd->info;
	
	F77_FUNC(dgesdd) (jobz, &m, &n, pa, &lda, ps, pu, &ldu, pvt, &ldvt, work, &lwork, iwork, info);
	assert(*info >= 0);
	return (*info == 0);
}


