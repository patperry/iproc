#include "port.h"
#include <assert.h>
#include "lapack-private.h"
#include "util.h"
#include "linalg.h"

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
	
	F77_FUNC(dsyevd)(jobz, uplo, &n1, a, &lda, w, &work, &inone, liwork,
			 &inone, &info);
	assert(info == 0);
	*lwork = (f77int)work;
}

void symeig_init(struct symeig *eig, ssize_t n, enum eig_job job)
{
	assert(eig);
	assert(0 <= n && n <= F77INT_MAX);
	
	f77int lwork, liwork;
	
	symeig_get_worksize(n, job, &lwork, &liwork);
	
	eig->n = n;
	eig->job = job;
	eig->work = xmalloc(lwork * sizeof(eig->work[0]));
	eig->lwork = lwork;
	eig->iwork = xmalloc(liwork * sizeof(eig->iwork[0]));
	eig->liwork = liwork;
	eig->info = -1;
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

	F77_FUNC(dsyevd)(jobz, suplo, &n, pa, &lda, pw, work, &lwork, iwork,
			 &liwork, &eig->info);
	assert(eig->info >= 0);

	return (eig->info == 0);
}
