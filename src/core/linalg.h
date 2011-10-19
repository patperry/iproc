#ifndef _LINALG_H
#define _LINALG_H

#include "vector.h"
#include "matrix.h"
#include "lapack.h"

struct symeig {
	ssize_t n;
	enum lapack_eigjob job;
	double *work;
	size_t lwork;
	ptrdiff_t *iwork;
	size_t liwork;
	ptrdiff_t info;
};

struct ldlfac {
	ssize_t n;
	double *work;
	size_t lwork;
	ptrdiff_t *ipiv;
	ptrdiff_t info;
};

struct svdfac {
	ssize_t m;
	ssize_t n;
	enum lapack_svdjob job;
	double *work;
	size_t lwork;
	ptrdiff_t *iwork;
	ptrdiff_t info;
};

// replace b with a \ b; replace a with its cholesky factor
// returns 0 on success, i > 0 if leading minor of order i is not pd
ssize_t chol_solve(enum blas_uplo uplo, struct matrix *a, struct matrix *b);

void ldlfac_init(struct ldlfac *fac, ssize_t n);
void ldlfac_reinit(struct ldlfac *fac, ssize_t n);
void ldlfac_deinit(struct ldlfac *fac);
static inline ssize_t ldlfac_dim(const struct ldlfac *fac)
{
	assert(fac);
	return fac->n;
}

ptrdiff_t ldlfac_solve(struct ldlfac *fac, enum blas_uplo uplo,
		       struct matrix *a, struct matrix *b);

void symeig_init(struct symeig *eig, ssize_t n, enum lapack_eigjob job);
void symeig_reinit(struct symeig *eig, ssize_t n, enum lapack_eigjob job);
void symeig_deinit(struct symeig *eig);

bool symeig_factor(struct symeig *eig, enum blas_uplo uplo, struct matrix *a, struct vector *w);	// destroys a

static inline ssize_t symeig_dim(const struct symeig *eig)
{
	assert(eig);
	return eig->n;
}

static inline enum lapack_eigjob symeig_job(const struct symeig *eig)
{
	assert(eig);
	return eig->job;
}

void svdfac_init(struct svdfac *svd, ssize_t m, ssize_t n, enum lapack_svdjob job);
void svdfac_reinit(struct svdfac *svd, ssize_t m, ssize_t n, enum lapack_svdjob job);
void svdfac_deinit(struct svdfac *svd);

bool svdfac_factor(struct svdfac *svd, struct matrix *a, struct vector *s, struct matrix *u, struct matrix *vt); // destroys a

static inline ssize_t svdfac_row_dim(const struct svdfac *svd)
{
	assert(svd);
	return svd->m;
}

static inline ssize_t svdfac_col_dim(const struct svdfac *svd)
{
	assert(svd);
	return svd->n;
}


static inline enum lapack_svdjob svdfac_job(const struct svdfac *svd)
{
	assert(svd);
	return svd->job;
}


#endif /* _LINALG_H */
