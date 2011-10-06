#ifndef F77LAPACK_H
#define F77LAPACK_H

#include "f77.h"

#define F77_SVDJOB(j) const char *j ## _ = (j == LA_SVD_ALL ? "A" \
					    : j == LA_SVD_SEP ? "S" \
					    : j == LA_SVD_OVER ? "O" : "N")
#define F77_EIGJOB(j) const char *j ## _ = (j == LA_EIG_VEC ? "V" : "N")


int F77_FUNC(dgesdd) (const char *jobz, const f77int *m, const f77int *n,
		      double *a, const f77int *lda, double *s, double *u,
		      const f77int *ldu, double *vt, const f77int *ldvt,
		      double *work, const f77int *lwork, f77int *iwork,
		      f77int *info);

int F77_FUNC(dposv) (const char *uplo, const f77int *n, const f77int *nrhs,
		     double *a, const f77int *lda, double *b, const f77int *ldb,
		     f77int *info);

int F77_FUNC(dsyevd) (const char *jobz, const char *uplo, const f77int *n,
		      double *a, const f77int *lda, double *w, double *work,
		      const f77int *lwork, f77int *iwork, const f77int *liwork,
		      f77int *info);

int F77_FUNC(dsysv) (const char *uplo, const f77int *n, const f77int *nrhs,
		     double *a, const f77int *lda, f77int *ipiv, double *b,
		     const f77int *ldb, double *work, const f77int *lwork,
		     f77int *info);

#endif /* F77LAPACK_H */
