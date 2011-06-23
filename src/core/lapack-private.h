#ifndef _LAPACK_PRIVATE_H
#define _LAPACK_PRIVATE_H

int F77_FUNC(dsyevd) (const char *jobz, const char *uplo, const f77int *n,
		      double *a, const f77int *lda, double *w, double *work,
		      const f77int *lwork, f77int *iwork, const f77int *liwork,
		      f77int *info);

#endif /* _LAPACK_PRIVATE_H */

