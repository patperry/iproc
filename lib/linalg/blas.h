#ifndef BLAS_H
#define BLAS_H

#include <stddef.h>

enum blas_trans {
	BLAS_NOTRANS = 111,
	BLAS_TRANS = 112,
	BLAS_CONJTRANS = 113
};

enum blas_uplo {
	BLAS_UPPER = 121,
	BLAS_LOWER = 122
};

enum blas_diag {
	BLAS_NONUNIT = 131,
	BLAS_UNIT = 132
};

enum blas_side {
	BLAS_LEFT = 141,
	BLAS_RIGHT = 142
};

/* Level 1 */

double blas_ddot(size_t n, const double *x, size_t incx, const double *y,
		 size_t incy);

double blas_dnrm2(size_t n, const double *x, size_t incx);

double blas_dasum(size_t n, const double *x, size_t incx);

size_t blas_idamax(size_t n, const double *x, size_t incx);

void blas_dswap(size_t n, double *x, size_t incx, double *y, size_t incy);

void blas_dcopy(size_t n, const double *x, size_t incx, double *y, size_t incy);

void blas_daxpy(size_t n, double alpha, const double *x, size_t incx,
		double *y, size_t incy);

void blas_drotg(double *a, double *b, double *c, double *s);

void blas_drotmg(double *d1, double *d2, double *b1, const double *b2,
		 double *p);

void blas_drot(size_t n, double *x, size_t incx, double *y, size_t incy,
	       double c, double s);

void blas_drotm(size_t n, double *x, size_t incx, double *y, size_t incy,
		const double *p);

void blas_dscal(size_t n, double alpha, double *x, size_t incx);

/* Level 2 */

void blas_dgemv(enum blas_trans trans, size_t m, size_t n, double alpha,
		const double *a, size_t lda, const double *x, size_t incx,
		double beta, double *y, size_t incy);

void blas_dgbmv(enum blas_trans trans, size_t m, size_t n, size_t kl, size_t ku,
		double alpha, const double *a, size_t lda,
		const double *x, size_t incx, double beta, double *y,
		size_t incy);

void blas_dsbmv(enum blas_uplo uplo, size_t n, size_t k, double alpha,
		const double *a, size_t lda, const double *x, size_t incx,
		double beta, double *y, size_t incy);

void blas_dspmv(enum blas_uplo uplo, size_t n, double alpha,
		const double *ap, const double *x, size_t incx,
		double beta, double *y, size_t incy);

void blas_dsymv(enum blas_uplo uplo, size_t n, double alpha,
		const double *a, size_t lda, const double *x, size_t incx,
		double beta, double *y, size_t incy);

void blas_dtbmv(enum blas_uplo uplo, enum blas_trans trans, enum blas_diag diag,
		size_t n, size_t k, const double *a, size_t lda, double *x,
		size_t incx);

void blas_dtpmv(enum blas_uplo uplo, enum blas_trans trans, enum blas_diag diag,
		size_t n, const double *ap, double *x, size_t incx);

void blas_dtrmv(enum blas_uplo uplo, enum blas_trans trans, enum blas_diag diag,
		size_t n, const double *a, size_t lda, double *x, size_t incx);

void blas_dtbsv(enum blas_uplo uplo, enum blas_trans trans, enum blas_diag diag,
		size_t n, size_t k, const double *a, size_t lda, double *x,
		size_t incx);

void blas_dtpsv(enum blas_uplo uplo, enum blas_trans trans, enum blas_diag diag,
		size_t n, const double *ap, double *x, size_t incx);

void blas_dtrsv(enum blas_uplo uplo, enum blas_trans trans, enum blas_diag diag,
		size_t n, const double *a, size_t lda, double *x, size_t incx);

void blas_dger(size_t m, size_t n, double alpha, const double *x,
	       size_t incx, const double *y, size_t incy, double *a,
	       size_t lda);

void blas_dsyr(enum blas_uplo uplo, size_t n, double alpha,
	       const double *x, size_t incx, double *a, size_t lda);

void blas_dspr(enum blas_uplo uplo, size_t n, double alpha,
	       const double *x, size_t incx, double *ap);

void blas_dsyr2(enum blas_uplo uplo, size_t n, double alpha,
		const double *x, size_t incx, const double *y, size_t incy,
		double *a, size_t lda);

void blas_dspr2(enum blas_uplo uplo, size_t n, double alpha,
		const double *x, size_t incx, const double *y, size_t incy,
		double *ap);

/* Level 3 */

void blas_dgemm(enum blas_trans transa, enum blas_trans transb, size_t m,
		size_t n, size_t k, double alpha, const double *a,
		size_t lda, const double *b, size_t ldb, double beta,
		double *c, size_t ldc);

void blas_dtrsm(enum blas_side side, enum blas_uplo uplo,
		enum blas_trans transa, enum blas_diag diag, size_t m, size_t n,
		double alpha, const double *a, size_t lda, double *b,
		size_t ldb);

void blas_dtrmm(enum blas_side side, enum blas_uplo uplo,
		enum blas_trans transa, enum blas_diag diag, size_t m, size_t n,
		double alpha, const double *a, size_t lda, double *b,
		size_t ldb);

void blas_dsymm(enum blas_side side, enum blas_uplo uplo, size_t m, size_t n,
		double alpha, const double *a, size_t lda,
		const double *b, size_t ldb, double beta, double *c,
		size_t ldc);

void blas_dsyrk(enum blas_uplo uplo, enum blas_trans trans, size_t n, size_t k,
		double alpha, const double *a, size_t lda,
		double beta, double *c, size_t ldc);

void blas_dsyr2k(enum blas_uplo uplo, enum blas_trans trans, size_t n, size_t k,
		 double alpha, const double *a, size_t lda,
		 const double *b, size_t ldb, double beta, double *c,
		 size_t ldc);

#endif /* BLAS_H */
