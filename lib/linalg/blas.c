#ifdef HAVE_CONFIG_H
# include <config.h>
#endif
#include "f77blas.h"
#include "blas.h"

/* Level 1 */

double blas_ddot(size_t n, const double *x, size_t incx, const double *y,
		 size_t incy)
{
	F77_INT(n);
	F77_INT(incx);
	F77_INT(incy);
	return F77_FUNC(ddot)(&n_, x, &incx_, y, &incy_);
}

double blas_dasum(size_t n, const double *x, size_t incx)
{
	F77_INT(n);
	F77_INT(incx);
	return F77_FUNC(dasum)(&n_, x, &incx_);
}

size_t blas_idamax(size_t n, const double *x, size_t incx)
{
	F77_INT(n);
	F77_INT(incx);
	f77int imax = F77_FUNC(idamax)(&n_, x, &incx_);
	return (size_t)imax;
}

void blas_dswap(size_t n, double *x, size_t incx, double *y, size_t incy)
{
	F77_INT(n);
	F77_INT(incx);
	F77_INT(incy);
	F77_FUNC(dswap)(&n_, x, &incx_, y, &incy_);
}

void blas_dcopy(size_t n, const double *x, size_t incx, double *y, size_t incy)
{
	F77_INT(n);
	F77_INT(incx);
	F77_INT(incy);
	F77_FUNC(dcopy)(&n_, x, &incx_, y, &incy_);
}

void blas_daxpy(size_t n, double alpha, const double *x, size_t incx,
		double *y, size_t incy)
{
	F77_INT(n);
	F77_INT(incx);
	F77_INT(incy);
	F77_FUNC(daxpy)(&n_, &alpha, x, &incx_, y, &incy_);
}

void blas_drotg(double *a, double *b, double *c, double *s)
{
	F77_FUNC(drotg)(a, b, c, s);
}

void blas_drotmg(double *d1, double *d2, double *b1, const double *b2,
		 double *p)
{
	F77_FUNC(drotmg)(d1, d2, b1, b2, p);
}

void blas_drot(size_t n, double *x, size_t incx, double *y, size_t incy,
	       double c, double s)
{
	F77_INT(n);
	F77_INT(incx);
	F77_INT(incy);
	F77_FUNC(drot)(&n_, x, &incx_, y, &incy_, &c, &s);
}

void blas_drotm(size_t n, double *x, size_t incx, double *y, size_t incy,
		const double *p)
{
	F77_INT(n);
	F77_INT(incx);
	F77_INT(incy);
	F77_FUNC(drotm)(&n_, x, &incx_, y, &incy_, p);
}

void blas_dscal(size_t n, double alpha, double *x, size_t incx)
{
	F77_INT(n);
	F77_INT(incx);
	F77_FUNC(dscal)(&n_, &alpha, x, &incx_);
}

/* Level 2 */

void blas_dgemv(enum blas_trans trans, size_t m, size_t n, double alpha,
		const double *a, size_t lda, const double *x, size_t incx,
		double beta, double *y, size_t incy)
{
	F77_TRANS(trans);
	F77_INT(m);
	F77_INT(n);
	F77_INT(lda);
	F77_INT(incx);
	F77_INT(incy);
	F77_FUNC(dgemv)(trans_, &m_, &n_, &alpha, a, &lda_, x, &incx_,
			&beta, y, &incy_);
}

void blas_dgbmv(enum blas_trans trans, size_t m, size_t n, size_t kl, size_t ku,
		double alpha, const double *a, size_t lda,
		const double *x, size_t incx, double beta, double *y,
		size_t incy)
{
	F77_TRANS(trans);
	F77_INT(m);
	F77_INT(n);
	F77_INT(kl);
	F77_INT(ku);
	F77_INT(lda);
	F77_INT(incx);
	F77_INT(incy);
	F77_FUNC(dgbmv)(trans_, &m_, &n_, &kl_, &ku_, &alpha, a, &lda_, x,
			&incx_, &beta, y, &incy_);
}

void blas_dsbmv(enum blas_uplo uplo, size_t n, size_t k, double alpha,
		const double *a, size_t lda, const double *x, size_t incx,
		double beta, double *y, size_t incy)
{
	F77_UPLO(uplo);
	F77_INT(n);
	F77_INT(k);
	F77_INT(lda);
	F77_INT(incx);
	F77_INT(incy);
	F77_FUNC(dsbmv)(uplo_, &n_, &k_, &alpha, a, &lda_, x, &incx_, &beta,
			y, &incy_);
}

void blas_dspmv(enum blas_uplo uplo, size_t n, double alpha,
		const double *ap, const double *x, size_t incx,
		double beta, double *y, size_t incy)
{
	F77_UPLO(uplo);
	F77_INT(n);
	F77_INT(incx);
	F77_INT(incy);
	F77_FUNC(dspmv)(uplo_, &n_, &alpha, ap, x, &incx_, &beta, y, &incy_);
}

void blas_dsymv(enum blas_uplo uplo, size_t n, double alpha,
		const double *a, size_t lda, const double *x, size_t incx,
		double beta, double *y, size_t incy)
{
	F77_UPLO(uplo);
	F77_INT(n);
	F77_INT(lda);
	F77_INT(incx);
	F77_INT(incy);
	F77_FUNC(dsymv)(uplo_, &n_, &alpha, a, &lda_, x, &incx_,
			&beta, y, &incy_);
}

void blas_dtbmv(enum blas_uplo uplo, enum blas_trans trans, enum blas_diag diag,
		size_t n, size_t k, const double *a, size_t lda, double *x,
		size_t incx)
{
	F77_UPLO(uplo);
	F77_TRANS(trans);
	F77_DIAG(diag);
	F77_INT(n);
	F77_INT(k);
	F77_INT(lda);
	F77_INT(incx);
	F77_FUNC(dtbmv)(uplo_, trans_, diag_, &n_, &k_, a, &lda_, x, &incx_);
}

void blas_dtpmv(enum blas_uplo uplo, enum blas_trans trans, enum blas_diag diag,
		size_t n, const double *ap, double *x, size_t incx)
{
	F77_UPLO(uplo);
	F77_TRANS(trans);
	F77_DIAG(diag);
	F77_INT(n);
	F77_INT(incx);
	F77_FUNC(dtpmv)(uplo_, trans_, diag_, &n_, ap, x, &incx_);
}

void blas_dtrmv(enum blas_uplo uplo, enum blas_trans trans, enum blas_diag diag,
		size_t n, const double *a, size_t lda, double *x, size_t incx)
{
	F77_UPLO(uplo);
	F77_TRANS(trans);
	F77_DIAG(diag);
	F77_INT(n);
	F77_INT(lda);
	F77_INT(incx);
	F77_FUNC(dtrmv)(uplo_, trans_, diag_, &n_, a, &lda_, x, &incx_);
}

void blas_dtbsv(enum blas_uplo uplo, enum blas_trans trans, enum blas_diag diag,
		size_t n, size_t k, const double *a, size_t lda, double *x,
		size_t incx)
{
	F77_UPLO(uplo);
	F77_TRANS(trans);
	F77_DIAG(diag);
	F77_INT(n);
	F77_INT(k);
	F77_INT(lda);
	F77_INT(incx);
	F77_FUNC(dtbsv)(uplo_, trans_, diag_, &n_, &k_, a, &lda_, x, &incx_);
}

void blas_dtpsv(enum blas_uplo uplo, enum blas_trans trans, enum blas_diag diag,
		size_t n, const double *ap, double *x, size_t incx)
{
	F77_UPLO(uplo);
	F77_TRANS(trans);
	F77_DIAG(diag);
	F77_INT(n);
	F77_INT(incx);
	F77_FUNC(dtpsv)(uplo_, trans_, diag_, &n_, ap, x, &incx_);
}

void blas_dtrsv(enum blas_uplo uplo, enum blas_trans trans, enum blas_diag diag,
		size_t n, const double *a, size_t lda, double *x, size_t incx)
{
	F77_UPLO(uplo);
	F77_TRANS(trans);
	F77_DIAG(diag);
	F77_INT(n);
	F77_INT(lda);
	F77_INT(incx);
	F77_FUNC(dtrsv)(uplo_, trans_, diag_, &n_, a, &lda_, x, &incx_);
}

void blas_dger(size_t m, size_t n, double alpha, const double *x,
	       size_t incx, const double *y, size_t incy, double *a,
	       size_t lda)
{
	F77_INT(m);
	F77_INT(n);
	F77_INT(incx);
	F77_INT(incy);
	F77_INT(lda);
	F77_FUNC(dger)(&m_, &n_, &alpha, x, &incx_, y, &incy_, a, &lda_);
}

void blas_dsyr(enum blas_uplo uplo, size_t n, double alpha,
	       const double *x, size_t incx, double *a, size_t lda)
{
	F77_UPLO(uplo);
	F77_INT(n);
	F77_INT(incx);
	F77_INT(lda);
	F77_FUNC(dsyr)(uplo_, &n_, &alpha, x, &incx_, a, &lda_);
}

void blas_dspr(enum blas_uplo uplo, size_t n, double alpha,
	       const double *x, size_t incx, double *ap)
{
	F77_UPLO(uplo);
	F77_INT(n);
	F77_INT(incx);
	F77_FUNC(dspr)(uplo_, &n_, &alpha, x, &incx_, ap);
}

void blas_dsyr2(enum blas_uplo uplo, size_t n, double alpha,
		const double *x, size_t incx, const double *y, size_t incy,
		double *a, size_t lda)
{
	F77_UPLO(uplo);
	F77_INT(n);
	F77_INT(incx);
	F77_INT(incy);
	F77_INT(lda);
	F77_FUNC(dsyr2)(uplo_, &n_, &alpha, x, &incx_, y, &incy_, a, &lda_);
}

void blas_dspr2(enum blas_uplo uplo, size_t n, double alpha,
		const double *x, size_t incx, const double *y, size_t incy,
		double *ap)
{
	F77_UPLO(uplo);
	F77_INT(n);
	F77_INT(incx);
	F77_INT(incy);
	F77_FUNC(dspr2)(uplo_, &n_, &alpha, x, &incx_, y, &incy_, ap);
}

/* Level 3 */

void blas_dgemm(enum blas_trans transa, enum blas_trans transb, size_t m,
		size_t n, size_t k, double alpha, const double *a,
		size_t lda, const double *b, size_t ldb, double beta,
		double *c, size_t ldc)
{
	F77_TRANS(transa);
	F77_TRANS(transb);
	F77_INT(m);
	F77_INT(n);
	F77_INT(k);
	F77_INT(lda);
	F77_INT(ldb);
	F77_INT(ldc);
	F77_FUNC(dgemm)(transa_, transb_, &m_, &n_, &k_, &alpha, a, &lda_, b,
			&ldb_, &beta, c, &ldc_);
}

void blas_dtrsm(enum blas_side side, enum blas_uplo uplo,
		enum blas_trans transa, enum blas_diag diag, size_t m, size_t n,
		double alpha, const double *a, size_t lda, double *b,
		size_t ldb)
{
	F77_SIDE(side);
	F77_UPLO(uplo);
	F77_TRANS(transa);
	F77_DIAG(diag);
	F77_INT(m);
	F77_INT(n);
	F77_INT(lda);
	F77_INT(ldb);
	F77_FUNC(dtrsm)(side_, uplo_, transa_, diag_, &m_, &n_, &alpha, a,
			&lda_, b, &ldb_);
}

void blas_dtrmm(enum blas_side side, enum blas_uplo uplo,
		enum blas_trans transa, enum blas_diag diag, size_t m, size_t n,
		double alpha, const double *a, size_t lda, double *b,
		size_t ldb)
{
	F77_SIDE(side);
	F77_UPLO(uplo);
	F77_TRANS(transa);
	F77_DIAG(diag);
	F77_INT(m);
	F77_INT(n);
	F77_INT(lda);
	F77_INT(ldb);
	F77_FUNC(dtrmm)(side_, uplo_, transa_, diag_, &m_, &n_, &alpha, a,
			&lda_, b, &ldb_);
}

void blas_dsymm(enum blas_side side, enum blas_uplo uplo, size_t m, size_t n,
		double alpha, const double *a, size_t lda,
		const double *b, size_t ldb, double beta, double *c,
		size_t ldc)
{
	F77_SIDE(side);
	F77_UPLO(uplo);
	F77_INT(m);
	F77_INT(n);
	F77_INT(lda);
	F77_INT(ldb);
	F77_INT(ldc);
	F77_FUNC(dsymm)(side_, uplo_, &m_, &n_, &alpha, a, &lda_, b, &ldb_,
			&beta, c, &ldc_);
}

void blas_dsyrk(enum blas_uplo uplo, enum blas_trans trans, size_t n, size_t k,
		double alpha, const double *a, size_t lda,
		double beta, double *c, size_t ldc)
{
	F77_UPLO(uplo);
	F77_TRANS(trans);
	F77_INT(n);
	F77_INT(k);
	F77_INT(lda);
	F77_INT(ldc);
	F77_FUNC(dsyrk)(uplo_, trans_, &n_, &k_, &alpha, a, &lda_,
			&beta, c, &ldc_);
}

void blas_dsyr2k(enum blas_uplo uplo, enum blas_trans trans, size_t n, size_t k,
		 double alpha, const double *a, size_t lda,
		 const double *b, size_t ldb, double beta, double *c,
		 size_t ldc)
{
	F77_UPLO(uplo);
	F77_TRANS(trans);
	F77_INT(n);
	F77_INT(k);
	F77_INT(lda);
	F77_INT(ldb);
	F77_INT(ldc);
	F77_FUNC(dsyr2k)(uplo_, trans_, &n_, &k_, &alpha, a, &lda_, b, &ldb_,
			&beta, c, &ldc_);
}

