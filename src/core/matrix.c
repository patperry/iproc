#include "port.h"

#include <assert.h>
#include <stdlib.h>
#include <stdio.h>
#include <strings.h>
#include "blas.h"
#include "xalloc.h"
#include "matrix.h"

void matrix_init(struct matrix *a, ssize_t nrow, ssize_t ncol)
{
	assert(a);
	assert(nrow >= 0);
	assert(ncol >= 0);
	assert(ncol == 0 || nrow <= SSIZE_MAX / ncol);

	a->data = xcalloc(nrow * ncol, sizeof(double));
	a->nrow = nrow;
	a->ncol = ncol;
	a->lda = MAX(1, nrow);
}

void matrix_reinit(struct matrix *a, ssize_t nrow, ssize_t ncol)
{
	assert(a);
	assert(nrow >= 0);
	assert(ncol >= 0);
	assert(ncol == 0 || nrow <= SSIZE_MAX / ncol);

	a->data = xrealloc(a->data, nrow * ncol * sizeof(double));
	a->nrow = nrow;
	a->ncol = ncol;
	a->lda = MAX(1, nrow);
}

struct matrix *matrix_alloc(ssize_t nrow, ssize_t ncol)
{
	struct matrix *a = xcalloc(1, sizeof(*a));
	matrix_init(a, nrow, ncol);
	return a;
}

void matrix_init_copy(struct matrix *a, enum blas_trans trans,
		      const struct matrix *src)
{
	assert(a);
	assert(src);

	ssize_t nrow, ncol;

	if (trans == BLAS_NOTRANS) {
		nrow = matrix_nrow(src);
		ncol = matrix_ncol(src);
	} else {
		nrow = matrix_ncol(src);
		ncol = matrix_nrow(src);
	}

	matrix_init(a, nrow, ncol);
	matrix_assign_copy(a, trans, src);
}

void matrix_assign_copy(struct matrix *a, enum blas_trans trans,
			const struct matrix *src)
{
	assert(a);
	assert(src);
	assert(trans != BLAS_NOTRANS || matrix_nrow(a) == matrix_nrow(src));
	assert(trans != BLAS_NOTRANS || matrix_ncol(a) == matrix_ncol(src));
	assert(trans == BLAS_NOTRANS || matrix_nrow(a) == matrix_ncol(src));
	assert(trans == BLAS_NOTRANS || matrix_ncol(a) == matrix_nrow(src));

	size_t m = matrix_nrow(a);
	size_t n = matrix_ncol(a);
	size_t i, j;
	double val;

	if (m == 0 || n == 0)
		return;

	if (trans == BLAS_NOTRANS) {
		if ((size_t)matrix_lda(a) == m && (size_t)matrix_lda(src) == m) {
			blas_dcopy(m * n, matrix_to_ptr(src), 1,
			      	   matrix_to_ptr(a), 1);
		} else {
			for (j = 0; j < n; j++) {
				for (i = 0; i < m; i++) {
					val = matrix_item(src, i, j);
					matrix_set_item(a, i, j, val);
				}
			}
		}
	} else {
		for (j = 0; j < n; j++) {
			for (i = 0; i < m; i++) {
				val = matrix_item(src, j, i);
				matrix_set_item(a, i, j, val);
			}
		}
	}
}

struct matrix *matrix_alloc_copy(enum blas_trans trans, const struct matrix *src)
{
	struct matrix *a = xcalloc(1, sizeof(*a));
	matrix_init_copy(a, trans, src);
	return a;
}

void matrix_deinit(struct matrix *a)
{
	free(a->data);
}

void matrix_free(struct matrix *a)
{
	if (a) {
		matrix_deinit(a);
		free(a);
	}
}

void matrix_fill(struct matrix *a, double value)
{
	assert(a);

	ssize_t m = matrix_nrow(a);
	ssize_t n = matrix_ncol(a);
	ssize_t i, j;

	for (j = 0; j < n; j++) {
		for (i = 0; i < m; i++) {
			matrix_set_item(a, i, j, value);
		}
	}
}

void matrix_assign_identity(struct matrix *a)
{
	assert(a);

	ssize_t m = matrix_nrow(a);
	ssize_t n = matrix_ncol(a);
	ssize_t mn = MIN(m, n);
	ssize_t i;

	matrix_fill(a, 0.0);
	for (i = 0; i < mn; i++) {
		matrix_set_item(a, i, i, 1.0);
	}

}

void matrix_fill_col(struct matrix *a, ssize_t j, double val)
{
	assert(a);
	assert(0 <= j && j < matrix_ncol(a));

	ssize_t i, m = matrix_nrow(a);

	for (i = 0; i < m; i++) {
		matrix_set_item(a, i, j, val);
	}
}

void matrix_fill_row(struct matrix *a, ssize_t i, double val)
{
	assert(a);
	assert(0 <= i && i < matrix_nrow(a));

	ssize_t j, n = matrix_ncol(a);

	for (j = 0; j < n; j++) {
		matrix_set_item(a, i, j, val);
	}
}

void matrix_set_row(struct matrix *a, ssize_t i, const double *src)
{
	assert(a);
	assert(0 <= i && i < matrix_nrow(a));
	assert(src || matrix_ncol(a) == 0);

	size_t n = (size_t)matrix_ncol(a);
	const double *x = src;
	size_t incx = 1;
	double *y = matrix_item_ptr(a, i, 0);
	size_t incy = (size_t)matrix_lda(a);

	blas_dcopy(n, x, incx, y, incy);
}

void matrix_get_row(const struct matrix *a, ssize_t i, double *dst)
{
	assert(a);
	assert(0 <= i && i < matrix_nrow(a));
	assert(dst || matrix_ncol(a) == 0);

	size_t n = (size_t)matrix_ncol(a);
	const double *x = matrix_item_ptr(a, i, 0);
	size_t incx = (size_t)matrix_lda(a);
	double *y = dst;
	size_t incy = 1;

	blas_dcopy(n, x, incx, y, incy);
}

void matrix_axpy_row(double alpha, const struct matrix *x, ssize_t i,
		     double *y)
{
	assert(x);
	assert(0 <= i && i < matrix_nrow(x));

	size_t n = (size_t)matrix_ncol(x);
	const double *px = matrix_item_ptr(x, i, 0);
	size_t incx = (size_t)matrix_lda(x);
	size_t incy = 1;

	blas_daxpy(n, alpha, px, incx, y, incy);
}

void matrix_axpy_col(double alpha, const struct matrix *x, ssize_t j,
		     double *y)
{
	assert(x);
	assert(0 <= j && j < matrix_ncol(x));

	size_t n = (size_t)matrix_nrow(x);
	const double *px = matrix_item_ptr(x, 0, j);
	size_t incx = 1;
	size_t incy = 1;

	blas_daxpy(n, alpha, px, incx, y, incy);
}

void matrix_fill_diag(struct matrix *a, ssize_t i, double val)
{
	assert(a);
	assert(0 <= i && i < matrix_nrow(a));

	ssize_t j, n = matrix_ncol(a);

	for (j = 0; j < n; j++) {
		matrix_set_item(a, i, j, val);
	}
}

void matrix_set_diag(struct matrix *a, ssize_t i, const double *src)
{
	assert(a);
	assert(-matrix_nrow(a) < i && i < matrix_ncol(a));
	assert(src);

	size_t n = (size_t)matrix_diag_dim(a, i);
	double *y = (i >= 0 ? matrix_item_ptr(a, 0, i)
		     : matrix_item_ptr(a, -i, 0));
	size_t incy = (size_t)(matrix_lda(a) + 1);
	const double *x = src;
	size_t incx = 1;

	blas_dcopy(n, x, incx, y, incy);
}

void matrix_get_diag(const struct matrix *a, ssize_t i, double *dst)
{
	assert(a);
	assert(dst);
	assert(-matrix_nrow(a) < i && i < matrix_ncol(a));

	size_t n = (size_t)matrix_diag_dim(a, i);
	const double *x = (i >= 0 ? matrix_item_ptr(a, 0, i)
			   : matrix_item_ptr(a, -i, 0));
	size_t incx = (size_t)(matrix_lda(a) + 1);
	double *y = dst;
	size_t incy = 1;

	blas_dcopy(n, x, incx, y, incy);
}

void matrix_axpy_diag(double alpha, const struct matrix *x, ssize_t i,
		      double *y)
{
	assert(x);
	assert(y);
	assert(-matrix_nrow(x) < i && i < matrix_ncol(x));

	size_t n = (size_t)matrix_diag_dim(x, i);
	const double *px = (i >= 0 ? matrix_item_ptr(x, 0, i)
			    : matrix_item_ptr(x, -i, 0));
	size_t incx = (size_t)(matrix_lda(x) + 1);
	size_t incy = 1;

	blas_daxpy(n, alpha, px, incx, y, incy);
}

void matrix_add(struct matrix *a, const struct matrix *src)
{
	matrix_axpy(+1.0, src, a);
}

void matrix_sub(struct matrix *a, const struct matrix *src)
{
	matrix_axpy(-1.0, src, a);
}

void matrix_axpy(double alpha, const struct matrix *x, struct matrix *y)
{
	assert(y);
	assert(x);
	assert(matrix_nrow(y) == matrix_nrow(x));
	assert(matrix_ncol(y) == matrix_ncol(x));

	size_t m = (size_t)matrix_nrow(y);
	size_t n = (size_t)matrix_ncol(y);
	size_t j;

	if (m == 0 || n == 0)
		return;

	if ((size_t)matrix_lda(x) == m && (size_t)matrix_lda(y) == m) {
		size_t mn = m * n;
		blas_daxpy(mn, alpha, matrix_to_ptr(x), 1,
			   matrix_to_ptr(y), 1);
	} else {
		for (j = 0; j < n; j++) {
			blas_daxpy(m, alpha, matrix_item_ptr(x, 0, j), 1,
				   matrix_item_ptr(y, 0, j), 1);
		}
	}

}

void matrix_scale(struct matrix *a, double scale)
{
	assert(a);

	size_t m = matrix_nrow(a);
	size_t n = matrix_ncol(a);
	double *col;
	size_t j;

	for (j = 0; j < n; j++) {
		col = matrix_col(a, j);
		blas_dscal(m, scale, col, 1);
	}
}

void matrix_scale_rows(struct matrix *a, const double *scale)
{
	double *col;
	size_t m = matrix_nrow(a);
	size_t n = matrix_ncol(a);
	size_t j;

	for (j = 0; j < n; j++) {
		col = matrix_col(a, j);
		blas_dtbmv(BLAS_UPPER, BLAS_NOTRANS, BLAS_NONUNIT,
			   m, 1, scale, 1, col, 1);
	}
}

void matrix_div_rows(struct matrix *a, const double *scale)
{
	double *col;
	size_t m = matrix_nrow(a);
	size_t n = matrix_ncol(a);
	size_t j;

	for (j = 0; j < n; j++) {
		col = matrix_col(a, j);
		blas_dtbsv(BLAS_UPPER, BLAS_NOTRANS, BLAS_NONUNIT,
			   m, 1, scale, 1, col, 1);
	}
}

void matrix_scale_cols(struct matrix *a, const double *scale)
{
	double *col;
	double alpha;
	size_t m = matrix_nrow(a);
	size_t n = matrix_ncol(a);
	size_t j;

	for (j = 0; j < n; j++) {
		col = matrix_col(a, j);
		alpha = scale[j];
		blas_dscal(m, alpha, col, 1);
	}
}

void matrix_div_cols(struct matrix *a, const double *scale)
{
	double *col;
	double alpha;
	size_t m = matrix_nrow(a);
	size_t n = matrix_ncol(a);
	size_t j;

	for (j = 0; j < n; j++) {
		col = matrix_col(a, j);
		alpha = scale[j];
		blas_dscal(m, 1.0 / alpha, col, 1);
	}
}

void matrix_mul(double alpha, enum blas_trans trans, const struct matrix *a,
		const double *x, double beta, double *y)
{
	size_t m = (size_t)matrix_nrow(a);
	size_t n = (size_t)matrix_ncol(a);
	void *pa = matrix_to_ptr(a);
	size_t lda = (size_t)matrix_lda(a);
	size_t incx = 1;
	size_t incy = 1;
	size_t nx = trans == BLAS_NOTRANS ? n : m;
	size_t ny = trans == BLAS_NOTRANS ? m : n;

	if (!nx) {
		blas_dscal(ny, beta, y, 1);
		return;
	}

	blas_dgemv(trans, m, n, alpha, pa, lda, x, incx, beta, y, incy);
}

void matrix_matmul(double alpha, enum blas_trans trans, const struct matrix *a,
		   const struct matrix *x, double beta, struct matrix *y)
{
	assert(a);
	assert(x);
	assert(y);
	assert(trans != BLAS_NOTRANS || matrix_nrow(x) == matrix_ncol(a));
	assert(trans != BLAS_NOTRANS || matrix_nrow(y) == matrix_nrow(a));
	assert(trans == BLAS_NOTRANS || matrix_nrow(x) == matrix_nrow(a));
	assert(trans == BLAS_NOTRANS || matrix_nrow(y) == matrix_ncol(a));
	assert(matrix_ncol(x) == matrix_ncol(y));

	char transb = BLAS_NOTRANS;
	size_t m = (size_t)matrix_nrow(y);
	size_t n = (size_t)matrix_ncol(y);
	size_t k = (size_t)((trans == BLAS_NOTRANS)
			    ? matrix_ncol(a)
			    : matrix_nrow(a));

	if (k == 0) {
		matrix_scale(y, beta);
		return;
	} else if (m == 0 || n == 0) {
		return;
	}

	void *pa = matrix_item_ptr(a, 0, 0);
	size_t lda = (size_t)matrix_lda(a);
	void *pb = matrix_item_ptr(x, 0, 0);
	size_t ldb = (size_t)matrix_lda(x);
	void *pc = matrix_item_ptr(y, 0, 0);
	size_t ldc = (size_t)matrix_lda(y);

	blas_dgemm(trans, transb, m, n, k, alpha, pa, lda, pb, ldb, beta,
		   pc, ldc);
}

void
matrix_update1(struct matrix *a,
	       double alpha, const double *x, const double *y)
{
	size_t m = (size_t)matrix_nrow(a);
	size_t n = (size_t)matrix_ncol(a);
	size_t incx = 1;
	size_t incy = 1;
	void *pa = matrix_item_ptr(a, 0, 0);
	size_t lda = (size_t)matrix_lda(a);

	blas_dger(m, n, alpha, x, incx, y, incy, pa, lda);
}

void matrix_sym_update1(enum blas_uplo uplo, struct matrix *a, double alpha,
			const double *x)
{
	size_t n = (size_t)matrix_ncol(a);
	size_t incx = 1;
	void *pa = matrix_to_ptr(a);
	size_t lda = (size_t)matrix_lda(a);	
	
	blas_dsyr(uplo, n, alpha, x, incx, pa, lda);
}

void matrix_sym_update2(enum blas_uplo uplo, struct matrix *a, double alpha,
			const double *x, const double *y)
{
	assert(matrix_nrow(a) == matrix_ncol(a));

	size_t n = (size_t)matrix_ncol(a);
	size_t incx = 1;
	size_t incy = 1;
	void *pa = matrix_to_ptr(a);
	size_t lda = (size_t)matrix_lda(a);	
	
	blas_dsyr2(uplo, n, alpha, x, incx, y, incy, pa, lda);
}

void matrix_printf(const struct matrix *a)
{
	ssize_t m = matrix_nrow(a);
	ssize_t n = matrix_ncol(a);
	ssize_t i, j;

	printf("a <- matrix(c(");
	for (j = 0; j < n; j++) {
		for (i = 0; i < m; i++) {
			double x = matrix_item(a, i, j);
			if (i == 0 && j == 0) {
				printf("%.8e", x);
			} else {
				printf(", %.8e", x);
			}
		}
	}
	printf("), %" SSIZE_FMT ", %" SSIZE_FMT ")\n", m, n);
}

void matrix_insert_col(struct matrix *a, ssize_t j)
{
	assert(a);
	assert(0 <= j && j <= matrix_ncol(a));
	
	ssize_t m = matrix_nrow(a);
	ssize_t n0 = matrix_ncol(a);
	ssize_t n1 = n0 + 1;
	
	matrix_reinit(a, m, n1);
	ssize_t j1;
	for (j1 = n1 - 1; j1 > j; j1--) {
		double *dst = matrix_col(a, j1);
		const double *src = matrix_col(a, j1 - 1);
		blas_dcopy(m, src, 1, dst, 1);
	}
	double *col = matrix_col(a, j);
	memset(col, 0, m * sizeof(*col));
}

void matrix_insert_row(struct matrix *a, ssize_t i)
{
	assert(a);
	assert(0 <= i && i <= matrix_nrow(a));
	
	ssize_t m0 = matrix_nrow(a);
	ssize_t m1 = m0 + 1;
	ssize_t n = matrix_ncol(a);
	
	matrix_reinit(a, m1, n);
	double *data = matrix_to_ptr(a);
	double *src = data + m0 * n;	
	double *dst = data + m1 * n;
	
	ssize_t j;
	for (j = 0; j < n; j++) {
		ssize_t k;
		for (k = 0; k < m0 - i; k++) {
			*(--dst) = *(--src);
		}
		*(--dst) = 0.0;
		for (k = 0; k < i; k++) {
			*(--dst) = *(--src);
		}
	}
	assert(dst == data);
	assert(src == data);
}

void matrix_insert_row_col(struct matrix *a, ssize_t i, ssize_t j)
{
	assert(a);
	assert(0 <= i && i <= matrix_nrow(a));
	assert(0 <= j && j <= matrix_ncol(a));
	
	size_t m0 = matrix_nrow(a);
	size_t m1 = m0 + 1;
	size_t n0 = matrix_ncol(a);
	size_t n1 = n0 + 1;
	
	matrix_reinit(a, m1, n1);
	double *data = matrix_to_ptr(a);
	double *src = data + m0 * n0;	
	double *dst = data + m1 * n1;
	
	ssize_t k, l;
	for (l = 0; l < n1; l++) {
		if (l == n0 - j) {
			for (k = 0; k < m1; k++) {
				*(--dst) = 0.0;
			}
			assert(dst >= src);
		} else {
			for (k = 0; k < m0 - i; k++) {
				*(--dst) = *(--src);
			}
			*(--dst) = 0.0;
			for (k = 0; k < i; k++) {
				*(--dst) = *(--src);
			}
		}
	}
	assert(dst == data);
	assert(src == data);
}



