#include "port.h"

#include <assert.h>
#include <stdlib.h>
#include <strings.h>
#include "blas-private.h"
#include "util.h"
#include "matrix.h"

bool matrix_init(struct matrix *a, ssize_t nrow, ssize_t ncol)
{
	assert(a);
	assert(nrow >= 0);
	assert(ncol >= 0);
	assert(ncol == 0 || nrow <= SSIZE_MAX / ncol);

	if (array_init(&a->array, nrow * ncol, sizeof(double))) {
		a->nrow = nrow;
		a->ncol = ncol;
		a->lda = MAX(1, nrow);
		return true;
	}
	return false;
}

struct matrix *matrix_alloc(ssize_t nrow, ssize_t ncol)
{
	struct matrix *a;

	if ((a = malloc(sizeof(*a)))) {
		if (matrix_init(a, nrow, ncol)) {
			return a;
		}
		free(a);
	}
	return NULL;
}

bool matrix_init_copy(struct matrix *a, const struct matrix *src)
{
	assert(a);
	assert(src);

	ssize_t nrow = matrix_nrow(src);
	ssize_t ncol = matrix_ncol(src);
	if (matrix_init(a, nrow, ncol)) {
		matrix_assign_copy(a, src);
		return true;
	}
	return false;
}

struct matrix *matrix_alloc_copy(const struct matrix *src)
{
	struct matrix *a;

	if ((a = malloc(sizeof(*a)))) {
		if (matrix_init_copy(a, src)) {
			return a;
		}
		free(a);
	}
	return NULL;
}

void matrix_deinit(struct matrix *a)
{
	array_deinit(&a->array);
}

void matrix_free(struct matrix *a)
{
	if (a) {
		matrix_deinit(a);
		free(a);
	}
}

void matrix_init_view(struct matrix *a, const double *ptr, ssize_t m, ssize_t n)
{
	matrix_init_view_with_lda(a, ptr, m, n, MAX(1, m));
}

void matrix_init_view_with_lda(struct matrix *a, const double *ptr, ssize_t m,
			       ssize_t n, ssize_t lda)
{
	assert(a);
	assert(ptr || m == 0 || n == 0);
	assert(m >= 0);
	assert(n >= 0);
	assert(lda >= MAX(1, m));
	assert(n <= SSIZE_MAX / lda);

	ssize_t array_size = lda * n;
	array_init_view(&a->array, ptr, array_size, sizeof(double));
	a->nrow = m;
	a->ncol = n;
	a->lda = lda;

}

void matrix_init_view_vector(struct matrix *a, const struct vector *v,
			     ssize_t m, ssize_t n)
{
	assert(a);
	assert(v);
	assert(m >= 0);
	assert(n >= 0);
	assert(m == 0 || n <= vector_dim(v) / m);

	const double *ptr = m > 0 && n > 0 ? vector_front(v) : NULL;
	matrix_init_view(a, ptr, m, n);
}

void matrix_init_slice(struct matrix *a, const struct matrix *parent,
		       ssize_t i, ssize_t j, ssize_t m, ssize_t n)
{
	assert(a);
	assert(parent);
	assert(i >= 0);
	assert(j >= 0);
	assert(m >= 0);
	assert(n >= 0);
	assert(i <= matrix_nrow(parent) - m);
	assert(j <= matrix_ncol(parent) - n);

	const double *ptr = (m > 0 && n > 0) ? matrix_at(parent, i, j) : NULL;
	ssize_t lda = matrix_lda(parent);

	matrix_init_view_with_lda(a, ptr, m, n, lda);
}

void matrix_init_slice_cols(struct matrix *a, const struct matrix *parent,
			    ssize_t j, ssize_t n)
{
	assert(a);
	assert(parent);
	assert(j >= 0);
	assert(n >= 0);
	assert(j <= matrix_ncol(parent) - n);

	matrix_init_slice(a, parent, 0, j, matrix_nrow(parent), n);
}

ssize_t matrix_nrow(const struct matrix *a)
{
	assert(a);
	return a->nrow;
}

ssize_t matrix_ncol(const struct matrix *a)
{
	assert(a);
	return a->ncol;
}

ssize_t matrix_lda(const struct matrix *a)
{
	assert(a);
	return a->lda;
}

ssize_t matrix_size(const struct matrix *a)
{
	return matrix_nrow(a) * matrix_ncol(a);
}

bool matrix_empty(const struct matrix *a)
{
	return matrix_nrow(a) == 0 || matrix_ncol(a) == 0;
}

void matrix_fill(struct matrix *a, double value)
{
	assert(a);

	ssize_t m = matrix_nrow(a);
	ssize_t n = matrix_ncol(a);
	ssize_t i, j;

	for (j = 0; j < n; j++) {
		for (i = 0; i < m; i++) {
			matrix_set(a, i, j, value);
		}
	}
}

void matrix_set_identity(struct matrix *a)
{
	assert(a);

	int64_t m = matrix_nrow(a);
	int64_t n = matrix_ncol(a);
	int64_t mn = MIN(m, n);
	int64_t i;

	matrix_fill(a, 0.0);
	for (i = 0; i < mn; i++) {
		matrix_set(a, i, i, 1.0);
	}

}

void matrix_assign_copy(struct matrix *a, const struct matrix *src)
{
	assert(a);
	assert(src);
	assert(matrix_nrow(a) == matrix_nrow(src));
	assert(matrix_ncol(a) == matrix_ncol(src));

	f77int m = (f77int) matrix_nrow(a);
	f77int n = (f77int) matrix_ncol(a);
	f77int i, j;

	if (m == 0 || n == 0)
		return;

	if (matrix_lda(a) == m && matrix_lda(src) == m) {
		f77int mn = m * n;
		f77int one = 1;

		F77_FUNC(dcopy) (&mn, matrix_at(src, 0, 0),
				 &one, matrix_at(a, 0, 0), &one);
	} else {
		double value;

		for (j = 0; j < n; j++) {
			for (i = 0; i < m; i++) {
				value = matrix_get(src, i, j);
				matrix_set(a, i, j, value);
			}
		}
	}
}

double matrix_get(const struct matrix *a, ssize_t i, ssize_t j)
{
	assert(a);
	assert(0 <= i && i < matrix_nrow(a));
	assert(0 <= j && j < matrix_ncol(a));

	double *ptr = matrix_at(a, i, j);
	return *ptr;
}

void matrix_set(struct matrix *a, ssize_t i, ssize_t j, double val)
{
	assert(a);
	assert(0 <= i && i < matrix_nrow(a));
	assert(0 <= j && j < matrix_ncol(a));

	double *ptr = matrix_at(a, i, j);
	*ptr = val;
}

double *matrix_at(const struct matrix *a, ssize_t i, ssize_t j)
{
	assert(a);
	assert(0 <= i && i < matrix_nrow(a));
	assert(0 <= j && j < matrix_ncol(a));

	return array_at(&a->array, i + j * matrix_lda(a));
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

	f77int m = (f77int) matrix_nrow(y);
	f77int n = (f77int) matrix_ncol(y);
	f77int one = 1;
	f77int j;

	if (m == 0 || n == 0)
		return;

	if (matrix_lda(x) == m && matrix_lda(y) == m) {
		f77int mn = m * n;
		F77_FUNC(daxpy) (&mn, &alpha,
				 matrix_at(x, 0, 0), &one,
				 matrix_at(y, 0, 0), &one);
	} else {
		for (j = 0; j < n; j++) {
			F77_FUNC(daxpy) (&m, &alpha,
					 matrix_at(x, 0, j), &one,
					 matrix_at(y, 0, j), &one);
		}
	}

}

void matrix_scale(struct matrix *a, double scale)
{
	assert(a);

	struct vector col;
	int64_t n = matrix_ncol(a);
	int64_t j;

	for (j = 0; j < n; j++) {
		vector_init_matrix_col(&col, a, j);
		vector_scale(&col, scale);
	}
}

void matrix_scale_rows(struct matrix *a, const struct vector *scale)
{
	assert(a);
	assert(scale);
	assert(matrix_nrow(a) == vector_dim(scale));

	struct vector col;
	int64_t n = matrix_ncol(a);
	int64_t j;

	for (j = 0; j < n; j++) {
		vector_init_matrix_col(&col, a, j);
		vector_mul(&col, scale);
	}
}

void vector_init_matrix_col(struct vector *v, const struct matrix *a, ssize_t j)
{
	assert(0 <= j && j < matrix_ncol(a));
	ssize_t m = matrix_nrow(a);
	double *ptr = m != 0 ? matrix_at(a, 0, j) : NULL;

	vector_init_view(v, ptr, m);

}

void matrix_mul(double alpha, enum trans_op trans, const struct matrix *a,
		const struct vector *x, double beta, struct vector *y)
{
	assert(a);
	assert(x);
	assert(y);
	assert(trans != TRANS_NOTRANS || vector_dim(x) == matrix_ncol(a));
	assert(trans != TRANS_NOTRANS || vector_dim(y) == matrix_nrow(a));
	assert(trans == TRANS_NOTRANS || vector_dim(x) == matrix_nrow(a));
	assert(trans == TRANS_NOTRANS || vector_dim(y) == matrix_ncol(a));

	if (vector_empty(x)) {
		vector_scale(y, beta);
		return;
	} else if (vector_empty(y)) {
		return;
	}

	char *ptrans = (trans == TRANS_NOTRANS) ? "N" : "T";
	f77int m = (f77int) matrix_nrow(a);
	f77int n = (f77int) matrix_ncol(a);
	void *pa = matrix_at(a, 0, 0);
	f77int lda = (f77int) matrix_lda(a);
	void *px = vector_front(x);
	f77int incx = 1;
	void *py = vector_front(y);
	f77int incy = 1;

	F77_FUNC(dgemv) (ptrans, &m, &n, &alpha, pa, &lda,
			 px, &incx, &beta, py, &incy);
}

void matrix_matmul(double alpha, enum trans_op trans, const struct matrix *a,
		   const struct matrix *x, double beta, struct matrix *y)
{
	assert(a);
	assert(x);
	assert(y);
	assert(trans != TRANS_NOTRANS || matrix_nrow(x) == matrix_ncol(a));
	assert(trans != TRANS_NOTRANS || matrix_nrow(y) == matrix_nrow(a));
	assert(trans == TRANS_NOTRANS || matrix_nrow(x) == matrix_nrow(a));
	assert(trans == TRANS_NOTRANS || matrix_nrow(y) == matrix_ncol(a));
	assert(matrix_ncol(x) == matrix_ncol(y));

	char *ptransa = (trans == TRANS_NOTRANS) ? "N" : "T";
	char *ptransb = "N";
	f77int m = (f77int) matrix_nrow(y);
	f77int n = (f77int) matrix_ncol(y);
	f77int k = (f77int) ((trans == TRANS_NOTRANS)
			     ? matrix_ncol(a)
			     : matrix_nrow(a));

	if (k == 0) {
		matrix_scale(y, beta);
		return;
	} else if (m == 0 || n == 0) {
		return;
	}

	void *pa = matrix_at(a, 0, 0);
	f77int lda = (f77int) matrix_lda(a);
	void *pb = matrix_at(x, 0, 0);
	f77int ldb = (f77int) matrix_lda(x);
	void *pc = matrix_at(y, 0, 0);
	f77int ldc = (f77int) matrix_lda(y);

	F77_FUNC(dgemm) (ptransa, ptransb, &m, &n, &k,
			 &alpha, pa, &lda, pb, &ldb, &beta, pc, &ldc);
}

void
matrix_update1(struct matrix *a,
	       double alpha, const struct vector *x, const struct vector *y)
{
	assert(a);
	assert(x);
	assert(y);
	assert(matrix_nrow(a) == vector_dim(x));
	assert(matrix_ncol(a) == vector_dim(y));

	if (vector_empty(x) || vector_empty(y))
		return;

	f77int m = (f77int) matrix_nrow(a);
	f77int n = (f77int) matrix_ncol(a);
	void *px = vector_front(x);
	f77int incx = 1;
	void *py = vector_front(y);
	f77int incy = 1;
	void *pa = matrix_at(a, 0, 0);
	f77int lda = (f77int) matrix_lda(a);

	F77_FUNC(dger) (&m, &n, &alpha, px, &incx, py, &incy, pa, &lda);
}
