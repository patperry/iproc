#ifndef _MATRIX_H
#define _MATRIX_H

#include "svector.h"
#include "vector.h"

struct matrix {
	struct vector data;
	ssize_t nrow;
	ssize_t ncol;
	ssize_t lda;
};

enum trans_op {
	TRANS_NOTRANS,
	TRANS_TRANS,
	TRANS_CONJTRANS
};

/* create, destroy */
void matrix_init(struct matrix *a, ssize_t nrow, ssize_t ncol);
void matrix_reinit(struct matrix *a, ssize_t nrow, ssize_t ncol);
void matrix_init_copy(struct matrix *a, const struct matrix *src);
void matrix_assign_copy(struct matrix *a, const struct matrix *src);
void matrix_deinit(struct matrix *a);

/* views */
static inline struct matrix matrix_make(const struct vector *v, ssize_t m,
					ssize_t n);
static inline struct matrix matrix_make_with_lda(const struct vector *v,
						 ssize_t m, ssize_t n,
						 ssize_t lda);
static inline struct matrix matrix_slice(struct matrix *a, ssize_t i, ssize_t j,
					 ssize_t m, ssize_t n);
static inline struct matrix matrix_slice_cols(struct matrix *a, ssize_t j,
					      ssize_t n);
static inline struct vector matrix_col(const struct matrix *a, ssize_t j);

/* properties */
static inline ssize_t matrix_nrow(const struct matrix *a);
static inline ssize_t matrix_ncol(const struct matrix *a);
static inline ssize_t matrix_lda(const struct matrix *a);
static inline ssize_t matrix_count(const struct matrix *a);
static inline double *matrix_to_ptr(const struct matrix *a);

static inline double matrix_item(const struct matrix *a, ssize_t i, ssize_t j);
static inline double *matrix_item_ptr(const struct matrix *a, ssize_t i,
				      ssize_t j);
static inline void matrix_set_item(struct matrix *a, ssize_t i, ssize_t j,
				   double val);

/* assign, copy, fill */
void matrix_fill(struct matrix *a, double value);
void matrix_assign_identity(struct matrix *a);

void matrix_fill_col(struct matrix *a, ssize_t j, double val);
void matrix_fill_row(struct matrix *a, ssize_t i, double val);
void matrix_set_row(struct matrix *a, ssize_t i, const double *src);
void matrix_get_row(const struct matrix *a, ssize_t i, double *dst);

/* arithmetic */
void matrix_scale(struct matrix *a, double scale);
void matrix_scale_rows(struct matrix *a, const struct vector *scale);
void matrix_add(struct matrix *a, const struct matrix *src);
void matrix_sub(struct matrix *a, const struct matrix *src);
void matrix_axpy(double alpha, const struct matrix *x, struct matrix *y);

void matrix_row_axpy(double alpha, const struct matrix *x, ssize_t i,
		     struct vector *y);
void matrix_col_axpy(double alpha, const struct matrix *x, ssize_t j,
		     struct vector *y);

/* linear algebra */
void matrix_mul(double alpha, enum trans_op trans, const struct matrix *a,
		const struct vector *x, double beta, struct vector *y);
void matrix_matmul(double alpha, enum trans_op trans, const struct matrix *a,
		   const struct matrix *x, double beta, struct matrix *y);
void matrix_muls(double alpha, enum trans_op trans, const struct matrix *a,
		 const struct svector *x, double beta, struct vector *y);
void matrix_update1(struct matrix *a,
		    double alpha, const struct vector *x,
		    const struct vector *y);

/* deprecated */
struct matrix *matrix_alloc(ssize_t nrow, ssize_t ncol);
struct matrix *matrix_alloc_copy(const struct matrix *src);
void matrix_free(struct matrix *a);

/* inline function definitions */
struct matrix matrix_make(const struct vector *v, ssize_t m, ssize_t n)
{
	return matrix_make_with_lda(v, m, n, MAX(1, m));
}

struct matrix matrix_make_with_lda(const struct vector *v, ssize_t m, ssize_t n,
				   ssize_t lda)
{
	assert(v);
	assert(m >= 0);
	assert(n >= 0);
	assert(lda >= MAX(1, m));
	assert(n <= SSIZE_MAX / lda);
	assert(vector_dim(v) == lda * n || (m == 0 && vector_dim(v) == 0));

	struct matrix a;
	a.data = *v;
	a.nrow = m;
	a.ncol = n;
	a.lda = lda;
	return a;
}

struct matrix matrix_slice(struct matrix *a, ssize_t i, ssize_t j, ssize_t m,
			   ssize_t n)
{
	assert(a);
	assert(i >= 0);
	assert(j >= 0);
	assert(m >= 0);
	assert(n >= 0);
	assert(i <= matrix_nrow(a) - m);
	assert(j <= matrix_ncol(a) - n);

	const double *ptr = (m > 0 && n > 0) ? matrix_item_ptr(a, i, j) : NULL;
	ssize_t lda = matrix_lda(a);
	struct vector v =
	    m > 0 ? vector_make(ptr, lda * n) : vector_make(ptr, 0);

	return matrix_make_with_lda(&v, m, n, lda);
}

struct matrix matrix_slice_cols(struct matrix *a, ssize_t j, ssize_t n)
{
	assert(a);
	assert(j >= 0);
	assert(n >= 0);
	assert(j <= matrix_ncol(a) - n);

	return matrix_slice(a, 0, j, matrix_nrow(a), n);
}

struct vector matrix_col(const struct matrix *a, ssize_t j)
{
	assert(0 <= j && j < matrix_ncol(a));

	ssize_t m = matrix_nrow(a);
	double *ptr = m != 0 ? matrix_item_ptr(a, 0, j) : NULL;

	return vector_make(ptr, m);
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

ssize_t matrix_count(const struct matrix *a)
{
	assert(a);
	return matrix_nrow(a) * matrix_ncol(a);
}

double *matrix_to_ptr(const struct matrix *a)
{
	assert(a);
	return vector_to_ptr(&a->data);
}

double matrix_item(const struct matrix *a, ssize_t i, ssize_t j)
{
	assert(a);
	assert(0 <= i && i < matrix_nrow(a));
	assert(0 <= j && j < matrix_ncol(a));

	double *ptr = matrix_item_ptr(a, i, j);
	return *ptr;
}

double *matrix_item_ptr(const struct matrix *a, ssize_t i, ssize_t j)
{
	assert(a);
	assert(0 <= i && i < matrix_nrow(a));
	assert(0 <= j && j < matrix_ncol(a));

	return vector_item_ptr(&a->data, i + j * matrix_lda(a));
}

void matrix_set_item(struct matrix *a, ssize_t i, ssize_t j, double val)
{
	assert(a);
	assert(0 <= i && i < matrix_nrow(a));
	assert(0 <= j && j < matrix_ncol(a));

	double *ptr = matrix_item_ptr(a, i, j);
	*ptr = val;
}

#endif /* _MATRIX_H */
