#ifndef _MATRIX_H
#define _MATRIX_H

#include "array.h"
#include "vector.h"

struct matrix {
	struct array array;
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
bool matrix_init(struct matrix *a, ssize_t nrow, ssize_t ncol);
bool matrix_init_copy(struct matrix *a, const struct matrix *src);
void matrix_deinit(struct matrix *a);

/* views */
void matrix_init_view(struct matrix *a, const double *ptr, ssize_t m,
		      ssize_t n);
void matrix_init_view_with_lda(struct matrix *a, const double *ptr, ssize_t m,
			       ssize_t n, ssize_t lda);
void matrix_init_view_vector(struct matrix *a, const struct vector *v,
			     ssize_t m, ssize_t n);
void matrix_init_slice(struct matrix *a, const struct matrix *parent, ssize_t i,
		       ssize_t j, ssize_t m, ssize_t n);
void matrix_init_slice_cols(struct matrix *a, const struct matrix *parent,
			    ssize_t j, ssize_t n);

void vector_init_matrix_col(struct vector *v, const struct matrix *a,
			    ssize_t j);

/* assign, copy, fill */
void matrix_assign_copy(struct matrix *a, const struct matrix *src);
void matrix_fill(struct matrix *a, double value);
void matrix_set_identity(struct matrix *a);

/* indexing */
double matrix_get(const struct matrix *a, ssize_t i, ssize_t j);
void matrix_set(struct matrix *a, ssize_t i, ssize_t j, double val);
double *matrix_at(const struct matrix *a, ssize_t i, ssize_t j);

/* informative */
ssize_t matrix_nrow(const struct matrix *a);
ssize_t matrix_ncol(const struct matrix *a);
ssize_t matrix_lda(const struct matrix *a);
ssize_t matrix_size(const struct matrix *a);
bool matrix_empty(const struct matrix *a);

/* arithmetic */
void matrix_scale(struct matrix *a, double scale);
void matrix_scale_rows(struct matrix *a, const struct vector *scale);
void matrix_add(struct matrix *a, const struct matrix *src);
void matrix_sub(struct matrix *a, const struct matrix *src);
void matrix_axpy(double alpha, const struct matrix *x, struct matrix *y);

/* linear algebra */
void matrix_mul(double alpha, enum trans_op trans, const struct matrix *a,
		const struct vector *x, double beta, struct vector *y);
void matrix_matmul(double alpha, enum trans_op trans, const struct matrix *a,
		   const struct matrix *x, double beta, struct matrix *y);
void matrix_update1(struct matrix *a,
		    double alpha, const struct vector *x,
		    const struct vector *y);

/* deprecated */
struct matrix *matrix_alloc(ssize_t nrow, ssize_t ncol);
struct matrix *matrix_alloc_copy(const struct matrix *src);
void matrix_free(struct matrix *a);

#endif /* _MATRIX_H */
