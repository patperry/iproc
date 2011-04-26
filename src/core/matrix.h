#ifndef _MATRIX_H
#define _MATRIX_H

#include "vector.h"

typedef struct _iproc_matrix_view iproc_matrix_view;

struct matrix {
	double *data;
	ssize_t nrow;
	ssize_t ncol;
	ssize_t lda;
};

struct _iproc_matrix_view {
	struct matrix matrix;
};

enum trans_op {
	TRANS_NOTRANS,
	TRANS_TRANS,
	TRANS_CONJTRANS
};

bool matrix_init(struct matrix *a, ssize_t nrow, ssize_t ncol);
bool matrix_init_copy(struct matrix *a, const struct matrix *src);
void matrix_deinit(struct matrix *a);

struct matrix *matrix_alloc(ssize_t nrow, ssize_t ncol);
struct matrix *matrix_alloc_copy(const struct matrix *src);
void matrix_free(struct matrix *a);

void matrix_assign_copy(struct matrix *a, const struct matrix *src);
void matrix_assign_identity(struct matrix *a);
void matrix_fill(struct matrix *a, double value);


ssize_t matrix_nrow(const struct matrix *a);
ssize_t matrix_ncol(const struct matrix *a);
ssize_t matrix_lda(const struct matrix *a);

double matrix_get(const struct matrix *a, ssize_t i, ssize_t j);
void matrix_set(struct matrix *a, ssize_t i, ssize_t j, double val);
double *matrix_at(const struct matrix *a, ssize_t i, ssize_t j);


void matrix_add(struct matrix *a, const struct matrix *src);
void matrix_sub(struct matrix *a, const struct matrix *src);
void matrix_axpy(double alpha, const struct matrix *x, struct matrix *y);
void matrix_scale(struct matrix *a, double scale);
void matrix_scale_rows(struct matrix *a, const struct vector *scale);

void vector_init_matrix_col(struct vector *v, const struct matrix *a, ssize_t j);

iproc_matrix_view iproc_matrix_cols(const struct matrix * matrix,
				    int64_t j, int64_t n);
iproc_matrix_view iproc_matrix_submatrix(const struct matrix * a,
					 int64_t i,
					 int64_t j, int64_t nrow, int64_t ncol);
iproc_matrix_view iproc_matrix_view_array(const double *data,
					  int64_t nrow, int64_t ncol);
iproc_matrix_view iproc_matrix_view_array_with_lda(const double *data,
						   int64_t nrow,
						   int64_t ncol, int64_t lda);
iproc_matrix_view iproc_matrix_view_vector(const struct vector *vector,
					   int64_t nrow, int64_t ncol);

void matrix_mul(double alpha, enum trans_op trans, const struct matrix *a,
		const struct vector *x, double beta, struct vector *y);
void matrix_matmul(double alpha, enum trans_op trans, const struct matrix *a,
		   const struct matrix *x, double beta, struct matrix *y);
void matrix_update1(struct matrix *a,
		    double alpha, const struct vector *x,
		    const struct vector *y);

#endif /* _MATRIX_H */
