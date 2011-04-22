#ifndef _IPROC_MATRIX_H
#define _IPROC_MATRIX_H

#include <stdint.h>
#include "refcount.h"
#include "vector.h"

typedef struct _iproc_matrix iproc_matrix;
typedef struct _iproc_matrix_view iproc_matrix_view;

struct _iproc_matrix {
	double *data;
	int64_t nrow;
	int64_t ncol;
	int64_t lda;
	struct refcount refcount;
};

struct _iproc_matrix_view {
	iproc_matrix matrix;
};

typedef enum _iproc_trans {
	IPROC_TRANS_NOTRANS,
	IPROC_TRANS_TRANS,
	IPROC_TRANS_CONJTRANS
} iproc_trans;

iproc_matrix *iproc_matrix_new(int64_t nrow, int64_t ncol);
iproc_matrix *iproc_matrix_new_copy(const iproc_matrix * matrix);
iproc_matrix *iproc_matrix_ref(iproc_matrix * matrix);
void iproc_matrix_unref(iproc_matrix * matrix);
void iproc_matrix_copy(iproc_matrix * dst, const iproc_matrix * src);
void iproc_matrix_set_all(iproc_matrix * matrix, double value);
void iproc_matrix_set_identity(iproc_matrix * matrix);

int64_t iproc_matrix_nrow(const iproc_matrix * a);
int64_t iproc_matrix_ncol(const iproc_matrix * a);
int64_t iproc_matrix_lda(const iproc_matrix * a);

double iproc_matrix_get(const iproc_matrix * a, int64_t i, int64_t j);
void iproc_matrix_set(iproc_matrix * a, int64_t i, int64_t j, double val);
double *iproc_matrix_ptr(const iproc_matrix * a, int64_t i, int64_t j);

void iproc_matrix_add(iproc_matrix * dst_matrix, const iproc_matrix * matrix);
void iproc_matrix_sub(iproc_matrix * dst_matrix, const iproc_matrix * matrix);
void iproc_matrix_acc(iproc_matrix * dst_matrix,
		      double scale, const iproc_matrix * matrix);
void iproc_matrix_scale(iproc_matrix * matrix, double scale);
void iproc_matrix_scale_rows(iproc_matrix * matrix, const struct vector *scale);

void vector_init_matrix_col(struct vector *v,
			    const iproc_matrix * a, ssize_t j);
iproc_matrix_view iproc_matrix_cols(const iproc_matrix * matrix,
				    int64_t j, int64_t n);
iproc_matrix_view iproc_matrix_submatrix(const iproc_matrix * a,
					 int64_t i,
					 int64_t j, int64_t nrow, int64_t ncol);
iproc_matrix_view iproc_matrix_view_array(const double *data,
					  int64_t nrow, int64_t ncol);
iproc_matrix_view iproc_matrix_view_array_with_lda(const double *data,
						   int64_t nrow,
						   int64_t ncol, int64_t lda);
iproc_matrix_view iproc_matrix_view_vector(const struct vector *vector,
					   int64_t nrow, int64_t ncol);

void iproc_matrix_mul(double alpha,
		      iproc_trans trans,
		      const iproc_matrix * matrix,
		      const struct vector *x, double beta, struct vector *y);
void iproc_matrix_matmul(double alpha,
			 iproc_trans trans,
			 const iproc_matrix * matrix,
			 const iproc_matrix * x, double beta, iproc_matrix * y);

void iproc_matrix_update1(iproc_matrix * matrix,
			  double alpha, const struct vector *x, const struct vector *y);

#endif /* _IPROC_MATRIX_H */
