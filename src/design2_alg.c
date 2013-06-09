#include "port.h"
#include <string.h>
#include "design2.h"



void design2_mul(double alpha, const struct design2 *d, size_t i,
		 const struct coefs2 *c, double beta, double *y)
{
	design2_traits_mul(alpha, d, i, c->traits, beta, y);
	design2_tvars_mul(alpha, d, i, c->tvars, 1.0, y);
}


void design2_tmul(double alpha, const struct design2 *d, size_t i,
		  const double *x, double beta, struct coefs2 *c)
{
	design2_traits_tmul(alpha, d, i, x, beta, c->traits);
	design2_tvars_tmul(alpha, d, i, x, beta, c->tvars);
}


void design2_axpy(double alpha, const struct design2 *d, size_t i, size_t j,
		  struct coefs2 *c)
{
	design2_traits_axpy(alpha, d, i, j, c->traits);
	design2_tvars_axpy(alpha, d, i, j, c->tvars);
}


void design2_traits_mul(double alpha, const struct design2 *d, size_t i,
			const double *x, double beta, double *y)
{
	assert(i < design2_count1(d));

	size_t n = design2_count2(d);
	size_t p = design2_trait_dim(d);
	const double *a = design2_trait_matrix(d, i);
	size_t lda = p;

	if (p) {
		blas_dgemv(BLAS_TRANS, p, n, alpha, a, lda, x, 1, beta, y, 1);
	} else if (beta == 0.0) {
		memset(y, 0, n * sizeof(*y));
	} else if (beta != 1.0) {
		blas_dscal(n, beta, y, 1);
	}
}


void design2_traits_tmul(double alpha, const struct design2 *d, size_t i, const double *x, double beta, double *y)
{
	assert(i < design2_count1(d));

	size_t n = design2_count2(d);
	size_t p = design2_trait_dim(d);
	const double *a = design2_trait_matrix(d, i);
	size_t lda = p;

	if (p) {
		blas_dgemv(BLAS_NOTRANS, p, n, alpha, a, lda, x, 1, beta, y, 1);
	}

}


void design2_traits_axpy(double alpha, const struct design2 *d, size_t i, size_t j, double *y)
{
	assert(i < design2_count1(d));
	assert(j < design2_count2(d));

	size_t p = design2_trait_dim(d);
	const double *x = design2_traits(d, i, j);

	blas_daxpy(p, alpha, x, 1, y, 1);
}


void design2_tvars_mul(double alpha, const struct design2 *d, size_t i,
		       const double *x, double beta, double *y)
{
	assert(i < design2_count1(d));

	const double *a;
	const size_t *j;
	size_t n = design2_count2(d);
	size_t dim = design2_tvar_dim(d);
	size_t nz;

	if (beta == 0.0) {
		memset(y, 0, n * sizeof(*y));
	} else if (beta != 1.0) {
		blas_dscal(n, beta, y, 1);
	}

	design2_get_tvar_matrix(d, i, &a, &j, &nz);
	for (; nz != 0; a += dim, j++, nz--) {
		y[*j] += alpha * blas_ddot(dim, a, 1, x, 1);
	}
}


void design2_tvars_tmul(double alpha, const struct design2 *d, size_t i, const double *x, double beta, double *y)
{
	assert(i < design2_count1(d));

	size_t dim = design2_tvar_dim(d);
	const double *a;
	const size_t *j;
	size_t nz;

	if (beta == 0.0) {
		memset(y, 0, dim * sizeof(*y));
	} else if (beta != 1.0) {
		blas_dscal(dim, beta, y, 1);
	}

	design2_get_tvar_matrix(d, i, &a, &j, &nz);
	for (; nz != 0; a += dim, j++, nz--) {
		blas_daxpy(dim, alpha * x[*j], a, 1, y, 1);
	}
}


void design2_tvars_axpy(double alpha, const struct design2 *d, size_t i, size_t j, double *y)
{
	assert(i < design2_count1(d));
	assert(j < design2_count2(d));

	const double *dx = design2_tvars(d, i, j);

	if (dx) {
		size_t dim = design2_tvar_dim(d);
		blas_daxpy(dim, alpha, dx, 1, y, 1);
	}
}
