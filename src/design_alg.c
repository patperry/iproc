#include "port.h"

#include <assert.h>
#include <inttypes.h>
#include <search.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>

#include "coreutil.h"
#include "lapack.h"
#include "sblas.h"
#include "strata.h"
#include "util.h"
#include "xalloc.h"

#include "design.h"


void design_traits_mul(double alpha, const struct design *d,
		       const double *x, double beta, double *y)
{
	size_t n = design_count(d);
	size_t p = design_trait_dim(d);
	const double *a = design_trait_matrix(d);
	size_t lda = p;

	if (p) {
		blas_dgemv(BLAS_TRANS, p, n, alpha, a, lda, x, 1, beta, y, 1);
	} else if (beta == 0) {
		memset(y, 0, n * sizeof(*y));
	} else if (beta != 1) {
		blas_dscal(n, beta, y, 1);
	}
}


void design_traits_tmul(double alpha, const struct design *d, const double *x, double beta, double *y)
{
	size_t n = design_count(d);
	size_t p = design_trait_dim(d);
	const double *a = design_trait_matrix(d);
	size_t lda = p;

	if (p) {
		blas_dgemv(BLAS_NOTRANS, p, n, alpha, a, lda, x, 1, beta, y, 1);
	}
}


void design_traits_axpy(double alpha, const struct design *d, size_t i, double *y)
{
	assert(i < design_count(d));

	const double *x = design_traits(d, i);
	size_t dim = design_trait_dim(d);

	blas_daxpy(dim, alpha, x, 1, y, 1);
}


void design_tvars_mul(double alpha, const struct design *d,
		      const double *x, double beta, double *y)
{
	size_t n = design_count(d);
	size_t dim = design_tvar_dim(d);
	const double *a;
	const size_t *k;
	size_t nz;

	if (beta == 0.0) {
		memset(y, 0, n * sizeof(*y));
	} else if (beta != 1.0) {
		blas_dscal(n, beta, y, 1);
	}

	design_get_tvar_matrix(d, &a, &k, &nz);
	for (; nz != 0; a += dim, k++, nz--) {
		y[*k] += alpha * blas_ddot(dim, a, 1, x, 1);
	}
}


void design_tvars_tmul(double alpha, const struct design *d, const double *x, double beta, double *y)
{
	size_t dim = design_tvar_dim(d);
	const double *a;
	const size_t *k;
	size_t nz;

	if (beta == 0.0) {
		memset(y, 0, dim * sizeof(*y));
	} else if (beta != 1.0) {
		blas_dscal(dim, beta, y, 1);
	}

	design_get_tvar_matrix(d, &a, &k, &nz);
	for (; nz != 0; a += dim, k++, nz--) {
		blas_daxpy(dim, alpha * x[*k], a, 1, y, 1);
	}
}


void design_tvars_axpy(double alpha, const struct design *d, size_t i, double *y)
{
	assert(i < design_count(d));

	const double *dx = design_tvars(d, i);

	if (dx) {
		size_t dim = design_tvar_dim(d);
		blas_daxpy(dim, alpha, dx, 1, y, 1);
	}
}


void design_mul(double alpha, const struct design *d,
		const struct coefs *c, double beta, double *y)
{
	design_traits_mul(alpha, d, c->traits, beta, y);
	design_tvars_mul(alpha, d, c->tvars, 1.0, y);
}


void design_tmul(double alpha, const struct design *d, const double *x,
		 double beta, struct coefs *c)
{
	design_traits_tmul(alpha, d, x, beta, c->traits);
	design_tvars_tmul(alpha, d, x, beta, c->tvars);
}


void design_axpy(double alpha, const struct design *d, size_t i, struct coefs *c)
{
	design_traits_axpy(alpha, d, i, c->traits);
	design_tvars_axpy(alpha, d, i, c->tvars);
}
