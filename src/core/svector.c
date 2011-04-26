#include "port.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include "blas-private.h"
#include "compare.h"
#include "svector.h"

DEFINE_COMPARE_AND_EQUALS_FN(int64_compare, int64_equals, int64_t)

int64_t iproc_svector_find_nz(const struct svector * svector, int64_t i)
{
	assert(svector);
	assert(0 <= i);
	assert(i < svector_dim(svector));
	ssize_t ix = darray_binary_search(&svector->index,
					  &i,
					  int64_compare);
	return ix;
}

struct svector *svector_new(ssize_t n)
{
	assert(n >= 0);

	struct svector *v;

	if ((v = malloc(sizeof(*v)))) {
		if (svector_init(v, n)) {
			return v;
		}
		free(v);
	}

	return NULL;
}

bool svector_init(struct svector *v, ssize_t n)
{
	assert(v);
	assert(n >= 0);

	if (darray_init(&v->index, sizeof(int64_t))) {
		if (darray_init(&v->value, sizeof(double))) {
			v->dim = n;
			return true;
		}
		darray_deinit(&v->index);
	}
	return false;
}

struct svector *svector_new_copy(const struct svector *src)
{
	assert(src);
	struct svector *v;

	if ((v = malloc(sizeof(*v)))) {
		if (svector_init_copy(v, src)) {
			return v;
		}
		free(v);
	}
	return NULL;
}

bool svector_init_copy(struct svector *v, const struct svector *src)
{
	assert(v);
	assert(src);

	if (svector_init(v, svector_dim(src))) {
		if (svector_assign_copy(v, src)) {
			return true;
		}
		svector_deinit(v);
	}
	return false;
}

void svector_free(struct svector *v)
{
	if (v) {
		svector_deinit(v);
		free(v);
	}
}

void svector_deinit(struct svector *v)
{
	darray_deinit(&v->value);
	darray_deinit(&v->index);
}

void svector_clear(struct svector *v)
{
	assert(v);
	darray_clear(&v->index);
	darray_clear(&v->value);
}

ssize_t svector_dim(const struct svector* v)
{
	assert(v);
	return v->dim;
}

double svector_get(const struct svector *v, ssize_t i)
{
	assert(v);
	assert(0 <= i);
	assert(i < svector_dim(v));

	double value = 0.0;
	int64_t ix = iproc_svector_find_nz(v, i);

	if (ix >= 0) {
		value = iproc_svector_nz_get(v, ix);
	}

	return value;
}

bool svector_set(struct svector* v, ssize_t i, double val)
{
	assert(v);
	assert(0 <= i);
	assert(i < svector_dim(v));

	ssize_t ix = iproc_svector_find_nz(v, i);

	if (ix < 0) {
		ix = ~ix;
		if (darray_insert(&v->index, ix, &i)) {
			if (darray_insert(&v->value, ix, &val)) {
				return true;
			}
			darray_erase(&v->index, ix);
		}
		return false;
		
	} else {
		*(double *)darray_at(&v->value, ix) = val;
		return true;
	}
}

double *svector_at(struct svector *v, ssize_t i)
{
	assert(v);
	assert(0 <= i && i < svector_dim(v));
	
	ssize_t ix = iproc_svector_find_nz(v, i);
	double zero = 0.0;
	
	if (ix < 0) {
		ix = ~ix;
		if (darray_insert(&v->index, ix, &i)) {
			if (darray_insert(&v->value, ix, &zero)) {
				return darray_at(&v->value, ix);
			}
			darray_erase(&v->index, ix);
		}
		return NULL;
	} else {
		return darray_at(&v->value, ix);
	}
}

void svector_scale(struct svector *v, double scale)
{
	assert(v);

	if (svector_size(v) == 0)
		return;

	f77int n = (f77int) svector_size(v);
	double alpha = scale;
	f77int incx = 1;
	double *px = darray_front(&v->value);

	F77_FUNC(dscal) (&n, &alpha, px, &incx);
}

ssize_t svector_size(const struct svector *v)
{
	assert(v);
	return darray_size(&v->index);
}

int64_t iproc_svector_nz(const struct svector * svector, int64_t i)
{
	assert(svector);
	assert(0 <= i);
	assert(i < svector_size(svector));
	return *(int64_t *)darray_at(&svector->index, i);
}

double iproc_svector_nz_get(const struct svector * svector, int64_t i)
{
	assert(svector);
	assert(0 <= i);
	assert(i < svector_size(svector));
	return *(double *)darray_at(&svector->value, i);
}

void iproc_svector_nz_set(struct svector * svector, int64_t i, double value)
{
	assert(svector);
	assert(0 <= i);
	assert(i < svector_size(svector));
	*(double *)darray_at(&svector->value, i) = value;
}

void iproc_svector_nz_inc(struct svector * svector, int64_t i, double inc)
{
	assert(svector);
	assert(0 <= i);
	assert(i < svector_size(svector));
	*(double *)darray_at(&svector->value, i) += inc;
}

iproc_vector_view iproc_svector_view_nz(const struct svector * svector)
{
	assert(svector);

	int64_t nnz = svector_size(svector);
	double *px = NULL;

	if (nnz > 0) {
		px = darray_front(&svector->value);
	}

	return iproc_vector_view_array(px, nnz);
}

double svector_dot(const struct svector *x, const struct vector *y)
{
	assert(x);
	assert(y);
	assert(vector_dim(y) == svector_dim(x));

	int64_t i, ix, n = svector_size(x);
	double e1, e2, dot = 0.0;

	for (i = 0; i < n; i++) {
		ix = iproc_svector_nz(x, i);
		e1 = *vector_at(y, ix);
		e2 = iproc_svector_nz_get(x, i);
		dot += e1 * e2;
	}
	return dot;
}

void svector_axpy(double scale, const struct svector *x, struct vector *y)
{
	assert(y);
	assert(x);
	assert(vector_dim(y) == svector_dim(x));

	int64_t i, ix, n = svector_size(x);
	double e;

	for (i = 0; i < n; i++) {
		e = iproc_svector_nz_get(x, i);
		ix = iproc_svector_nz(x, i);
		*vector_at(y, ix) += scale * e;
	}
}

double svector_dots(const struct svector *v1, const struct svector *v2)
{
	assert(v1);
	assert(v2);
	assert(svector_dim(v1) == svector_dim(v2));

	int64_t n1 = svector_size(v1);
	int64_t n2 = svector_size(v2);

	if (n1 == 0 || n2 == 0)
		return 0.0;

	int64_t i1 = 0;
	int64_t i2 = 0;
	int64_t ix1, ix2;
	double dot = 0.0, val1, val2;

	while (i1 < n1 && i2 < n2) {
		ix1 = iproc_svector_nz(v1, i1);
		ix2 = iproc_svector_nz(v2, i2);

		if (ix1 == ix2) {
			val1 = iproc_svector_nz_get(v1, i1);
			val2 = iproc_svector_nz_get(v2, i2);
			dot += val1 * val2;
			i1++;
			i2++;
		} else if (ix1 < ix2) {
			i1++;
		} else {	/* ix1 > ix2 */
			i2++;
		}
	}

	return dot;
}

void svector_axpys(double scale, const struct svector* x, struct svector *y)
{
	assert(y);
	assert(x);
	assert(svector_dim(y) == svector_dim(x));

	int64_t inz, nnz = svector_size(x);
	for (inz = 0; inz < nnz; inz++) {
		int64_t i = iproc_svector_nz(x, inz);
		double val = iproc_svector_nz_get(x, inz);
		*svector_at(y, i) += scale * val;
	}
}

void svector_printf(const struct svector * svector)
{
	printf("\nsvector {");
	printf("\n  dim: %" SSIZE_FMT "", svector_dim(svector));
	printf("\n   nz: {");

	int64_t i, n = svector_size(svector);
	for (i = 0; i < n; i++) {
		printf("\n         %" PRId64 ", %.8f",
		       iproc_svector_nz(svector, i),
		       iproc_svector_nz_get(svector, i));
	}
	printf("\n       }");
	printf("\n}\n");
}

bool svector_assign_copy(struct svector *dst, const struct svector *src)
{
	assert(dst);
	assert(src);
	assert(svector_dim(dst) == svector_dim(src));

	int64_t nnz = svector_size(src);

	darray_resize(&dst->index, nnz);
	if (nnz > 0) {
		int64_t *idst = darray_front(&dst->index);
		int64_t *isrc = darray_front(&src->index);
		memcpy(idst, isrc, nnz * sizeof(int64_t));
	}

	darray_resize(&dst->value, nnz);
	if (nnz > 0) {
		double *vdst = darray_front(&dst->value);
		double *vsrc = darray_front(&src->value);
		memcpy(vdst, vsrc, nnz * sizeof(double));
	}

	return true;
}
