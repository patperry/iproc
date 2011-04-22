#include "port.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include <inttypes.h>
#include "blas-private.h"
#include "compare.h"
#include "svector.h"

DEFINE_COMPARE_AND_EQUALS_FN(int64_compare, int64_equals, int64_t)


int64_t iproc_svector_find_nz(const iproc_svector * svector, int64_t i)
{
	assert(svector);
	assert(0 <= i);
	assert(i < iproc_svector_dim(svector));
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

	if (svector_init(v, iproc_svector_dim(src))) {
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




void svector_clear(struct svector * v)
{
	assert(v);
	darray_clear(&v->index);
	darray_clear(&v->value);
}

int64_t iproc_svector_dim(const iproc_svector * svector)
{
	assert(svector);
	return svector->dim;
}

double iproc_svector_get(const iproc_svector * svector, int64_t i)
{
	assert(svector);
	assert(0 <= i);
	assert(i < iproc_svector_dim(svector));

	double value = 0.0;
	int64_t ix = iproc_svector_find_nz(svector, i);

	if (ix >= 0) {
		value = iproc_svector_nz_get(svector, ix);
	}

	return value;
}

void iproc_svector_set(iproc_svector * svector, int64_t i, double value)
{
	assert(svector);
	assert(0 <= i);
	assert(i < iproc_svector_dim(svector));

	int64_t ix = iproc_svector_find_nz(svector, i);

	if (ix < 0) {
		ix = ~ix;
		darray_insert(&svector->index, ix, &i);
		darray_insert(&svector->value, ix, &value);
	} else {
		*(double *)darray_at(&svector->value, ix) = value;
	}
}

void iproc_svector_inc(iproc_svector * svector, int64_t i, double value)
{
	assert(svector);
	assert(0 <= i);
	assert(i < iproc_svector_dim(svector));

	int64_t ix = iproc_svector_find_nz(svector, i);

	if (ix < 0) {
		ix = ~ix;
		darray_insert(&svector->index, ix, &i);
		darray_insert(&svector->value, ix, &value);
	} else {
		*(double *)darray_at(&svector->value, ix) += value;
	}
}

void iproc_svector_scale(iproc_svector * svector, double scale)
{
	assert(svector);

	if (iproc_svector_nnz(svector) == 0)
		return;

	f77int n = (f77int) iproc_svector_nnz(svector);
	double alpha = scale;
	f77int incx = 1;
	double *px = darray_front(&svector->value);

	F77_FUNC(dscal) (&n, &alpha, px, &incx);
}

int64_t iproc_svector_nnz(const iproc_svector * svector)
{
	assert(svector);
	return darray_size(&svector->index);
}

int64_t iproc_svector_nz(const iproc_svector * svector, int64_t i)
{
	assert(svector);
	assert(0 <= i);
	assert(i < iproc_svector_nnz(svector));
	return *(int64_t *)darray_at(&svector->index, i);
}

double iproc_svector_nz_get(const iproc_svector * svector, int64_t i)
{
	assert(svector);
	assert(0 <= i);
	assert(i < iproc_svector_nnz(svector));
	return *(double *)darray_at(&svector->value, i);
}

void iproc_svector_nz_set(iproc_svector * svector, int64_t i, double value)
{
	assert(svector);
	assert(0 <= i);
	assert(i < iproc_svector_nnz(svector));
	*(double *)darray_at(&svector->value, i) = value;
}

void iproc_svector_nz_inc(iproc_svector * svector, int64_t i, double inc)
{
	assert(svector);
	assert(0 <= i);
	assert(i < iproc_svector_nnz(svector));
	*(double *)darray_at(&svector->value, i) += inc;
}

iproc_vector_view iproc_svector_view_nz(const iproc_svector * svector)
{
	assert(svector);

	int64_t nnz = iproc_svector_nnz(svector);
	double *px = NULL;

	if (nnz > 0) {
		px = darray_front(&svector->value);
	}

	return iproc_vector_view_array(px, nnz);
}

double vector_dots(const struct vector *vector,
			 const iproc_svector * svector)
{
	assert(vector);
	assert(svector);
	assert(vector_size(vector) == iproc_svector_dim(svector));

	int64_t i, ix, n = iproc_svector_nnz(svector);
	double e1, e2, dot = 0.0;

	for (i = 0; i < n; i++) {
		ix = iproc_svector_nz(svector, i);
		e1 = *vector_at(vector, ix);
		e2 = iproc_svector_nz_get(svector, i);
		dot += e1 * e2;
	}
	return dot;
}

void
iproc_vector_sacc(struct vector *dst_vector,
		  double scale, const iproc_svector * svector)
{
	assert(dst_vector);
	assert(svector);
	assert(vector_size(dst_vector) == iproc_svector_dim(svector));

	int64_t i, ix, n = iproc_svector_nnz(svector);
	double e;

	for (i = 0; i < n; i++) {
		e = iproc_svector_nz_get(svector, i);
		ix = iproc_svector_nz(svector, i);
		*vector_at(dst_vector, ix) += scale * e;
	}
}

double iproc_svector_sdot(const iproc_svector * svector1,
			  const iproc_svector * svector2)
{
	assert(svector1);
	assert(svector2);
	assert(iproc_svector_dim(svector1) == iproc_svector_dim(svector2));

	int64_t n1 = iproc_svector_nnz(svector1);
	int64_t n2 = iproc_svector_nnz(svector2);

	if (n1 == 0 || n2 == 0)
		return 0.0;

	int64_t i1 = 0;
	int64_t i2 = 0;
	int64_t ix1, ix2;
	double dot = 0.0, val1, val2;

	while (i1 < n1 && i2 < n2) {
		ix1 = iproc_svector_nz(svector1, i1);
		ix2 = iproc_svector_nz(svector2, i2);

		if (ix1 == ix2) {
			val1 = iproc_svector_nz_get(svector1, i1);
			val2 = iproc_svector_nz_get(svector2, i2);
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

void
iproc_svector_sacc(iproc_svector * dst_svector,
		   double scale, const iproc_svector * svector)
{
	assert(dst_svector);
	assert(svector);
	assert(iproc_svector_dim(dst_svector) == iproc_svector_dim(svector));

	int64_t inz, nnz = iproc_svector_nnz(svector);
	for (inz = 0; inz < nnz; inz++) {
		int64_t i = iproc_svector_nz(svector, inz);
		double val = iproc_svector_nz_get(svector, inz);
		iproc_svector_inc(dst_svector, i, scale * val);
	}
}

void iproc_svector_printf(const iproc_svector * svector)
{
	printf("\nsvector {");
	printf("\n  dim: %" PRId64 "", iproc_svector_dim(svector));
	printf("\n   nz: {");

	int64_t i, n = iproc_svector_nnz(svector);
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
	assert(iproc_svector_dim(dst) ==
	       iproc_svector_dim(src));

	int64_t nnz = iproc_svector_nnz(src);

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
