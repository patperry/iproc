#include "port.h"

#include <assert.h>
#include <stdio.h>
#include <stdlib.h>
#include "blas-private.h"
#include "compare.h"
#include "svector.h"

DEFINE_COMPARE_FN(ssize_compare, ssize_t)

static ssize_t svector_find_index(const struct svector *v, ssize_t i)
{
	assert(v);
	assert(0 <= i);
	assert(i < svector_dim(v));
	ssize_t ix = darray_binary_search(&v->index,
					  &i,
					  ssize_compare);
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

	if (darray_init(&v->index, sizeof(ssize_t))) {
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
	ssize_t ix = svector_find_index(v, i);

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

	ssize_t ix = svector_find_index(v, i);

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
	
	ssize_t ix = svector_find_index(v, i);
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

double svector_dot(const struct svector *x, const struct vector *y)
{
	assert(x);
	assert(y);
	assert(vector_dim(y) == svector_dim(x));

	ssize_t i;
	double dot = 0.0, valx, valy;
	struct svector_iter itx;
	
	svector_iter_init(x, &itx);
	while (svector_iter_advance(x, &itx)) {
		i = svector_iter_current_index(x, &itx);
		valx = *svector_iter_current(x, &itx);
		valy = *vector_at(y, i);
		dot += valx * valy;
	}
	svector_iter_deinit(x, &itx);
	
	return dot;
}

void svector_axpy(double scale, const struct svector *x, struct vector *y)
{
	assert(y);
	assert(x);
	assert(vector_dim(y) == svector_dim(x));

	struct svector_iter itx;
	ssize_t i;
	double val;
	
	svector_iter_init(x, &itx);
	while (svector_iter_advance(x, &itx)) {
		i = svector_iter_current_index(x, &itx);
		val = *svector_iter_current(x, &itx);
		*vector_at(y, i) += scale * val;
	}
	svector_iter_deinit(x, &itx);
}

double svector_dots(const struct svector *v1, const struct svector *v2)
{
	assert(v1);
	assert(v2);
	assert(svector_dim(v1) == svector_dim(v2));

	ssize_t n1 = svector_size(v1);
	ssize_t n2 = svector_size(v2);

	if (n1 == 0 || n2 == 0)
		return 0.0;

	const struct svector *x, *y;
	if (n1 <= n2) {
		x = v1;
		y = v2;
	} else {
		x = v2;
		y = v1;
	}
	
	ssize_t i;
	double dot = 0.0, valx, valy;
	struct svector_iter itx;

	svector_iter_init(x, &itx);
	while (svector_iter_advance(x, &itx)) {
		i = svector_iter_current_index(x, &itx);
		valx = *svector_iter_current(x, &itx);
		valy = svector_get(y, i);
		dot += valx * valy;
	}
	svector_iter_deinit(x, &itx);

	return dot;
}

void svector_axpys(double scale, const struct svector* x, struct svector *y)
{
	assert(y);
	assert(x);
	assert(svector_dim(y) == svector_dim(x));

	struct svector_iter itx;
	ssize_t i;
	double val;

	svector_iter_init(x, &itx);
	while (svector_iter_advance(x, &itx)) {
		i = svector_iter_current_index(x, &itx);
		val = *svector_iter_current(x, &itx);
		*svector_at(y, i) += scale * val;
	}
	svector_iter_deinit(x, &itx);
}

void svector_printf(const struct svector *v)
{
	struct svector_iter it;
	ssize_t i;
	double val;
	
	printf("\nsvector {");
	printf("\n  dim: %" SSIZE_FMT "", svector_dim(v));
	printf("\n   nz: {");

	svector_iter_init(v, &it);
	while (svector_iter_advance(v, &it)) {
		i = svector_iter_current_index(v, &it);
		val = *svector_iter_current(v, &it);
		printf("\n         %" SSIZE_FMT ", %.8f", i, val);

	}
	svector_iter_deinit(v, &it);

	printf("\n       }");
	printf("\n}\n");
}

bool svector_assign_copy(struct svector *dst, const struct svector *src)
{
	assert(dst);
	assert(src);
	assert(svector_dim(dst) == svector_dim(src));

	ssize_t nnz = svector_size(src);

	darray_resize(&dst->index, nnz);
	if (nnz > 0) {
		ssize_t *idst = darray_front(&dst->index);
		ssize_t *isrc = darray_front(&src->index);
		memcpy(idst, isrc, nnz * sizeof(ssize_t));
	}

	darray_resize(&dst->value, nnz);
	if (nnz > 0) {
		double *vdst = darray_front(&dst->value);
		double *vsrc = darray_front(&src->value);
		memcpy(vdst, vsrc, nnz * sizeof(double));
	}

	return true;
}

void svector_iter_init(const struct svector *v, struct svector_iter *it)
{
	it->pat_index = -1;
}

bool svector_iter_advance(const struct svector *v, struct svector_iter *it)
{
	it->pat_index++;
	return it->pat_index < svector_size(v);
}

double *svector_iter_current(const struct svector *v, struct svector_iter *it)
{
	return darray_at(&v->value, it->pat_index);
}

ssize_t svector_iter_current_index(const struct svector *v, struct svector_iter *it)
{
	return *(ssize_t *)darray_at(&v->index, it->pat_index);
}

void svector_iter_deinit(const struct svector *v, struct svector_iter *it)
{
}


/* DEPRECATED */
ssize_t iproc_svector_nz(const struct svector * svector, ssize_t i)
{
	assert(svector);
	assert(0 <= i);
	assert(i < svector_size(svector));
	return *(ssize_t *)darray_at(&svector->index, i);
}

double iproc_svector_nz_get(const struct svector * svector, ssize_t i)
{
	assert(svector);
	assert(0 <= i);
	assert(i < svector_size(svector));
	return *(double *)darray_at(&svector->value, i);
}

void iproc_svector_nz_set(struct svector * svector, ssize_t i, double value)
{
	assert(svector);
	assert(0 <= i);
	assert(i < svector_size(svector));
	*(double *)darray_at(&svector->value, i) = value;
}

void iproc_svector_nz_inc(struct svector * svector, ssize_t i, double inc)
{
	assert(svector);
	assert(0 <= i);
	assert(i < svector_size(svector));
	*(double *)darray_at(&svector->value, i) += inc;
}

iproc_vector_view iproc_svector_view_nz(const struct svector * svector)
{
	assert(svector);
	
	ssize_t nnz = svector_size(svector);
	double *px = NULL;
	
	if (nnz > 0) {
		px = darray_front(&svector->value);
	}
	
	return iproc_vector_view_array(px, nnz);
}
