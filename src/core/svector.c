#include "port.h"

#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "blas-private.h"
#include "hash.h"
#include "compare.h"
#include "svector.h"

struct svector *svector_alloc(ssize_t n)
{
	assert(n >= 0);

	struct svector *v = xcalloc(1, sizeof(*v));
	svector_init(v, n);
	return v;
}

void svector_init(struct svector *v, ssize_t n)
{
	assert(v);
	assert(n >= 0);
	assert(n <= INTPTR_MAX);

	intmap_init(&v->map, sizeof(double), alignof(double));
	v->dim = n;
}

struct svector *svector_alloc_copy(const struct svector *src)
{
	assert(src);
	struct svector *v = xcalloc(1, sizeof(*v));

	svector_init_copy(v, src);
	return v;
}

void svector_init_copy(struct svector *v, const struct svector *src)
{
	assert(v);
	assert(src);

	intmap_init_copy(&v->map, &src->map);
	v->dim = src->dim;
}

void svector_assign_copy(struct svector *dst, const struct svector *src)
{
	assert(dst);
	assert(src);
	assert(svector_dim(dst) == svector_dim(src));

	intmap_assign_copy(&dst->map, &src->map);
}

void svector_free(struct svector *v)
{
	if (v) {
		svector_deinit(v);
		xfree(v);
	}
}

void svector_deinit(struct svector *v)
{
	intmap_deinit(&v->map);
}

void svector_clear(struct svector *v)
{
	assert(v);
	intmap_clear(&v->map);
}

void svector_set_basis(struct svector *v, ssize_t i)
{
	assert(v);
	assert(0 <= i && i < svector_dim(v));

	svector_clear(v);
	svector_set_item(v, i, 1.0);
}

double svector_item(const struct svector *v, ssize_t i)
{
	assert(v);
	assert(0 <= i);
	assert(i < svector_dim(v));

	double *ptr = intmap_item(&v->map, i);
	return ptr ? *ptr : 0.0;
}

double *svector_item_ptr(struct svector *v, ssize_t i)
{
	assert(v);
	assert(0 <= i && i < svector_dim(v));

	struct intmap_pos pos;
	double zero = 0.0;
	double *val;

	if ((val = intmap_find(&v->map, i, &pos))) {
		return val;
	} else {
		return intmap_insert(&v->map, &pos, &zero);
	}
}

double *svector_set_item(struct svector *v, ssize_t i, double val)
{
	assert(v);
	assert(0 <= i);
	assert(i < svector_dim(v));

	return intmap_set_item(&v->map, i, &val);
}

void svector_remove(struct svector *v, ssize_t i)
{
	assert(v);
	assert(0 <= i && i < svector_dim(v));

	struct svector_pos pos;

	if ((svector_find(v, i, &pos))) {
		svector_remove_at(v, &pos);
	}
}

void svector_scale(struct svector *v, double scale)
{
	assert(v);

	struct svector_iter it;
	double *ptr;

	SVECTOR_FOREACH(it, v) {
		ptr = SVECTOR_PTR(it);
		*ptr *= scale;
	}
}

double svector_max(const struct svector *v)
{
	assert(v);

	struct svector_iter it;
	double val;
	double max = NAN;

	SVECTOR_FOREACH(it, v) {
		val = SVECTOR_VAL(it);
		if (isnan(max)) {
			max = val;
		} else if (val > max) {
			max = val;
		}
	}

	return max;
}

double svector_dot(const struct svector *x, const struct vector *y)
{
	assert(x);
	assert(y);
	assert(vector_dim(y) == svector_dim(x));

	ssize_t i;
	double dot = 0.0, valx, valy;
	struct svector_iter itx;

	SVECTOR_FOREACH(itx, x) {
		i = SVECTOR_IDX(itx);
		valx = SVECTOR_VAL(itx);
		valy = *vector_item_ptr(y, i);
		dot += valx * valy;
	}

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

	SVECTOR_FOREACH(itx, x) {
		i = SVECTOR_IDX(itx);
		val = SVECTOR_VAL(itx);
		*vector_item_ptr(y, i) += scale * val;
	}
}

double svector_dots(const struct svector *v1, const struct svector *v2)
{
	assert(v1);
	assert(v2);
	assert(svector_dim(v1) == svector_dim(v2));

	ssize_t n1 = svector_count(v1);
	ssize_t n2 = svector_count(v2);

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

	SVECTOR_FOREACH(itx, x) {
		i = SVECTOR_IDX(itx);
		valx = SVECTOR_VAL(itx);
		valy = svector_item(y, i);
		dot += valx * valy;
	}

	return dot;
}

void svector_axpys(double scale, const struct svector *x, struct svector *y)
{
	assert(y);
	assert(x);
	assert(svector_dim(y) == svector_dim(x));

	struct svector_iter itx;
	ssize_t i;
	double val;
	double *yi;

	SVECTOR_FOREACH(itx, x) {
		i = SVECTOR_IDX(itx);
		val = SVECTOR_VAL(itx);
		yi = svector_item_ptr(y, i);
		*yi += scale * val;
	}
}

void svector_printf(const struct svector *v)
{
	struct svector_iter it;
	ssize_t i;
	double val;

	printf("\nsvector {");
	printf("\n  dim: %" SSIZE_FMT "", svector_dim(v));
	printf("\n   nz: {");

	SVECTOR_FOREACH(it, v) {
		i = SVECTOR_IDX(it);
		val = SVECTOR_VAL(it);
		printf("\n         %" SSIZE_FMT ", %.8f", i, val);
	}

	printf("\n       }");
	printf("\n}\n");
}

double *svector_find(const struct svector *v, ssize_t i,
		     struct svector_pos *pos)
{
	assert(v);
	assert(0 <= i && i < svector_dim(v));
	assert(pos);

	return intmap_find(&v->map, i, &pos->map_pos);
}

double *svector_insert(struct svector *v, struct svector_pos *pos, double val)
{
	assert(v);
	assert(pos);

	return intmap_insert(&v->map, &pos->map_pos, &val);
}

void svector_remove_at(struct svector *v, struct svector_pos *pos)
{
	assert(v);
	assert(pos);

	intmap_remove_at(&v->map, &pos->map_pos);
}

struct svector_iter svector_iter_make(const struct svector *v)
{
	assert(v);
	struct svector_iter it;
	it.map_it = intmap_iter_make(&v->map);
	return it;
}

bool svector_iter_advance(struct svector_iter *it)
{
	assert(it);
	return intmap_iter_advance(&it->map_it);
}

void svector_iter_reset(struct svector_iter *it)
{
	assert(it);
	intmap_iter_reset(&it->map_it);
}
