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
	assert(n <= INTPTR_MAX);

	if (intmap_init(&v->map, sizeof(double), alignof(double))) {
		v->dim = n;
		return true;
	}
	return false;
}

struct svector *svector_alloc_copy(const struct svector *src)
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

	if (intmap_init_copy(&v->map, &src->map)) {
		v->dim = src->dim;
		return true;
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
	intmap_deinit(&v->map);
}

void svector_clear(struct svector *v)
{
	assert(v);
	intmap_clear(&v->map);
}

ssize_t svector_dim(const struct svector *v)
{
	assert(v);
	return v->dim;
}

double svector_get(const struct svector *v, ssize_t i)
{
	assert(v);
	assert(0 <= i);
	assert(i < svector_dim(v));

	double val0 = 0.0;
	return *(double *)intmap_lookup_with(&v->map, i, &val0);
}

bool svector_set(struct svector *v, ssize_t i, double val)
{
	assert(v);
	assert(0 <= i);
	assert(i < svector_dim(v));

	return intmap_add(&v->map, i, &val);
}

double *svector_at(struct svector *v, ssize_t i)
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

void svector_scale(struct svector *v, double scale)
{
	assert(v);

	struct svector_iter it;
	svector_iter_init(v, &it);
	while (svector_iter_advance(v, &it)) {
		*svector_iter_current(v, &it) *= scale;
	}
	svector_iter_deinit(v, &it);
}

ssize_t svector_size(const struct svector *v)
{
	assert(v);
	return intmap_size(&v->map);
}

double svector_max(const struct svector *v)
{
	assert(v);

	struct svector_iter it;
	double val;
	double max = NAN;

	svector_iter_init(v, &it);
	while (svector_iter_advance(v, &it)) {
		val = *svector_iter_current(v, &it);
		if (isnan(max)) {
			max = val;
		} else if (val > max) {
			max = val;
		}
	}
	svector_iter_deinit(v, &it);

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

bool svector_axpys(double scale, const struct svector *x, struct svector *y)
{
	assert(y);
	assert(x);
	assert(svector_dim(y) == svector_dim(x));

	struct svector_iter itx;
	ssize_t i;
	double val;
	double *yi;
	bool ok = true;

	svector_iter_init(x, &itx);
	while (svector_iter_advance(x, &itx)) {
		i = svector_iter_current_index(x, &itx);
		val = *svector_iter_current(x, &itx);
		yi = svector_at(y, i);
		
		if (!yi) {
			ok = false;
			break;
		}
		
		*yi += scale * val;
	}
	svector_iter_deinit(x, &itx);
	return ok;
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

	return intmap_assign_copy(&dst->map, &src->map);
}

double *svector_find(const struct svector *v, ssize_t i, struct svector_pos *pos)
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

void svector_erase(struct svector *v, struct svector_pos *pos)
{
	assert(v);
	assert(pos);
	
	intmap_erase(&v->map, &pos->map_pos);
}

void svector_iter_init(const struct svector *v, struct svector_iter *it)
{
	intmap_iter_init(&v->map, &it->map_it);
}

bool svector_iter_advance(const struct svector *v, struct svector_iter *it)
{
	return intmap_iter_advance(&v->map, &it->map_it);
}

double *svector_iter_current(const struct svector *v, struct svector_iter *it)
{
	return intmap_iter_current(&v->map, &it->map_it);
}

ssize_t svector_iter_current_index(const struct svector *v,
				   struct svector_iter *it)
{
	return (ssize_t)intmap_iter_current_key(&v->map, &it->map_it);
}

void svector_iter_deinit(const struct svector *v, struct svector_iter *it)
{
	intmap_iter_deinit(&v->map, &it->map_it);

}
