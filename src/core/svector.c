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
	assert(n <= SSIZE_MAX);

	v->data = NULL;
	v->index = NULL;
	v->dim = n;
	v->nnz = 0;
	v->nnzmax = 0;
	v->is_view = false;
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

	svector_init(v, svector_dim(src));
	svector_assign_copy(v, src);
}

void svector_assign_copy(struct svector *dst, const struct svector *src)
{
	assert(dst);
	assert(src);
	assert(svector_dim(dst) == svector_dim(src));

	ssize_t nnz = src->nnz;
	if (dst->nnz < nnz) {
		dst->data = xrealloc(dst->data, nnz * sizeof(dst->data[0]));
		dst->index = xrealloc(dst->index, nnz * sizeof(dst->index[0]));
		dst->nnzmax = nnz;
	}
	memcpy(dst->data, src->data, nnz * sizeof(dst->data[0]));
	memcpy(dst->index, src->index, nnz * sizeof(dst->index[0]));
	dst->nnz = nnz;
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
	if (svector_owner(v)) {
		xfree(v->index);
		xfree(v->data);
	}
}

void svector_clear(struct svector *v)
{
	assert(v);
	assert(svector_owner(v));
	v->nnz = 0;
}

void svector_set_basis(struct svector *v, ssize_t i)
{
	assert(v);
	assert(0 <= i && i < svector_dim(v));
	assert(svector_owner(v));

	if (!v->nnzmax) {
		v->data = xmalloc(sizeof(v->data[0]));
		v->index = xmalloc(sizeof(v->index[0]));
		v->nnzmax = 1;
	}
	v->data[0] = 1.0;
	v->index[0] = i;
	v->nnz = 1;
}

double svector_item(const struct svector *v, ssize_t i)
{
	assert(v);
	assert(0 <= i);
	assert(i < svector_dim(v));

	struct svector_pos pos;
	double *ptr;
	
	if ((ptr = svector_find(v, i, &pos))) {
		return *ptr;
	} else {
		return 0.0;
	}
}

double *svector_item_ptr(struct svector *v, ssize_t i)
{
	assert(v);
	assert(0 <= i && i < svector_dim(v));

	struct svector_pos pos;
	double *ptr;
	
	if (!(ptr = svector_find(v, i, &pos))) {
		ptr = svector_insert(v, &pos, 0.0);
	}
	return ptr;
}

void svector_set_item(struct svector *v, ssize_t i, double val)
{
	assert(v);
	assert(0 <= i);
	assert(i < svector_dim(v));

	struct svector_pos pos;
	double *ptr;
	
	if ((ptr = svector_find(v, i, &pos))) {
		*ptr = val;
	} else {
		svector_insert(v, &pos, val);
	}
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
		valy = vector_item(y, i);
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
	
	const ssize_t *base = v->index;
	ssize_t len = v->nnz;
	const ssize_t *begin = base, *end = base + len, *ptr;
	
	while (begin < end) {
		ptr = begin + ((end - begin) >> 1);
		if (*ptr < i) {
			begin = ptr + 1;			
		} else if (*ptr > i) {
			end = ptr;
		} else {
			pos->inz = ptr - v->index;
			pos->i = i;
			return &v->data[pos->inz];
		}
	}
	assert(begin == end);
	
	pos->inz = begin - v->index;
	pos->i = i;
	return NULL;
}

double *svector_insert(struct svector *v, struct svector_pos *pos, double val)
{
	assert(v);
	assert(svector_owner(v));
	assert(pos);
	assert(pos->inz == v->nnz || v->index[pos->inz] > pos->i);

	if (v->nnz == v->nnzmax) {
		ssize_t nnzmax1 = (v->nnzmax + 1) * 2;
		v->data = xrealloc(v->data, nnzmax1 * sizeof(v->data[0]));
		v->index = xrealloc(v->index, nnzmax1 * sizeof(v->index[0]));
		v->nnzmax = nnzmax1;
	}
	assert(v->nnz < v->nnzmax);
	
	double *data = v->data;
	ssize_t *index = v->index;
	ssize_t inz = pos->inz;
	ssize_t nnz = v->nnz;
	ssize_t ntail = nnz - inz;

	// Switch is just an optimization for memmove in default case
	switch (ntail) {
	case 2:
		data[inz + 2] = data[inz + 1];
		index[inz + 2] = index[inz + 1];
	case 1:
		data[inz + 1] = data[inz];
		index[inz + 1] = index[inz];
		break;
	default:
		memmove(data + inz + 1, data + inz, ntail * sizeof(data[0]));
		memmove(index + inz + 1, index + inz, ntail * sizeof(index[0]));
	}
	
	data[inz] = val;		
	index[inz] = pos->i;
	v->nnz = nnz + 1;	

	return &data[inz];
}

void svector_remove_at(struct svector *v, struct svector_pos *pos)
{
	assert(v);
	assert(svector_owner(v));
	assert(pos);
	assert(0 <= pos->inz && pos->inz < v->nnz);
	
	double *data = v->data;
	ssize_t *index = v->index;
	ssize_t inz = pos->inz;
	ssize_t nnz = v->nnz;
	ssize_t ntail = nnz - inz - 1;
	
	memmove(data + inz, data + inz + 1, ntail * sizeof(data[0]));
	memmove(index + inz, index + inz + 1, ntail * sizeof(index[0]));
	v->nnz = nnz - 1;
}

