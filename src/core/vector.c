#include "port.h"

#include <assert.h>
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "xalloc.h"

#include "blas-private.h"
#include "compare.h"
#include "hash.h"
#include "ieee754.h"
#include "vector.h"

void vector_init(struct vector *v, ssize_t n)
{
	assert(v);
	assert(n >= 0);
	assert(n <= F77INT_MAX);
	assert(n <= (ssize_t)(SSIZE_MAX / sizeof(double)));

	v->data = xcalloc(n, sizeof(v->data[0]));
	v->dim = n;
	v->is_view = false;
}

void vector_reinit(struct vector *v, ssize_t n)
{
	assert(v);
	assert(n >= 0);
	assert(n <= F77INT_MAX);
	assert(n <= (ssize_t)(SSIZE_MAX / sizeof(double)));
	assert(vector_owner(v));

	ssize_t nold = vector_dim(v);

	v->data = xrealloc(v->data, n * sizeof(v->data[0]));
	v->dim = n;

	if (nold < n) {
		struct vector vnew = vector_slice(v, nold, n - nold);
		vector_fill(&vnew, 0.0);
	}
}

void vector_init_copy(struct vector *v, const struct vector *src)
{
	assert(v);
	assert(src);

	vector_init(v, vector_dim(src));
	vector_assign_copy(v, src);
}

void vector_assign_copy(struct vector *v, const struct vector *src)
{
	assert(v);
	assert(src);
	assert(vector_dim(v) == vector_dim(src));

	if (!vector_dim(v))
		return;

	vector_copy_to(src, vector_to_ptr(v));
}

void vector_deinit(struct vector *v)
{
	assert(v);
	if (vector_owner(v))
		free(v->data);
}

struct vector *vector_alloc(ssize_t n)
{
	struct vector *v = xcalloc(1, sizeof(*v));
	vector_init(v, n);
	return v;
}

struct vector *vector_alloc_copy(const struct vector *src)
{
	struct vector *v = xcalloc(1, sizeof(*v));

	vector_init_copy(v, src);
	return v;
}

void vector_free(struct vector *v)
{
	if (v) {
		vector_deinit(v);
		free(v);
	}
}

void vector_copy_to(const struct vector *v, double *dst)
{
	assert(v);
	assert(dst || vector_dim(v) == 0);

	memcpy(dst, v->data, vector_dim(v) * sizeof(v->data[0]));
}

void vector_fill(struct vector *vector, double value)
{
	assert(vector);

	if (!vector_dim(vector))
		return;

	ssize_t n = vector_dim(vector);
	double *ptr = vector_item_ptr(vector, 0);
	double *end = ptr + n;

	while (ptr != end) {
		*ptr++ = value;
	}
}

void vector_set_basis(struct vector *v, ssize_t i)
{
	assert(v);
	assert(vector_dim(v) > 0);
	assert(0 <= i && i < vector_dim(v));

	vector_fill(v, 0.0);
	*vector_item_ptr(v, i) = 1.0;
}

void vector_scale(struct vector *vector, double scale)
{
	assert(vector);

	if (!vector_dim(vector))
		return;

	f77int n = (f77int)vector_dim(vector);
	double alpha = scale;
	void *px = vector_to_ptr(vector);
	f77int incx = 1;

	F77_FUNC(dscal) (&n, &alpha, px, &incx);
}

void vector_shift(struct vector *vector, double shift)
{
	assert(vector);

	if (!vector_dim(vector))
		return;

	ssize_t n = vector_dim(vector);
	double *ptr = vector_to_ptr(vector);
	double *end = ptr + n;

	while (ptr != end) {
		*ptr++ += shift;
	}
}

void vector_add(struct vector *dst_vector, const struct vector *vector)
{
	assert(dst_vector);
	assert(vector);
	assert(vector_dim(dst_vector) == vector_dim(vector));

	vector_axpy(1.0, vector, dst_vector);
}

void vector_sub(struct vector *dst_vector, const struct vector *vector)
{
	assert(dst_vector);
	assert(vector);
	assert(vector_dim(dst_vector) == vector_dim(vector));

	vector_axpy(-1.0, vector, dst_vector);
}

void vector_mul(struct vector *dst_vector, const struct vector *vector)
{
	assert(dst_vector);
	assert(vector);
	assert(vector_dim(dst_vector) == vector_dim(vector));

	if (!vector_dim(dst_vector))
		return;

	f77int n = (f77int)vector_dim(dst_vector);
	f77int k = 0;
	double *px = vector_to_ptr(vector);
	f77int incx = 1;
	double *py = vector_to_ptr(dst_vector);
	f77int incy = 1;

	F77_FUNC(dtbmv) ("U", "N", "N", &n, &k, px, &incx, py, &incy);
}

void vector_div(struct vector *dst_vector, const struct vector *vector)
{
	assert(dst_vector);
	assert(vector);
	assert(vector_dim(dst_vector) == vector_dim(vector));

	if (!vector_dim(dst_vector))
		return;

	f77int n = (f77int)vector_dim(dst_vector);
	f77int k = 0;
	double *px = vector_to_ptr(vector);
	f77int incx = 1;
	double *py = vector_to_ptr(dst_vector);
	f77int incy = 1;

	F77_FUNC(dtbsv) ("U", "N", "N", &n, &k, px, &incx, py, &incy);
}

void vector_axpy(double alpha, const struct vector *x, struct vector *y)
{
	assert(x);
	assert(y);
	assert(vector_dim(x) == vector_dim(y));

	if (!vector_dim(x))
		return;

	f77int n = (f77int)vector_dim(y);
	double *px = vector_to_ptr(x);
	f77int incx = 1;
	double *py = vector_to_ptr(y);
	f77int incy = 1;

	F77_FUNC(daxpy) (&n, &alpha, px, &incx, py, &incy);
}

double vector_dot(const struct vector *vector1, const struct vector *vector2)
{
	assert(vector1);
	assert(vector2);
	assert(vector_dim(vector1) == vector_dim(vector2));

	if (!vector_dim(vector1))
		return 0.0;

	f77int n = (f77int)vector_dim(vector1);
	double *px = vector_to_ptr(vector1);
	f77int incx = 1;
	double *py = vector_to_ptr(vector2);
	f77int incy = 1;

	double dot = F77_FUNC(ddot) (&n, px, &incx, py, &incy);
	return dot;
}

double vector_dist(const struct vector *v1, const struct vector *v2)
{
	assert(v1);
	assert(v2);
	assert(vector_dim(v1) == vector_dim(v2));

	ssize_t i, n = vector_dim(v1);
	double dist, dist2 = 0.0;
	double delta;

	for (i = 0; i < n; i++) {
		delta = vector_item(v1, i) - vector_item(v2, i);
		dist2 += delta * delta;
	}

	dist = sqrt(dist2);
	return dist;
}

double vector_norm(const struct vector *vector)
{
	assert(vector);

	if (!vector_dim(vector))
		return 0.0;

	f77int n = (f77int)vector_dim(vector);
	void *px = vector_to_ptr(vector);
	f77int incx = 1;

	double norm = F77_FUNC(dnrm2) (&n, px, &incx);
	return norm;
}

double vector_norm2(const struct vector *v)
{
	assert(v);
	return vector_dot(v, v);
}

double vector_sum_abs(const struct vector *vector)
{
	assert(vector);

	if (!vector_dim(vector))
		return 0.0;

	f77int n = (f77int)vector_dim(vector);
	void *px = vector_to_ptr(vector);
	f77int incx = 1;

	double sum_abs = F77_FUNC(dasum) (&n, px, &incx);
	return sum_abs;
}

double vector_max_abs(const struct vector *vector)
{
	assert(vector);

	if (!vector_dim(vector))
		return 0.0;

	ssize_t i = vector_max_abs_index(vector);
	double e = *vector_item_ptr(vector, i);
	double max_abs = fabs(e);

	return max_abs;
}

ssize_t vector_max_abs_index(const struct vector *vector)
{
	assert(vector);
	assert(vector_dim(vector));

	f77int n = (f77int)vector_dim(vector);
	void *px = vector_to_ptr(vector);
	f77int incx = 1;

	f77int index1 = F77_FUNC(idamax) (&n, px, &incx);

	ssize_t index = index1 - 1;
	return index;
}

ssize_t vector_max_index(const struct vector *vector)
{
	assert(vector);
	assert(vector_dim(vector));

	ssize_t n = vector_dim(vector);
	ssize_t i, imax;
	double x, max = NAN;

	/* Find the first non-NaN entry of the vector */
	imax = 0;
	while (imax < n && isnan(vector_item(vector, imax))) {
		imax++;
	}

	/* If all of the entries are NaN, define imax as 0. */
	if (imax == n)
		return 0;

	/* Otherwise, search for the largest entry in the tail of the vector */
	max = *vector_item_ptr(vector, imax);

	for (i = imax + 1; i < n; i++) {
		x = *vector_item_ptr(vector, i);
		if (x > max) {
			max = x;
			imax = i;
		}
	}

	return imax;
}

double vector_max(const struct vector *vector)
{
	assert(vector);

	if (!vector_dim(vector)) {
		return -INFINITY;
	} else {
		ssize_t i = vector_max_index(vector);
		return *vector_item_ptr(vector, i);
	}
}

void vector_exp(struct vector *vector)
{
	assert(vector);

	ssize_t n = vector_dim(vector);
	ssize_t i;
	double x;

	for (i = 0; i < n; i++) {
		x = *vector_item_ptr(vector, i);
		*vector_item_ptr(vector, i) = exp(x);
	}
}

void vector_log(struct vector *vector)
{
	assert(vector);

	ssize_t n = vector_dim(vector);
	ssize_t i;
	double x;

	for (i = 0; i < n; i++) {
		x = *vector_item_ptr(vector, i);
		*vector_item_ptr(vector, i) = log(x);
	}
}

void vector_sqrt(struct vector *v)
{
	ssize_t n = vector_dim(v);
	ssize_t i;
	double x;

	for (i = 0; i < n; i++) {
		x = *vector_item_ptr(v, i);
		*vector_item_ptr(v, i) = sqrt(x);
	}
}

double vector_log_sum_exp(const struct vector *vector)
{
	assert(vector);

	ssize_t n = vector_dim(vector);

	if (n == 0)
		return -INFINITY;

	ssize_t imax = vector_max_index(vector);
	double max = *vector_item_ptr(vector, imax);
	double summ1 = 0.0;
	ssize_t i;

	for (i = 0; i < n; i++) {
		if (i == imax)
			continue;

		summ1 += exp(*vector_item_ptr(vector, i) - max);
	}

	return max + log1p(summ1);
}

void vector_printf(const struct vector *v)
{
	printf("x <- c(");

	ssize_t i, n = vector_dim(v);

	if (n > 0) {
		printf("%.22f", vector_item(v, 0));
	}

	for (i = 1; i < n; i++) {
		printf(", %.22f", vector_item(v, i));
	}
	printf(")\n");
}

uint32_t vector_hash(const void *v)
{
	const struct vector *vector = v;
	uint32_t seed = 0;
	ssize_t i, n = vector_dim(vector);

	for (i = 0; i < n; i++) {
		double x = *vector_item_ptr(vector, i);
		uint32_t hash_value = double_hash(&x);
		hash_combine(&seed, hash_value);
	}

	hash_finalize(&seed);

	return seed;
}

bool vector_equals(const void *v1, const void *v2)
{
	const struct vector *vector1 = v1;
	const struct vector *vector2 = v2;

	if (vector1 == vector2)
		return true;

	ssize_t n = vector_dim(vector1);

	if (vector_dim(vector2) != n)
		return false;

	ssize_t i;
	for (i = 0; i < n; i++) {
		double x1 = *vector_item_ptr(vector1, i);
		double x2 = *vector_item_ptr(vector2, i);

		if (!double_identical(x1, x2))
			return false;
	}

	return true;
}

int vector_compare(const void *x1, const void *x2)
{
	const struct vector *vector1 = x1;
	const struct vector *vector2 = x2;
	ssize_t n1 = vector_dim(vector1);
	ssize_t n2 = vector_dim(vector2);

	if (n1 < n2) {
		return -1;
	} else if (n1 > n2) {
		return +1;
	}

	if (n1 == 0)
		return 0;

	double *p1 = vector_to_ptr(vector1);
	double *p2 = vector_to_ptr(vector2);
	ssize_t i, n = n1;

	for (i = 0; i < n; i++) {
		int cmp = double_compare(p1 + i, p2 + i);
		if (cmp != 0)
			return cmp;
	}

	return 0;
}

int vector_ptr_compare(const void *px1, const void *px2)
{
	struct vector *const *pvector1 = px1;
	struct vector *const *pvector2 = px2;

	if (!pvector1)
		return pvector2 ? 0 : -1;
	if (!pvector2)
		return +1;

	return vector_compare(*pvector1, *pvector2);
}

void vector_insert(struct vector *v, ssize_t i)
{
	assert(v);
	assert(0 <= i && i <= vector_dim(v));
	assert(vector_owner(v));
	
	ssize_t n0 = vector_dim(v);
	ssize_t n1 = n0 + 1;
	
	vector_reinit(v, n1);
	double *ptr = vector_to_ptr(v);
	ssize_t j;
	for (j = n1 - 1; j > i; j--) {
		ptr[j] = ptr[j - 1];
	}
	ptr[i] = 0.0;
}

