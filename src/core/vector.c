#include "port.h"

#include <assert.h>
#include <inttypes.h>
#include <stddef.h>
#include <stdio.h>
#include <math.h>
#include "blas-private.h"
#include "compare.h"
#include "hash.h"
#include "ieee754.h"
#include "memory.h"
#include "vector.h"

struct vector *vector_init(struct vector *v, ssize_t n)
{
	assert(v);
	assert(n >= 0);
	assert(n <= F77_INT_MAX);
	assert(n <= SSIZE_MAX / sizeof(double));

	if (array_init(&v->array, double, n)) {
		return v;
	}

	return NULL;
}

struct vector *vector_init_view(struct vector *v, double *ptr, ssize_t n)
{
	assert(v);
	assert(ptr || n == 0);
	assert(n >= 0);

	if (array_init_view(&v->array, double, ptr, n)) {
		return v;
	}

	return NULL;
}

struct vector *vector_init_slice(struct vector *v,
				 const struct vector *parent,
				 ssize_t i, ssize_t n)
{
	assert(v);
	assert(parent);
	assert(n >= 0);
	assert(0 <= i && i <= vector_size(parent) - n);

	return vector_init_view(v, vector_ptr(parent, i), n);
}

struct vector *vector_init_copy(struct vector *v, const struct vector *src)
{
	assert(v);
	assert(src);

	if (vector_init(v, vector_size(src))) {
		if (vector_assign_copy(v, src)) {
			return v;
		}
		vector_deinit(v);
	}

	return NULL;
}

struct vector *vector_new(ssize_t n)
{
	struct vector *v = malloc(sizeof(*v));

	if (vector_init(v, n)) {
		return v;
	}

	free(v);
	return NULL;
}

struct vector *vector_new_copy(const struct vector *src)
{
	struct vector *v = malloc(sizeof(*v));

	if (vector_init_copy(v, src)) {
		return v;
	}

	free(v);
	return NULL;
}

void vector_deinit(struct vector *v)
{
	assert(v);
	array_deinit(&v->array);
}

void vector_free(struct vector *v)
{
	if (v) {
		vector_deinit(v);
		free(v);
	}
}

struct vector *vector_assign_copy(struct vector *v, const struct vector *src)
{
	assert(v);
	assert(src);
	assert(vector_size(v) == vector_size(src));

	vector_copy_to(src, vector_begin(v));
	return v;
}

double *vector_copy_to(const struct vector *v, double *dst)
{
	assert(v);
	assert(dst || vector_size(v) == 0);

	return array_copy_to(&v->array, dst);
}

void vector_fill(struct vector *vector, double value)
{
	assert(vector);

	ssize_t n = vector_size(vector);
	double *ptr = vector_ptr(vector, 0);
	double *end = vector_ptr(vector, n);

	while (ptr != end) {
		*ptr++ = value;
	}
}

void vector_set_basis(struct vector *v, ssize_t i)
{
	assert(v);
	assert(vector_size(v) > 0);
	assert(0 <= i && i < vector_size(v));

	vector_fill(v, 0.0);
	vector_index(v, i) = 1.0;
}

iproc_vector_view vector_slice(struct vector *v, ssize_t i, ssize_t n)
{
	iproc_vector_view view;
	vector_init_slice(&view.vector, v, i, n);
	return view;
}

iproc_vector_view iproc_vector_view_array(double *ptr, ssize_t n)
{
	iproc_vector_view view;
	vector_init_view(&view.vector, ptr, n);
	return view;
}

void vector_scale(struct vector *vector, double scale)
{
	assert(vector);

	f77int n = (f77int) vector_size(vector);
	double alpha = scale;
	void *px = vector_ptr(vector, 0);
	f77int incx = 1;

	F77_FUNC(dscal) (&n, &alpha, px, &incx);
}

void vector_shift(struct vector *vector, double shift)
{
	assert(vector);

	ssize_t n = vector_size(vector);
	double *ptr = vector_ptr(vector, 0);
	double *end = vector_ptr(vector, n);

	while (ptr != end) {
		*ptr++ += shift;
	}
}

void vector_add(struct vector *dst_vector, const struct vector *vector)
{
	assert(dst_vector);
	assert(vector);
	assert(vector_size(dst_vector) == vector_size(vector));

	vector_axpy(1.0, vector, dst_vector);
}

void vector_sub(struct vector *dst_vector, const struct vector *vector)
{
	assert(dst_vector);
	assert(vector);
	assert(vector_size(dst_vector) == vector_size(vector));

	vector_axpy(-1.0, vector, dst_vector);
}

void vector_mul(struct vector *dst_vector, const struct vector *vector)
{
	assert(dst_vector);
	assert(vector);
	assert(vector_size(dst_vector) == vector_size(vector));

	f77int n = (f77int) vector_size(dst_vector);
	f77int k = 0;
	double *px = vector_ptr(vector, 0);
	f77int incx = 1;
	double *py = vector_ptr(dst_vector, 0);
	f77int incy = 1;

	F77_FUNC(dtbmv) ("U", "N", "N", &n, &k, px, &incx, py, &incy);
}

void vector_div(struct vector *dst_vector, const struct vector *vector)
{
	assert(dst_vector);
	assert(vector);
	assert(vector_size(dst_vector) == vector_size(vector));

	f77int n = (f77int) vector_size(dst_vector);
	f77int k = 0;
	double *px = vector_ptr(vector, 0);
	f77int incx = 1;
	double *py = vector_ptr(dst_vector, 0);
	f77int incy = 1;

	F77_FUNC(dtbsv) ("U", "N", "N", &n, &k, px, &incx, py, &incy);
}

void vector_axpy(double alpha, const struct vector *x, struct vector *y)
{
	assert(x);
	assert(y);
	assert(vector_size(x) == vector_size(y));

	f77int n = (f77int) vector_size(y);
	double *px = vector_begin(x);
	f77int incx = 1;
	double *py = vector_begin(y);
	f77int incy = 1;

	F77_FUNC(daxpy) (&n, &alpha, px, &incx, py, &incy);
}

double vector_dot(const struct vector *vector1, const struct vector *vector2)
{
	assert(vector1);
	assert(vector2);
	assert(vector_size(vector1) == vector_size(vector2));

	f77int n = (f77int) vector_size(vector1);
	double *px = vector_ptr(vector1, 0);
	f77int incx = 1;
	double *py = vector_ptr(vector2, 0);
	f77int incy = 1;

	double dot = F77_FUNC(ddot) (&n, px, &incx, py, &incy);
	return dot;
}

double vector_norm(const struct vector *vector)
{
	assert(vector);

	f77int n = (f77int) vector_size(vector);
	void *px = vector_ptr(vector, 0);
	f77int incx = 1;

	double norm = F77_FUNC(dnrm2) (&n, px, &incx);
	return norm;
}

double vector_norm1(const struct vector *vector)
{
	assert(vector);

	f77int n = (f77int) vector_size(vector);
	void *px = vector_ptr(vector, 0);
	f77int incx = 1;

	double sum_abs = F77_FUNC(dasum) (&n, px, &incx);
	return sum_abs;
}

double vector_max_abs(const struct vector *vector)
{
	assert(vector);

	double max_abs = 0;
	double e;
	ssize_t i;

	if (vector_size(vector) > 0) {
		i = vector_max_abs_index(vector);
		e = vector_index(vector, i);
		max_abs = fabs(e);
	}

	return max_abs;
}

ssize_t vector_max_abs_index(const struct vector * vector)
{
	assert(vector);
	if (!(vector_size(vector) > 0))
		return -1;

	f77int n = (f77int) vector_size(vector);
	void *px = vector_ptr(vector, 0);
	f77int incx = 1;

	f77int index1 = F77_FUNC(idamax) (&n, px, &incx);

	ssize_t index = index1 - 1;
	return index;
}

ssize_t vector_max_index(const struct vector * vector)
{
	assert(vector);
	assert(vector_size(vector) > 0);

	ssize_t n = vector_size(vector);
	ssize_t i, imax;
	double x, max = NAN;

	/* Find the first non-NaN entry of the vector */
	for (imax = 0; imax < n && isnan(max); imax++) {
		max = vector_index(vector, imax);
	}

	/* If all of the entries are NaN, define imax as 0. */
	if (imax == n)
		return 0;

	/* Otherwise, search for the largest entry in the tail of the vector */
	for (i = imax + 1; i < n; i++) {
		x = vector_index(vector, i);
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
	if (vector_size(vector) == 0) {
		return -INFINITY;
	} else {
		ssize_t i = vector_max_index(vector);
		return vector_index(vector, i);
	}
}

void vector_exp(struct vector *vector)
{
	assert(vector);
	ssize_t n = vector_size(vector);
	ssize_t i;
	double x;

	for (i = 0; i < n; i++) {
		x = vector_index(vector, i);
		vector_index(vector, i) = exp(x);
	}
}

double vector_log_sum_exp(const struct vector *vector)
{
	assert(vector);
	ssize_t n = vector_size(vector);

	if (n == 0)
		return -INFINITY;

	ssize_t imax = vector_max_index(vector);
	double max = vector_index(vector, imax);
	double summ1 = 0.0;
	ssize_t i;

	for (i = 0; i < n; i++) {
		if (i == imax)
			continue;

		summ1 += exp(vector_index(vector, i) - max);
	}

	return max + log1p(summ1);
}

void vector_printf(const struct vector *v)
{
	printf("\nvector {");
	printf("\n  dim: %" SSIZE_FMT "", vector_size(v));
	printf("\n   nz: {");

	ssize_t i, n = vector_size(v);
	for (i = 0; i < n; i++) {
		if (vector_index(v, i) == 0.0)
			continue;

		printf("\n         %" SSIZE_FMT ", %.8f", i,
		       vector_index(v, i));
	}
	printf("\n       }");
	printf("\n}\n");
}

uint32_t vector_hash(const void *v)
{
	const struct vector *vector = v;
	if (!vector)
		return 0;

	uint32_t seed = 0;
	ssize_t i, n = vector_size(vector);

	for (i = 0; i < n; i++) {
		double x = vector_index(vector, i);
		uint32_t hash_value = double_hash(&x);
		hash_combine(&seed, hash_value);
	}

	hash_finalize(&seed);

	return seed;
}

bool
vector_identical(const struct vector * vector1, const struct vector * vector2)
{
	if (vector1 == vector2)
		return 1;

	ssize_t n = vector_size(vector1);

	if (vector_size(vector2) != n)
		return 0;

	ssize_t i;
	for (i = 0; i < n; i++) {
		double x1 = vector_index(vector1, i);
		double x2 = vector_index(vector2, i);

		if (!iproc_identical(x1, x2))
			return 0;
	}

	return 1;
}

int vector_compare(const void *x1, const void *x2)
{
	const struct vector *vector1 = x1;
	const struct vector *vector2 = x2;
	ssize_t n1 = vector_size(vector1);
	ssize_t n2 = vector_size(vector2);

	if (n1 < n2) {
		return -1;
	} else if (n1 > n2) {
		return +1;
	}

	double *p1 = vector_ptr(vector1, 0);
	double *p2 = vector_ptr(vector2, 0);
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
