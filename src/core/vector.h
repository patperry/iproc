#ifndef _VECTOR_H
#define _VECTOR_H

#include <stddef.h>

struct vector {
	double *data;
	ssize_t dim;
	bool owner;
};

/* create, destroy */
void vector_init(struct vector *v, ssize_t n);
void vector_reinit(struct vector *v, ssize_t n);
void vector_init_copy(struct vector *v, const struct vector *src);
void vector_assign_copy(struct vector *v, const struct vector *src);
void vector_deinit(struct vector *v);

/* views */
static inline struct vector vector_make(const double *ptr, ssize_t n);
static inline struct vector vector_slice(const struct vector *v, ssize_t i,
					 ssize_t n);

/* properties */
static inline ssize_t vector_dim(const struct vector *v);
static inline double *vector_to_ptr(const struct vector *v);

/* index */
static inline double vector_item(const struct vector *v, ssize_t i);
static inline double *vector_item_ptr(const struct vector *v, ssize_t i);
static inline void vector_set_item(struct vector *v, ssize_t i, double val);

/* copy, fill */
void vector_copy_to(const struct vector *v, double *dst);
void vector_fill(struct vector *v, double val);
void vector_set_basis(struct vector *v, ssize_t i);

/* hash, compare */
uint32_t vector_hash(const void *v);
bool vector_equals(const void *v1, const void *v2);
int vector_compare(const void *v1, const void *v2);
int vector_ptr_compare(const void *pv1, const void *pv2);

/* arithmetic operations */
void vector_scale(struct vector *v, double scale);
void vector_shift(struct vector *v, double shift);
void vector_add(struct vector *v1, const struct vector *v2);
void vector_sub(struct vector *v1, const struct vector *v2);
void vector_mul(struct vector *v1, const struct vector *v2);
void vector_div(struct vector *v1, const struct vector *v2);

/* special functions */
void vector_exp(struct vector *v);
void vector_log(struct vector *v);
void vector_sqrt(struct vector *v);

/* other arithmetic operations */
double vector_log_sum_exp(const struct vector *vector);

/* linear algebra */
double vector_dot(const struct vector *v1, const struct vector *v2);
double vector_norm(const struct vector *v);	// Eucludean norm
double vector_norm2(const struct vector *v);	// Square of Euclidean
double vector_sum_abs(const struct vector *v);
void vector_axpy(double alpha, const struct vector *x, struct vector *y);
double vector_dist(const struct vector *v1, const struct vector *v2);

/* min and max */
double vector_max(const struct vector *vector);
ssize_t vector_max_index(const struct vector *vector);
double vector_max_abs(const struct vector *vector);
ssize_t vector_max_abs_index(const struct vector *vector);

/* DEPRECATED */
struct vector *vector_alloc(ssize_t n);
struct vector *vector_alloc_copy(const struct vector *v);
void vector_free(struct vector *v);

/* inline function definitions */
struct vector vector_make(const double *ptr, ssize_t n)
{
	assert(ptr || n == 0);
	assert(n >= 0);

	struct vector v;
	v.data = (double *)ptr;
	v.dim = n;
	v.owner = false;
	return v;
}

struct vector vector_slice(const struct vector *v, ssize_t i, ssize_t n)
{
	assert(v);
	assert(n >= 0);
	assert(0 <= i && i <= vector_dim(v) - n);

	return vector_make(v->data + i, n);
}

ssize_t vector_dim(const struct vector *v)
{
	return v->dim;
}

double *vector_to_ptr(const struct vector *v)
{
	return v->data;
}

double *vector_item_ptr(const struct vector *v, ssize_t i)
{
	assert(0 <= i && i < vector_dim(v));
	return &v->data[i];
}

double vector_item(const struct vector *v, ssize_t i)
{
	return *vector_item_ptr(v, i);
}

void vector_set_item(struct vector *v, ssize_t i, double val)
{
	*vector_item_ptr(v, i) = val;
}

#endif /* _VECTOR_H */
