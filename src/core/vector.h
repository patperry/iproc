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
void vector_init_copy(struct vector *v, const struct vector *src);
void vector_assign_copy(struct vector *v, const struct vector *src);
void vector_deinit(struct vector *v);

/* views */
void vector_init_view(struct vector *v, const double *ptr, ssize_t n);
void vector_init_slice(struct vector *v,
		       const struct vector *parent, ssize_t i, ssize_t n);

/* assign, copy, fill */
void vector_assign_array(struct vector *v, const double *src);

void vector_copy_to(const struct vector *v, double *dst);

void vector_fill(struct vector *v, double val);
void vector_set_basis(struct vector *v, ssize_t i);

/* hash, compare */
uint32_t vector_hash(const void *v);
bool vector_equals(const void *v1, const void *v2);
int vector_compare(const void *v1, const void *v2);
int vector_ptr_compare(const void *pv1, const void *pv2);

/* index */
#define vector_front(v)          vector_at(v, 0)
#define vector_get_front(v)      vector_get(v, 0)
#define vector_set_front(v, val) vector_set(v, 0, val)
#define vector_back(v)           vector_at(v, vector_dim(v) - 1)
#define vector_get_back(v)       vector_get(v, vector_dim(v) - 1)
#define vector_set_back(v)       vector_set(v, vector_dim(v) - 1, val)
static inline double *vector_at(const struct vector *v, ssize_t i);
static inline double vector_get(const struct vector *v, ssize_t i);
static inline void vector_set(struct vector *v, ssize_t i, double val);

/* informative */
static inline ssize_t vector_dim(const struct vector *v);

/* arithmetic operations */
void vector_scale(struct vector *v, double scale);
void vector_shift(struct vector *v, double shift);
void vector_add(struct vector *v1, const struct vector *v2);
void vector_sub(struct vector *v1, const struct vector *v2);
void vector_mul(struct vector *v1, const struct vector *v2);
void vector_div(struct vector *v1, const struct vector *v2);

/* special functions */
void vector_exp(struct vector *v);

/* other arithmetic operations */
double vector_log_sum_exp(const struct vector *vector);

/* linear algebra */
double vector_dot(const struct vector *v1, const struct vector *v2);
double vector_norm(const struct vector *v);
double vector_norm1(const struct vector *v);
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

struct vector vector_slice(const struct vector *vector,
			   ssize_t index, ssize_t dim);

/* inline function definitions */
ssize_t vector_dim(const struct vector *v)
{
	return v->dim;
}

double *vector_at(const struct vector *v, ssize_t i)
{
	assert(0 <= i && i < vector_dim(v));
	return &v->data[i];
}

double vector_get(const struct vector *v, ssize_t i)
{
	return *vector_at(v, i);
}

void vector_set(struct vector *v, ssize_t i, double val)
{
	*vector_at(v, i) = val;
}

#endif /* _VECTOR_H */
