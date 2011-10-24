#ifndef _SVECTOR_H
#define _SVECTOR_H

#include "sblas.h"
#include "vector.h"

struct svector {
	double *data;
	struct vpattern pattern;
	size_t dim;
	bool is_view;
};

struct svector_pos {
	size_t i;
	ssize_t inz;
};

struct svector_iter {
	const struct svector *x;
	ssize_t inz;
};

#define SVECTOR_IDX(it) ((it).x->pattern.indx[(it).inz])
#define SVECTOR_PTR(it) (&(it).x->data[(it).inz])
#define SVECTOR_VAL(it) (*(const double *)SVECTOR_PTR(it))

#define SVECTOR_FOREACH(it, v) \
	for ((it) = svector_iter_make(v); svector_iter_advance(&(it));)

/* init/assign/deinit */
void svector_init(struct svector *v, ssize_t n);
void svector_init_copy(struct svector *v, const struct svector *src);
void svector_assign_copy(struct svector *v, const struct svector *src);
void svector_deinit(struct svector *v);

/* views */
static inline struct svector svector_make(const double *data, const struct vpattern *pat, size_t n);


/* alloc/dealloc */
struct svector *svector_alloc(ssize_t n);
struct svector *svector_alloc_copy(const struct svector *v);
void svector_free(struct svector *v);

/* properties */
static inline size_t svector_dim(const struct svector *v);
static inline size_t svector_count(const struct svector *v);
static inline bool svector_owner(const struct svector *v);

double svector_item(const struct svector *v, size_t i);
double *svector_item_ptr(struct svector *v, size_t i);
void svector_set_item(struct svector *v, size_t i, double val);

double svector_max(const struct svector *v);

/* methods */
void svector_clear(struct svector *v);
void svector_remove(struct svector *v, size_t i);
void svector_set_basis(struct svector *v, size_t i);

/* linear algebra */
void svector_scale(struct svector *v, double scale);
double svector_dot(const struct svector *x, const struct vector *y);
double svector_dots(const struct svector *v1, const struct svector *v2);
void svector_axpy(double scale, const struct svector *x, struct vector *y);
void svector_axpys(double scale, const struct svector *x, struct svector *y);

/* sparse operations */
void svector_scatter(const struct svector *src, struct vector *dst);
void svector_gather(struct svector *dst, const struct vector *src);

/* position-based operations */
double *svector_find(const struct svector *v, size_t i,
		     struct svector_pos *pos);
double *svector_insert(struct svector *v, struct svector_pos *pos, double val);
void svector_remove_at(struct svector *v, struct svector_pos *pos);

/* iteration */
static inline struct svector_iter svector_iter_make(const struct svector *v);
static inline bool svector_iter_advance(struct svector_iter *it);
static inline void svector_iter_reset(struct svector_iter *it);

/* low-level access */
static inline double *svector_data_ptr(const struct svector *v);
static inline size_t *svector_index_ptr(const struct svector *v);

/* debug */
void svector_printf(const struct svector *v);

/* inline function definitions */
struct svector svector_make(const double *data, const struct vpattern *pat,
			    size_t n)
{
	struct svector v;
	v.data = (double *)data;
	v.pattern = *pat;
	v.dim = n;
	v.is_view = true;
	return v;
}

size_t svector_dim(const struct svector *v)
{
	assert(v);
	return v->dim;
}

size_t svector_count(const struct svector *v)
{
	assert(v);
	return v->pattern.nz;
}

bool svector_owner(const struct svector *v)
{
	assert(v);
	return !v->is_view;
}

double *svector_data_ptr(const struct svector *v)
{
	assert(v);
	return v->data;
}

size_t *svector_index_ptr(const struct svector *v)
{
	assert(v);
	return (size_t *)v->pattern.indx;
}

struct svector_iter svector_iter_make(const struct svector *v)
{
	assert(v);
	struct svector_iter it;
	it.x = v;
	it.inz = -1;
	return it;
}

bool svector_iter_advance(struct svector_iter *it)
{
	assert(it);
	it->inz++;
	return (size_t)it->inz < it->x->pattern.nz;
}

void svector_iter_reset(struct svector_iter *it)
{
	assert(it);
	it->inz = -1;
}


#endif /* _SVECTOR_H */
