#ifndef _SVECTOR_H
#define _SVECTOR_H

#include "darray.h"
#include "vector.h"

struct svector {
	ssize_t dim;
	struct darray index;
	struct darray value;
};

struct svector_iter {
	ssize_t pat_index;
};

bool svector_init(struct svector *v, ssize_t n);
bool svector_init_copy(struct svector *v, const struct svector *src);
void svector_deinit(struct svector *v);

struct svector *svector_new(ssize_t n);
struct svector *svector_new_copy(const struct svector *v);
void svector_free(struct svector *v);

bool svector_assign_copy(struct svector *v, const struct svector *src);
void svector_clear(struct svector *v);

ssize_t svector_dim(const struct svector* v);
ssize_t svector_size(const struct svector *v);
double svector_get(const struct svector *v, ssize_t i);
bool svector_set(struct svector* v, ssize_t i, double val);
double *svector_at(struct svector *v, ssize_t i);

void svector_scale(struct svector *v, double scale);
double svector_dot(const struct svector *x, const struct vector *y);
double svector_dots(const struct svector *v1, const struct svector *v2);
void svector_axpy(double scale, const struct svector *x, struct vector *y);
void svector_axpys(double scale, const struct svector* x, struct svector *y);

/* iteration */
void svector_iter_init(const struct svector *v, struct svector_iter *it);
bool svector_iter_advance(const struct svector *v, struct svector_iter *it);
double *svector_iter_current(const struct svector *v, struct svector_iter *it);
ssize_t svector_iter_current_index(const struct svector *v, struct svector_iter *it);
void svector_iter_deinit(const struct svector *v, struct svector_iter *it);

/* deprecated */
ssize_t iproc_svector_nz(const struct svector * svector, ssize_t inz);
double iproc_svector_nz_get(const struct svector * svector, ssize_t inz);
void iproc_svector_nz_set(struct svector * svector, ssize_t inz, double value);
void iproc_svector_nz_inc(struct svector * svector, ssize_t inz, double inc);
iproc_vector_view iproc_svector_view_nz(const struct svector * svector);


#endif /* _SVECTOR_H */
