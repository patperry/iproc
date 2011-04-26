#ifndef _SVECTOR_H
#define _SVECTOR_H

#include "intmap.h"
#include "vector.h"

struct svector {
	struct intmap map;
	ssize_t dim;
};

struct svector_iter {
	struct intmap_iter map_it;
};

bool svector_init(struct svector *v, ssize_t n);
bool svector_init_copy(struct svector *v, const struct svector *src);
void svector_deinit(struct svector *v);

struct svector *svector_alloc(ssize_t n);
struct svector *svector_alloc_copy(const struct svector *v);
void svector_free(struct svector *v);

bool svector_assign_copy(struct svector *v, const struct svector *src);
void svector_clear(struct svector *v);

ssize_t svector_dim(const struct svector* v);
ssize_t svector_size(const struct svector *v);
double svector_get(const struct svector *v, ssize_t i);
bool svector_set(struct svector* v, ssize_t i, double val);
double *svector_at(struct svector *v, ssize_t i);

double svector_max(const struct svector *v);

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

#endif /* _SVECTOR_H */
