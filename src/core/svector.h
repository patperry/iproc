#ifndef _IPROC_SVECTOR_H
#define _IPROC_SVECTOR_H

#include <stdint.h>
#include "darray.h"
#include "refcount.h"
#include "vector.h"

typedef struct svector iproc_svector;

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
void svector_free(iproc_svector * svector);

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


void iproc_svector_inc(iproc_svector * svector, int64_t i, double value);


int64_t iproc_svector_nz(const iproc_svector * svector, int64_t inz);
double iproc_svector_nz_get(const iproc_svector * svector, int64_t inz);
void iproc_svector_nz_set(iproc_svector * svector, int64_t inz, double value);
void iproc_svector_nz_inc(iproc_svector * svector, int64_t inz, double inc);

iproc_vector_view iproc_svector_view_nz(const iproc_svector * svector);
int64_t iproc_svector_find_nz(const iproc_svector * svector, int64_t i);


#endif /* _IPROC_SVECTOR_H */
