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

bool svector_init(struct svector *v, ssize_t n);
bool svector_init_copy(struct svector *v, const struct svector *src);
void svector_deinit(struct svector *v);

struct svector *svector_new(ssize_t n);
struct svector *svector_new_copy(const struct svector *v);
void svector_free(iproc_svector * svector);


bool svector_assign_copy(struct svector *v, const struct svector *src);
void svector_clear(struct svector *v);

int64_t iproc_svector_dim(const iproc_svector * svector);
double iproc_svector_get(const iproc_svector * svector, int64_t i);
void iproc_svector_set(iproc_svector * svector, int64_t i, double value);
void iproc_svector_inc(iproc_svector * svector, int64_t i, double value);
void iproc_svector_scale(iproc_svector * svector, double scale);

int64_t iproc_svector_nnz(const iproc_svector * svector);
int64_t iproc_svector_nz(const iproc_svector * svector, int64_t inz);
double iproc_svector_nz_get(const iproc_svector * svector, int64_t inz);
void iproc_svector_nz_set(iproc_svector * svector, int64_t inz, double value);
void iproc_svector_nz_inc(iproc_svector * svector, int64_t inz, double inc);

iproc_vector_view iproc_svector_view_nz(const iproc_svector * svector);
int64_t iproc_svector_find_nz(const iproc_svector * svector, int64_t i);

double vector_dots(const struct vector *vector,
			 const iproc_svector * svector);
void iproc_vector_sacc(struct vector *dst_vector,
		       double scale, const iproc_svector * svector);

double iproc_svector_sdot(const iproc_svector * svector1,
			  const iproc_svector * svector2);
void iproc_svector_sacc(iproc_svector * dst_svector,
			double scale, const iproc_svector * svector);

void iproc_svector_printf(const iproc_svector * svector);

#endif /* _IPROC_SVECTOR_H */
