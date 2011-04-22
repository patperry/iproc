#ifndef _IPROC_LOGLIK_H
#define _IPROC_LOGLIK_H

#include <stdbool.h>

#include "darray.h"
#include "messages.h"
#include "model.h"
#include "refcount.h"
#include "vector.h"

typedef struct _iproc_loglik iproc_loglik;

struct _iproc_loglik {
	iproc_model *model;
	struct darray sloglik_array;
	struct vector *grad;
	bool grad_cached;
	int64_t nsend;
	int64_t nrecv;
	struct refcount refcount;
};

iproc_loglik *iproc_loglik_new(iproc_model * model, iproc_messages * messages);
iproc_loglik *iproc_loglik_ref(iproc_loglik * loglik);
void iproc_loglik_unref(iproc_loglik * loglik);

void iproc_loglik_insert(iproc_loglik * loglik,
			 iproc_history * history, int64_t from, int64_t to);
void iproc_loglik_insertm(iproc_loglik * loglik,
			  iproc_history * history,
			  int64_t from, int64_t *to, int64_t nto);

double iproc_loglik_value(iproc_loglik * loglik);
struct vector *iproc_loglik_grad(iproc_loglik * loglik);

#endif /* _IPROC_LOGLIK_H */
