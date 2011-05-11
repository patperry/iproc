#ifndef _IPROC_LOGLIK_H
#define _IPROC_LOGLIK_H

#include "array.h"
#include "frame.h"
#include "messages.h"
#include "model.h"
#include "refcount.h"
#include "vector.h"

typedef struct _iproc_loglik iproc_loglik;

struct _iproc_loglik {
	iproc_model *model;
	struct array sloglik_array;
	struct vector *grad;
	bool grad_cached;
	ssize_t nsend;
	ssize_t nrecv;
	struct refcount refcount;
};

iproc_loglik *iproc_loglik_new(iproc_model * model, struct messages *messages);
iproc_loglik *iproc_loglik_ref(iproc_loglik * loglik);
void iproc_loglik_unref(iproc_loglik * loglik);

void iproc_loglik_insert(iproc_loglik * loglik,
			 const struct frame *f, const struct message *msg);

double iproc_loglik_value(iproc_loglik * loglik);
struct vector *iproc_loglik_grad(iproc_loglik * loglik);

#endif /* _IPROC_LOGLIK_H */
