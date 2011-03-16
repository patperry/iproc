#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <assert.h>
#include <math.h>
#include "memory.h"
#include "sloglik.h"
#include "loglik.h"

static void
iproc_loglik_free (iproc_loglik *loglik)
{
    if (loglik) {
        iproc_array *array = loglik->sloglik_array;
        int64_t i, n = iproc_array_size(array);
        
        for (i = 0; i < n; i++) {
            iproc_sloglik *sll = iproc_array_index(array, iproc_sloglik *, i);
            iproc_sloglik_unref(sll);
        }

        iproc_array_unref(array);
        iproc_model_unref(loglik->model);
        iproc_vector_unref(loglik->grad);
        iproc_free(loglik);
    }
}

static iproc_loglik *
iproc_loglik_new_empty (iproc_model *model)
{
    iproc_loglik *loglik = iproc_malloc(sizeof(*loglik));
    iproc_design *design = iproc_model_design(model);
    int64_t nsender = iproc_design_nsender(design);
    
    if (!loglik)
        return NULL;
    
    loglik->sloglik_array = iproc_array_new(sizeof(iproc_sloglik *));
    loglik->model = iproc_model_ref(model);
    loglik->grad = iproc_vector_new(iproc_design_dim(design));
    loglik->grad_cached = 0;
    iproc_refcount_init(&loglik->refcount);
    
    if (!(loglik->sloglik_array && loglik->grad)) {
        iproc_loglik_free(loglik);
        loglik = NULL;
    }
    
    iproc_array_set_size(loglik->sloglik_array, nsender);
    
    return loglik;
}


iproc_loglik *
iproc_loglik_new (iproc_model    *model,
                  iproc_messages *messages)
{
    assert(model);
    assert(!messages || iproc_model_nsender(model) <= iproc_messages_max_from(messages));
    assert(!messages || iproc_model_nreceiver(model) <= iproc_messages_max_to(messages));
    
    iproc_loglik *loglik = iproc_loglik_new_empty(model);
    
    if (!messages)
        return loglik;
    
    iproc_message_iter *it = iproc_message_iter_new(messages);
    
    while (iproc_message_iter_next(it)) {
        iproc_history *history = iproc_message_iter_history(it);
        int64_t tie, ntie = iproc_message_iter_ntie(it);
        
        for (tie = 0; tie < ntie; tie++) {
            iproc_message_iter_select(it, tie);
            
            int64_t from = iproc_message_iter_from(it);
            int64_t *to = iproc_message_iter_to(it);
            int64_t nto = iproc_message_iter_nto(it);
            
            iproc_loglik_insertm(loglik, history, from, to, nto);
        }
    }
    
    iproc_message_iter_unref(it);
    return loglik;
}


iproc_loglik *
iproc_loglik_ref (iproc_loglik *loglik)
{
    if (loglik) {
        iproc_refcount_get(&loglik->refcount);
    }

    return loglik;
}

static void
iproc_loglik_release (iproc_refcount *refcount)
{
    iproc_loglik *loglik = container_of(refcount, iproc_loglik, refcount);
    iproc_loglik_free(loglik);
}

void
iproc_loglik_unref (iproc_loglik *loglik)
{
    if (loglik) {
        iproc_refcount_put(&loglik->refcount, iproc_loglik_release);
    }
}

static iproc_sloglik *
iproc_loglik_sloglik (iproc_loglik *loglik,
                      int64_t       isend)
{
    iproc_array *array = loglik->sloglik_array;
    iproc_sloglik *sll = iproc_array_index(array, iproc_sloglik *, isend);

    if (!sll) {
        iproc_model *model = loglik->model;
        sll = iproc_sloglik_new(model, isend);
        iproc_array_index(array, iproc_sloglik *, isend) = sll;
    }

    return sll;
}

void
iproc_loglik_insert (iproc_loglik  *loglik,
                     iproc_history *history,
                     int64_t        from,
                     int64_t        to)

{
    iproc_sloglik *sll = iproc_loglik_sloglik(loglik, from);
    iproc_sloglik_insert(sll, history, to);
}

void
iproc_loglik_insertm  (iproc_loglik  *loglik,
                       iproc_history *history,
                       int64_t        from,
                       int64_t       *to,
                       int64_t        nto)
{
    iproc_sloglik *sll = iproc_loglik_sloglik(loglik, from);
    iproc_sloglik_insertm(sll, history, to, nto);
    loglik->grad_cached = 0;
}

double
iproc_loglik_value (iproc_loglik *loglik)
{
    if (!loglik)
        return 0.0;

    iproc_array *array = loglik->sloglik_array;
    int64_t i, n = iproc_array_size(array);
    iproc_sloglik *sll;
    double value = 0.0;
    
    for (i = 0; i < n; i++) {
        sll = iproc_array_index(array, iproc_sloglik *, i);
        value += iproc_sloglik_value(sll);
    }
    
    return value;
}


static void
iproc_vector_acc_loglik_grad_nocache (iproc_vector *dst_vector,
                                      double        scale,
                                      iproc_loglik *loglik)
{
    iproc_array *array = loglik->sloglik_array;
    int64_t nsend = iproc_array_size(array);
    int64_t i;
    iproc_sloglik *sll;
    
    for (i = 0; i < nsend; i++) {
        sll = iproc_array_index(array, iproc_sloglik *, i);
        if (sll) {
            iproc_vector *g = iproc_sloglik_grad(sll);
            iproc_vector_acc(dst_vector, scale, g);
        }
    }
}


iproc_vector *
iproc_loglik_grad (iproc_loglik *loglik)
{
    if (!loglik)
        return NULL;

    if (!loglik->grad_cached) {
        iproc_vector_set_all(loglik->grad, 0.0);
        iproc_vector_acc_loglik_grad_nocache(loglik->grad, 1.0, loglik);
        loglik->grad_cached = 1;
    }
    return loglik->grad;
}
