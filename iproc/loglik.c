#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <assert.h>
#include <math.h>
#include <iproc/memory.h>
#include <iproc/sloglik.h>
#include <iproc/loglik.h>

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
        iproc_free(loglik);
    }
}

iproc_loglik *
iproc_loglik_new (iproc_model *model)
{
    iproc_loglik *loglik = iproc_malloc(sizeof(*loglik));
    iproc_vars *vars = model->vars;
    int64_t nsend = iproc_vars_nsend(vars);

    if (!loglik)
        return NULL;

    loglik->sloglik_array = iproc_array_new(sizeof(iproc_sloglik *));
    loglik->model = iproc_model_ref(model);
    iproc_refcount_init(&loglik->refcount);

    if (!loglik->sloglik_array) {
        iproc_loglik_free(loglik);
        loglik = NULL;
    }

    iproc_array_set_size(loglik->sloglik_array, nsend);

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

void
iproc_loglik_insert (iproc_loglik  *loglik,
                     int64_t        from,
                     int64_t        to,
                     iproc_history *history)
{
    iproc_array *array = loglik->sloglik_array;
    iproc_sloglik *sll = iproc_array_index(array, iproc_sloglik *, from);

    if (!sll) {
        iproc_model *model = loglik->model;
        sll = iproc_sloglik_new(model, from);
        iproc_array_index(array, iproc_sloglik *, from) = sll;
    }

    iproc_sloglik_insert(sll, to, history);
}

double
iproc_vector_acc_loglik_grad (iproc_vector *dst_vector,
                              double        scale,
                              iproc_loglik *loglik)
{
    iproc_array *array = loglik->sloglik_array;
    int64_t nsend = iproc_array_size(array);
    int64_t i;
    iproc_sloglik *sll;
    double value = 0.0;
    
    for (i = 0; i < nsend; i++) {
        sll = iproc_array_index(array, iproc_sloglik *, i);
        if (sll) {
            value += iproc_vector_acc_sloglik_grad(dst_vector, scale, sll);
        }
    }
    
    return value;
}
