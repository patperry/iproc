#ifdef HAVE_CONFIG_H
#include <config.h>
#endif

#include <assert.h>
#include <math.h>
#include <iproc/memory.h>
#include <iproc/loglik.h>

static void
iproc_sloglik_free (iproc_sloglik *sll)
{
    if (sll) {
        iproc_svector_unref(sll->ovarsdiff);
        iproc_svector_unref(sll->newprob);
        iproc_svector_unref(sll->evarsdiff);
        iproc_free(sll);
    }
}

iproc_sloglik *
iproc_sloglik_new (iproc_model *model,
                   int64_t      isend)
{
    iproc_sloglik *sll = iproc_malloc(sizeof(*sll));
    if (!sll)
        return NULL;

    iproc_vars *vars = model->vars;
    int64_t n = iproc_vars_nrecv(vars);
    int64_t p = iproc_vars_dim(vars);

    sll->model = iproc_model_ref(model);
    sll->isend = isend;
    sll->nsend = 0;
    sll->nrecv = 0;
    sll->ovarsdiff = iproc_svector_new(p);
    sll->value = 0;
    sll->suminvwt = 0;
    sll->suminvwt_scale = -INFINITY;
    sll->newprob = iproc_svector_new(n);
    sll->evarsdiff = iproc_svector_new(p);
    iproc_refcount_init(&sll->refcount);

    if (!(sll->ovarsdiff && sll->newprob && sll->evarsdiff)) {
        iproc_sloglik_free(sll);
        sll = NULL;
    }

    return sll;
}

iproc_sloglik *
iproc_sloglik_ref (iproc_sloglik *sll)
{
    if (sll) {
        iproc_refcount_get(&sll->refcount);
    }
    return sll;
}

static void
iproc_sloglik_release (iproc_refcount *refcount)
{
    iproc_sloglik *sll = container_of(refcount, iproc_sloglik, refcount);
    iproc_sloglik_free(sll);
}

void
iproc_sloglik_unref (iproc_sloglik *sll)
{
    if (sll) {
        iproc_refcount_put(&sll->refcount, iproc_sloglik_release);
    }
}

void
iproc_sloglik_insert (iproc_sloglik *sll,
                      int64_t        jrecv,
                      iproc_history *history)
{

}

double
iproc_vector_acc_sloglik_grad (iproc_vector  *dst_vector,
                               double         scale,
                               iproc_sloglik *sll)
{
    return 0.0;
}
