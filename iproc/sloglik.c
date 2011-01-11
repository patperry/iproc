#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <assert.h>
#include <math.h>
#include <iproc/memory.h>
#include <iproc/sloglik.h>

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
    int64_t isend = sll->isend;
    iproc_model *model = sll->model;
    iproc_vars *vars = model->vars;
    iproc_vars_ctx *ctx = iproc_vars_ctx_new(vars, history, isend);
    int64_t nrecv = iproc_vars_nrecv(vars);
    iproc_svector *logprobs = iproc_svector_new(nrecv);
    double logprob0_shift = -INFINITY;

    iproc_model_get_new_logprobs(model, ctx, &logprob0_shift, logprobs);
    int64_t jnz = iproc_svector_find_nz(logprobs, jrecv);

    if (jrecv == isend && !model->has_loops) {
        sll->value += -INFINITY;
    } else if (jnz >= 0) {
        sll->value += iproc_svector_nz(logprobs, jnz);
    } else {
        sll->value += (iproc_model_logprob0(model, isend, jrecv) + logprob0_shift);
    }
    
    /* update gradient vars */

    iproc_svector_unref(logprobs);
    iproc_vars_ctx_unref(ctx);
}

double
iproc_vector_acc_sloglik_grad (iproc_vector  *dst_vector,
                               double         scale,
                               iproc_sloglik *sll)
{
    return sll->value;
}
