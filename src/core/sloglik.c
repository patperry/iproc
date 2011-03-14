#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <assert.h>
#include <inttypes.h>
#include <math.h>
#include <stdio.h>
#include "memory.h"
#include "sloglik.h"

static void
iproc_sloglik_free (iproc_sloglik *sll)
{
    if (sll) {
        iproc_model_unref(sll->model);
        iproc_svector_unref(sll->sum_obs_var_diff);
        iproc_svector_unref(sll->sum_active_probs);
        iproc_svector_unref(sll->sum_mean_var_diff);
        iproc_vector_unref(sll->grad);
        iproc_svector_unref(sll->nrecv);
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

    int64_t n = iproc_model_nreceiver(model);
    int64_t p = iproc_model_dim(model);

    sll->model = iproc_model_ref(model);
    sll->isend = isend;
    sll->nsend = 0;
    sll->nrecv = iproc_svector_new(n);
    sll->grad = iproc_vector_new(p);
    sll->grad_cached = 0;
    sll->sum_obs_var_diff = iproc_svector_new(p);
    sll->value = 0.0;
    sll->suminvwt = 0.0;
    sll->sum_active_probs = iproc_svector_new(n);
    sll->sum_mean_var_diff = iproc_svector_new(p);
    iproc_refcount_init(&sll->refcount);

    if (!(sll->grad && sll->sum_obs_var_diff && sll->sum_active_probs
          && sll->sum_mean_var_diff)) {
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
                      iproc_history *history,
                      int64_t        jrecv)
{
    iproc_sloglik_insertm(sll, history, &jrecv, 1);
}

/* update active probs;
 * If we add a new nonzero, we need to retroactively make the entry active;
 * if a prob is no longer active, we need to add the old prob.
 */
static void
insertm_update_active_probs (iproc_sloglik   *sll,
                             int64_t          n,
                             iproc_model_ctx *ctx)
{
    int64_t isend = sll->isend;
    iproc_model *model = sll->model;
    iproc_vector *logprobs0 = iproc_model_logprobs0(model, isend);
    double logsuminvwt = log(sll->suminvwt);
    double logsumweight_diff = iproc_model_ctx_logsumweight_diff(ctx);

    iproc_svector *s = sll->sum_active_probs;
    iproc_svector *x = iproc_model_ctx_active_probs(ctx);
    int64_t inz_s = 0, nnz_s = iproc_svector_nnz(s);
    int64_t inz_x = 0, nnz_x = iproc_svector_nnz(x);

    while (inz_s < nnz_s && inz_x < nnz_x) {
        int64_t jrecv_s = iproc_svector_nz(s, inz_s);
        int64_t jrecv_x = iproc_svector_nz(x, inz_x);

        /* Update old nonzeros in s up to jrecv_x (these are not active in x) */
        if (jrecv_s < jrecv_x) {
            double lp0 = iproc_vector_get(logprobs0, jrecv_s);
            double p1 = exp(lp0 - logsumweight_diff);
            iproc_svector_nz_inc(s, inz_s, n * p1);
            inz_s++;

        /* If we are at an index which is active in both, update it */
        } else if (jrecv_s == jrecv_x) {
            double p1 = iproc_svector_nz_get(x, inz_x);
            iproc_svector_nz_inc(s, inz_s, n * p1);

            inz_s++;
            inz_x++;

        /* Otherwise, we are at a new nonzero; retroactively make it active */
        } else {
            double lp0 = iproc_vector_get(logprobs0, jrecv_x);
            double np0 = exp(lp0 + logsuminvwt);
            double p1 = iproc_svector_nz_get(x, inz_x);
            double value = np0 + n * p1;

            iproc_array_insert(s->index, inz_s, &jrecv_x);
            iproc_array_insert(s->value, inz_s, &value);
            inz_s++;
            nnz_s++;
            inz_x++;
        }
    }

    /* Leftovers in s are old nonzeros */
    while (inz_s < nnz_s) {
        int64_t jrecv = iproc_svector_nz(s, inz_s);
        double lp0 = iproc_vector_get(logprobs0, jrecv);
        double p1 = exp(lp0 - logsumweight_diff);
        iproc_svector_nz_inc(s, inz_s, n * p1);

        inz_s++;
    }

    /* Leftovers in x are new nonzeros */
    while (inz_x < nnz_x) {
        int64_t jrecv = iproc_svector_nz(x, inz_x);
        double lp0 = iproc_vector_get(logprobs0, jrecv);
        double np0 = exp(lp0 + logsuminvwt);
        double p1 = iproc_svector_nz_get(x, inz_x);
        double value = np0 + n * p1;

        iproc_array_append(s->index, &jrecv);
        iproc_array_append(s->value, &value);

        inz_x++;
    }


    sll->suminvwt += n * iproc_model_ctx_invsumweight_ratio(ctx);
    // iproc_svector_sacc(sll->sum_active_probs, n, active_probs);
}

static void
insertm_update_sum_mean_var_diff (iproc_sloglik   *sll,
                                  int64_t          n,
                                  iproc_model_ctx *ctx)
{
    iproc_svector *active_probs = iproc_model_ctx_active_probs(ctx);
    iproc_design_ctx_diff_muls(n, IPROC_TRANS_TRANS, ctx->design_ctx, active_probs,
                             1.0, sll->sum_mean_var_diff);
}

void
iproc_sloglik_insertm (iproc_sloglik *sll,
                       iproc_history *history,
                       int64_t       *jrecv,
                       int64_t        n)
{
    int64_t isend = sll->isend;
    iproc_model *model = sll->model;
    int64_t nreceiver = iproc_model_nreceiver(model);
    iproc_model_ctx *ctx = iproc_model_ctx_new(model, isend, history);
    int64_t i;
    iproc_svector *wt = iproc_svector_new(nreceiver);

    sll->grad_cached = 0;

    for (i = 0; i < n; i++) {
        assert(jrecv[i] >= 0);
        assert(jrecv[i] < nreceiver);

        /* update log likelihood */
        sll->value += iproc_model_ctx_logprob(ctx, jrecv[i]);
        iproc_svector_inc(wt, jrecv[i], 1.0);

        /* update number of receives */
        iproc_svector_inc(sll->nrecv, jrecv[i], 1.0);
    }

    /* update observed variables */
    iproc_design_ctx_diff_muls(1.0, IPROC_TRANS_TRANS, ctx->design_ctx, wt,
                             1.0, sll->sum_obs_var_diff);


    // printf("\nnrecv");
    // printf("\n-----");
    // iproc_svector_printf(sll->nrecv);

    // printf("\ninvsumweight_ratio");
    // printf("\n------------------");
    // printf("\n%.8f\n", iproc_model_ctx_invsumweight_ratio(ctx));

    // printf("\nsuminvwt");
    // printf("\n--------");
    // printf("\n%.8f\n", sll->suminvwt);

    insertm_update_active_probs(sll, n, ctx);

    // printf("\nactive_probs");
    // printf("\n------------");
    // iproc_svector_printf(active_probs);

    // printf("\nsum_active_probs");
    // printf("\n----------------");
    // iproc_svector_printf(sll->sum_active_probs);

    insertm_update_sum_mean_var_diff(sll, n, ctx);

    /* update number of sends */
    sll->nsend += n;

    /* update sum of weights
     * IMPORTANT: this must get done *after* active_probs gets updated
     */

    // printf("\nsum_obs_var_diff");
    // printf("\n----------------");
    // iproc_svector_printf(sll->sum_obs_var_diff);

    // printf("\nsum_mean_var_diff");
    // printf("\n-----------------");
    // iproc_svector_printf(sll->sum_mean_var_diff);


    iproc_svector_unref(wt);
    iproc_model_ctx_unref(ctx);
}

double
iproc_sloglik_value (iproc_sloglik *sll)
{
    if (!sll)
        return 0.0;

    return sll->value;
}

static void
iproc_vector_acc_sloglik_grad_nocache (iproc_vector  *dst_vector,
                                       double         scale,
                                       iproc_sloglik *sll)
{
    if (!sll)
        return;
    
    double suminvwt = sll->suminvwt;

    // printf("\nisend");
    // printf("\n-----");
    // printf("\n%"PRId64"\n", sll->isend);

    // printf("\nsuminvwt");
    // printf("\n--------");
    // printf("\n%.8f\n", suminvwt);

    /* compute relative difference in active and initial probabilities */
    iproc_vector *probs0 = iproc_model_probs0(sll->model, sll->isend);
    iproc_svector *dp = iproc_svector_new_copy(sll->sum_active_probs);
    int64_t inz, nnz = iproc_svector_nnz(dp);
    for (inz = 0; inz < nnz; inz++) {
        int64_t j = iproc_svector_nz(dp, inz);
        double p0 = iproc_vector_get(probs0, j);
        iproc_svector_inc(dp, j, -suminvwt * p0);
    }

    // printf("\ndp");
    // printf("\n--");
    // iproc_svector_printf(dp);

    iproc_vector *mean0 = iproc_model_mean0(sll->model, sll->isend);

    // printf("\nmean0");
    // printf("\n-----");
    // iproc_vector_printf(mean0);

    /* sum of observed variables */
    iproc_design_sender0_muls(scale, IPROC_TRANS_TRANS, sll->model->design, sll->isend,
                            sll->nrecv, 1.0, dst_vector);
    iproc_vector_sacc(dst_vector, scale, sll->sum_obs_var_diff);

    // printf("\nsum of observed");
    // printf("\n---------------");
    // iproc_vector_printf(dst_vector);

    /* sum of expected variables */
    iproc_vector_acc(dst_vector, -scale * suminvwt, mean0);

    // printf("\nsum of expected (I)");
    // printf("\n-------------------");
    // iproc_vector_printf(dst_vector);

    iproc_design_sender0_muls(-scale, IPROC_TRANS_TRANS,
                            sll->model->design, sll->isend, dp,
                            1.0, dst_vector);

    // printf("\nsum of expected (II)");
    // printf("\n--------------------");
    // iproc_vector_printf(dst_vector);

    iproc_vector_sacc(dst_vector, -scale, sll->sum_mean_var_diff);

    // printf("\nsum of expected (III)");
    // printf("\n---------------------");
    // iproc_vector_printf(dst_vector);

    iproc_svector_unref(dp);
}


iproc_vector *
iproc_sloglik_grad (iproc_sloglik *sll)
{
    if (!sll)
        return NULL;

    if (!sll->grad_cached) {
        iproc_vector_set_all(sll->grad, 0.0);
        iproc_vector_acc_sloglik_grad_nocache(sll->grad, 1.0, sll);
        sll->grad_cached = 1;
    }
    return sll->grad;
}

