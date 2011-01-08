#ifndef _IPROC_MODEL_H
#define _IPROC_MODEL_H

#include <iproc/refcount.h>
#include <iproc/vars.h>
#include <iproc/vector.h>


typedef struct _iproc_model iproc_model;

struct _iproc_model {
    iproc_vars    *vars;
    iproc_vector  *coef;
    int            has_loops;
    iproc_refcount refcount;
};

iproc_model * iproc_model_new              (iproc_vars     *vars,
                                            iproc_vector   *coef,
                                            int             has_loops);
iproc_model * iproc_model_ref              (iproc_model    *model);
void          iproc_model_unref            (iproc_model    *model);


double        iproc_model_logprob0         (iproc_model    *model,
                                            int64_t         i,
                                            int64_t         j);
double        iproc_model_prob0            (iproc_model    *model,
                                            int64_t         i,
                                            int64_t         j);

void          iproc_model_get_probs        (iproc_model    *model,
                                            iproc_vars_ctx *ctx,
                                            iproc_vector   *probs);
void          iproc_model_get_logprobs     (iproc_model    *model,
                                            iproc_vars_ctx *ctx,
                                            iproc_vector   *logprobs);

/* The output of this function (logprob0_shift and logprobs) can be used
 * to compute the probabilites for the given context.  To compute log(p[j])
 * do the following:
 *
 *   1. If j == i and no self-loops are allowed, return -INFINITY
 *   2. If j is a sparse index of logprobs, return logprobs[j]
 *   3. Otherwise return iproc_model_logprob0(i,j) + logprob0_shift.
 *
 * Note the special handling for self loops.  An invalid self-loop may or may
 * not be have an index in 'logprobs'.  We only gaurantee that if it *does*
 * have an index in logprobs, then the corresponding value will be -INFINITY
 * (thus, steps 1 and 2 commute).
 */
void          iproc_model_get_new_logprobs (iproc_model    *model,
                                            iproc_vars_ctx *ctx,
                                            double         *plogprob0_shift,
                                            iproc_svector  *logprobs);


#endif /* _IPROC_MODEL_H */
