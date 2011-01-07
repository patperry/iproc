#ifndef _IPROC_MODEL_H
#define _IPROC_MODEL_H

#include <iproc/vars.h>
#include <iproc/vector.h>

typedef struct _iproc_model iproc_model;

struct _iproc_model {
    iproc_vars   *vars;
    iproc_vector *coef;
    int           has_loops;
    int           refcount;
};

iproc_model * iproc_model_new              (iproc_vars     *vars,
                                            iproc_vector   *coef,
                                            int             has_loops);
iproc_model * iproc_model_ref              (iproc_model    *model);
void          iproc_model_unref            (iproc_model    *model);


double        iproc_model_logprob0         (iproc_model    *model,
                                            int64_t         j);
double        iproc_model_prob0            (iproc_model    *model,
                                            int64_t         j);

void          iproc_model_get_probs        (iproc_model    *model,
                                            iproc_vars_ctx *ctx,
                                            iproc_vector   *probs);
void          iproc_model_get_logprobs     (iproc_model    *model,
                                            iproc_vars_ctx *ctx,
                                            iproc_vector   *logprobs);
void          iproc_model_get_new_logprobs (iproc_model    *model,
                                            iproc_vars_ctx *ctx,
                                            double         *plogprob0_shift,
                                            iproc_svector  *logprobs);


#endif /* _IPROC_MODEL_H */
