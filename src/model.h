#ifndef _IPROC_MODEL_H
#define _IPROC_MODEL_H

#include "refcount.h"
#include "vars.h"
#include "vector.h"


typedef struct _iproc_model     iproc_model;
typedef struct _iproc_model_ctx iproc_model_ctx;

struct _iproc_model {
    iproc_vars    *vars;
    iproc_vector  *coefs;
    int            has_loops;
    iproc_array   *group_logprobs0;
    iproc_refcount refcount;
};

struct _iproc_model_ctx {
    iproc_model    *model;
    int64_t         isend;
    iproc_vector   *logprobs0;
    iproc_svector  *logprobs;
    double          logprob0_shift;
    iproc_refcount  refcount;
};

iproc_model *     iproc_model_new                (iproc_vars   *vars,
                                                  iproc_vector *coefs,
                                                  int           has_loops);
iproc_model *     iproc_model_ref                (iproc_model  *model);
void              iproc_model_unref              (iproc_model  *model);

iproc_vars *      iproc_model_vars               (iproc_model *model);
iproc_vector *    iproc_model_coefs              (iproc_model *model);
int               iproc_model_has_loops          (iproc_model *model);
int64_t           iproc_model_nsender            (iproc_model *model);
int64_t           iproc_model_nreceiver          (iproc_model *model);
int64_t           iproc_model_dim                (iproc_model *model);

/* Initial probability, and expectations, without adjustment for self-loops. */
double            iproc_model_sumweight0         (iproc_model *model,
                                                  int64_t      isend);
iproc_vector *    iproc_model_logprobs0          (iproc_model *model,
                                                  int64_t      isend);
iproc_vector *    iproc_model_probs0             (iproc_model *model,
                                                  int64_t      isend);


iproc_model_ctx * iproc_model_ctx_new            (iproc_model     *model,
                                                  int64_t          isend,
                                                  iproc_history   *h);
iproc_model_ctx * iproc_model_ctx_ref            (iproc_model_ctx *ctx);
void              iproc_model_ctx_unref          (iproc_model_ctx *ctx);

int64_t           iproc_model_ctx_nreceiver      (iproc_model_ctx *ctx);
double            iproc_model_ctx_prob           (iproc_model_ctx *ctx,
                                                  int64_t          jrecv);
double            iproc_model_ctx_logprob        (iproc_model_ctx *ctx,
                                                  int64_t          jrecv);
void              iproc_model_ctx_get_probs      (iproc_model_ctx *ctx,
                                                  iproc_vector    *probs);
void              iproc_model_ctx_get_logprobs   (iproc_model_ctx *ctx,
                                                  iproc_vector    *logprobs);

double            iproc_model_sumweight          (iproc_model_ctx *ctx);
double            iproc_model_sumweight_ratio    (iproc_model_ctx *ctx);
double            iproc_model_sumweight_logratio (iproc_model_ctx *ctx);
iproc_svector *   iproc_model_active_probs       (iproc_model_ctx *ctx);
iproc_svector *   iproc_model_active_logprobs    (iproc_model_ctx *ctx);


/*
void              iproc_vector_acc_diffmean (iproc_vector    *dst_vector,
                                             double           scale,
                                             iproc_model_ctx *ctx);
*/

#endif /* _IPROC_MODEL_H */
