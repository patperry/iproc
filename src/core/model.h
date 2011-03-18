#ifndef _IPROC_MODEL_H
#define _IPROC_MODEL_H

#include <stdbool.h>
#include "refcount.h"
#include "array.h"
#include "design.h"
#include "vector.h"


typedef struct _iproc_group_model iproc_group_model;
typedef struct _iproc_model       iproc_model;
typedef struct _iproc_model_ctx   iproc_model_ctx;

struct _iproc_group_model {
    double        invsumweight0;
    double        logsumweight0;
    iproc_vector *logprobs0;
    iproc_vector *probs0;
    iproc_vector *mean0;
};

struct _iproc_model {
    iproc_design  *design;
    iproc_vector  *coefs;
    bool           has_loops;
    iproc_array   *group_models;
    iproc_array   *ctxs;
    iproc_refcount refcount;
};

struct _iproc_model_ctx {
    iproc_model       *model;
    iproc_design_ctx  *design_ctx;
    iproc_group_model *group;
    iproc_svector     *active_logprobs;
    iproc_svector     *active_probs;
    double             invsumweight_ratio;
    double             logsumweight_diff;
    iproc_refcount     refcount;
};

iproc_model *       iproc_model_new                (iproc_design *design,
                                                    iproc_vector *coefs,
                                                    bool          has_loops);
iproc_model *       iproc_model_ref                (iproc_model  *model);
void                iproc_model_unref              (iproc_model  *model);

iproc_design *      iproc_model_design             (iproc_model *model);
iproc_vector *      iproc_model_coefs              (iproc_model *model);
bool                iproc_model_has_loops          (iproc_model *model);
int64_t             iproc_model_nsender            (iproc_model *model);
int64_t             iproc_model_nreceiver          (iproc_model *model);
int64_t             iproc_model_dim                (iproc_model *model);

/* Initial probability, and expectations, without adjustment for self-loops. */
iproc_group_model * iproc_model_send_group         (iproc_model *model,
                                                    int64_t      isend);
double              iproc_model_invsumweight0      (iproc_model *model,
                                                    int64_t      isend);
double              iproc_model_logsumweight0      (iproc_model *model,
                                                    int64_t      isend);
iproc_vector *      iproc_model_logprobs0          (iproc_model *model,
                                                    int64_t      isend);
iproc_vector *      iproc_model_probs0             (iproc_model *model,
                                                    int64_t      isend);
iproc_vector *      iproc_model_mean0              (iproc_model *model,
                                                    int64_t      isend);


iproc_model_ctx *   iproc_model_ctx_new            (iproc_model     *model,
                                                    int64_t          isend,
                                                    iproc_history   *h);
iproc_model_ctx *   iproc_model_ctx_ref            (iproc_model_ctx *ctx);
void                iproc_model_ctx_unref          (iproc_model_ctx *ctx);

void                iproc_model_ctx_set            (iproc_model_ctx *ctx,
                                                    int64_t          isend,
                                                    iproc_history   *h);

int64_t             iproc_model_ctx_nreceiver      (iproc_model_ctx *ctx);
double              iproc_model_ctx_prob           (iproc_model_ctx *ctx,
                                                    int64_t          jrecv);
double              iproc_model_ctx_logprob        (iproc_model_ctx *ctx,
                                                    int64_t          jrecv);
void                iproc_model_ctx_get_probs      (iproc_model_ctx *ctx,
                                                    iproc_vector    *probs);
void                iproc_model_ctx_get_logprobs   (iproc_model_ctx *ctx,
                                                    iproc_vector    *logprobs);
void                iproc_model_ctx_get_mean       (iproc_model_ctx *ctx,
                                                    iproc_vector    *mean);

double              iproc_model_ctx_invsumweight       (iproc_model_ctx *ctx);
double              iproc_model_ctx_logsumweight       (iproc_model_ctx *ctx);
double              iproc_model_ctx_invsumweight_ratio (iproc_model_ctx *ctx);
double              iproc_model_ctx_logsumweight_diff  (iproc_model_ctx *ctx);
iproc_svector *     iproc_model_ctx_active_probs       (iproc_model_ctx *ctx);
iproc_svector *     iproc_model_ctx_active_logprobs    (iproc_model_ctx *ctx);


/*
void                iproc_vector_acc_diffmean (iproc_vector    *dst_vector,
                                               double           scale,
                                               iproc_model_ctx *ctx);
*/

#endif /* _IPROC_MODEL_H */
