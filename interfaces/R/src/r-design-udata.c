
#include <Rdefines.h>
#include "memory.h"
#include "r-design.h"

Riproc_design_udata *
Riproc_design_udata_new (SEXP Rrecip_intervals)
{
    Riproc_design_udata *udata = iproc_malloc(sizeof(*udata));
    int64_t n;
    double *intvls;

    if (Rrecip_intervals == NULL_USER_OBJECT) {
        n = 0;
        intvls = NULL;
    } else {
        n = GET_LENGTH(Rrecip_intervals);
        intvls = NUMERIC_POINTER(Rrecip_intervals);
    }

    udata->recip = iproc_vrecip_new(intvls, n);
    return udata;
}

void
Riproc_design_udata_free (void *design_udata)
{
    Riproc_design_udata *udata = design_udata;
    if (udata) {
        iproc_vrecip_unref(udata->recip);
        iproc_free(udata);
    }
}

int64_t
Riproc_design_udata_dim (Riproc_design_udata *udata)
{
    if (!udata)
        return 0;

    return iproc_vrecip_dim(udata->recip);
}

void
Riproc_design_udata_get_sdesign_vars (iproc_design_ctx *ctx)
{
    if (iproc_design_ndynamic(ctx->design) == 0)
        return;

    Riproc_design_udata *udata = ctx->design->user_data;
    iproc_history *history = ctx->history;
    int64_t isend = ctx->isend;
    iproc_array *dst = ctx->dxs;
    int64_t offset = iproc_design_idynamic(ctx->design, 0);
    iproc_vrecip_get(ctx->design, udata->recip, history, isend, dst, offset);
}
