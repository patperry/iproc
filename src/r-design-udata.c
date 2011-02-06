
#include <Rdefines.h>
#include "memory.h"
#include "r-design.h"

Riproc_design_udata *
Riproc_design_udata_new (SEXP Rreceive_intervals)
{
    Riproc_design_udata *udata = iproc_malloc(sizeof(*udata));
    int64_t n;
    double *intvls;

    if (Rreceive_intervals == NULL_USER_OBJECT) {
        n = 0;
        intvls = NULL;
    } else {
        n = GET_LENGTH(Rreceive_intervals);
        intvls = NUMERIC_POINTER(Rreceive_intervals);
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
Riproc_design_udata_get_sender_design (iproc_design_ctx *ctx)
{
    Riproc_design_udata *udata = ctx->design->user_data;
    iproc_history *history = ctx->history;
    int64_t isend = ctx->isend;
    iproc_array *dst = ctx->sender_design;
    int64_t offset = 0;
    int64_t parent_dim = Riproc_design_udata_dim(udata);

    iproc_vrecip_get(udata->recip, history, isend, dst, offset, parent_dim);
}
