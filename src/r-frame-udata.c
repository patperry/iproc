
#include <Rdefines.h>
#include "memory.h"
#include "r-frame.h"

Riproc_frame_udata *
Riproc_frame_udata_new (SEXP Rreceive_intervals)
{
    Riproc_frame_udata *udata = iproc_malloc(sizeof(*udata));
    int64_t i, n;
    int64_t *intvls;

    if (Rreceive_intervals == NULL_USER_OBJECT) {
        n = 0;
        intvls = NULL;
    } else {
        int *intvls_src = INTEGER_POINTER(Rreceive_intervals);
        n = GET_LENGTH(Rreceive_intervals);

        intvls = (int64_t *)R_alloc(n, sizeof(int64_t)); /* reclaimed on exit */
        for (i = 0; i < n; i++) {
            intvls[i] = intvls_src[i];
        }
    }
    
    udata->recip = iproc_vrecip_new(intvls, n);
    return udata;
}

void
Riproc_frame_udata_free (void *frame_udata)
{
    Riproc_frame_udata *udata = frame_udata;
    if (udata) {
        iproc_vrecip_unref(udata->recip);
        iproc_free(udata);
    }
}

int64_t
Riproc_frame_udata_dim (Riproc_frame_udata *udata)
{
    if (!udata)
        return 0;

    return iproc_vrecip_dim(udata->recip);
}

void
Riproc_frame_udata_get_sender_frame (iproc_frame_ctx *ctx)
{
    Riproc_frame_udata *udata = ctx->frame->user_data;
    iproc_history *history = ctx->history;
    int64_t isend = ctx->isend;
    iproc_array *dst = ctx->sender_frame;
    int64_t offset = 0;
    int64_t parent_dim = Riproc_frame_udata_dim(udata);

    iproc_vrecip_get(udata->recip, history, isend, dst, offset, parent_dim);
}
