
#include <Rdefines.h>
#include "memory.h"
#include "r-vars.h"

Riproc_vars_udata *
Riproc_vars_udata_new (SEXP Rreceive_intervals)
{
    Riproc_vars_udata *udata = iproc_malloc(sizeof(*udata));
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
    
    udata->recip = iproc_v_recip_new(intvls, n);
    return udata;
}

void
Riproc_vars_udata_free (Riproc_vars_udata *udata)
{
    if (udata) {
        iproc_v_recip_unref(udata->recip);
        iproc_free(udata);
    }
}

int64_t
Riproc_vars_udata_dim (Riproc_vars_udata *udata)
{
    if (!udata)
        return 0;

    return iproc_v_recip_dim(udata->recip);
}

void
Riproc_vars_udata_get_sender_vars (iproc_vars_ctx *ctx)
{
    Riproc_vars_udata *udata = ctx->vars->user_data;
    iproc_history *history = ctx->history;
    int64_t isend = ctx->isend;
    iproc_array *dst = ctx->sender_vars;
    int64_t offset = 0;
    int64_t parent_dim = Riproc_vars_udata_dim(udata);

    iproc_v_recip_get (udata->recip, history, isend, dst, offset, parent_dim);
}
