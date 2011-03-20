#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <assert.h>
#include <stddef.h>
#include "memory.h"
#include "design.h"
#include "utils.h"
#include "vrecip.h"


static int
compare_double (void *px,
                void *py)
{
    double x = *((double *)px);
    double y = *((double *)py);
    
    if (x < y) {
        return -1;
    } else if (x > y) {
        return +1;
    } else {
        return 0;
    }
}

static void
iproc_vrecip_free (iproc_vrecip *v)
{
    if (v) {
        iproc_array_unref(v->intvls);
        iproc_free(v);
    }
}

static void
design_var_free (iproc_design_var *var)
{
    iproc_vrecip *v = container_of(var, iproc_vrecip, var);
    iproc_vrecip_unref(v);
}

static void
design_var_get_dxs (iproc_design_var *var,
                    iproc_design_ctx *ctx,
                    int64_t           offset)
{
    iproc_vrecip *v = container_of(var, iproc_vrecip, var);
    iproc_array *intvls = v->intvls;
    int64_t nintvl = iproc_array_size(intvls);
    iproc_history *history = ctx->history;
    int64_t isend = ctx->isend;

    if (!history || nintvl == 0)
        return;
    
    double tcur = iproc_history_tcur(history);
    iproc_trace *trace = iproc_history_recv(history, isend);
    
    int64_t i, n = iproc_trace_size(trace);
    for (i = 0; i < n; i++) {
        iproc_events *events = iproc_trace_get(trace, i);
        iproc_event_meta *meta = iproc_events_last(events);
        
        if (!meta)
            continue;
        
        int64_t jsend = iproc_events_id(events);
        double t = meta->time;
        double dt = tcur - t;
        int64_t pos = iproc_array_bsearch(intvls, &dt, compare_double);
        
        if (pos < 0)
            pos = ~pos;

        if (pos < nintvl) {
            /* (jsend, [(pos, +1.0)]) */
            iproc_svector *dx = iproc_design_ctx_dx(ctx, jsend, false);
            iproc_svector_inc(dx, offset + pos, 1.0);
        }
    }
}

iproc_vrecip *
iproc_vrecip_new (double       *intvls,
                  int64_t       n)
{
    assert(n >= 0);
    assert(n == 0 || intvls);

    iproc_vrecip *v = iproc_malloc(sizeof(*v));

    if (!v)
        return NULL;
    
    iproc_design_var_init(&v->var, n, design_var_get_dxs, design_var_free);
    v->intvls = iproc_array_new(sizeof(double));
    iproc_refcount_init(&v->refcount);

    if (!v->intvls) {
        iproc_vrecip_free(v);
        v = NULL;
    } else {
        int64_t i;
        for (i = 0; i < n; i++) {
            assert(intvls[i] > 0.0);
            assert(i == 0 || intvls[i] > intvls[i-1]);

            iproc_array_append(v->intvls, intvls + i);
        }
    }

    return v;
}

iproc_vrecip *
iproc_vrecip_ref (iproc_vrecip *v)
{
    if (v) {
        iproc_refcount_get(&v->refcount);
    }

    return v;
}

static void
iproc_vrecip_release (iproc_refcount *refcount)
{
    iproc_vrecip *v = container_of(refcount, iproc_vrecip, refcount);
    iproc_vrecip_free(v);
}

void
iproc_vrecip_unref (iproc_vrecip *v)
{
    if (v) {
        iproc_refcount_put(&v->refcount, iproc_vrecip_release);
    }
}
