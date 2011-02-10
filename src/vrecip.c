#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <assert.h>
#include <stddef.h>
#include "memory.h"
#include "design.h"
#include "vrecip.h"

static void
iproc_vrecip_free (iproc_vrecip *v)
{
    if (v) {
        iproc_array_unref(v->intvls);
        iproc_free(v);
    }
}

iproc_vrecip *
iproc_vrecip_new (double       *intvls,
                  int64_t       n)
{
    iproc_vrecip *v = iproc_malloc(sizeof(*v));

    if (!v)
        return NULL;
    
    v->intvls = iproc_array_new(sizeof(double));
    iproc_refcount_init(&v->refcount);

    if (!v->intvls) {
        iproc_vrecip_free(v);
        v = NULL;
    } else {
        int64_t i;
        for (i = 0; i < n; i++) {
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

int64_t
iproc_vrecip_dim (iproc_vrecip *v)
{
    if (v)
        return iproc_array_size(v->intvls);

    return 0;
}

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


static int
compare_sdesign_var_jrecv (void *px,
                           void *py)
{
    int64_t x = ((iproc_sdesign_var *)px)->jrecv;
    int64_t y = ((iproc_sdesign_var *)py)->jrecv;

    if (x < y) {
        return -1;
    } else if (x > y) {
        return +1;
    } else {
        return 0;
    }
}

void
iproc_vrecip_get (iproc_design  *design,
                  iproc_vrecip  *v,
                  iproc_history *history,
                  int64_t        isend,
                  iproc_array   *dst,
                  int64_t        offset)
{
    assert(v);
    assert(dst);
    assert(offset >= 0);
    assert(iproc_design_dim(design) >= iproc_vrecip_dim(v) + offset);

    if (!history)
        return;

    iproc_array *intvls = v->intvls;
    int64_t nintvl = iproc_array_size(intvls);

    if (nintvl == 0)
        return;

    iproc_events *events = iproc_history_recv(history, isend);

    int64_t i, n = iproc_events_npast(events);
    for (i = 0; i < n; i++) {
        int64_t jsend = iproc_events_past(events, i);
        double dt     = iproc_events_past_dt(events, i);
        int64_t pos   = iproc_array_bsearch(intvls, &dt, compare_double);

        if (pos < 0)
            pos = ~pos;
        
        /* (jsend, [(pos, +1.0)]) */
        iproc_svector *jdiff;
        if (pos < nintvl) {
            int64_t k = iproc_array_bsearch(dst, &jsend, compare_sdesign_var_jrecv);
            
            if (k < 0) {
                jdiff = iproc_sdesign_var_new(design);
                iproc_sdesign_var new_sv = { jsend, jdiff };
                iproc_array_insert(dst, ~k, &new_sv);
            } else {
                jdiff = (iproc_array_index(dst, iproc_sdesign_var, k)).jdiff;
            }

            iproc_svector_inc(jdiff, offset + pos, 1.0);
        }
    }

}

