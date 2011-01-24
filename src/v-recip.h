#ifndef _IPROC_V_RECIP
#define _IPROC_V_RECIP

#include <stdint.h>
#include "array.h"
#include "history.h"
#include "refcount.h"

typedef struct _iproc_v_recip iproc_v_recip;

struct _iproc_v_recip {
    iproc_array    *intvls;
    iproc_refcount  refcount;
};

iproc_v_recip * iproc_v_recip_new   (int64_t       *intvls,
                                     int64_t        n);
iproc_v_recip * iproc_v_recip_ref   (iproc_v_recip *v);
void            iproc_v_recip_unref (iproc_v_recip *v);
int64_t         iproc_v_recip_dim   (iproc_v_recip *v);

void            iproc_v_recip_get   (iproc_v_recip *v,
                                     iproc_history *history,
                                     int64_t        isend,
                                     iproc_array   *dst,
                                     int64_t        offset,
                                     int64_t        parent_dim);
                                     



#endif /* _IPROC_V_RECIP */
