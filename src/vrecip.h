#ifndef _IPROC_VRECIP
#define _IPROC_VRECIP

#include <stdint.h>
#include "array.h"
#include "design.h"
#include "history.h"
#include "refcount.h"

typedef struct _iproc_vrecip iproc_vrecip;

struct _iproc_vrecip {
    iproc_array    *intvls;
    iproc_refcount  refcount;
};

iproc_vrecip * iproc_vrecip_new   (double        *intvls,
                                   int64_t        n);
iproc_vrecip * iproc_vrecip_ref   (iproc_vrecip *v);
void           iproc_vrecip_unref (iproc_vrecip *v);
int64_t        iproc_vrecip_dim   (iproc_vrecip *v);

void           iproc_vrecip_get   (iproc_design  *design,
                                   iproc_vrecip  *v,
                                   iproc_history *history,
                                   int64_t        isend,
                                   iproc_array   *dst,
                                   int64_t        offset);

                                     



#endif /* _IPROC_VRECIP */
