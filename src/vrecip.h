#ifndef _IPROC_VRECIP
#define _IPROC_VRECIP

#include <stdint.h>
#include "darray.h"
#include "design.h"
#include "refcount.h"

typedef struct _iproc_vrecip iproc_vrecip;

struct _iproc_vrecip {
    iproc_design_var var;
    struct darray   *intvls;
    iproc_refcount   refcount;
};

iproc_vrecip * iproc_vrecip_new   (double        *intvls,
                                   int64_t        n);
iproc_vrecip * iproc_vrecip_ref   (iproc_vrecip *v);
void           iproc_vrecip_unref (iproc_vrecip *v);


#endif /* _IPROC_VRECIP */
