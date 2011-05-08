#ifndef _IPROC_VRECIP
#define _IPROC_VRECIP

#include <stdint.h>
#include "darray.h"
#include "design.h"
#include "hashset.h"
#include "refcount.h"

struct vrecv_active {
	struct dyad dyad;
	ssize_t id;
	ssize_t intvl;
};

struct vrecv {
	struct dyad_var dyad_var;
	struct hashset active;
};

bool vrecv_init(struct vrecv *v, const struct design *d);
void vrecv_deinit(struct vrecv *v);



/* DEPRECATED */
typedef struct _iproc_vrecip iproc_vrecip;

struct _iproc_vrecip {
	iproc_design_var var;
	struct darray intvls;
	struct refcount refcount;
};

iproc_vrecip *iproc_vrecip_new(double *intvls, ssize_t n);
iproc_vrecip *iproc_vrecip_ref(iproc_vrecip * v);
void iproc_vrecip_unref(iproc_vrecip * v);

#endif /* _IPROC_VRECIP */
