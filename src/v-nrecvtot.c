#include "port.h"
#include <assert.h>
#include <string.h>
#include "frame.h"
#include "vars.h"

static void nrecvtot_message_add(void *udata, struct frame *f,
			      const struct message *msg)
{
	if (!frame_interval_count(f))
		return;

	const struct tvar *tv = udata;
	const struct var *v = &tv->var;
	struct design *d = v->design;

	size_t jrecv = msg->from;
	double dx_data[1] = { +1.0 };
	size_t dx_index[1] = { 0 };
	size_t nz = 1;
	
	size_t ito, nto = msg->nto;
	for (ito = 0; ito < nto; ito++) {
		design_update(d, v, jrecv, dx_data, dx_index, nz);
	}
}


static void nrecvtot_message_advance(void *udata, struct frame *f,
				  const struct message *msg, size_t intvl)
{
	const struct tvar *tv = udata;
	const struct var *v = &tv->var;
	struct design *d = v->design;
	size_t jrecv = msg->from;

	double dx_data[2] = { -1.0, +1.0 };
	size_t dx_index[2];
	size_t nz = 2;
	
	size_t ito, nto = msg->nto;
	for (ito = 0; ito < nto; ito++) {
		size_t ix1 = intvl;
		size_t ix0 = ix1 - 1;
		
		dx_index[0] = ix0;
		dx_index[1] = ix1;
		
		design_update(d, v, jrecv, dx_data, dx_index, nz);
	}
}


static struct frame_callbacks nrecvtot_frame_callbacks = {
	nrecvtot_message_add,
	nrecvtot_message_advance,
	NULL
};


static void nrecvtot_init(struct tvar *tv, struct design *d, va_list ap)
{
	(void)ap;		// unused;

	struct frame *f = design_frame(d);
	size_t n = frame_interval_count(f);

	tv->var.dim = n;
	tv->udata = NULL;
	
	frame_add_observer(f, tv, &nrecvtot_frame_callbacks);
}


static void nrecvtot_deinit(struct tvar *tv, struct design *d)
{
	struct frame *f = design_frame(d);
	frame_remove_observer(f, tv);
}


static struct tvar_type VAR_NRECVTOT_REP = {
	nrecvtot_init,
	nrecvtot_deinit,
	
};


const struct tvar_type *VAR_NRECVTOT = &VAR_NRECVTOT_REP;
