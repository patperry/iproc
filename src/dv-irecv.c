#include "port.h"
#include <assert.h>
#include <string.h>

#include "frame.h"
#include "vars.h"


static void irecv_message_add(void *udata, struct frame *f,
			      const struct message *msg)
{
	const struct tvar *tv = udata;
	const struct var *v = &tv->var;
	struct design *d = frame_dyad_design(f);	
	size_t index = v->index;
	double one = 1.0;	
	size_t jrecv = msg->from;
	
	size_t ito, nto = msg->nto;
	for (ito = 0; ito < nto; ito++) {
		if (msg->from == msg->to[ito])
			continue;
		
		size_t isend = msg->to[ito];
		size_t ix = frame_dyad_ix(f, isend, jrecv);
		const double *dx = design_tvars(d, ix);
		if (!dx || dx[index] == 0.0) {
			design_update(d, v, ix, &one, NULL);
		}
	}
}


static void irecv_clear(void *udata, struct frame *f)
{
	const struct tvar *tv = udata;
	const struct var *v = &tv->var;
	struct design *d = frame_dyad_design(f);	
	design_clear(d, v);
}


static struct frame_callbacks irecv_frame_callbacks = {
	irecv_message_add,
	NULL,			// message_advance,
	NULL,			// dyad_update
	irecv_clear
};


static void irecv_init(struct tvar *tv, const struct design *d, va_list ap)
{
	(void)ap; // unused
	
	struct frame *f = design_frame(d);
	
	tv->var.dim = 1;
	tv->udata = NULL;

	frame_add_observer(f, tv, &irecv_frame_callbacks);
}


static void irecv_deinit(struct tvar *tv, const struct design *d)
{
	struct frame *f = design_frame(d);
	frame_remove_observer(f, tv);
}


static struct tvar_type DYAD_VAR_IRECV_REP = {
	irecv_init,
	irecv_deinit
};


const struct tvar_type *DYAD_VAR_IRECV = &DYAD_VAR_IRECV_REP;
