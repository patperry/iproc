#include "port.h"
#include <assert.h>
#include <string.h>

#include "frame.h"
#include "vars.h"


static void irecvtot_message_add(void *udata, struct frame *f,
			      const struct message *msg)
{
	const struct tvar *tv = udata;
	const struct var *v = &tv->var;
	struct design *d = v->design;
	const double one = 1.0;	
	const size_t izero = 0;	

	size_t ito, nto = msg->nto;
	for (ito = 0; ito < nto; ito++) {
		size_t jrecv = msg->to[ito];
		const double *dx = design_tvar(d, v, jrecv);
		if (!dx || *dx == 0.0) {
			design_update(d, v, jrecv, &one, &izero, 1);
		}
	}
}


static struct frame_callbacks irecvtot_frame_callbacks = {
	irecvtot_message_add,
	NULL,			// message_advance,
	NULL
};


static void irecvtot_init(struct tvar *tv, struct design *d, va_list ap)
{
	(void)ap; // unused
	
	struct frame *f = design_frame(d);
	
	tv->var.dim = 1;
	tv->udata = NULL;

	frame_add_observer(f, tv, &irecvtot_frame_callbacks);
}


static void irecvtot_deinit(struct tvar *tv, struct design *d)
{
	struct frame *f = design_frame(d);
	frame_remove_observer(f, tv);
}


static struct tvar_type VAR_IRECVTOT_REP = {
	irecvtot_init,
	irecvtot_deinit
};


const struct tvar_type *VAR_IRECVTOT = &VAR_IRECVTOT_REP;
