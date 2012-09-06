#include "port.h"
#include <assert.h>
#include <string.h>

#include "frame.h"
#include "vars.h"


static void irecv_message_add(void *udata, struct frame *f,
			      const struct message *msg)
{
	const struct tvar2 *tv = udata;
	const struct var2 *v = &tv->var;
	struct design2 *d = frame_dyad_design(f);	
	const double one = 1.0;	
	const size_t izero = 0;	
	size_t jrecv = msg->from;
	
	size_t ito, nto = msg->nto;
	for (ito = 0; ito < nto; ito++) {
		if (msg->from == msg->to[ito])
			continue;
		
		size_t isend = msg->to[ito];
		const double *dx = design2_tvar(d, v, isend, jrecv);
		if (!dx || *dx == 0.0) {
			design2_update(d, v, isend, jrecv, &one, &izero, 1);
		}
	}
}


static struct frame_callbacks irecv_frame_callbacks = {
	irecv_message_add,
	NULL,			// message_advance,
	NULL
};


static void irecv_init(struct tvar2 *tv, const struct design2 *d, va_list ap)
{
	(void)ap; // unused
	
	struct frame *f = design2_frame(d);
	
	tv->var.dim = 1;
	tv->udata = NULL;

	frame_add_observer(f, tv, &irecv_frame_callbacks);
}


static void irecv_deinit(struct tvar2 *tv, const struct design2 *d)
{
	struct frame *f = design2_frame(d);
	frame_remove_observer(f, tv);
}


static struct tvar2_type DYAD_VAR_IRECV_REP = {
	irecv_init,
	irecv_deinit
};


const struct tvar2_type *DYAD_VAR_IRECV = &DYAD_VAR_IRECV_REP;
