#include "port.h"
#include <assert.h>
#include <string.h>

#include "frame.h"
#include "vars.h"


static void isend_message_add(void *udata, struct frame *f,
			      const struct message *msg)
{
	const struct tvar *tv = udata;
	const struct var *v = &tv->var;
	struct design *d = frame_dyad_design(f);	
	size_t index = v->index;
	double one = 1.0;	
	size_t isend = msg->from;
		
	size_t ito, nto = msg->nto;
	for (ito = 0; ito < nto; ito++) {
		if (msg->to[ito] == msg->from)
			continue;
		
		size_t jrecv = msg->to[ito];
		size_t ix = frame_dyad_ix(f, isend, jrecv);
		const double *dx = design_tvars(d, ix);
		if (!dx || dx[index] == 0.0) {
			design_update(d, v, ix, &one, NULL);
		}
	}
}


static struct frame_callbacks isend_frame_callbacks = {
	isend_message_add,
	NULL,			// message_advance,
	NULL
};


static void isend_init(struct tvar *tv, const struct design *d, va_list ap)
{
	(void)ap;		// unused;

	struct frame *f = design_frame(d);
	
	tv->var.dim = 1;
	tv->udata = NULL;

	frame_add_observer(f, tv, &isend_frame_callbacks);
}


static void isend_deinit(struct tvar *tv, const struct design *d)
{
	struct frame *f = design_frame(d);
	frame_remove_observer(f, tv);
}


static struct tvar_type DYAD_VAR_ISEND_REP = {
	isend_init,
	isend_deinit
};


const struct tvar_type *DYAD_VAR_ISEND = &DYAD_VAR_ISEND_REP;
