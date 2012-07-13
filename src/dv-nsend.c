#include "port.h"
#include <assert.h>
#include <string.h>
#include "frame.h"
#include "vars.h"


static void nsend_message_add(void *udata, struct frame *f,
			      const struct message *msg)
{
	const struct tvar2 *tv = udata;
	const struct var2 *v = &tv->var;
	struct design2 *d = frame_dyad_design(f);
	size_t isend = msg->from;
	
	double dx_data[1] = { +1.0 };
	size_t dx_index[1] = { 0 };
	struct vpattern pat = vpattern_make(dx_index, 1);
	
	size_t ito, nto = msg->nto;
	for (ito = 0; ito < nto; ito++) {
		if (msg->to[ito] == msg->from)
			continue;
		
		size_t jrecv = msg->to[ito];
		
		design2_update(d, v, isend, jrecv, dx_data, &pat);
	}
}


static void nsend_message_advance(void *udata, struct frame *f,
				  const struct message *msg, size_t intvl)
{
	const struct tvar2 *tv = udata;
	const struct var2 *v = &tv->var;
	struct design2 *d = frame_dyad_design(f);	

	size_t isend = msg->from;
		
	double dx_data[2] = { -1.0, +1.0 };
	size_t dx_index[2] = { 0, 1 }; // values are unused
	struct vpattern pat = vpattern_make(dx_index, 2);
	
	size_t ito, nto = msg->nto;
	for (ito = 0; ito < nto; ito++) {
		if (msg->to[ito] == msg->from)
			continue;
		
		size_t jrecv = msg->to[ito];
		
		size_t ix1 = intvl;
		size_t ix0 = ix1 - 1;
		
		dx_index[0] = ix0;
		dx_index[1] = ix1;
		
		design2_update(d, v, isend, jrecv, dx_data, &pat);
	}
}


struct frame_callbacks nsend_frame_callbacks = {
	nsend_message_add,
	nsend_message_advance,
	NULL
};


static void nsend_init(struct tvar2 *tv, const struct design2 *d, va_list ap)
{
	(void)ap;		// unused;

	struct frame *f = design2_frame(d);
	size_t n = frame_interval_count(f);

	tv->var.dim = n + 1;
	tv->udata = NULL;
	
	frame_add_observer(f, tv, &nsend_frame_callbacks);
}

static void nsend_deinit(struct tvar2 *tv, const struct design2 *d)
{
	struct frame *f = design2_frame(d);
	frame_remove_observer(f, tv);
}


static struct tvar2_type DYAD_VAR_NSEND_REP = {
	nsend_init,
	nsend_deinit,
	};

const struct tvar2_type *DYAD_VAR_NSEND = &DYAD_VAR_NSEND_REP;
