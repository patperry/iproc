#include "port.h"
#include <assert.h>
#include <string.h>

#include "frame.h"
#include "vars.h"


static void isendtot_message_add(void *udata, struct frame *f,
				 const struct message *msg)
{
	const struct tvar *tv = udata;
	const struct var *v = &tv->var;
	struct design *d = v->design;
	const double one = 1.0;	
	const size_t izero = 0;	

	size_t isend = msg->from;
	const double *dx = design_tvar(d, v, isend);
	if (!dx || *dx == 0.0) {
		design_update(d, v, isend, &one, &izero, 1);
	}
}


static struct frame_callbacks isendtot_frame_callbacks = {
	isendtot_message_add,
	NULL,			// message_advance,
	NULL
};


static void isendtot_init(struct tvar *tv, struct design *d, va_list ap)
{
	(void)ap; // unused
	
	struct frame *f = design_frame(d);
	
	tv->var.dim = 1;
	tv->udata = NULL;

	frame_add_observer(f, tv, &isendtot_frame_callbacks);
}


static void isendtot_deinit(struct tvar *tv, struct design *d)
{
	struct frame *f = design_frame(d);
	frame_remove_observer(f, tv);
}


static struct tvar_type VAR_ISENDTOT_REP = {
	isendtot_init,
	isendtot_deinit
};


const struct tvar_type *VAR_ISENDTOT = &VAR_ISENDTOT_REP;
