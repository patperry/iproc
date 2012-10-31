#include "port.h"
#include <assert.h>
#include <string.h>

#include "frame.h"
#include "vars.h"


static void isendtot_message_add(void *udata, struct history *h,
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


static struct history_callbacks isendtot_history_callbacks = {
	isendtot_message_add,
	NULL,			// message_advance,
	NULL
};


static void isendtot_init(struct tvar *tv, struct history *h, va_list ap)
{
	(void)ap; // unused
	
	tv->var.rank = 0;
	tv->udata = NULL;

	history_add_observer(h, tv, &isendtot_history_callbacks);
}


static void isendtot_deinit(struct tvar *tv, struct history *h)
{
	history_remove_observer(h, tv);
}


static struct tvar_type VAR_ISENDTOT_REP = {
	isendtot_init,
	isendtot_deinit
};


const struct tvar_type *VAR_ISENDTOT = &VAR_ISENDTOT_REP;
