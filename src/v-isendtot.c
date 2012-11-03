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


static void isendtot_init(struct tvar *tv, const char *name, struct history *h, va_list ap)
{
	(void)ap; // unused

	var_meta_init(&tv->var.meta, name, VAR_TYPE_TVAR, NULL, 0);
	tv->udata = NULL;

	history_add_observer(h, tv, &isendtot_history_callbacks);
}


static void isendtot_deinit(struct tvar *tv, struct history *h)
{
	history_remove_observer(h, tv);
	var_meta_deinit(&tv->var.meta);
}


static struct tvar_type VAR_ISENDTOT_REP = {
	isendtot_init,
	isendtot_deinit
};


const struct tvar_type *VAR_ISENDTOT = &VAR_ISENDTOT_REP;
