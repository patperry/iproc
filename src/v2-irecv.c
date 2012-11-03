#include "port.h"
#include <assert.h>
#include <string.h>

#include "frame.h"
#include "vars.h"


static void irecv_message_add(void *udata, struct history *h,
			      const struct message *msg)
{
	const struct tvar2 *tv = udata;
	const struct var2 *v = &tv->var;
	struct design2 *d = v->design;
	const double one = 1.0;	
	const size_t izero = 0;	
	size_t jrecv = msg->from;
	
	size_t ito, nto = msg->nto;
	for (ito = 0; ito < nto; ito++) {
		size_t isend = msg->to[ito];
		const double *dx = design2_tvar(d, v, isend, jrecv);
		if (!dx || *dx == 0.0) {
			design2_update(d, v, isend, jrecv, &one, &izero, 1);
		}
	}
}


static struct history_callbacks irecv_history_callbacks = {
	irecv_message_add,
	NULL,			// message_advance,
	NULL
};


static void irecv_init(struct tvar2 *tv, const char *name, struct history *h, va_list ap)
{
	(void)ap; // unused

	var_meta_init(&tv->var.meta, name, VAR_TYPE_TVAR, NULL, 0);
	tv->udata = NULL;

	history_add_observer(h, tv, &irecv_history_callbacks);
}


static void irecv_deinit(struct tvar2 *tv, struct history *h)
{
	history_remove_observer(h, tv);
	var_meta_deinit(&tv->var.meta);
}


static struct tvar2_type VAR2_IRECV_REP = {
	irecv_init,
	irecv_deinit
};


const struct tvar2_type *VAR2_IRECV = &VAR2_IRECV_REP;
