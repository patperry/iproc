#include "port.h"
#include <assert.h>
#include "design2.h"
#include "var.h"


/*
 *   /----> i
 *  k
 *   \----> j
 */

static void isib_message_add(void *udata, struct history *h,
			     const struct message *msg)
{
	const struct tvar2 *tv = udata;
	const struct var2 *v = &tv->var;
	struct design2 *d = v->design;
	const double one = 1.0;
	const size_t izero = 0;
	size_t hsend = msg->from;
	size_t ito, nto = msg->nto;
	const size_t *indx;
	size_t iz, nz;

	history_get_send_active(h, hsend, &indx, &nz);
	
	for (iz = 0; iz < nz; iz++) {
		size_t jrecv = indx[iz];
		size_t coisend = jrecv;
		
		for (ito = 0; ito < nto; ito++) {
			size_t isend = msg->to[ito];
				
			if (hsend == isend || hsend == jrecv || isend == jrecv)
				continue;
				
			const double *dx = design2_tvar(d, v, isend, jrecv);
			if (!dx || *dx == 0.0) {
				design2_update(d, v, isend, jrecv, &one, &izero, 1);
			}
		}
			
		for (ito = 0; ito < nto; ito++) {
			size_t cojrecv = msg->to[ito];
				
			if (hsend == coisend || hsend == cojrecv || coisend == cojrecv)
				continue;
			
			const double *dx = design2_tvar(d, v, coisend, cojrecv);
			if (!dx || *dx == 0.0) {
				design2_update(d, v, coisend, cojrecv, &one, &izero, 1);
			}
		}
	}
}


static struct history_callbacks isib_history_callbacks = {
	isib_message_add,
	NULL,			// message_advance,
	NULL
};

static void isib_init(struct tvar2 *tv, const char *name, struct history *h, va_list ap)
{
	(void)ap; // unused

	var_meta_init(&tv->var.meta, name, VAR_TYPE_TVAR, NULL, 0);
	tv->udata = NULL;

	history_add_observer(h, tv, &isib_history_callbacks);
}

static void isib_deinit(struct tvar2 *tv, struct history *h)
{
	history_remove_observer(h, tv);
	var_meta_deinit(&tv->var.meta);
}


static struct tvar2_type VAR2_ISIB_REP = {
	isib_init,
	isib_deinit
};


const struct tvar2_type *VAR2_ISIB = &VAR2_ISIB_REP;

