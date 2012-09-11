#include "port.h"
#include <assert.h>

#include "frame.h"
#include "vars.h"

/*
 *  i ----> h ----> j
 */


static void isend2_message_add(void *udata, struct frame *f,
			       const struct message *msg)
{
	const struct tvar2 *tv = udata;
	const struct var2 *v = &tv->var;
	struct design2 *d = frame_dyad_design(f);
	const double one = 1.0;
	const size_t izero = 0;
	const struct history *h = frame_history(f);
	size_t ito, nto = msg->nto;

	const size_t *indx;
	size_t iz, nz;
	size_t isend = msg->from;
	
	for (ito = 0; ito < nto; ito++) {
		size_t hsend = msg->to[ito];
		history_get_send_active(h, hsend, &indx, &nz);
		
		for (iz = 0; iz < nz; iz++) {
			size_t jrecv = indx[iz];
				
			if (hsend == isend || hsend == jrecv || isend == jrecv)
				continue;
			
			const double *dx = design2_tvar(d, v, isend, jrecv);
			
			if (!dx || *dx == 0.0) {
				design2_update(d, v, isend, jrecv, &one, &izero, 1);
			}
		}
	}
	
	size_t cohrecv = isend;
	history_get_recv_active(h, cohrecv, &indx, &nz);
	
	for (iz = 0; iz < nz; iz++) {
		size_t coisend = indx[iz];
		
		for (ito = 0; ito < nto; ito++) {
			size_t cojrecv = msg->to[ito];
				
			if (cohrecv == coisend || cohrecv == cojrecv || coisend == cojrecv)
				continue;
			
			const double *dx = design2_tvar(d, v, coisend, cojrecv);
			
			if (!dx || *dx == 0.0) {
				design2_update(d, v, coisend, cojrecv, &one, &izero, 1);
			}
		}
	}
}


static struct frame_callbacks isend2_frame_callbacks = {
	isend2_message_add,
	NULL,			// message_advance,
	NULL
};

static void isend2_init(struct tvar2 *tv, struct design2 *d, va_list ap)
{
	(void)ap; // unused

	struct frame *f = design2_frame(d);

	tv->var.dim = 1;
	tv->udata = NULL;

	frame_add_observer(f, tv, &isend2_frame_callbacks);
}

static void isend2_deinit(struct tvar2 *tv, struct design2 *d)
{
	struct frame *f = design2_frame(d);
	frame_remove_observer(f, tv);
}


static struct tvar2_type DYAD_VAR_ISEND2_REP = {
	isend2_init,
	isend2_deinit
};


const struct tvar2_type *DYAD_VAR_ISEND2 = &DYAD_VAR_ISEND2_REP;


