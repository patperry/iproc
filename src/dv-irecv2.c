#include "port.h"
#include <assert.h>

#include "frame.h"
#include "vars.h"

/*
 *  i <---- k <---- j
 */


static void irecv2_message_add(void *udata, struct frame *f,
			       const struct message *msg)
{
	const struct tvar2 *tv = udata;
	const struct var2 *v = &tv->var;
	struct design2 *d = frame_dyad_design(f);
	size_t index = v->index;
	const double one = 1.0;
	const size_t izero = 0;
	const struct history *h = frame_history(f);
	size_t ito, nto = msg->nto;
	size_t hrecv = msg->from;
	size_t iz, nz;
	const size_t *indx;

	history_get_recv_active(h, hrecv, &indx, &nz);
	
	for (iz = 0; iz < nz; iz++) {
		size_t jrecv = indx[iz];
		
		for (ito = 0; ito < nto; ito++) {
			size_t isend = msg->to[ito];
				
			if (hrecv == isend || hrecv == jrecv || isend == jrecv)
				continue;
			
			const double *dx = design2_tvars(d, isend, jrecv);
			
			if (!dx || dx[index] == 0.0) {
				design2_update(d, v, isend, jrecv, &one, &izero, 1);
			}
		}
	}
	
	size_t cojrecv = msg->from;
	
	for (ito = 0; ito < nto; ito++) {
		size_t cohsend = msg->to[ito];
		history_get_send_active(h, cohsend, &indx, &nz);
		
		for (iz = 0; iz < nz; iz++) {
			size_t coisend = indx[iz];
				
			if (cohsend == coisend || cohsend == cojrecv || coisend == cojrecv)
				continue;
			
			const double *dx = design2_tvars(d, coisend, cojrecv);
			
			if (!dx || dx[index] == 0.0) {
				design2_update(d, v, coisend, cojrecv, &one, &izero, 1);
			}
		}
	}

}


static struct frame_callbacks irecv2_frame_callbacks = {
	irecv2_message_add,
	NULL,			// message_advance,
	NULL
};

static void irecv2_init(struct tvar2 *tv, const struct design2 *d, va_list ap)
{
	(void)ap; // unused

	struct frame *f = design2_frame(d);

	tv->var.dim = 1;
	tv->udata = NULL;

	frame_add_observer(f, tv, &irecv2_frame_callbacks);
}

static void irecv2_deinit(struct tvar2 *tv, const struct design2 *d)
{
	struct frame *f = design2_frame(d);
	frame_remove_observer(f, tv);
}


static struct tvar2_type DYAD_VAR_IRECV2_REP = {
	irecv2_init,
	irecv2_deinit
};


const struct tvar2_type *DYAD_VAR_IRECV2 = &DYAD_VAR_IRECV2_REP;
