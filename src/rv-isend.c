#include "port.h"
#include <assert.h>

#include "frame.h"
#include "vars.h"

struct isend_active {
	struct intset jrecv;
};

struct isend_udata {
	struct array active;
};

static void isend_init(struct design_var *dv, const struct design *d,
		       void *params)
{
	assert(dv);
	assert(d);
	assert(!params);

	dv->dim = 1;
}

static void isend_deinit(struct design_var *dv)
{
	assert(dv);
}

static void isend_frame_init(struct frame_var *fv, struct frame *f)
{
	assert(fv);
	assert(f);

	const struct design *d = frame_design(f);
	struct isend_udata *udata = xcalloc(1, sizeof(*udata));
	array_init(&udata->active, sizeof(struct isend_active));

	ssize_t isend, nsend = design_send_count(d);
	for (isend = 0; isend < nsend; isend++) {
		struct isend_active *active = array_add(&udata->active, NULL);
		intset_init(&active->jrecv);
	}

	fv->udata = udata;
}

static void isend_frame_deinit(struct frame_var *fv)
{
	assert(fv);
	assert(fv->udata);

	struct isend_udata *udata = fv->udata;
	struct isend_active *active;

	ARRAY_FOREACH(active, &udata->active) {
		intset_deinit(&active->jrecv);
	}

	array_deinit(&udata->active);
	xfree(udata);
}

static void handle_clear(void *udata, struct frame *f)
{
	struct frame_var *fv = udata;
	struct isend_udata *fv_udata = fv->udata;

	assert(fv);
	assert(fv->udata);

	struct isend_active *active;

	ARRAY_FOREACH(active, &fv_udata->active) {
		intset_clear(&active->jrecv);
	}
}

static void handle_message_add(void *udata, struct frame *f,
			       const struct message *msg)
{
	struct frame_var *fv = udata;
	struct isend_udata *fv_udata = fv->udata;

	assert(fv);
	assert(f);
	assert(msg);
	assert(fv->design);
	assert(fv->design->dyn_index >= 0);
	assert(fv->design->dyn_index + fv->design->dim
	       <= design_recv_dyn_dim(f->design));

	ssize_t isend = msg->from;
	ssize_t dyn_index = fv->design->dyn_index;

	ssize_t ito, nto = msg->nto;
	for (ito = 0; ito < nto; ito++) {
		ssize_t jrecv = msg->to[ito];

		struct isend_active *active =
		    array_item(&fv_udata->active, isend);
		if (intset_add(&active->jrecv, jrecv)) {
			frame_recv_update(f, isend, jrecv, dyn_index, 1.0);
		}
	}
}

static struct var_type RECV_VAR_ISEND_REP = {
	VAR_RECV_VAR,
	isend_init,
	isend_deinit,
	isend_frame_init,
	isend_frame_deinit,
	{
	 handle_message_add,
	 NULL,			// message_advance,
	 NULL,			// recv_update
	 NULL,			// send_update
	 handle_clear}
};

const struct var_type *RECV_VAR_ISEND = &RECV_VAR_ISEND_REP;
