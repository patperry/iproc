#include "port.h"
#include <assert.h>

#include "frame.h"
#include "vnrecv.h"

static bool vnrecv_init(struct design_var *dv, const struct design *d)
{
	assert(dv);
	assert(d);

	ssize_t n = vector_dim(design_intervals(d));
	dv->dim = n + 1;
	return true;
}

static void vnrecv_deinit(struct design_var *dv)
{
	assert(dv);
}

static bool vnrecv_handle_dyad(struct frame_var *fv, const struct dyad_event *e,
			       struct frame *f)
{
	assert(fv);
	assert(e);
	assert(f);
	assert(fv->design);
	assert(fv->design->index >= 0);
	assert(fv->design->index <= design_dim(f->design) - fv->design->dim);

	ssize_t index = fv->design->index;
	struct svector *dx = frame_dx(f, e->dyad.jrecv, e->dyad.isend);
	struct svector_pos pos;
	double *val;

	if (!dx)
		return false;

	if (e->type == DYAD_EVENT_MOVE) {
		val = svector_find(dx, index + e->intvl - 1, &pos);
		assert(val);
		*val -= 1.0;
		if (*val == 0.0)
			svector_erase(dx, &pos);
	}

	val = svector_at(dx, index + e->intvl);
	if (!val)
		return false;

	*val += 1.0;
	return true;
}

static struct var_type VAR_TYPE_NRECV_REP = {
	DYAD_EVENT_INIT | DYAD_EVENT_MOVE,
	vnrecv_init,
	vnrecv_deinit,
	NULL,			// frame_init,
	NULL,			// frame_deinit,
	NULL,			// frame_clear,
	vnrecv_handle_dyad
};

const struct var_type *VAR_TYPE_NRECV = &VAR_TYPE_NRECV_REP;
