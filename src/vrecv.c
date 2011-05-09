#include "port.h"

#include <assert.h>
#include <stddef.h>
#include <stdlib.h>
#include "compare.h"
#include "ieee754.h"
#include "design.h"
#include "frame.h"
#include "hashset.h"
#include "intmap.h"
#include "util.h"
#include "vrecv.h"

struct vrecv_active {
	struct dyad dyad;
	ssize_t id;
	ssize_t intvl;
};

struct vrecv_udata {
	struct hashset active;
};


static uint32_t vrecv_active_hash(const void *x)
{
	const struct vrecv_active *active = x;
	return memory_hash(&active->dyad, sizeof(active->dyad));
}

static bool vrecv_active_equals(const void *x, const void *y)
{
	const struct vrecv_active *a = x, *b = y;
	return a->dyad.isend == b->dyad.isend && a->dyad.jrecv == b->dyad.jrecv;
}

static bool vrecv_init(struct design_var *dv, const struct design *d)
{
	assert(dv);
	assert(d);

	ssize_t n = vector_dim(design_intervals(d));
	dv->dim = n + 1;
	return true;
}

static void vrecv_deinit(struct design_var *dv)
{
	assert(dv);
}

static bool vrecv_frame_init(struct frame_var *fv, struct frame *f)
{
	assert(fv);
	assert(f);
	
	struct vrecv_udata *udata;
	
	if (!(udata = malloc(sizeof(*udata))))
		goto fail_malloc;

	if (!hashset_init(&udata->active, vrecv_active_hash, vrecv_active_equals,
			  sizeof(struct vrecv_active)))
		goto fail_active;
	
	fv->udata = udata;
	return true;
	
	hashset_deinit(&udata->active);
fail_active:
	free(udata);
fail_malloc:
	return false;
}

static void vrecv_frame_deinit(struct frame_var *fv)
{
	assert(fv);
	assert(fv->udata);
	
	struct vrecv_udata *udata = fv->udata;
	hashset_deinit(&udata->active);
	free(udata);
}

static void vrecv_frame_clear(struct frame_var *fv)
{
	assert(fv);
	assert(fv->udata);
	
	struct vrecv_udata *udata = fv->udata;
	hashset_clear(&udata->active);
}

static bool vrecv_handle_dyad (struct frame_var *fv, const struct dyad_event *e,
			       struct frame *f)
{
	assert(fv);
	assert(e);
	assert(f);
	assert(fv->design);
	assert(fv->design->index >= 0);
	assert(fv->design->index <= design_dim(f->design) - fv->design->dim);
	assert(fv->udata);

	ssize_t index = fv->design->index;
	struct vrecv_udata *udata = fv->udata;
	const struct vrecv_active *key =
	    container_of(&e->dyad, const struct vrecv_active, dyad);
	struct vrecv_active *active;
	struct svector *dx = frame_dx(f, e->dyad.jrecv, e->dyad.isend);

	if (!dx)
		return false;

	if (e->type == DYAD_EVENT_INIT) {
		struct hashset_pos pos;

		if ((active = hashset_find(&udata->active, key, &pos)))
			goto move;

		if ((active = hashset_insert(&udata->active, &pos, key)))
			goto init;

		return false;
	} else {		// e->type == DYAD_EVENT_MOVE
		active = hashset_lookup(&udata->active, key);
		assert(active);

		if (active->id == e->id) {
			assert(active->intvl == e->intvl - 1);
			goto move;
		}

		return true;
	}

move:
	svector_clear_at(dx, index + active->intvl);
init:
	if (!svector_set(dx, index + e->intvl, 1.0))
		return false;
	active->id = e->id;
	active->intvl = e->intvl;
	return true;
}

static struct var_type VAR_TYPE_RECV_REP = {
	DYAD_EVENT_INIT | DYAD_EVENT_MOVE,
	vrecv_init,
	vrecv_deinit,
	vrecv_frame_init,
	vrecv_frame_deinit,
	vrecv_frame_clear,
	vrecv_handle_dyad
};

const struct var_type *VAR_TYPE_RECV = &VAR_TYPE_RECV_REP;
