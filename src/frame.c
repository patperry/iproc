#include "port.h"
#include <assert.h>
#include <stdlib.h>
#include <math.h>
#include "ieee754.h"
#include "frame.h"

static bool send_frame_init(struct send_frame *sf, struct frame *f)
{
	assert(sf);

	sf->frame = f;
	return true;
}


static void send_frame_deinit(struct send_frame *sf)
{
	assert(sf);
}

static struct svector *frame_dx_alloc(struct frame *f)
{
	assert(f);
	assert(0 <= f->next_dx && f->next_dx <= darray_size(&f->dxs));
	
	struct svector *dx;
	
	if (f->next_dx == darray_size(&f->dxs)) {
		if (!((dx = malloc(sizeof(*dx)))))
			goto fail_malloc;
		
		if (!svector_init(dx, design_dim(f->design)))
			goto fail_init;
		
		if (!darray_push_back(&f->dxs, &dx))
			goto fail_push_back;
		
		goto success;
		
	fail_push_back:
		svector_deinit(dx);
	fail_init:
		free(dx);
	fail_malloc:
		return NULL;
	}
success:
	assert(f->next_dx < darray_size(&f->dxs));

	dx = *(struct svector **)darray_at(&f->dxs, f->next_dx++);
	svector_clear(dx);

	assert(svector_dim(dx) == design_dim(f->design));
	return dx;
}

bool frame_init(struct frame *f, struct design *design)
{
	assert(f);
	assert(design);
	
	if (!history_init(&f->history))
	    goto fail_history;
	
	if (!intmap_init(&f->send_frames, sizeof(struct send_frame),
			 alignof(struct send_frame)))
		goto fail_send_frames;
	
	if (!intmap_init(&f->jrecv_dxs, sizeof(struct svector *),
			 alignof(struct svector *)))
		goto fail_jrecv_dxs;

	if (!darray_init(&f->dxs, sizeof(struct svector *)))
		goto fail_dxs;
	
	f->design = design;
	f->time = -INFINITY;
	f->isend = -1;
	f->next_dx = 0;
	return true;
	
fail_dxs:
	intmap_deinit(&f->jrecv_dxs);
fail_jrecv_dxs:
	intmap_deinit(&f->send_frames);
fail_send_frames:
	history_deinit(&f->history);
fail_history:
	return false;
}


void frame_deinit(struct frame *f)
{
	assert(f);
	
	struct intmap_iter it;
	struct send_frame *sf;
	struct svector *dx;

	ssize_t i, n = darray_size(&f->dxs);
	for (i = 0; i < n; i++) {
		dx = *(struct svector **)darray_at(&f->dxs, i);
		svector_deinit(dx);
		free(dx);
	}
	darray_deinit(&f->dxs);
	intmap_deinit(&f->jrecv_dxs);
	
	intmap_iter_init(&f->send_frames, &it);
	while (intmap_iter_advance(&f->send_frames, &it)) {
		sf = intmap_iter_current(&f->send_frames, &it);
		send_frame_deinit(sf);
	}
	intmap_iter_deinit(&f->send_frames, &it);
	
	intmap_deinit(&f->send_frames);
	history_deinit(&f->history);
}

static void send_frame_clear(struct send_frame *sf)
{
	assert(sf);
}

void frame_clear(struct frame *f)
{
	assert(f);
	
	struct intmap_iter it;
	struct send_frame *sf;
	
	intmap_iter_init(&f->send_frames, &it);
	while (intmap_iter_advance(&f->send_frames, &it)) {
		sf = intmap_iter_current(&f->send_frames, &it);
		send_frame_clear(sf);
	}
	intmap_iter_deinit(&f->send_frames, &it);
	
	f->time = -INFINITY;
	f->isend = -1;
	
	intmap_clear(&f->jrecv_dxs);
	f->next_dx = 0;
}

double frame_time(const struct frame *f)
{
	assert(f);
	return f->time;
}

const struct history *frame_history(const struct frame *f)
{
	assert(f);
	return &f->history;
}

bool frame_insert(struct frame *f, const struct message *msg)
{
	assert(f);
	assert(msg);
	assert(msg->time == frame_time(f));
	
	ssize_t i, n = darray_size(&f->design->design_dyad_vars);
	const struct design_dyad_var *design_dyad_var;
	
	for (i = 0; i < n; i++) {
		design_dyad_var = darray_at(&f->design->design_dyad_vars, i);
	}
	
	return history_insert(&f->history, msg->from, msg->to, msg->nto, msg->attr);
}


static struct send_frame *frame_get_send(struct frame *f, ssize_t isend)
{
	assert(f);
	assert(0 <= isend && isend < design_nsender(f->design));
	
	struct intmap_pos pos;
	struct send_frame *sf;

	if ((sf = intmap_find(&f->send_frames, isend, &pos)))
		return sf;
		
	if ((sf = intmap_insert(&f->send_frames, &pos, NULL))) {
		if (send_frame_init(sf, f))
			return sf;
		
		intmap_erase(&f->send_frames, &pos);
	}

	return NULL;
}


struct svector *frame_dx(struct frame *f, ssize_t jrecv)
{
	assert(f);
	assert(0 <= jrecv && jrecv < design_nreceiver(f->design));
	
	struct intmap_pos pos;
	struct svector **pdx;	
	struct svector *dx;
	
	if ((pdx = intmap_find(&f->jrecv_dxs, jrecv, &pos))) {
		return *pdx;
	}
	
	if ((dx = frame_dx_alloc(f))
	    && intmap_insert(&f->jrecv_dxs, &pos, &dx)) {
		return dx;
	}
	
	return NULL;
}


bool frame_advance_to(struct frame *f, double t)
{
	assert(f);
	assert(t >= frame_time(f));
	
	if (t == frame_time(f))
		return true;

	if (!history_advance_to(&f->history, t))
		return false;
	
	// TODO: update var sums
	
	
	intmap_clear(&f->jrecv_dxs);
	f->next_dx = 0;
	f->isend = -1;
	f->time = t;
	
	return true;
}

ssize_t frame_sender(const struct frame *f)
{
	assert(f);
	return f->isend;
}

bool frame_set_sender(struct frame *f, ssize_t isend)
{
	assert(f);
	assert(0 <= isend && isend < design_nsender(f->design));

	// set sender and clear old jrecv_dxs
	f->isend = isend;
	intmap_clear(&f->jrecv_dxs);
	f->next_dx = 0;

	// get new jrecv_dxs
	const struct darray *design_dyad_vars = &f->design->design_dyad_vars;
	struct design_dyad_var *v;
	ssize_t i, n = darray_size(design_dyad_vars);
	bool ok;
	
	for (i = 0; i < n; i++) {
		v = darray_at(design_dyad_vars, i);
		ok = v->var->get_jrecv_dxs((struct dyad_var *)v->var, f, v->index);
		if (!ok)
			return false;
	}
	
	return true;
}


bool frame_mul(double alpha, enum trans_op trans,
	       const struct frame *f,
	       const struct vector *x, double beta, struct vector *y)
{
	assert(f);
	assert(f->design);
	assert(f->isend >= 0);
	assert(x);
	assert(y);
	assert(trans != TRANS_NOTRANS
	       || vector_dim(x) == design_dim(f->design));
	assert(trans != TRANS_NOTRANS
	       || vector_dim(y) == design_nreceiver(f->design));
	assert(trans == TRANS_NOTRANS
	       || vector_dim(x) == design_nreceiver(f->design));
	assert(trans == TRANS_NOTRANS
	       || vector_dim(y) == design_dim(f->design));
	
	
	struct svector diffprod;
	bool ok = false;
	
	if (svector_init(&diffprod, vector_dim(y))) {
		if (frame_dmul(alpha, trans, f, x, 0.0, &diffprod)) {
			design_mul0(alpha, trans, f->design, f->isend, x, beta, y);
			svector_axpy(1.0, &diffprod, y);
			ok = true;
		}
		svector_deinit(&diffprod);
	}
	
	return ok;
}


bool frame_muls(double alpha, enum trans_op trans,
		const struct frame *f,
		const struct svector *x, double beta, struct vector *y)
{
	assert(f);
	assert(f->design);
	assert(f->isend >= 0);
	assert(x);
	assert(y);
	assert(trans != TRANS_NOTRANS
	       || svector_dim(x) == design_dim(f->design));
	assert(trans != TRANS_NOTRANS
	       || vector_dim(y) == design_nreceiver(f->design));
	assert(trans == TRANS_NOTRANS
	       || svector_dim(x) == design_nreceiver(f->design));
	assert(trans == TRANS_NOTRANS
	       || vector_dim(y) == design_dim(f->design));
	
	
	struct svector diffprod;
	bool ok = false;
	
	if (svector_init(&diffprod, vector_dim(y))) {
		if (frame_dmuls(alpha, trans, f, x, 0.0, &diffprod)) {
			design_muls0(alpha, trans, f->design, f->isend, x, beta, y);
			svector_axpy(1.0, &diffprod, y);
			ok = true;
		}
		
		svector_deinit(&diffprod);
	}
	
	return ok;
}


bool frame_dmul(double alpha, enum trans_op trans,
		const struct frame *f,
		const struct vector *x, double beta, struct svector *y)
{
	assert(f);
	assert(f->design);
	assert(f->isend >= 0);
	assert(x);
	assert(y);
	assert(trans != TRANS_NOTRANS
	       || vector_dim(x) == design_dim(f->design));
	assert(trans != TRANS_NOTRANS
	       || svector_dim(y) == design_nreceiver(f->design));
	assert(trans == TRANS_NOTRANS
	       || vector_dim(x) == design_nreceiver(f->design));
	assert(trans == TRANS_NOTRANS
	       || svector_dim(y) == design_dim(f->design));
	
	/* y := beta y */
	if (beta == 0.0) {
		svector_clear(y);
	} else if (beta != 1.0) {
		svector_scale(y, beta);
	}
	
	struct design *design = f->design;
	ssize_t ndynamic = design->ndynamic;
	
	if (ndynamic == 0)
		return true;

	struct intmap_iter it;
	ssize_t jrecv;
	const struct svector *dx;
	bool ok = true;
	
	intmap_iter_init(&f->jrecv_dxs, &it);
	
	while (intmap_iter_advance(&f->jrecv_dxs, &it)) {
		jrecv = intmap_iter_current_key(&f->jrecv_dxs, &it);
		dx = *(struct svector **)intmap_iter_current(&f->jrecv_dxs, &it);
		if (trans == TRANS_NOTRANS) {
			double dot = svector_dot(dx, x);
			double *yjrecv = svector_at(y, jrecv);
			if (!yjrecv) {
				ok = false;
				break;
			}
			
			*yjrecv += alpha * dot;
		} else {
			double xjrecv = *vector_at(x, jrecv);
			
			if (xjrecv == 0.0)
				continue;
			
			/* y := y + alpha * x[j] * dx[j] */
			double jscale = alpha * xjrecv;
			if (!svector_axpys(jscale, dx, y)) {
				ok = false;
				break;
			}
		}
	}

	intmap_iter_deinit(&f->jrecv_dxs, &it);
	return ok;
}


bool frame_dmuls(double alpha, enum trans_op trans,
		 const struct frame *f,
		 const struct svector *x, double beta, struct svector *y)
{
	assert(f);
	assert(f->design);
	assert(f->isend >= 0);	
	assert(x);
	assert(y);
	assert(trans != TRANS_NOTRANS
	       || svector_dim(x) == design_dim(f->design));
	assert(trans != TRANS_NOTRANS
	       || svector_dim(y) == design_nreceiver(f->design));
	assert(trans == TRANS_NOTRANS
	       || svector_dim(x) == design_nreceiver(f->design));
	assert(trans == TRANS_NOTRANS
	       || svector_dim(y) == design_dim(f->design));
	
	/* y := beta y */
	if (beta == 0.0) {
		svector_clear(y);
	} else if (beta != 1.0) {
		svector_scale(y, beta);
	}
	
	struct design *design = f->design;
	ssize_t ndynamic = design->ndynamic;
	
	if (ndynamic == 0)
		return true;
	
	struct intmap_iter it;
	ssize_t jrecv;
	const struct svector *dx;
	bool ok = true;
	
	intmap_iter_init(&f->jrecv_dxs, &it);
	
	while (intmap_iter_advance(&f->jrecv_dxs, &it)) {
		jrecv = intmap_iter_current_key(&f->jrecv_dxs, &it);
		dx = *(struct svector **)intmap_iter_current(&f->jrecv_dxs, &it);

		if (trans == TRANS_NOTRANS) {
			double dot = svector_dots(dx, x);
			double *yjrecv = svector_at(y, jrecv);
			if (!yjrecv) {
				ok = false;
				break;
			}
			
			*yjrecv += alpha * dot;
		} else {
			double xjrecv = svector_get(x, jrecv);
			
			if (xjrecv == 0.0)
				continue;
			
			/* y := y + alpha * x[j] * dx[j] */
			double jscale = alpha * xjrecv;
			if (!svector_axpys(jscale, dx, y)) {
				ok = false;
				break;
			}
		}
	}
	
	return ok;
}

