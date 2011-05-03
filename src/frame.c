#include "port.h"
#include <assert.h>
#include <math.h>
#include "ieee754.h"
#include "frame.h"

static bool send_frame_init(struct send_frame *sf, struct frame *f)
{
	assert(sf);

	if (!intmap_init(&sf->jrecv_dxs, sizeof(struct svector),
			 alignof(struct svector)))
		goto fail_jrecv_dxs;
	
	sf->frame = f;
	return true;
fail_jrecv_dxs:
	return false;
}


static void send_frame_deinit(struct send_frame *sf)
{
	assert(sf);
	
	struct intmap_iter it;
	struct svector *dx;
	
	intmap_iter_init(&sf->jrecv_dxs, &it);
	while (intmap_iter_advance(&sf->jrecv_dxs, &it)) {
		dx = intmap_iter_current(&sf->jrecv_dxs, &it);
		svector_deinit(dx);
	}
	intmap_iter_deinit(&sf->jrecv_dxs, &it);

	intmap_deinit(&sf->jrecv_dxs);
}

static int dyad_var_diff_rcompare(const void *px, const void *py)
{
	return double_compare(&((struct dyad_var_diff *)py)->time,
			      &((struct dyad_var_diff *)px)->time);
}


bool frame_init(struct frame *f, struct design *design)
{
	assert(f);
	assert(design);
	
	if (!pqueue_init(&f->dyad_var_diffs, dyad_var_diff_rcompare,
			 sizeof(struct dyad_var_diff)))
		goto fail_dyad_var_diffs;

	if (!intmap_init(&f->send_frames, sizeof(struct send_frame),
			 alignof(struct send_frame)))
		goto fail_send_frames;
	
	f->design = design;
	f->time = -INFINITY;
	return true;
	
fail_send_frames:
	pqueue_deinit(&f->dyad_var_diffs);
fail_dyad_var_diffs:
	return false;
}


void frame_deinit(struct frame *f)
{
	assert(f);
	
	struct intmap_iter it;
	struct send_frame *sf;
	
	intmap_iter_init(&f->send_frames, &it);
	while (intmap_iter_advance(&f->send_frames, &it)) {
		sf = intmap_iter_current(&f->send_frames, &it);
		send_frame_deinit(sf);
	}
	intmap_iter_deinit(&sf->jrecv_dxs, &it);
	
	intmap_deinit(&f->send_frames);
	pqueue_deinit(&f->dyad_var_diffs);
}

static void send_frame_clear(struct send_frame *sf)
{
	assert(sf);

	struct intmap_iter it;
	struct svector *dx;
	
	intmap_iter_init(&sf->jrecv_dxs, &it);
	while (intmap_iter_advance(&sf->jrecv_dxs, &it)) {
		dx = intmap_iter_current(&sf->jrecv_dxs, &it);
		svector_clear(dx);
	}
	intmap_iter_deinit(&sf->jrecv_dxs, &it);
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
	
	pqueue_clear(&f->dyad_var_diffs);
	f->time = -INFINITY;
}

double frame_time(const struct frame *f)
{
	assert(f);
	return f->time;
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
		if (!design_dyad_var->var->insert(design_dyad_var->var,
						  msg,
						  f,
						  design_dyad_var->index))
			return false;
	}
	
	return true;
}


double frame_next_update(const struct frame *f)
{
	assert(f);
	
	if (pqueue_empty(&f->dyad_var_diffs))
		return INFINITY;
	
	return ((struct dyad_var_diff *)pqueue_top(&f->dyad_var_diffs))->time;
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

static struct svector *send_frame_get_recv(struct send_frame *sf, ssize_t jrecv)
{
	assert(sf);
	assert(0 <= jrecv && jrecv < design_nreceiver(sf->frame->design));
	
	struct intmap_pos pos;
	struct svector *dx;
	
	if ((dx = intmap_find(&sf->jrecv_dxs, jrecv, &pos))) {
		return dx;
	}
	
	if ((dx = intmap_insert(&sf->jrecv_dxs, &pos, NULL))) {
		if (svector_init(dx, design_dim(sf->frame->design)))
			return dx;
		
		intmap_erase(&sf->jrecv_dxs, &pos);
	}
	
	return NULL;
}


static bool frame_apply_dyad_var_diff(struct frame *f,
				      const struct dyad_var_diff *diff)
{
	assert(f);
	assert(diff);
	assert(diff->time >= frame_time(f));
	assert(f->design->idynamic <= diff->index);
	assert(diff->index < f->design->idynamic + f->design->ndynamic);
	       
	struct send_frame *sf;
	struct svector *dx;
	double *val;
	
	if ((sf = frame_get_send(f, diff->isend))
	    && (dx = send_frame_get_recv(sf, diff->jrecv))
	    && (val = svector_at(dx, diff->index))) {
		*val += diff->delta; // TODO: clear zeros?
		f->time = diff->time;
		return true;
	}
	return false;
}


bool frame_advance_to(struct frame *f, double t, struct frame_diff *diff)
{
	assert(f);
	assert(t >= frame_time(f));

	if (t == frame_time(f))
		return true;
	
	struct dyad_var_diff *dyad_var_diff;
	
	while (!pqueue_empty(&f->dyad_var_diffs) && frame_next_update(f) <= t) {
		dyad_var_diff = pqueue_top(&f->dyad_var_diffs);
		
		if (diff && !darray_push_back(&diff->dyad_var_diffs,
					      dyad_var_diff)) {
			return false;
		}
		
		if (!frame_apply_dyad_var_diff(f, dyad_var_diff)) {
			if (diff)
				darray_pop_back(&diff->dyad_var_diffs);
			return false;
		}
		pqueue_pop(&f->dyad_var_diffs);
	}
	f->time = t;
	return true;
}


bool frame_add_dyad_event(struct frame *f, const struct dyad_var_diff *delta)
{
	assert(f);
	assert(delta);
	assert(!(delta->time <= frame_time(f)));
	
	return pqueue_push(&f->dyad_var_diffs, delta);
}


bool frame_reserve_dyad_events(struct frame *f, ssize_t nadd)
{
	assert(f);
	assert(nadd >= 0);

	return pqueue_reserve_push(&f->dyad_var_diffs, nadd);
}


bool frame_mul(double alpha, enum trans_op trans,
	       const struct frame *f, ssize_t isend,
	       const struct vector *x, double beta, struct vector *y)
{
	assert(f);
	assert(f->design);
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
		if (frame_dmul(alpha, trans, f, isend, x, 0.0, &diffprod)) {
			design_mul0(alpha, trans, f->design, isend, x, beta, y);
			svector_axpy(1.0, &diffprod, y);
			ok = true;
		}
		svector_deinit(&diffprod);
	}
	
	return ok;
}


bool frame_muls(double alpha, enum trans_op trans,
		const struct frame *f, ssize_t isend,
		const struct svector *x, double beta, struct vector *y)
{
	assert(f);
	assert(f->design);
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
		if (frame_dmuls(alpha, trans, f, isend, x, 0.0, &diffprod)) {
			design_muls0(alpha, trans, f->design, isend, x, beta, y);
			svector_axpy(1.0, &diffprod, y);
			ok = true;
		}
		
		svector_deinit(&diffprod);
	}
	
	return ok;
}


bool frame_dmul(double alpha, enum trans_op trans,
		const struct frame *f, ssize_t isend,
		const struct vector *x, double beta, struct svector *y)
{
	assert(f);
	assert(f->design);
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

	const struct send_frame *sf = intmap_lookup(&f->send_frames, isend);
	
	if (!sf)
		return true;
	
	struct intmap_iter it;
	ssize_t jrecv;
	const struct svector *dx;
	bool ok = true;
	
	intmap_iter_init(&sf->jrecv_dxs, &it);
	
	while (intmap_iter_advance(&sf->jrecv_dxs, &it)) {
		jrecv = intmap_iter_current_key(&sf->jrecv_dxs, &it);
		dx = intmap_iter_current(&sf->jrecv_dxs, &it);
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

	intmap_iter_deinit(&sf->jrecv_dxs, &it);
	return ok;
}


bool frame_dmuls(double alpha, enum trans_op trans,
		 const struct frame *f, ssize_t isend,
		 const struct svector *x, double beta, struct svector *y)
{
	assert(f);
	assert(f->design);
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
	
	const struct send_frame *sf = intmap_lookup(&f->send_frames, isend);
	
	if (!sf)
		return true;
	
	struct intmap_iter it;
	ssize_t jrecv;
	const struct svector *dx;
	bool ok = true;
	
	intmap_iter_init(&sf->jrecv_dxs, &it);
	
	while (intmap_iter_advance(&sf->jrecv_dxs, &it)) {
		jrecv = intmap_iter_current_key(&sf->jrecv_dxs, &it);
		dx = intmap_iter_current(&sf->jrecv_dxs, &it);

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

