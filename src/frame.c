#include "port.h"
#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include "sblas.h"
#include "coreutil.h"
#include "xalloc.h"
#include "ieee754.h"
#include "vars.h"
#include "frame.h"

static void frame_history_message_add (void *udata, struct history * h, const struct message * msg)
{
	struct frame *f = udata;
	size_t i, n = f->nobs;
	const struct frame_observer *obs;
	for (i = 0; i < n; i++) {
		obs = &f->observers[i];
		if (obs->callbacks.message_add) {
			obs->callbacks.message_add(obs->udata, f, msg);
		}
	}
}

static void frame_history_message_advance (void *udata, struct history * h, const struct message * msg, size_t intvl)
{
	struct frame *f = udata;
	size_t i, n = f->nobs;
	const struct frame_observer *obs;
	for (i = 0; i < n; i++) {
		obs = &f->observers[i];
		if (obs->callbacks.message_advance) {
			obs->callbacks.message_advance(obs->udata, f, msg, intvl);
		}
	}
}

static void frame_history_clear(void *udata, struct history * h)
{
	struct frame *f = udata;
	size_t i, n = f->nobs;
	const struct frame_observer *obs;
	for (i = 0; i < n; i++) {
		obs = &f->observers[i];
		if (obs->callbacks.clear) {
			obs->callbacks.clear(obs->udata, f);
		}
	}
}

static struct history_callbacks frame_history_callbacks = {
	frame_history_message_add,
	frame_history_message_advance,
	frame_history_clear
};

static void recv_frame_init(struct recv_frame *rf, struct frame *f)
{
	(void)f;		// unused;
	assert(rf);

	vpattern_init(&rf->active);
	rf->dx = NULL;
}

static void recv_frame_clear(struct recv_frame *rf)
{
	assert(rf);
	vpattern_clear(&rf->active);
}

static void recv_frame_deinit(struct recv_frame *rf)
{
	assert(rf);

	recv_frame_clear(rf);
	free(rf->dx);
	vpattern_deinit(&rf->active);
}

static double *recv_frame_dx(struct recv_frame *rf, size_t jrecv,
			     size_t dyn_dim)
{
	assert(rf);

	int ins;
	size_t nzmax = rf->active.nzmax;
	size_t ix = vpattern_search(&rf->active, jrecv, &ins);

	if (ins) {
		if (nzmax != rf->active.nzmax) {
			nzmax = rf->active.nzmax;
			rf->dx = xrealloc(rf->dx, nzmax * dyn_dim * sizeof(rf->dx[0]));
		}

		memmove(rf->dx + (ix + 1) * dyn_dim,
			rf->dx + ix * dyn_dim,
			(rf->active.nz - 1 - ix) * dyn_dim * sizeof(rf->dx[0]));
		memset(rf->dx + ix * dyn_dim, 0, dyn_dim * sizeof(rf->dx[0]));
	}

	return rf->dx + ix * dyn_dim;
}

static void recv_frames_init(struct frame *f)
{
	assert(f);

	size_t isend, nsend = frame_send_count(f);

	struct recv_frame *rfs = xcalloc(nsend, sizeof(struct recv_frame));

	for (isend = 0; isend < nsend; isend++) {
		recv_frame_init(&rfs[isend], f);
	}
	f->recv_frames = rfs;
}

static void recv_frames_deinit(struct frame *f)
{
	assert(f);
	size_t isend, nsend = frame_send_count(f);

	struct recv_frame *rfs = f->recv_frames;

	for (isend = 0; isend < nsend; isend++) {
		recv_frame_deinit(&rfs[isend]);
	}
	free(rfs);
}

static void recv_frames_clear(struct frame *f)
{
	assert(f);
	size_t isend, nsend = frame_send_count(f);
	struct recv_frame *rfs = f->recv_frames;

	for (isend = 0; isend < nsend; isend++) {
		recv_frame_clear(&rfs[isend]);
	}
}

static struct recv_frame *recv_frames_item(const struct frame *f, size_t isend)
{
	assert(f);
	assert(isend < frame_send_count(f));

	struct recv_frame *rf = &f->recv_frames[isend];
	return rf;
}

void frame_init(struct frame *f, size_t nsend, size_t nrecv, int has_loops,
		const double *intvls, size_t nintvl)
{
	assert(f);
	assert(intvls || !nintvl);

	history_init(&f->history, nsend, nrecv, intvls, nintvl);
	history_add_observer(&f->history, f, &frame_history_callbacks);
	design_init(&f->send_design, f, nsend);
	design_init(&f->recv_design, f, nrecv);
	f->has_loops = has_loops;

	f->observers = NULL;
	f->nobs = 0;
	f->nobs_max = 0;
	recv_frames_init(f);
	frame_clear(f);
}

void frame_deinit(struct frame *f)
{
	assert(f);

	size_t nsend = frame_send_count(f);
	size_t nrecv = frame_recv_count(f);

	recv_frames_deinit(f);
	free(f->observers);
	design_deinit(&f->recv_design);
	design_deinit(&f->send_design);
	history_remove_observer(&f->history, f);
	history_deinit(&f->history);
}

void frame_clear(struct frame *f)
{
	assert(f);

	size_t nsend = frame_send_count(f);
	size_t nrecv = frame_recv_count(f);

	recv_frames_clear(f);
	history_clear(&f->history);

}

static void frame_observers_grow(struct frame *f, size_t delta)
{
	size_t nmax = array_grow(f->nobs, f->nobs_max, delta, SIZE_MAX);
	if (nmax > f->nobs_max) {
		f->observers = xrealloc(f->observers, nmax * sizeof(f->observers[0]));
		f->nobs_max = nmax;
	}
}

void frame_add_observer(struct frame *f, void *udata,
			const struct frame_callbacks *callbacks)
{
	assert(f);
	assert(udata);
	assert(callbacks);

	frame_observers_grow(f, 1);
	struct frame_observer *obs = &f->observers[f->nobs++];
	obs->udata = udata;
	obs->callbacks = *callbacks;
}

void frame_remove_observer(struct frame *f, void *udata)
{
	assert(f);
	assert(udata);

	size_t i, n = f->nobs;
	for (i = n; i > 0; i--) {
		if (f->observers[i-1].udata == udata) {
			memmove(f->observers + i - 1, f->observers + i,
				(n - i) * sizeof(f->observers[0]));
			f->nobs = n - 1;
			break;
		}
	}
}

void frame_add(struct frame *f, const struct message *msg)
{
	assert(f);
	assert(msg);
	assert(msg->time == frame_time(f));

	history_add(&f->history, msg);
}


void frame_advance(struct frame *f, double time)
{
	assert(f);
	history_advance(&f->history, time);
}

void frame_recv_update(struct frame *f, size_t isend, size_t jrecv,
		       const double *delta, const struct vpattern *pat)
{
	assert(isend < frame_send_count(f));
	assert(jrecv < frame_recv_count(f));

	double *dx = (double *)frame_recv_dx(f, isend, jrecv);
	sblas_daxpyi(1.0, delta, pat, dx);

	size_t i, n = f->nobs;
	const struct frame_observer *obs;
	for (i = 0; i < n; i++) {
		obs = &f->observers[i];
		if (obs->callbacks.recv_update) {
			obs->callbacks.recv_update(obs->udata, f, isend,
						   jrecv, delta, pat);
		}
	}
}

const double *frame_recv_dx(const struct frame *f, size_t isend, size_t jrecv)
{
	assert(f);
	assert(isend < frame_send_count(f));
	assert(jrecv < frame_recv_count(f));

	struct recv_frame *rf = recv_frames_item((struct frame *)f, isend);
	const struct design *d = frame_recv_design(f);
	size_t dyn_dim = design_dvars_dim(d);
	return recv_frame_dx(rf, jrecv, dyn_dim);
}

void frame_recv_get_dx(const struct frame *f, size_t isend,
		       double **dxp, size_t **activep,
		       size_t *nactivep)
{
	assert(isend < frame_send_count(f));

	struct recv_frame *rf = recv_frames_item((struct frame *)f, isend);

	*dxp = rf->dx;
	*activep = rf->active.indx;
	*nactivep = rf->active.nz;
}

void frame_recv_mul(double alpha, enum blas_trans trans,
		    const struct frame *f, size_t isend,
		    const double *x, double beta, double *y)
{
	assert(isend < frame_send_count(f));

	const struct design *d = frame_recv_design(f);
	size_t off = design_dvars_index(d);

	design_mul0(alpha, trans, d, x, beta, y);

	if (trans == BLAS_NOTRANS) {
		frame_recv_dmul(alpha, trans, f, isend, x + off, 1.0, y);
	} else {
		frame_recv_dmul(alpha, trans, f, isend, x, 1.0, y + off);
	}
}

void frame_recv_muls(double alpha, enum blas_trans trans,
		     const struct frame *f, size_t isend,
		     const double *x, const struct vpattern *pat,
		     double beta, double *y)
{
	assert(isend < frame_send_count(f));

	const struct design *d = frame_recv_design(f);
	size_t off = design_dvars_index(d);
	size_t dim = design_dvars_dim(d);

	design_muls0(alpha, trans, d, x, pat, beta, y);

	if (trans == BLAS_NOTRANS) {
		ptrdiff_t jx0 = vpattern_find(pat, off);
		ptrdiff_t jx1 = vpattern_find(pat, off + dim);
		size_t jz0 = (jx0 < 0 ? ~jx0 : jx0);
		size_t jz1 = (jx1 < 0 ? ~jx1 : jx1);
		size_t jz;

		const struct recv_frame *rf = recv_frames_item(f, isend);
		size_t iz, nz = rf->active.nz;

		for (iz = 0; iz < nz; iz++) {
			size_t jrecv = rf->active.indx[iz];
			const double *dx = rf->dx + iz * dim;

			double dot = 0;
			for (jz = jz0; jz < jz1; jz++) {
				size_t i = pat->indx[jz] - off;
				dot += x[jz] * dx[i];
			}

			y[jrecv] += alpha * dot;
		}
	} else {
		frame_recv_dmuls(alpha, trans, f, isend, x, pat, 1.0, y + off);
	}
}

void frame_recv_dmul(double alpha, enum blas_trans trans,
		     const struct frame *f, size_t isend,
		     const double *x, double beta, double *y)
{
	assert(isend < frame_send_count(f));

	const struct design *d = frame_recv_design(f);
	size_t n = design_count(d);
	size_t dim = design_dvars_dim(d);
	size_t ny = (trans == BLAS_NOTRANS ? n : dim);

	/* y := beta y */
	if (beta == 0.0) {
		memset(y, 0, ny * sizeof(y[0]));
	} else if (beta != 1.0) {
		blas_dscal(ny, beta, y, 1);
	}

	if (dim == 0)
		return;

	const struct recv_frame *rf = recv_frames_item(f, isend);
	size_t iz, nz = rf->active.nz;
	size_t jrecv;
	const double *dx;

	if (trans == BLAS_NOTRANS) {
		for (iz = 0; iz < nz; iz++) {
			jrecv = rf->active.indx[iz];
			dx = rf->dx + iz * dim;

			double dot = blas_ddot(dim, dx, 1, x, 1);
			y[jrecv] += alpha * dot;
		}
	} else {
		for (iz = 0; iz < nz; iz++) {
			jrecv = rf->active.indx[iz];
			dx = rf->dx + iz * dim;

			if (x[jrecv] == 0.0)
				continue;

			blas_daxpy(dim, alpha * x[jrecv], dx, 1, y, 1);
		}
	}
}

void frame_recv_dmuls(double alpha, enum blas_trans trans,
		      const struct frame *f, size_t isend,
		      const double *x, const struct vpattern *pat,
		      double beta, double *y)
{
	assert(isend < frame_send_count(f));

	const struct design *d = frame_recv_design(f);
	size_t n = design_count(d);
	size_t dim = design_dvars_dim(d);
	size_t ny = (trans == BLAS_NOTRANS ? n : dim);

	/* y := beta y */
	if (beta == 0.0) {
		memset(y, 0, ny * sizeof(y[0]));
	} else if (beta != 1.0) {
		blas_dscal(ny, beta, y, 1);
	}

	if (dim == 0)
		return;

	const struct recv_frame *rf = recv_frames_item(f, isend);
	size_t iz, nz = rf->active.nz;
	size_t jrecv;
	const double *dx;

	if (trans == BLAS_NOTRANS) {

		for (iz = 0; iz < nz; iz++) {
			jrecv = rf->active.indx[iz];
			dx = rf->dx + iz * dim;

			double dot = sblas_ddoti(x, pat, dx);
			y[jrecv] += alpha * dot;
		}
	} else {
		size_t jz, mz = pat->nz;
		for (jz = 0; jz < mz; jz++) {
			size_t jrecv = pat->indx[jz];
			ptrdiff_t ix = vpattern_find(&rf->active, jrecv);
			if (ix < 0)
				continue;

			dx = rf->dx + ix * dim;
			/* y := y + alpha * x[j] * dx[j] */
			double jscale = alpha * x[jz];
			blas_daxpy(dim, jscale, dx, 1, y, 1);
		}
	}
}
