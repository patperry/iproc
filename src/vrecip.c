#include "port.h"

#include <assert.h>
#include <stddef.h>
#include <stdlib.h>
#include "compare.h"
#include "ieee754.h"
#include "design.h"
#include "util.h"
#include "vrecip.h"

static void iproc_vrecip_free(iproc_vrecip * v)
{
	if (v) {
		darray_deinit(&v->intvls);
		free(v);
	}
}

static void design_var_free(iproc_design_var * var)
{
	iproc_vrecip *v = container_of(var, iproc_vrecip, var);
	iproc_vrecip_unref(v);
}

static void
design_var_get_dxs(iproc_design_var * var,
		   iproc_design_ctx * ctx, ssize_t offset)
{
	iproc_vrecip *v = container_of(var, iproc_vrecip, var);
	struct darray *intvls = &v->intvls;
	ssize_t nintvl = darray_size(intvls);
	iproc_history *history = ctx->history;
	ssize_t isend = ctx->isend;

	if (!history || nintvl == 0)
		return;

	double tcur = iproc_history_tcur(history);
	iproc_trace *trace = iproc_history_recv(history, isend);

	ssize_t i, n = iproc_trace_size(trace);
	for (i = 0; i < n; i++) {
		iproc_events *events = iproc_trace_get(trace, i);
		iproc_event_meta *meta = iproc_events_last(events);

		if (!meta)
			continue;

		ssize_t jsend = iproc_events_id(events);
		double t = meta->time;
		double dt = tcur - t;
		ssize_t pos = darray_binary_search(intvls,
						   &dt, double_compare);

		if (pos < 0)
			pos = ~pos;

		if (pos < nintvl) {
			/* (jsend, [(pos, +1.0)]) */
			struct svector *dx =
			    iproc_design_ctx_dx(ctx, jsend, false);
			*svector_at(dx, offset + pos) += 1.0;
		}
	}
}

iproc_vrecip *iproc_vrecip_new(double *intvls, ssize_t n)
{
	assert(n >= 0);
	assert(n == 0 || intvls);

	iproc_vrecip *v = calloc(1, sizeof(*v));

	if (!v)
		return NULL;

	iproc_design_var_init(&v->var, n, design_var_get_dxs, design_var_free);
	refcount_init(&v->refcount);

	if (!darray_init(&v->intvls, sizeof(double))) {
		iproc_vrecip_free(v);
		v = NULL;
	} else {
		ssize_t i;
		for (i = 0; i < n; i++) {
			assert(intvls[i] > 0.0);
			assert(i == 0 || intvls[i] > intvls[i - 1]);

			darray_push_back(&v->intvls, intvls + i);
		}
	}

	return v;
}

iproc_vrecip *iproc_vrecip_ref(iproc_vrecip * v)
{
	if (v) {
		refcount_get(&v->refcount);
	}

	return v;
}

static void iproc_vrecip_release(struct refcount *refcount)
{
	iproc_vrecip *v = container_of(refcount, iproc_vrecip, refcount);
	iproc_vrecip_free(v);
}

void iproc_vrecip_unref(iproc_vrecip * v)
{
	if (v) {
		refcount_put(&v->refcount, iproc_vrecip_release);
	}
}
