#include "port.h"
#include <assert.h>
#include <stdlib.h>
#include "ieee754.h"
#include "pqueue.h"
#include "xalloc.h"
#include "design.h"


struct update {
	double t;
	size_t i;
};

static int update_rcompare(const struct pqueue *q, const void *x1, const void *x2)
{
	(void)q;
	const struct update *u1 = x1;
	const struct update *u2 = x2;
	return double_rcompare(&u1->t, &u2->t);
}

struct irecvtot_thunk {
	struct pqueue updates;
	double window;
};


static void irecvtot_init(struct var_meta *meta, void **thunk, const char *name,
			  struct design *d, va_list ap)
{
	double window = va_arg(ap, double);
	struct irecvtot_thunk *irecvtot = xmalloc(sizeof(*irecvtot));

	assert(window > 0);

	var_meta_init(meta, VAR_TYPE_TVAR, name, NULL, 0);
	pqueue_init(&irecvtot->updates, sizeof(struct update), update_rcompare);
	irecvtot->window = window;
	*thunk = irecvtot;
}


static void irecvtot_deinit(struct var_meta *meta, void *thunk,
			    struct design *d)
{
	struct irecvtot_thunk *irecvtot = (struct irecvtot_thunk *)thunk;
	pqueue_deinit(&irecvtot->updates);
	free(irecvtot);
	var_meta_deinit(meta);
}


static void process_updates(struct tvar *tv, double t0, double t)
{
	struct design *d = tv->var.design;
	struct irecvtot_thunk *irecvtot = (struct irecvtot_thunk *)tv->thunk;
	struct pqueue *updates = &irecvtot->updates;
	struct update *u;
	double *x;

	if (t0 == -INFINITY)
		pqueue_clear(updates);

	/* process all elapsed updates */
	while (pqueue_count(updates)) {
		u = pqueue_top(updates);
		if (u->t >= t)
			break;

		x = design_make_active(d, tv, u->i);
		assert(x[0] == 1.0);
		x[0] = 0.0;
		design_update(d, tv, u->i, u->t);
		pqueue_pop(updates);
	}
}

static void irecvtot_update(struct tvar *tv, double t0, const struct history *h)
{
	struct design *d = tv->var.design;
	struct irecvtot_thunk *irecvtot = (struct irecvtot_thunk *)tv->thunk;
	struct pqueue *updates = &irecvtot->updates;
	struct update u;
	double t = history_time(h);
	double w = irecvtot->window;
	double tmin = t - irecvtot->window;
	double *x;

	
	process_updates(tv, t0, t);

	/* process all recent messages */
	size_t imsg, nmsg = history_count(h);
	size_t ito, nto;
	const struct message *msg;

	for (imsg = nmsg; imsg > 0; imsg--) {
		msg = history_item(h, imsg - 1);

		if (msg->time < tmin)
			break;

		nto = msg->nto;
		for (ito = 0; ito < nto; ito++) {
			size_t i = msg->to[ito];
			x = design_make_active(d, tv, i);
			if (x[0] == 0.0) {
				x[1] = 1.0;
				u.t = msg->time + w;
				u.i = i;
				design_update(d, tv, i, msg->time);
				pqueue_push(updates, &u);
			}
		}
	}
}


static struct tvar_type VAR_IRECVTOT_REP = {
	irecvtot_init,
	irecvtot_deinit,
	irecvtot_update
};


const struct tvar_type *VAR_IRECVTOT = &VAR_IRECVTOT_REP;
