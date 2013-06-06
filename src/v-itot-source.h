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

struct itot_thunk {
	struct pqueue updates;
	double window;
};


static void itot_init(struct var_meta *meta, void **thunk, const char *name,
		struct design *d, va_list ap)
{
	double window = va_arg(ap, double);
	struct itot_thunk *itot = xmalloc(sizeof(*itot));

	assert(window > 0);

	var_meta_init(meta, VAR_TYPE_TVAR, name, NULL, 0);
	pqueue_init(&itot->updates, sizeof(struct update), update_rcompare);
	itot->window = window;
	*thunk = itot;
}


static void itot_deinit(struct var_meta *meta, void *thunk, struct design *d)
{
	struct itot_thunk *itot = (struct itot_thunk *)thunk;
	pqueue_deinit(&itot->updates);
	free(itot);
	var_meta_deinit(meta);
}


static void process_updates(struct tvar *tv, double t0, double t)
{
	struct design *d = tv->var.design;
	struct itot_thunk *itot = (struct itot_thunk *)tv->thunk;
	struct pqueue *updates = &itot->updates;
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


static void process_messages(struct tvar *tv, double t0, const struct history *h)
{
	struct design *d = tv->var.design;
	struct itot_thunk *itot = (struct itot_thunk *)tv->thunk;
	struct pqueue *updates = &itot->updates;
	struct update u;
	double t = history_time(h);
	double w = itot->window;
	double tmin = t - itot->window;
	double *x;
	size_t imsg, nmsg = history_count(h);
	const struct message *msg;
	size_t k, i;


	for (imsg = nmsg; imsg > 0; imsg--) {
		msg = history_item(h, imsg - 1);

		if (msg->time < tmin)
			break;

		FOREACH_ACTOR(k, i, msg) {
			x = design_make_active(d, tv, i);
			if (x[0] == 0.0) {
				x[0] = 1.0;
				u.t = msg->time + w;
				u.i = i;
				design_update(d, tv, i, msg->time);
				pqueue_push(updates, &u);
			}
		}
	}

}


static double itot_update(struct tvar *tv, double t0, const struct history *h)
{
	struct itot_thunk *itot = (struct itot_thunk *)tv->thunk;
	double t = history_time(h);
	double tnext = INFINITY;

	process_updates(tv, t0, t);
	process_messages(tv, t0, h);

	if (pqueue_count(&itot->updates)) {
		const struct update *u = pqueue_top(&itot->updates);
		tnext = u->t;
	}

	return tnext;
}


static struct tvar_type VAR_ITOT_REP = {
	itot_init,
	itot_deinit,
	itot_update
};

