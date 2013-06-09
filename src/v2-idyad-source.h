#include "port.h"
#include <assert.h>
#include <stdlib.h>
#include "ieee754.h"
#include "pqueue.h"
#include "xalloc.h"
#include "design2.h"


struct update {
	double t;
	size_t j;
};

static int update_rcompare(const struct pqueue *q, const void *x1, const void *x2)
{
	(void)q;
	const struct update *u1 = x1;
	const struct update *u2 = x2;
	return double_rcompare(&u1->t, &u2->t);
}

struct idyad_thunk {
	struct pqueue *updates;
	double window;
};


static void idyad_init(struct var_meta *meta, void **thunk, const char *name,
		       struct design2 *d, va_list ap)
{
	double window = va_arg(ap, double);
	struct idyad_thunk *idyad = xmalloc(sizeof(*idyad));
	size_t i, m = design2_count1(d);

	assert(window > 0);

	var_meta_init(meta, VAR_TYPE_TVAR, name, NULL, 0);
	idyad->updates = xmalloc(m * sizeof(struct pqueue));
	for (i = 0; i < m; i++) {
		pqueue_init(&idyad->updates[i], sizeof(struct update),
			    update_rcompare);
	}
	idyad->window = window;
	*thunk = idyad;
}


static void idyad_deinit(struct var_meta *meta, void *thunk, struct design2 *d)
{
	struct idyad_thunk *idyad = thunk;
	size_t i, m = design2_count1(d);

	for (i = 0; i < m; i++) {
		pqueue_deinit(&idyad->updates[i]);
	}
	free(idyad->updates);
	free(idyad);
	var_meta_deinit(meta);
}


static void process_updates(struct tvar2 *tv, size_t i, double t0, double t)
{
	struct design2 *d = tv->var.design;
	struct idyad_thunk *idyad = tv->thunk;
	struct pqueue *updates = &idyad->updates[i];
	struct update *u;
	double *x;

	if (t0 == -INFINITY)
		pqueue_clear(updates);

	/* process all elapsed updates */
	while (pqueue_count(updates)) {
		u = pqueue_top(updates);
		if (u->t >= t)
			break;

		x = design2_make_active(d, tv, i, u->j);
		assert(x[0] == 1.0);
		x[0] = 0.0;
		design2_update(d, tv, i, u->j, u->t);
		pqueue_pop(updates);
	}
}


static void process_messages(struct tvar2 *tv, size_t i, double t0,
			     const struct history *h)
{
	struct design2 *d = tv->var.design;
	struct idyad_thunk *idyad = tv->thunk;
	struct pqueue *updates = &idyad->updates[i];
	struct update u;
	double t = history_time(h);
	double w = idyad->window;
	double tmin = t - idyad->window;
	double *x;
	size_t imsg, nmsg = history_count(h);
	const struct message *msg;
	size_t it, j;


	for (imsg = nmsg; imsg > 0; imsg--) {
		msg = history_item(h, imsg - 1);

		if (msg->time < tmin)
			break;

		FOREACH_ACTOR(it, j, msg) {
			x = design2_make_active(d, tv, i, j);
			if (x[0] == 0.0) {
				x[0] = 1.0;
				u.t = msg->time + w;
				u.j = j;
				design2_update(d, tv, i, j, msg->time);
				pqueue_push(updates, &u);
			}
		}
	}
}


static double idyad_update(struct tvar2 *tv, size_t i, double t0,
			   const struct history *h)
{
	struct idyad_thunk *idyad = tv->thunk;
	double t = history_time(h);
	double tnext = INFINITY;

	process_updates(tv, i, t0, t);
	process_messages(tv, i, t0, h);

	if (pqueue_count(&idyad->updates[i])) {
		const struct update *u = pqueue_top(&idyad->updates[i]);
		tnext = u->t;
	}

	return tnext;
}


static struct tvar2_type VAR2_IDYAD_REP = {
	idyad_init,
	idyad_deinit,
	idyad_update
};

