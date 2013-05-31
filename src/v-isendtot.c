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

struct isendtot_thunk {
	struct pqueue updates;
	double window;
};


static void isendtot_init(struct var_meta *meta, void **thunk, const char *name,
			  struct design *d, va_list ap)
{
	double window = va_arg(ap, double);
	struct isendtot_thunk *isendtot = xmalloc(sizeof(*isendtot));

	assert(window > 0);

	var_meta_init(meta, VAR_TYPE_TVAR, name, NULL, 0);
	pqueue_init(&isendtot->updates, sizeof(struct update), update_rcompare);
	isendtot->window = window;
	*thunk = isendtot;
}


static void isendtot_deinit(struct var_meta *meta, void *thunk,
			    struct design *d)
{
	struct isendtot_thunk *isendtot = (struct isendtot_thunk *)thunk;
	pqueue_deinit(&isendtot->updates);
	free(isendtot);
	var_meta_deinit(meta);
}


static void process_updates(struct tvar *tv, double t0, double t)
{
	struct design *d = tv->var.design;
	struct isendtot_thunk *isendtot = (struct isendtot_thunk *)tv->thunk;
	struct pqueue *updates = &isendtot->updates;
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

static void isendtot_update(struct tvar *tv, double t0, const struct history *h)
{
	struct design *d = tv->var.design;
	struct isendtot_thunk *isendtot = (struct isendtot_thunk *)tv->thunk;
	struct pqueue *updates = &isendtot->updates;
	struct update u;
	double t = history_time(h);
	double w = isendtot->window;
	double tmin = t - isendtot->window;
	double *x;

	process_updates(tv, t0, t);

	/* process all recent messages */
	size_t imsg, nmsg = history_count(h);
	const struct message *msg;
	size_t i;

	for (imsg = nmsg; imsg > 0; imsg--) {
		msg = history_item(h, imsg - 1);

		if (msg->time < tmin)
			break;

		i = msg->from;
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


static struct tvar_type VAR_ISENDTOT_REP = {
	isendtot_init,
	isendtot_deinit,
	isendtot_update
};


const struct tvar_type *VAR_ISENDTOT = &VAR_ISENDTOT_REP;
