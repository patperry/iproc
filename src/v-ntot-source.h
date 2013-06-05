#include "port.h"
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include "ieee754.h"
#include "pqueue.h"
#include "xalloc.h"
#include "design.h"


struct update {
	double t;
	size_t i;
	size_t intvl;
	double wt;
	size_t t0;
};

static int update_rcompare(const struct pqueue *q, const void *x1, const void *x2)
{
	(void)q;
	const struct update *u1 = x1;
	const struct update *u2 = x2;
	return double_rcompare(&u1->t, &u2->t);
}

struct ntot_thunk {
	struct pqueue updates;
	double *intvls;
	size_t nintvl;
};



static void ntot_init(struct var_meta *meta, void **thunk, const char *name,
		      struct design *d, va_list ap)
{
	const double *intvls = va_arg(ap, double *);
	size_t nintvl = va_arg(ap, size_t);
	struct ntot_thunk *ntot = xmalloc(sizeof(*ntot));
	size_t i;

	if (nintvl) {
		assert(intvls[0] > 0);
	}

	for (i = 1; i < nintvl; i++) {
		assert(intvls[i-1] < intvls[i]);
	}

	var_meta_init(meta, VAR_TYPE_TVAR, name, &nintvl, 1);
	pqueue_init(&ntot->updates, sizeof(struct update), update_rcompare);

	ntot->intvls = xmalloc(nintvl * sizeof(double));
	memcpy(ntot->intvls, intvls, nintvl * sizeof(double));
	ntot->nintvl = nintvl;

	*thunk = ntot;
}



static void ntot_deinit(struct var_meta *meta, void *thunk, struct design *d)
{
	struct ntot_thunk *ntot = (struct ntot_thunk *)thunk;
	free(ntot->intvls);
	pqueue_deinit(&ntot->updates);
	free(ntot);
	var_meta_deinit(meta);
}


static size_t find_intvl(double dt, const double *intvls, size_t nintvl, size_t kmin)
{
	size_t k;
	assert(dt > 0);
	assert(kmin <= nintvl);

	for (k = kmin; k < nintvl; k++) {
		if (dt <= intvls[k])
			break;
	}

	return k;
}


static void process_updates(struct tvar *tv, double t0, double t)
{
	struct design *d = tv->var.design;
	struct ntot_thunk *ntot = (struct ntot_thunk *)tv->thunk;
	struct pqueue *updates = &ntot->updates;
	const double *intvls = ntot->intvls;
	size_t nintvl = ntot->nintvl;
	struct update *u;
	double *x;
	size_t k;

	if (t0 == -INFINITY)
		pqueue_clear(updates);

	/* process all elapsed updates */
	while (pqueue_count(updates)) {
		u = pqueue_top(updates);
		if (u->t >= t)
			break;

		x = design_make_active(d, tv, u->i);
		x[u->intvl] -= u->wt;

		k = find_intvl(t - u->t0, intvls, nintvl, u->intvl + 1);
		design_update(d, tv, u->i, u->t0 + intvls[k - 1]);

		if (k < nintvl) {
			x[k] += u->wt;
			u->t = u->t0 + intvls[k];
			if (u->t < INFINITY) {
				u->intvl = k;
				pqueue_update_top(updates);

			} else {
				pqueue_pop(updates);
			}
		} else {
			pqueue_pop(updates);
		}
	}
}


static void process_messages(struct tvar *tv, double t0, const struct history *h)
{
	struct design *d = tv->var.design;
	struct ntot_thunk *ntot = (struct ntot_thunk *)tv->thunk;
	const double *intvls = ntot->intvls;
	size_t nintvl = ntot->nintvl;
	struct pqueue *updates = &ntot->updates;
	struct update u;
	double t = history_time(h);

	size_t imsg, nmsg = history_count(h);
	const struct message *msg;
	size_t it, i;
	size_t k;
	double *x;
	double wt;
	double t1;

	for (imsg = nmsg; imsg > 0; imsg--) {
		msg = history_item(h, imsg - 1);

		if (msg->time < t0)
			break;

		k = find_intvl(t - msg->time, intvls, nintvl, 0);
		if (k == nintvl)
			break;

		wt = MSG_WEIGHT(msg);
		t1 = msg->time + (k ? intvls[k - 1] : 0);

		FOREACH_ACTOR(it, i, msg) {
			x = design_make_active(d, tv, i);
			x[k] += wt;

			u.t = msg->time + intvls[k];
			if (u.t < INFINITY) {
				u.i = i;
				u.intvl = k;
				u.wt = wt;
				u.t0 = msg->time;
				pqueue_push(updates, &u);
			}

			design_update(d, tv, i, t1);
		}
	}
}


static void ntot_update(struct tvar *tv, double t0, const struct history *h)
{
	double t = history_time(h);
	process_updates(tv, t0, t);
	process_messages(tv, t0, h);
}


static struct tvar_type VAR_NTOT_REP = {
	ntot_init,
	ntot_deinit,
	ntot_update
};

