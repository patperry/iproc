#include "port.h"
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include "ieee754.h"
#include "pqueue.h"
#include "xalloc.h"
#include "design2.h"


struct update {
	double t;
	size_t j;
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

struct ndyad_thunk {
	struct pqueue *updates;
	double *intvls;
	size_t nintvl;
};



static void ndyad_init(struct var_meta *meta, void **thunk, const char *name,
		       struct design2 *d, va_list ap)
{
	const double *intvls = va_arg(ap, double *);
	size_t nintvl = va_arg(ap, size_t);
	struct ndyad_thunk *ndyad = xmalloc(sizeof(*ndyad));
	size_t i, m = design2_count1(d);

#ifndef NDEBUG
	size_t k;
	if (nintvl) {
		assert(intvls[0] > 0);
	}
	for (k = 1; k < nintvl; k++) {
		assert(intvls[k-1] < intvls[k]);
	}
#endif

	var_meta_init(meta, VAR_TYPE_TVAR, name, &nintvl, 1);

	ndyad->updates = xmalloc(m * sizeof(struct pqueue));
	for (i = 0; i < m; i++) {
		pqueue_init(&ndyad->updates[i], sizeof(struct update),
			    update_rcompare);
	}

	ndyad->intvls = xmalloc(nintvl * sizeof(double));
	memcpy(ndyad->intvls, intvls, nintvl * sizeof(double));
	ndyad->nintvl = nintvl;

	*thunk = ndyad;
}



static void ndyad_deinit(struct var_meta *meta, void *thunk, struct design2 *d)
{
	struct ndyad_thunk *ndyad = thunk;
	size_t i, m = design2_count1(d);

	free(ndyad->intvls);

	for (i = 0; i < m; i++) {
		pqueue_deinit(&ndyad->updates[i]);
	}
	free(ndyad->updates);

	free(ndyad);
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


static void process_updates(struct tvar2 *tv, size_t i, double t0, double t)
{
	struct design2 *d = tv->var.design;
	struct ndyad_thunk *ndyad = tv->thunk;
	struct pqueue *updates = &ndyad->updates[i];
	const double *intvls = ndyad->intvls;
	size_t nintvl = ndyad->nintvl;
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

		x = design2_make_active(d, tv, i, u->j);
		x[u->intvl] -= u->wt;

		k = find_intvl(t - u->t0, intvls, nintvl, u->intvl + 1);
		design2_update(d, tv, i, u->j, u->t0 + intvls[k - 1]);

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


static void process_messages(struct tvar2 *tv, size_t i, double t0, const struct history *h)
{
	struct design2 *d = tv->var.design;
	struct ndyad_thunk *ndyad = tv->thunk;
	const double *intvls = ndyad->intvls;
	size_t nintvl = ndyad->nintvl;
	struct pqueue *updates = &ndyad->updates[i];
	struct update u;
	double t = history_time(h);
	size_t k;
	double *x;
	double wt;
	double t1;

	const struct history_actor *ha = HISTORY_ACTOR(h, i);
	const size_t *ind;
	size_t iz, nz;

	history_actor_get_msgs(ha, &ind, &nz);

	for (iz = nz; iz > 0; iz--) {
		size_t imsg = ind[iz - 1];
		const struct message *msg = history_item(h, imsg);
		size_t it, j;

		if (msg->time < t0)
			break;

		k = find_intvl(t - msg->time, intvls, nintvl, 0);
		if (k == nintvl)
			break;

		wt = MSG_WEIGHT(msg);
		t1 = msg->time + (k ? intvls[k - 1] : 0);

		FOREACH_ACTOR(it, j, msg) {
			x = design2_make_active(d, tv, i, j);
			x[k] += wt;

			u.t = msg->time + intvls[k];
			if (u.t < INFINITY) {
				u.j = j;
				u.intvl = k;
				u.wt = wt;
				u.t0 = msg->time;
				pqueue_push(updates, &u);
			}

			design2_update(d, tv, i, j, t1);
		}
	}
}


static double ndyad_update(struct tvar2 *tv, size_t i, double t0, const struct history *h)
{
	struct ndyad_thunk *ndyad = tv->thunk;
	double t = history_time(h);
	double tnext = INFINITY;

	process_updates(tv, i, t0, t);
	process_messages(tv, i, t0, h);

	if (pqueue_count(&ndyad->updates[i])) {
		const struct update *u = pqueue_top(&ndyad->updates[i]);
		tnext = u->t;
	}

	return tnext;
}


static struct tvar2_type VAR2_NDYAD_REP = {
	ndyad_init,
	ndyad_deinit,
	ndyad_update
};

