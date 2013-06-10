#include "port.h"
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include "coreutil.h"
#include "ieee754.h"
#include "pqueue.h"
#include "xalloc.h"
#include "design2.h"


struct update {
	double t;
	size_t j;
	size_t intvl1, intvl2;
	double wt;
	size_t tmsg1, tmsg2;
};

static int update_rcompare(const struct pqueue *q, const void *x1, const void *x2)
{
	(void)q;
	const struct update *u1 = x1;
	const struct update *u2 = x2;
	return double_rcompare(&u1->t, &u2->t);
}

struct ntriad_thunk {
	struct pqueue *updates;
	double *intvls1, *intvls2;
	size_t nintvl1, nintvl2;
};



static void ntriad_init(struct var_meta *meta, void **thunk, const char *name,
		       struct design2 *d, va_list ap)
{
	const double *intvls1 = va_arg(ap, double *);
	size_t nintvl1 = va_arg(ap, size_t);
	const double *intvls2 = va_arg(ap, double *);
	size_t nintvl2 = va_arg(ap, size_t);
	struct ntriad_thunk *ntriad = xmalloc(sizeof(*ntriad));
	size_t i, m = design2_count1(d);

#ifndef NDEBUG
	size_t k;
	if (nintvl1) {
		assert(intvls1[0] > 0);
	}
	for (k = 1; k < nintvl1; k++) {
		assert(intvls1[k-1] < intvls1[k]);
	}
	if (nintvl2) {
		assert(intvls2[0] > 0);
	}
	for (k = 1; k < nintvl2; k++) {
		assert(intvls2[k-1] < intvls2[k]);
	}
#endif

	size_t dims[2] = { nintvl1, nintvl2 };
	var_meta_init(meta, VAR_TYPE_TVAR, name, dims, 2);

	ntriad->updates = xmalloc(m * sizeof(struct pqueue));
	for (i = 0; i < m; i++) {
		pqueue_init(&ntriad->updates[i], sizeof(struct update),
			    update_rcompare);
	}

	ntriad->intvls1 = xmalloc(nintvl1 * sizeof(double));
	memcpy(ntriad->intvls1, intvls1, nintvl1 * sizeof(double));
	ntriad->nintvl1 = nintvl1;

	ntriad->intvls2 = xmalloc(nintvl2 * sizeof(double));
	memcpy(ntriad->intvls2, intvls2, nintvl2 * sizeof(double));
	ntriad->nintvl2 = nintvl2;

	*thunk = ntriad;
}


static void ntriad_deinit(struct var_meta *meta, void *thunk, struct design2 *d)
{
	struct ntriad_thunk *ntriad = thunk;
	size_t i, m = design2_count1(d);

	free(ntriad->intvls2);
	free(ntriad->intvls1);

	for (i = 0; i < m; i++) {
		pqueue_deinit(&ntriad->updates[i]);
	}
	free(ntriad->updates);

	free(ntriad);
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
	struct ntriad_thunk *ntriad = tv->thunk;
	struct pqueue *updates = &ntriad->updates[i];
	const double *intvls1 = ntriad->intvls1;
	const double *intvls2 = ntriad->intvls2;
	size_t nintvl1 = ntriad->nintvl1;
	size_t nintvl2 = ntriad->nintvl2;
	struct update *u;
	double *x;
	size_t k, l;

	if (t0 == -INFINITY)
		pqueue_clear(updates);

	/* process all elapsed updates */
	while (pqueue_count(updates)) {
		u = pqueue_top(updates);
		if (u->t >= t)
			break;

		x = design2_make_active(d, tv, i, u->j);
		x[u->intvl1 * nintvl2 + u->intvl2] -= u->wt;

		k = find_intvl(t - u->tmsg1, intvls1, nintvl1, u->intvl1 + 1);
		l = find_intvl(t - u->tmsg2, intvls2, nintvl2, u->intvl2 + 1);
		design2_update(d, tv, i, u->j,
			       MAX(u->tmsg1 + intvls1[k - 1],
				   u->tmsg2 + intvls2[l - 1]));

		if (k < nintvl1 && l < nintvl2) {
			x[k * nintvl2 + l] += u->wt;
			u->t = MIN(u->tmsg1 + intvls1[k],
				   u->tmsg2 + intvls2[l]);
			if (u->t < INFINITY) {
				u->intvl1 = k;
				u->intvl2 = l;
				pqueue_update_top(updates);

			} else {
				pqueue_pop(updates);
			}
		} else {
			pqueue_pop(updates);
		}
	}
}


static void process_messages(struct tvar2 *tv, size_t i, double t0,
			     const struct history *history)
{
	struct design2 *d = tv->var.design;
	struct ntriad_thunk *ntriad = tv->thunk;
	const double *intvls1 = ntriad->intvls1;
	const double *intvls2 = ntriad->intvls1;
	size_t nintvl1 = ntriad->nintvl1;
	size_t nintvl2 = ntriad->nintvl2;
	struct pqueue *updates = &ntriad->updates[i];
	struct update u;
	double t = history_time(history);
	size_t k, l;
	double *x;
	double wt1, wt2, wt;
	double tnext1, tnext2, tnext;

	const struct history_actor *ha1 = HISTORY_ACTOR1(history, i);
	const size_t *ind1;
	size_t iz1, nz1;
	history_actor_get_msgs(ha1, &ind1, &nz1);

	for (iz1 = nz1; iz1 > 0; iz1--) {
		size_t imsg1 = ind1[iz1 - 1];
		const struct message *msg1 = history_item(history, imsg1);
		size_t it1, h;

		if (msg1->time < t0)
			break;

		k = find_intvl(t - msg1->time, intvls1, nintvl1, 0);
		if (k == nintvl1)
			break;

		wt1 = MSG_WEIGHT1(msg1);
		tnext1 = msg1->time + (k ? intvls1[k - 1] : 0);

		FOREACH_ACTOR1(it1, h, msg1) {
			const struct history_actor *ha2 = HISTORY_ACTOR2(history, h);
			const size_t *ind2;
			size_t iz2, nz2;
			history_actor_get_msgs(ha2, &ind2, &nz2);

			for (iz2 = nz2; iz2 > 0; iz2--) {
				size_t imsg2 = ind2[iz2 - 1];
				const struct message *msg2 = history_item(history, imsg2);
				size_t it2, j;

				if (msg2->time < t0)
					break;

				l = find_intvl(t - msg2->time, intvls2, nintvl2, 0);
				if (l == nintvl2)
					break;

				wt2 = MSG_WEIGHT2(msg2);
				tnext2 = msg2->time + (l ? intvls2[l - 1] : 0);

				wt = wt1 * wt2;
				tnext = MAX(tnext1, tnext2);

				FOREACH_ACTOR2(it2, j, msg2) {
					x = design2_make_active(d, tv, i, j);
					x[k * nintvl2 + l] += wt;

					u.t = MIN(msg1->time + intvls1[k],
						  msg2->time + intvls2[l]);
					if (u.t < INFINITY) {
						u.j = j;
						u.intvl1 = k;
						u.intvl2 = l;
						u.wt = wt;
						u.tmsg1 = msg1->time;
						u.tmsg2 = msg2->time;
						pqueue_push(updates, &u);
					}

					design2_update(d, tv, i, j, tnext);
				}
			}

		}
	}
}


static double ntriad_update(struct tvar2 *tv, size_t i, double t0, const struct history *h)
{
	struct ntriad_thunk *ntriad = tv->thunk;
	double t = history_time(h);
	double tnext = INFINITY;

	process_updates(tv, i, t0, t);
	process_messages(tv, i, t0, h);

	if (pqueue_count(&ntriad->updates[i])) {
		const struct update *u = pqueue_top(&ntriad->updates[i]);
		tnext = u->t;
	}

	return tnext;
}


static struct tvar2_type VAR2_NTRIAD_REP = {
	ntriad_init,
	ntriad_deinit,
	ntriad_update
};

