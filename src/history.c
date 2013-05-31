#include "port.h"
#include "coreutil.h"
#include "xalloc.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include "history.h"


static void actors_init(struct history_actor **has, size_t n)
{
	struct history_actor *ha;
	size_t i;

	*has = xmalloc(n * sizeof(struct history_actor));
	for (i = 0; i < n; i++) {
		ha = &(*has)[i];
		uintset_init(&ha->msgs);
		uintset_init(&ha->alters);
		ha->wts = NULL;
	}
}


static void actors_deinit(struct history_actor *has, size_t n)
{
	struct history_actor *ha;
	size_t i;

	for (i = 0; i < n; i++) {
		ha = &has[i];
		free(ha->wts);
		uintset_deinit(&ha->alters);
		uintset_deinit(&ha->msgs);
	}
	free(has);
}


static void actors_clear(struct history_actor *has, size_t n)
{
	struct history_actor *ha;
	size_t i;

	for (i = 0; i < n; i++) {
		ha = &has[i];
		uintset_clear(&ha->msgs);
		uintset_clear(&ha->alters);
	}
}


static void actor_add_msg(struct history_actor *ha, size_t i)
{
	size_t n = uintset_count(&ha->msgs);

	/* we always add messages in increasing order of index */
	assert(n == 0 || ha->msgs.vals[n - 1] < i);

	uintset_insert(&ha->msgs, n, i);
}


static void actor_add_alter(struct history_actor *ha, size_t j, double wt)
{
	size_t index;

	assert(!isnan(wt));

	/* find the alter index */
	if (!uintset_find(&ha->alters, j, &index)) {
		size_t n0, n1, ntail;

		/* insert the new alter  */
		n0 = uintset_capacity(&ha->alters);
		ntail = uintset_insert(&ha->alters, index, j);
		n1 = uintset_capacity(&ha->alters);

		/* expand the weights array if necessary */
		if (n1 != n0) {
			ha->wts = xrealloc(ha->wts, n1 * sizeof(double));
		}

		/* initialize the new weight to zero */
		memmove(ha->wts + index + 1, ha->wts + index,
			ntail * sizeof(double));
		ha->wts[index] = 0;
	}

	/* increment the alter weight */
	ha->wts[index] += wt;
}


void history_init(struct history *h, size_t nsend, size_t nrecv)
{
	h->nsend = nsend;
	h->nrecv = nrecv;
	h->msgs = NULL;
	h->nmsg = 0;
	h->nmsg_max = 0;
	h->to = NULL;
	h->nto = 0;
	h->nto_max = 0;
	actors_init(&h->sends, nsend);
	actors_init(&h->recvs, nrecv);
	h->ncur = 0;
	h->tcur = -INFINITY;
}


void history_deinit(struct history *h)
{
	actors_deinit(h->recvs, history_nrecv(h));
	actors_deinit(h->sends, history_nsend(h));
	free(h->to);
	free(h->msgs);
}


void history_clear(struct history *h)
{
	h->nmsg = 0;
	h->nto = 0;
	history_reset(h);
}


void history_reset(struct history *h)
{
	actors_clear(h->sends, history_nsend(h));
	actors_clear(h->recvs, history_nrecv(h));
	h->ncur = 0;
	h->tcur = -INFINITY;
}


void history_advance(struct history *h, double time)
{
	size_t imsg, nmsg = h->nmsg;
	const struct message *msgs = h->msgs;
	const struct message *msg;
	struct history_actor *s, *r;
	size_t ito, nto;
	size_t from, to;
	double wt;

	assert(time >= history_time(h));

	for (imsg = h->ncur, msg = &msgs[imsg];
			imsg < nmsg && msg->time < time; imsg++, msg++) {
		nto = msg->nto;
		wt = 1 / (double)nto;
		from = msg->from;
		s = history_send(h, from);
		actor_add_msg(s, imsg);

		for (ito = 0; ito < nto; ito++) {
			to = msg->to[ito];
			r = history_recv(h, to);
			actor_add_msg(r, imsg);
			actor_add_alter(r, from, wt);
			actor_add_alter(s, to, wt);
		}
	}
	assert((ptrdiff_t)imsg == msg - msgs);
	h->ncur = imsg;
	h->tcur = time;
}


void history_add(struct history *h, size_t from, size_t *to, size_t nto,
		 intptr_t attr)
{
	size_t imsg;
	size_t *to_buf;
	struct message *msg;
	double time = h->tcur;

	assert(history_time(h) < INFINITY);

	/* make space for the new message */
	if (needs_grow(h->nmsg + 1, &h->nmsg_max)) {
		h->msgs = xrealloc(h->msgs,
				   h->nmsg_max * sizeof(struct message));
	}

	/* make space for the new 'to' list */
	if (needs_grow(h->nto + nto, &h->nto_max)) {
		h->to = xrealloc(h->to, h->nto_max * sizeof(size_t));

		to_buf = h->to;
		for (imsg = 0; imsg < h->ncur; imsg++) {
			h->msgs[imsg].to = to_buf;
			to_buf += h->msgs[imsg].nto;
		}
	} else {
		if (h->ncur == h->nmsg) {
			to_buf = h->to + h->nto;
		} else if (h->ncur) {
			struct message *prev = &h->msgs[h->ncur - 1];
			to_buf = prev->to + prev->nto;
		} else {
			to_buf = h->to;
		}
	}

	/* find a space for the new message, update to_buf for current messages */
	for (imsg = h->ncur; imsg < h->nmsg && h->msgs[imsg].time == time; imsg++) {
		h->msgs[imsg].to = to_buf;
		to_buf += h->msgs[imsg].nto;
	}

	/* insert the new message */
	msg = h->msgs + imsg;
	memmove(msg + 1, msg, (h->nmsg - imsg) * sizeof(struct message));
	msg->time = time;
	msg->from = from;
	msg->to = to_buf;
	msg->nto = nto;
	msg->attr = attr;
	h->nmsg++;

	/* insert the new 'to' field */
	memmove(to_buf + nto, to_buf,
		(h->nto - (to_buf - h->to)) * sizeof(size_t));
	memcpy(to_buf, to, nto * sizeof(size_t));
	to_buf += nto;
	h->nto += nto;

	/* update the pointers for the subsequent 'to' fields */
	for (imsg++; imsg < h->nmsg; imsg++) {
		h->msgs[imsg].to = to_buf;
		to_buf += h->msgs[imsg].nto;
	}
	assert(to_buf == h->to + h->nto);
}


