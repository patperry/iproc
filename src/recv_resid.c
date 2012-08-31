#include "port.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include "coreutil.h"
#include "matrixutil.h"
#include "xalloc.h"
#include "recv_resid.h"

static void recv_resid_count_init(struct recv_resid_count *c,
				  size_t nsend, size_t nrecv)
{
	c->dyad = xcalloc(nsend * nrecv, sizeof(double));
	c->send = xcalloc(nsend, sizeof(double));
	c->recv = xcalloc(nrecv, sizeof(double));
	c->tot = 0.0;
}

static void recv_resid_count_deinit(struct recv_resid_count *c)
{
	free(c->recv);
	free(c->send);
	free(c->dyad);
}

static void recv_resid_clear(struct recv_resid_count *c, size_t nsend,
			     size_t nrecv)
{
	memset(c->dyad, 0, nsend * nrecv * sizeof(double));
	memset(c->send, 0, nsend * sizeof(double));
	memset(c->recv, 0, nrecv * sizeof(double));
	c->tot = 0;
}

static void axpy_probs(double alpha, const struct recv_model *m, size_t isend,
		       double *y)
{
	const struct catdist1 *dist = recv_model_dist(m, isend);
	const struct frame *f = recv_model_frame(m);
	size_t j, n = frame_recv_count(f);

	for (j = 0; j < n; j++) {
		double p = catdist1_cached_prob(dist, j);
		y[j] += alpha * p;
	}
}

static void update_obs(struct recv_resid_count *obs, const struct message *msg, size_t nrecv)
{
	size_t isend = msg->from;
	size_t ito, nto = msg->nto;

	obs->send[isend] += (double)nto;

	for (ito = 0; ito < nto; ito++) {
		size_t jrecv = msg->to[ito];

		obs->dyad[isend * nrecv + jrecv] += 1.0;
		obs->recv[jrecv] += 1.0;
	}

	obs->tot += nto;
}

static void update_exp(struct recv_resid_count *exp, const struct message *msg,
		       const struct recv_model *model)
{
	const struct frame *f = recv_model_frame(model);
	size_t nrecv = frame_recv_count(f);
	size_t isend = msg->from;
	size_t nto = msg->nto;

	exp->tot += nto;

	exp->send[isend] += (double)nto;

	axpy_probs(nto, model, isend, exp->recv);

	double *row = exp->dyad + isend * nrecv;
	axpy_probs(nto, model, isend, row);
}

static void recv_resid_set(struct recv_resid *resid,
		           struct frame *f,
			   const struct messages *msgs,
			   const struct recv_coefs *coefs)
{
	size_t nsend = frame_send_count(f);
	size_t nrecv = frame_recv_count(f);

	recv_model_set_coefs(&resid->model, coefs);
	recv_resid_clear(&resid->obs, nsend, nrecv);
	recv_resid_clear(&resid->exp, nsend, nrecv);

	struct messages_iter it;
	MESSAGES_FOREACH(it, msgs) {
		size_t i, n = MESSAGES_COUNT(it);
		double t = MESSAGES_TIME(it);

		frame_advance(f, t);
		for (i = 0; i < n; i++) {
			struct message *msg = MESSAGES_VAL(it, i);
			frame_add(f, msg);

			update_obs(&resid->obs, msg, nrecv);
			update_exp(&resid->exp, msg, &resid->model);
		}
	}
}

void recv_resid_init(struct recv_resid *resid,
		     struct frame *f,
		     const struct messages *msgs,
		     const struct recv_coefs *coefs)
{
	size_t nsend = frame_send_count(f);
	size_t nrecv = frame_recv_count(f);

	recv_model_init(&resid->model, f, coefs);
	recv_resid_count_init(&resid->obs, nsend, nrecv);
	recv_resid_count_init(&resid->exp, nsend, nrecv);

	recv_resid_set(resid, f, msgs, coefs);
}

void recv_resid_deinit(struct recv_resid *resid)
{
	recv_resid_count_deinit(&resid->exp);
	recv_resid_count_deinit(&resid->obs);
	recv_model_deinit(&resid->model);
}
