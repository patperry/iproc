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
	size_t j, n = recv_model_count(m);

	for (j = 0; j < n; j++) {
		double p = catdist1_prob(dist, j);
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


static void update_fit(struct recv_resid_count *fit, const struct message *msg,
		       const struct recv_model *model)
{
	size_t nrecv = recv_model_count(model);
	size_t isend = msg->from;
	size_t nto = msg->nto;

	fit->tot += nto;

	fit->send[isend] += (double)nto;

	axpy_probs(nto, model, isend, fit->recv);

	double *row = fit->dyad + isend * nrecv;
	axpy_probs(nto, model, isend, row);
}


void recv_resid_set(struct recv_resid *resid,
		    struct recv_model *model,
		    const struct message *msgs,
		    size_t nmsg)
{
	size_t nsend = recv_model_send_count(model);
	size_t nrecv = recv_model_count(model);
	const struct design *r = recv_model_design(model);
	const struct design2 *d = recv_model_design2(model);
	struct history *hr = design_history(r);
	struct history *hd = design2_history(d);
	size_t imsg;

	recv_resid_clear(&resid->obs, nsend, nrecv);
	recv_resid_clear(&resid->fit, nsend, nrecv);

	for (imsg = 0; imsg < nmsg; imsg++) {
		const struct message *msg = &msgs[imsg];
		double t = msg->time;

		if (history_time(hr) < t) {
			history_reset(hr);
		}

		history_advance(hr, t);
		if (history_time(hd) < t) {
			history_reset(hd);
		}
		history_advance(hd, t);

		update_obs(&resid->obs, msg, nrecv);
		update_fit(&resid->fit, msg, model);
	}
}


void recv_resid_init(struct recv_resid *resid,
		     struct recv_model *model,
		     const struct message *msgs,
		     size_t nmsg)
{
	size_t nsend = recv_model_send_count(model);
	size_t nrecv = recv_model_count(model);

	recv_resid_count_init(&resid->obs, nsend, nrecv);
	recv_resid_count_init(&resid->fit, nsend, nrecv);

	recv_resid_set(resid, model, msgs, nmsg);
}

void recv_resid_deinit(struct recv_resid *resid)
{
	recv_resid_count_deinit(&resid->fit);
	recv_resid_count_deinit(&resid->obs);
}
