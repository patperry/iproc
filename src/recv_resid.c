#include "port.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include "matrixutil.h"
#include "xalloc.h"
#include "recv_resid.h"

static void recv_resid_count_init(struct recv_resid_count *c,
				  size_t nsend, size_t nrecv)
{
	c->dyad = (struct dmatrix) { xcalloc(nsend * nrecv, sizeof(double)), MAX(1, nsend) };
	c->dyad_trans = (struct dmatrix) { xcalloc(nrecv * nsend, sizeof(double)), MAX(1, nrecv) };
	c->send = xcalloc(nsend, sizeof(double));
	c->recv = xcalloc(nrecv, sizeof(double));
	c->tot = 0.0;
	c->dyad_cached = true;
}

static void recv_resid_count_deinit(struct recv_resid_count *c)
{
	free(c->recv);
	free(c->send);
	free(c->dyad_trans.data);
	free(c->dyad.data);
}

static void recv_resid_clear(struct recv_resid_count *c, size_t nsend,
			     size_t nrecv)
{
	matrix_dzero(nsend, nrecv, &c->dyad);
	matrix_dzero(nrecv, nsend, &c->dyad_trans);
	memset(c->send, 0, nsend * sizeof(double));
	memset(c->recv, 0, nrecv * sizeof(double));
	c->tot = 0;
	c->dyad_cached = true;
}

static void recv_resid_count_compute_dyad(struct recv_resid_count *c,
					  size_t nsend, size_t nrecv)
{
	if (!c->dyad_cached) {
		matrix_dtrans(nrecv, nsend, &c->dyad_trans, &c->dyad);
		c->dyad_cached = true;
	}
}

static void update_obs(struct recv_resid_count *obs, const struct message *msg)
{
	size_t isend = msg->from;
	size_t ito, nto = msg->nto;

	obs->send[isend] += (double)nto;

	for (ito = 0; ito < nto; ito++) {
		size_t jrecv = msg->to[ito];

		MATRIX_ITEM(&obs->dyad_trans, jrecv, isend) += 1.0;
		obs->recv[jrecv] += 1.0;
	}

	obs->tot += nto;
	obs->dyad_cached = false;
}

static void update_exp(struct recv_resid_count *exp, const struct message *msg,
		       const struct recv_model *model)
{
	size_t isend = msg->from;
	size_t nto = msg->nto;

	exp->tot += nto;

	exp->send[isend] += (double)nto;

	recv_model_axpy_probs(nto, model, isend, exp->recv);

	double *col = MATRIX_COL(&exp->dyad_trans, isend);
	recv_model_axpy_probs(nto, model, isend, col);

	exp->dyad_cached = false;
}

static void recv_resid_set(struct recv_resid *resid,
		           struct frame *f,
			   const struct messages *msgs,
			   const struct dmatrix *coefs)
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

			update_obs(&resid->obs, msg);
			update_exp(&resid->exp, msg, &resid->model);
		}
	}

	recv_resid_count_compute_dyad(&resid->obs, nsend, nrecv);
	recv_resid_count_compute_dyad(&resid->exp, nsend, nrecv);
}

void recv_resid_init(struct recv_resid *resid,
		     struct frame *f,
		     const struct messages *msgs,
		     size_t ncohort,
		     const size_t *cohorts, const struct dmatrix *coefs)
{
	size_t nsend = frame_send_count(f);
	size_t nrecv = frame_recv_count(f);

	recv_model_init(&resid->model, f, ncohort, cohorts, coefs);
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
