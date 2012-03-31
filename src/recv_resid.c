#include "port.h"

#include <assert.h>
#include <stdlib.h>
#include "recv_resid.h"

static void recv_resid_count_init(struct recv_resid_count *c,
				  size_t nsend, size_t nrecv)
{
	matrix_init(&c->dyad, nsend, nrecv);
	matrix_init(&c->dyad_trans, nrecv, nsend);
	vector_init(&c->send, nsend);
	vector_init(&c->recv, nrecv);
	c->tot = 0.0;
	c->dyad_cached = true;
}

static void recv_resid_count_deinit(struct recv_resid_count *c)
{
	vector_deinit(&c->recv);
	vector_deinit(&c->send);
	matrix_deinit(&c->dyad_trans);
	matrix_deinit(&c->dyad);
}

static void recv_resid_clear(struct recv_resid_count *c)
{
	matrix_fill(&c->dyad, 0.0);
	matrix_fill(&c->dyad_trans, 0.0);
	vector_fill(&c->send, 0.0);
	vector_fill(&c->recv, 0.0);
	c->tot = 0;
	c->dyad_cached = true;
}

static void recv_resid_count_compute_dyad(struct recv_resid_count *c)
{
	if (!c->dyad_cached) {
		matrix_assign_copy(&c->dyad, BLAS_TRANS, &c->dyad_trans);
		c->dyad_cached = true;
	}
}

static void update_obs(struct recv_resid_count *obs, const struct message *msg)
{
	size_t isend = msg->from;
	size_t ito, nto = msg->nto;

	double *ps = vector_item_ptr(&obs->send, isend);
	*ps += nto;

	for (ito = 0; ito < nto; ito++) {
		size_t jrecv = msg->to[ito];

		double *pd = matrix_item_ptr(&obs->dyad_trans, jrecv, isend);
		*pd += 1.0;

		double *pr = vector_item_ptr(&obs->recv, jrecv);
		*pr += 1.0;
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

	double *ps = vector_item_ptr(&exp->send, isend);
	*ps += nto;

	recv_model_axpy_probs(nto, model, isend, exp->recv.data);

	double *col = matrix_col(&exp->dyad_trans, isend);
	recv_model_axpy_probs(nto, model, isend, col);

	exp->dyad_cached = false;
}

static void recv_resid_set(struct recv_resid *resid,
		           struct frame *f,
			   const struct messages *msgs,
			   const struct matrix *coefs)
{
	recv_model_set_coefs(&resid->model, coefs);
	recv_resid_clear(&resid->obs);
	recv_resid_clear(&resid->exp);

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

	recv_resid_count_compute_dyad(&resid->obs);
	recv_resid_count_compute_dyad(&resid->exp);
}

void recv_resid_init(struct recv_resid *resid,
		     struct frame *f,
		     const struct messages *msgs,
		     size_t ncohort,
		     const size_t *cohorts, const struct matrix *coefs)
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
