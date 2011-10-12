#include "port.h"

#include <assert.h>
#include <stdlib.h>
#include "recv_resid.h"

static void recv_resid_count_init(struct recv_resid_count *c,
				  ssize_t nsend, ssize_t nrecv)
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
	ssize_t isend = msg->from;
	ssize_t ito, nto = msg->nto;

	double *ps = vector_item_ptr(&obs->send, isend);
	*ps += nto;

	for (ito = 0; ito < nto; ito++) {
		ssize_t jrecv = msg->to[ito];

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
	ssize_t isend = msg->from;
	ssize_t nto = msg->nto;

	exp->tot += nto;

	double *ps = vector_item_ptr(&exp->send, isend);
	*ps += nto;

	recv_model_axpy_probs(nto, model, isend, &exp->recv);

	struct vector col = matrix_col(&exp->dyad_trans, isend);
	recv_model_axpy_probs(nto, model, isend, &col);

	exp->dyad_cached = false;
}

static void recv_resid_set(struct recv_resid *resid,
			   const struct messages *msgs,
			   const struct matrix *coefs)
{
	frame_clear(&resid->frame);
	recv_model_set_coefs(&resid->model, coefs);
	recv_resid_clear(&resid->obs);
	recv_resid_clear(&resid->exp);

	struct messages_iter it;
	MESSAGES_FOREACH(it, msgs) {
		ssize_t i, n = MESSAGES_COUNT(it);
		double t = MESSAGES_TIME(it);

		frame_advance(&resid->frame, t);
		for (i = 0; i < n; i++) {
			struct message *msg = MESSAGES_VAL(it, i);
			frame_add(&resid->frame, msg);

			update_obs(&resid->obs, msg);
			update_exp(&resid->exp, msg, &resid->model);
		}
	}

	recv_resid_count_compute_dyad(&resid->obs);
	recv_resid_count_compute_dyad(&resid->exp);
}

void recv_resid_init(struct recv_resid *resid,
		     const struct messages *msgs,
		     const struct design *design,
		     size_t ncohort,
		     const ptrdiff_t *cohorts, const struct matrix *coefs)
{
	frame_init(&resid->frame, design);
	recv_model_init(&resid->model, &resid->frame, ncohort, cohorts, coefs);

	ssize_t nsend = design_send_count(design);
	ssize_t nrecv = design_recv_count(design);

	recv_resid_count_init(&resid->obs, nsend, nrecv);
	recv_resid_count_init(&resid->exp, nsend, nrecv);

	recv_resid_set(resid, msgs, coefs);
}

void recv_resid_deinit(struct recv_resid *resid)
{
	recv_resid_count_deinit(&resid->exp);
	recv_resid_count_deinit(&resid->obs);
	recv_model_deinit(&resid->model);
	frame_deinit(&resid->frame);
}
