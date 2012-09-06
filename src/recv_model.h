#ifndef RECV_MODEL_H
#define RECV_MODEL_H

#include <stddef.h>
#include "design.h"
#include "frame.h"
#include "mlogit.h"
#include "mlogitaug.h"

/* Two senders, i1 and i2, are in the same cohort if and only if their
 * covariates agree at time 0, i.e.
 *
 *   x[0,i1,j] = x[0,i2,j] for all j.
 *
 * Owing to the convention that J[0,i] includes all receivers,
 * all senders in the same group share the same values of
 * w[0,i,j], W[0,i,j], p[0,i,j], and xbar[0,i].
 */

struct recv_coefs {
	double *all;
	size_t dim;
	struct coefs recv;
	struct coefs2 dyad;
};

struct recv_model_cohort {
	struct mlogit mlogit;
};

struct recv_model_sender {
	struct mlogitaug mlogitaug;
};

struct recv_model {
	struct frame *frame;
	struct recv_coefs coefs;
	struct recv_model_cohort *cohort_models;
	struct recv_model_sender *sender_models;
	size_t *ind_buf;
	struct mlogitaug_work work;
	int moments;
};


void recv_coefs_init(struct recv_coefs *c, const struct frame *f);
void recv_coefs_deinit(struct recv_coefs *c);

void recv_model_init(struct recv_model *model,
		     struct frame *f,
		     const struct recv_coefs *coefs);
void recv_model_deinit(struct recv_model *model);

struct frame *recv_model_frame(const struct recv_model *model);
const size_t *recv_model_cohorts(const struct recv_model *model);
const struct recv_coefs *recv_model_coefs(const struct recv_model *model);
size_t recv_model_send_count(const struct recv_model *model);
size_t recv_model_cohort_count(const struct recv_model *model);
size_t recv_model_cohort(const struct recv_model *model, size_t isend);
size_t recv_model_count(const struct recv_model *model);
size_t recv_model_dim(const struct recv_model *model);


int recv_model_moments(const struct recv_model *m);
void recv_model_set_moments(struct recv_model *m, int k);


void recv_model_set_coefs(struct recv_model *m, const struct recv_coefs *coefs);

struct catdist1 *recv_model_dist(const struct recv_model *m, size_t isend);
struct mlogitaug *recv_model_mlogit(const struct recv_model *m, size_t isend);


#endif /* RECV_MODEL_H */
