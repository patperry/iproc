#ifndef RECV_MODEL_H
#define RECV_MODEL_H

#include <stddef.h>
#include "design.h"
#include "design2.h"
#include "mlogit.h"
#include "mlogitaug.h"


struct recv_params {
	struct coefs recv;
	struct coefs2 dyad;
};

struct recv_model_cohort {
	struct mlogit mlogit;
	struct version_watch version;
	double tcur;
};

struct recv_model_sender {
	struct mlogitaug mlogitaug;
	struct version_watch version;
	double tcur;
};

struct recv_model {
	struct design *recv;
	struct design2 *dyad;
	int exclude_loops;
	struct recv_params params;
	struct recv_model_cohort *cohort_models;
	struct recv_model_sender *sender_models;
	struct mlogit_work work;
	struct mlogitaug_work augwork;
	int moments;
};


void recv_params_init(struct recv_params *p, const struct design *r, const struct design2 *d);
void recv_params_deinit(struct recv_params *p);
void recv_params_set(struct recv_params *dst, const struct recv_params *src, const struct design *r, const struct design2 *d);
void recv_params_axpy(double alpha, const struct recv_params *x, struct recv_params *y, const struct design *r, const struct design2 *d);


void recv_model_init(struct recv_model *m, const struct recv_params *p,
		     struct design *r, struct design2 *d);
void recv_model_deinit(struct recv_model *m);

struct design *recv_model_design(const struct recv_model *m);
struct design2 *recv_model_design2(const struct recv_model *m);
size_t recv_model_send_count(const struct recv_model *model);
size_t recv_model_dim(const struct recv_model *model);
size_t recv_model_count(const struct recv_model *model);

static inline int recv_model_exclude_loops(const struct recv_model *m)
{
	return m->exclude_loops;
}

void recv_model_set_exclude_loops(struct recv_model *m, int exclude_loops);


int recv_model_moments(const struct recv_model *m);
void recv_model_set_moments(struct recv_model *m, int k);

const struct recv_params *recv_model_params(const struct recv_model *m);
void recv_model_set_params(struct recv_model *m, const struct recv_params *p);

struct catdist1 *recv_model_dist(const struct recv_model *m, size_t i);
struct mlogitaug *recv_model_mlogit(const struct recv_model *m, size_t i);


#endif /* RECV_MODEL_H */
