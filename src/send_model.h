#ifndef SEND_MODEL_H
#define SEND_MODEL_H

#include <stddef.h>
#include "design.h"
#include "mlogit.h"
#include "mlogitaug.h"


struct send_params {
	double *all;
	size_t dim;
	struct coefs coefs;
	int owner;
};

struct send_model {
	struct design *design;
	struct send_params params;
	struct mlogit stat;
	struct mlogitaug dyn;
	struct mlogit_work work;
	struct mlogitaug_work augwork;
	int moments;
};


size_t send_params_dim(const struct design *d);
void send_params_init(struct send_params *p, const struct design *d);
void send_params_init_view(struct send_params *p, const struct design *d,
			   const double *data);
void send_params_deinit(struct send_params *p);


void send_model_init(struct send_model *m, struct design *d,
		     const struct send_params *p);
void send_model_deinit(struct send_model *model);

struct design *send_model_design(const struct send_model *model);
const struct send_params *send_model_params(const struct send_model *model);
size_t send_model_dim(const struct send_model *model);


int send_model_moments(const struct send_model *m);
void send_model_set_moments(struct send_model *m, int k);


void send_model_set_params(struct send_model *m, const struct send_params *p);

struct catdist1 *send_model_dist(const struct send_model *m, size_t i);
struct mlogitaug *send_model_mlogit(const struct send_model *m, size_t i);


#endif /* SEND_MODEL_H */
