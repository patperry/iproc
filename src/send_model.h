#ifndef SEND_MODEL_H
#define SEND_MODEL_H

#include <stddef.h>
#include "design.h"
#include "mlogit.h"
#include "mlogitaug.h"


struct send_params {
	struct coefs coefs;
};

struct send_model {
	struct design *design;
	struct version_watch history_version;
	struct send_params params;
	struct mlogit_work base_work;
	struct mlogit base;
	struct mlogitaug_work aug_work;
	struct mlogitaug aug;
	double tcur;
};




void send_model_init(struct send_model *m, const struct send_params *p,
		     struct design *d);
void send_model_deinit(struct send_model *m);

struct design *send_model_design(const struct send_model *m);
const struct send_params *send_model_params(const struct send_model *m);
size_t send_model_dim(const struct send_model *m);
size_t send_model_count(const struct send_model *m);


int send_model_moments(const struct send_model *m);
void send_model_set_moments(struct send_model *m, int k);


void send_model_set_params(struct send_model *m, const struct send_params *p);

struct catdist1 *send_model_dist(const struct send_model *m);
struct mlogitaug *send_model_mlogit(const struct send_model *m);


#endif /* SEND_MODEL_H */
