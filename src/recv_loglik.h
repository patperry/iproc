#ifndef RECV_LOGLIK_H
#define RECV_LOGLIK_H

#include <stddef.h>
#include "messages.h"
#include "recv_model.h"


struct recv_loglik_cohort {
};

struct recv_loglik_sender {
	size_t count;
	double dev;
	struct recv_coefs mean;
	struct recv_coefs score;
	double *cov;
};

struct recv_loglik_update {
	size_t count;
	double dev;
	struct recv_coefs mean;
	struct recv_coefs score;
	double *cov;
};

struct recv_loglik {
	struct recv_model *model;
	struct recv_loglik_cohort *cohorts;
	struct recv_loglik_sender *senders;
	struct recv_loglik_update last;
};

void recv_loglik_init(struct recv_loglik *ll, struct recv_model *m);
void recv_loglik_deinit(struct recv_loglik *ll);

static inline struct recv_model *recv_loglik_model(const struct recv_loglik *ll);

void recv_loglik_add(struct recv_loglik *ll,
		     const struct frame *f, const struct message *msg);
void recv_loglik_add_all(struct recv_loglik *ll,
			 struct frame *f, const struct messages *msgs);
void recv_loglik_clear(struct recv_loglik *ll);


size_t recv_loglik_count(const struct recv_loglik *ll);
double recv_loglik_dev(const struct recv_loglik *ll);
void recv_loglik_axpy_mean(double alpha, const struct recv_loglik *ll, struct recv_coefs *y);
void recv_loglik_axpy_score(double alpha, const struct recv_loglik *ll, struct recv_coefs *y);
void recv_loglik_axpy_imat(double alpha, const struct recv_loglik *ll, double *y);


size_t recv_loglik_last_count(const struct recv_loglik *ll);
double recv_loglik_last_dev(const struct recv_loglik *ll);
void recv_loglik_axpy_last_mean(double alpha, const struct recv_loglik *ll, struct recv_coefs *y);
void recv_loglik_axpy_last_score(double alpha, const struct recv_loglik *ll, struct recv_coefs *y);
void recv_loglik_axpy_last_imat(double alpha, const struct recv_loglik *ll, double *y);


/* inline funciton definitions */
struct recv_model *recv_loglik_model(const struct recv_loglik *ll)
{
	return ll->model;
}


#endif /* RECV_LOGLIK_H */
