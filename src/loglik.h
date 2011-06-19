#ifndef _RECV_LOGLIK_H
#define _RECV_LOGLIK_H

#include "array.h"
#include "frame.h"
#include "messages.h"
#include "model.h"
#include "vector.h"

struct recv_sloglik;

struct recv_loglik {
	struct model *model;
	struct array slogliks;
	struct vector grad;
	bool grad_cached;
	ssize_t nsend;
	ssize_t nrecv;
	
	struct recv_sloglik *last;
};

void recv_loglik_init(struct recv_loglik *ll,
		      struct model *m);
void recv_loglik_deinit(struct recv_loglik *ll);

void recv_loglik_update(struct model *m, const struct frame *f);

struct recv_loglik *recv_loglik_alloc(struct model *m, struct messages *messages);
void recv_loglik_free(struct recv_loglik *ll);

void recv_loglik_add(struct recv_loglik *ll,
		     const struct frame *f,
		     const struct message *msg);

void recv_loglik_add_all(struct recv_loglik *ll,
			 struct frame *f,
			 const struct messages *msgs);


double recv_loglik_value(const struct recv_loglik *ll);
void recv_loglik_axpy_grad(double alpha,
			   const struct recv_loglik * ll,
			   struct vector *y);

ssize_t recv_loglik_count(const struct recv_loglik *ll);

double recv_loglik_avg_dev(const struct recv_loglik *ll);
double recv_loglik_last_dev(const struct recv_loglik *ll);

void recv_loglik_axpy_avg_mean(double alpha, const struct recv_loglik *ll, struct vector *y);
void recv_loglik_axpy_last_mean(double alpha, const struct recv_loglik *ll, struct vector *y);

void recv_loglik_axpy_avg_imat(double alpha, const struct recv_loglik *ll, struct matrix *y);
void recv_loglik_axpy_last_imat(double alpha, const struct recv_loglik *ll, struct matrix *y);


#endif /* _RECV_LOGLIK_H */
