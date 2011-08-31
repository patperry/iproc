#include "port.h"

#include <assert.h>
#include <stdlib.h>
#include "vector.h"
#include "recv_boot.h"


static ssize_t sample1(const double *probs, ssize_t n, dsfmt_t *dsfmt)
{
	ssize_t i;
	double u = dsfmt_genrand_close_open(dsfmt);
	double sum = 0;
	
	for (i = 0; i < n - 1; i++) {
		sum += probs[i];
		if (sum >= u)
			break;
	}

	return i;
}

static bool unique(ssize_t *vals, ssize_t len)
{
	ssize_t i, j;
	
	if (len == 1)
		return true;
	if (len == 2)
		return vals[0] != vals[1];
	
	for (i = 0; i < len; i++) {
		ssize_t val = vals[i];
		
		for (j = i + 1; j < len; j++) {
			if (val == vals[j])
				return false;
		}
	}
	return true;
}


static bool sample_subset(const double *probs, ssize_t n, dsfmt_t *dsfmt,
			  ssize_t maxntry,
			  ssize_t *out, ssize_t nout)
{
	ssize_t itry;
	ssize_t i;
	
	for (itry = 0; itry < maxntry; itry++) {
		for (i = 0; i < nout; i++) {
			out[i] = sample1(probs, n, dsfmt);
		}
		if (unique(out, nout))
			return true;
	}
	return false;
}


void recv_boot_init(struct recv_boot *boot,
		    const struct messages *msgs,
		    const struct design *design,
		    const struct actors *senders,
		    const struct matrix *coefs,
		    dsfmt_t *dsfmt)
{
	messages_init(&boot->messages);	
	frame_init(&boot->frame, design);
	recv_model_init(&boot->model, &boot->frame, senders, coefs);
	ssize_t nrecv = design_recv_count(design);
	
	ssize_t max_nto = messages_max_nto(msgs);
	ssize_t *to = xcalloc(max_nto, sizeof(to[0]));
	
	struct vector probs;

	struct messages_iter it;
	const struct message *msg;
	double t;
	ssize_t i, n;
	ssize_t from, nto;
	
	vector_init(&probs, nrecv);
	double *p = vector_to_ptr(&probs);	
	
	ssize_t maxntry = 100000000;
	
	MESSAGES_FOREACH(it, msgs) {
		t = MESSAGES_TIME(it);
		n = MESSAGES_COUNT(it);
		
		frame_advance(&boot->frame, t);
		
		for (i = 0; i < n; i++) {
			msg = MESSAGES_VAL(it, i);
			from = msg->from;
			nto = msg->nto;
			
			vector_fill(&probs, 0.0);			
			recv_model_axpy_probs(1.0, &boot->model, from, &probs);			
			
			if (!sample_subset(p, nrecv, dsfmt, maxntry, to, nto)) {
				fprintf(stdout, "Failed to sample subset of size %"SSIZE_FMT"\n", nto);
				fprintf(stderr, "Failed to sample subset of size %"SSIZE_FMT"\n", nto);				
				
				printf("probs:\n");
				for (i = 0; i < vector_dim(&probs); i++) {
					printf("%.10e, ", vector_item(&probs, i));
				}
				exit(1);
			}

			messages_add(&boot->messages, t, from, to, nto, msg->attr);
		}
		
		for (i = 0; i < n; i++) {
			msg = MESSAGES_VAL(it, i);
			frame_add(&boot->frame, msg);
		}
	}
	
	vector_deinit(&probs);
	xfree(to);
}


void recv_boot_deinit(struct recv_boot *boot)
{
	recv_model_deinit(&boot->model);
	frame_deinit(&boot->frame);
	messages_deinit(&boot->messages);	
}

