#include "port.h"

#include <assert.h>
#include <stdlib.h>
#include "xalloc.h"
#include "vector.h"
#include "recv_boot.h"

static size_t sample1(const double *probs, size_t n, dsfmt_t * dsfmt)
{
	size_t i;
	double u = dsfmt_genrand_close_open(dsfmt);
	double sum = 0;

	for (i = 0; i < n - 1; i++) {
		sum += probs[i];
		if (sum >= u)
			break;
	}

	return i;
}

static bool unique(size_t *vals, size_t len)
{
	size_t i, j;

	if (len == 1)
		return true;
	if (len == 2)
		return vals[0] != vals[1];

	for (i = 0; i < len; i++) {
		size_t val = vals[i];

		for (j = i + 1; j < len; j++) {
			if (val == vals[j])
				return false;
		}
	}
	return true;
}

static bool sample_subset(const double *probs, size_t n, dsfmt_t * dsfmt,
			  size_t maxntry, size_t *out, size_t nout)
{
	size_t itry;
	size_t i;

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
		    size_t ncohort,
		    const size_t *cohorts,
		    const struct matrix *coefs, dsfmt_t * dsfmt)
{
	messages_init(&boot->messages);
	frame_init(&boot->frame, design);
	recv_model_init(&boot->model, &boot->frame, ncohort, cohorts, coefs);
	size_t nrecv = design_recv_count(design);

	size_t max_nto = messages_max_nto(msgs);
	size_t *to = xcalloc(max_nto, sizeof(to[0]));

	struct vector probs;

	struct messages_iter it;
	const struct message *msg;
	double t;
	size_t i, n;
	size_t from, nto;

	vector_init(&probs, nrecv);
	double *p = vector_to_ptr(&probs);

	size_t maxntry = 100000000;

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
				fprintf(stdout,
					"Failed to sample subset of size %"
					SSIZE_FMT "\n", nto);
				fprintf(stderr,
					"Failed to sample subset of size %"
					SSIZE_FMT "\n", nto);

				printf("probs:\n");
				for (i = 0; i < (size_t)vector_dim(&probs); i++) {
					printf("%.10e, ",
					       vector_item(&probs, i));
				}
				exit(1);
			}

			messages_add(&boot->messages, t, from, to, nto,
				     msg->attr);
		}

		for (i = 0; i < n; i++) {
			msg = MESSAGES_VAL(it, i);
			frame_add(&boot->frame, msg);
		}
	}

	vector_deinit(&probs);
	free(to);
}

void recv_boot_deinit(struct recv_boot *boot)
{
	recv_model_deinit(&boot->model);
	frame_deinit(&boot->frame);
	messages_deinit(&boot->messages);
}
