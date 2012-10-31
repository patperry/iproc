#include "port.h"

#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include "xalloc.h"
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

static int unique(size_t *vals, size_t len)
{
	size_t i, j;

	if (len == 1)
		return 1;
	if (len == 2)
		return vals[0] != vals[1];

	for (i = 0; i < len; i++) {
		size_t val = vals[i];

		for (j = i + 1; j < len; j++) {
			if (val == vals[j])
				return 0;
		}
	}
	return 1;
}

static int sample_subset(const double *probs, size_t n, dsfmt_t * dsfmt,
			  size_t maxntry, size_t *out, size_t nout)
{
	size_t itry;
	size_t i;

	for (itry = 0; itry < maxntry; itry++) {
		for (i = 0; i < nout; i++) {
			out[i] = sample1(probs, n, dsfmt);
		}
		if (unique(out, nout))
			return 1;
	}
	return 0;
}


static void axpy_probs(double alpha, const struct recv_model *m, size_t isend,
		       double *y)
{
	const struct catdist1 *dist = recv_model_dist(m, isend);
	const struct frame *f = recv_model_frame(m);
	size_t j, n = frame_recv_count(f);

	for (j = 0; j < n; j++) {
		double p = catdist1_cached_prob(dist, j);
		y[j] += alpha * p;
	}
}

void recv_boot_init(struct recv_boot *boot,
		    struct frame *f,
		    const struct messages *msgs,
		    const struct recv_coefs *coefs, dsfmt_t * dsfmt)
{
	messages_init(&boot->messages);
	recv_model_init(&boot->model, f, coefs);

	struct history *h = frame_history(&boot->frame);
	size_t nrecv = frame_recv_count(f);
	size_t max_nto = messages_max_nto(msgs);
	size_t *to = xcalloc(max_nto, sizeof(to[0]));
	double *p = xmalloc(nrecv * sizeof(p[0]));
	size_t maxntry = 100000000;

	struct messages_iter it;
	const struct message *msg;
	double t;
	size_t i, n;
	size_t from, nto;

	MESSAGES_FOREACH(it, msgs) {
		t = MESSAGES_TIME(it);
		n = MESSAGES_COUNT(it);

		history_advance(h, t);

		for (i = 0; i < n; i++) {
			msg = MESSAGES_VAL(it, i);
			from = msg->from;
			nto = msg->nto;

			memset(p, 0, nrecv * sizeof(p[0]));
			axpy_probs(1.0, &boot->model, from, p);

			if (!sample_subset(p, nrecv, dsfmt, maxntry, to, nto)) {
				fprintf(stdout,
					"Failed to sample subset of size %zu\n",
					nto);
				fprintf(stderr,
					"Failed to sample subset of size %zu\n",
					nto);

				printf("probs:\n");
				for (i = 0; i < nrecv; i++) {
					printf("%.10e, ", p[i]);
				}
				exit(1);
			}

			messages_add(&boot->messages, t, from, to, nto,
				     msg->attr);
		}

		for (i = 0; i < n; i++) {
			msg = MESSAGES_VAL(it, i);
			history_add(h, msg);
		}
	}

	free(p);
	free(to);
}

void recv_boot_deinit(struct recv_boot *boot)
{
	recv_model_deinit(&boot->model);
	frame_deinit(&boot->frame);
	messages_deinit(&boot->messages);
}
