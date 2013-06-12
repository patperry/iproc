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
	size_t j, n = recv_model_count(m);

	for (j = 0; j < n; j++) {
		double p = catdist1_prob(dist, j);
		y[j] += alpha * p;
	}
}

void recv_boot_init(struct recv_boot *boot,
		    const struct recv_model *m,
		    const struct message *msgs,
		    size_t nmsg,
		    dsfmt_t * dsfmt)

{
	size_t nsend = recv_model_send_count(m);
	size_t j, nrecv = recv_model_count(m);
	const struct design *r = recv_model_design(m);
	const struct design2 *d = recv_model_design2(m);
	struct history *hr = design_history(r);
	struct history *hd = design2_history(d);
	size_t imsg;

	size_t max_nto = nrecv;
	size_t *to = xcalloc(max_nto, sizeof(to[0]));
	double *p = xmalloc(nrecv * sizeof(p[0]));
	size_t maxntry = 100000000;

	history_init(&boot->history, nsend, nrecv);

	for (imsg = 0; imsg < nmsg; imsg++) {
		const struct message *msg = &msgs[imsg];
		double t = msg->time;
		size_t from = msg->from;
		size_t nto = msg->nto;

		if (t < history_time(hr)) {
			history_reset(hr);
		}
		history_advance(hr, t);

		if (t < history_time(hd)) {
			history_reset(hd);
		}
		history_advance(hd, t);


		memset(p, 0, nrecv * sizeof(p[0]));
		axpy_probs(1.0, m, from, p);

		if (!sample_subset(p, nrecv, dsfmt, maxntry, to, nto)) {
			fprintf(stdout,
				"Failed to sample subset of size %zu\n",
				nto);
			fprintf(stderr,
				"Failed to sample subset of size %zu\n",
				nto);

			printf("probs:\n");
			for (j = 0; j < nrecv; j++) {
				printf("%.10e, ", p[j]);
			}
			exit(1);
		}

		if (t < history_time(&boot->history)) {
			history_reset(&boot->history);
		}
		history_advance(&boot->history, t);
		history_add(&boot->history, from, to, nto, msg->attr);
	}

	free(p);
	free(to);
}

void recv_boot_deinit(struct recv_boot *boot)
{
	history_deinit(&boot->history);
}
