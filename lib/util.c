
#include "logsumexp.h"
#include "util.h"
#include "xalloc.h"
#include <math.h>
#include <stdlib.h>



void free2(void **ptrs, size_t len)
{
	size_t i;
	
	if (!ptrs)
		return;
	
	for (i = len; i > 0; i--) {
		free(ptrs[i-1]);
	}
	free(ptrs);
}


char **xstrdup2(const char *const *strs, size_t len)
{
	if (!strs)
		return NULL;

	char **res = xmalloc(len * sizeof(res[0]));
	size_t i;

	for (i = 0; i < len; i++) {
		res[i] = xstrdup(strs[i]);
	}
	return res;
}


void vector_shift(size_t n, double alpha, double *x)
{
	size_t i;
	for (i = 0; i < n; i++) {
		x[i] += alpha;
	}
}


void vector_exp(size_t n, double *x)
{
	size_t i;
	for (i = 0; i < n; i++) {
		x[i] = exp(x[i]);
	}
}


double vector_logsumexp(size_t n, double *x)
{
	struct logsumexp lse;
	size_t i;

	logsumexp_init(&lse);
	for (i = 0; i < n; i++) {
		logsumexp_insert(&lse, x[i]);
	}

	return logsumexp_value(&lse);
}


double vector_max(size_t n, double *x)
{
	size_t i;
	double max = -INFINITY;

	for (i = 0; i < n; i++) {
		if (x[i] > max)
			max = x[i];
	}

	return max;
}

