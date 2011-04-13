#include "port.h"

#include <assert.h>
#include <math.h>
#include "logsumexp.h"

/* Algorithm logsumexp
 *   INPUT: a[1], ..., a[N]
 *   OUTPUT: log(sum(exp(ai)))
 *
 *   INVARIANT: S[n] = sum(exp(a[i] - M[n])) - 1
 *              M[n] = max(a[i])
 *
 *   S := -1
 *   M := -INFINITY
 *
 *   for i = 1 to N do
 *       if a[i] <= M
 *       then
 *           S := S + exp(a[i] - M)
 *       else
 *           S := (S + 1) * exp(M - a[i])
 *           M := a[i]
 *       end
 *   end
 *
 *   return M + log(S + 1)
 */
void iproc_logsumexp_init(iproc_logsumexp * lse)
{
	assert(lse);
	lse->max = -INFINITY;
	lse->sumexpm1 = -1.0;
}

void iproc_logsumexp_insert(iproc_logsumexp * lse, double val)
{
	assert(lse);

	if (isinf(val)) {
		if (val > 0) {
			lse->max = val;
			lse->sumexpm1 = val;
		}
		return;
	}

	double sumexpm1 = lse->sumexpm1;
	double max = lse->max;

	if (!(val > max)) {
		lse->sumexpm1 += exp(val - max);
	} else {
		double scale = exp(max - val);
		lse->sumexpm1 = sumexpm1 * scale + scale;
		lse->max = val;
	}
}

double iproc_logsumexp_value(iproc_logsumexp * lse)
{
	assert(lse);

	double sumexpm1 = lse->sumexpm1;
	double max = lse->max;

	return max + log1p(sumexpm1);
}

double iproc_logsumexp_max(iproc_logsumexp * lse)
{
	assert(lse);
	return lse->max;
}
