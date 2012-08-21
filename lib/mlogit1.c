#include "port.h"
#include <math.h>
#include <stdlib.h>
#include "mlogit1.h"


void mlogit1_init(struct mlogit1 *m1, const struct mlogit *parent)
{
	assert(parent);

	m1->parent = parent;
	m1->eta1 = NULL;
	m1->ind = NULL;
	m1->nz = 0;
	m1->nzmax = 0;

	mlogit1_clear(m1);
}


void mlogit1_deinit(struct mlogit1 *m1)
{
	free(m1->ind);
	free(m1->eta1);
}

void mlogit1_clear(struct mlogit1 *m1)
{
	m1->nz = 0;
}

double mlogit1_eta(const struct mlogit1 *m1, size_t i)
{
	assert(i < mlogit1_ncat(m1));
	return mlogit_eta(m1->parent, i);
}

double mlogit1_prob(const struct mlogit1 *m1, size_t i)
{
	assert(i < mlogit1_ncat(m1));
	double lp = mlogit1_lprob(m1, i);
	return exp(lp);
}

double mlogit1_lprob(const struct mlogit1 *m1, size_t i)
{
	assert(i < mlogit1_ncat(m1));
	return mlogit_lprob(m1->parent, i);
}


void mlogit1_set_eta(struct mlogit1 *m1, size_t i, double eta)
{
	assert(i < mlogit1_ncat(m1));
}

int _mlogit1_check(const struct mlogit1 *m)
{
	int fail = 0;
	return fail;
}
