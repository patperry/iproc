#include "port.h"
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include "xalloc.h"
#include "../../src/history.h"
#include "../../src/design.h"
#include "../../src/design2.h"
#include "enron/actors.h"
#include "history.h"
#include "design2.h"


void design2_fixture_setup(struct design2_fixture *f, size_t count1, size_t count2)
{
	memset(f, 0, sizeof(*f));
	f->count1 = count1;
	f->count2 = count2;
}


void design2_fixture_setup_enron(struct design2_fixture *f)
{
	design2_fixture_setup(f, ENRON_ACTOR_COUNT, ENRON_ACTOR_COUNT);
}


void design2_fixture_teardown(struct design2_fixture *f)
{
	if (f->nsend.exists)
		free((void *)f->nsend.intvls);
	if (f->nrecv.exists)
		free((void *)f->nrecv.intvls);
	free((void *)f->trait_x);
}


void design2_fixture_add_irecv(struct design2_fixture *f, double window)
{
	f->irecv.exists = 1;
	f->irecv.name = "IRecv";
	f->irecv.window = window;
	f->tvar_dim += 1;
}


void design2_fixture_add_isend(struct design2_fixture *f, double window)
{
	f->isend.exists = 1;
	f->isend.name = "ISend";
	f->isend.window = window;
	f->tvar_dim += 1;
}


void design2_fixture_add_nrecv(struct design2_fixture *f, const double *intvls, size_t nintvl)
{
	f->nrecv.exists = 1;
	f->nrecv.name = "NRecv";
	f->nrecv.intvls = xmalloc(nintvl * sizeof(double));
	memcpy((void *)f->nrecv.intvls, intvls, nintvl * sizeof(double));
	f->nrecv.nintvl = nintvl;
	f->tvar_dim += nintvl;
}


void design2_fixture_add_nsend(struct design2_fixture *f, const double *intvls, size_t nintvl)
{
	f->nsend.exists = 1;
	f->nsend.name = "NSend";
	f->nsend.intvls = xmalloc(nintvl * sizeof(double));
	memcpy((void *)f->nsend.intvls, intvls, nintvl * sizeof(double));
	f->nsend.nintvl = nintvl;
	f->tvar_dim += nintvl;
}


void design2_test_setup(struct design2 *d, struct history *h, const struct design2_fixture *f)
{
	design2_init(d, h, f->count1, f->count2);
	if (f->irecv.exists)
		design2_add_tvar(d, f->irecv.name, VAR2_IRECV, f->irecv.window);
	if (f->isend.exists)
		design2_add_tvar(d, f->isend.name, VAR2_ISEND, f->isend.window);
	//if (f->nrecv.exists)
	//	design2_add_tvar(d, f->nrecv.name, VAR2_NRECV, f->nrecv.intvls, f->nrecv.nintvl);
	//if (f->nsend.exists)
	//	design2_add_tvar(d, f->nsend.name, VAR2_NSEND, f->nsend.intvls, f->nsend.nintvl);
	assert(design2_tvar_dim(d) == f->tvar_dim);
}

void design2_test_teardown(struct design2 *d)
{
	design2_deinit(d);
}
