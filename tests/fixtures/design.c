#include "port.h"
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include "xalloc.h"
#include "../../src/history.h"
#include "../../src/design.h"
#include "enron/actors.h"
#include "history.h"
#include "design.h"


void design_fixture_setup(struct design_fixture *f, size_t count)
{
	memset(f, 0, sizeof(*f));
	f->count = count;
}


void design_fixture_setup_enron(struct design_fixture *f)
{
	design_fixture_setup(f, ENRON_ACTOR_COUNT);
}

void design_fixture_add_traits_enron(struct design_fixture *f)
{
	enron_actors_init(ENRON_TERMS_MAX, &f->count, (double **)&f->trait_x,
			  &f->trait_names, &f->trait_dim);
}


void design_fixture_teardown(struct design_fixture *f)
{
	if (f->nsendtot.exists)
		free((void *)f->nsendtot.intvls);
	if (f->nrecvtot.exists)
		free((void *)f->nrecvtot.intvls);
	free((void *)f->trait_x);
}


void design_fixture_add_irecvtot(struct design_fixture *f, double window)
{
	f->irecvtot.exists = 1;
	f->irecvtot.name = "IRecvTot";
	f->irecvtot.window = window;
	f->tvar_dim += 1;
}


void design_fixture_add_isendtot(struct design_fixture *f, double window)
{
	f->isendtot.exists = 1;
	f->isendtot.name = "ISendTot";
	f->isendtot.window = window;
	f->tvar_dim += 1;
}


void design_fixture_add_nrecvtot(struct design_fixture *f, const double *intvls, size_t nintvl)
{
	f->nrecvtot.exists = 1;
	f->nrecvtot.name = "NRecvTot";
	f->nrecvtot.intvls = xmalloc(nintvl * sizeof(double));
	memcpy((void *)f->nrecvtot.intvls, intvls, nintvl * sizeof(double));
	f->nrecvtot.nintvl = nintvl;
	f->tvar_dim += nintvl;
}


void design_fixture_add_nsendtot(struct design_fixture *f, const double *intvls, size_t nintvl)
{
	f->nsendtot.exists = 1;
	f->nsendtot.name = "NSendTot";
	f->nsendtot.intvls = xmalloc(nintvl * sizeof(double));
	memcpy((void *)f->nsendtot.intvls, intvls, nintvl * sizeof(double));
	f->nsendtot.nintvl = nintvl;
	f->tvar_dim += nintvl;
}


void design_test_setup(struct design *d, struct history *h, const struct design_fixture *f)
{
	design_init(d, h, f->count);
	design_add_traits(d, f->trait_names, f->trait_x, f->trait_dim);
	if (f->irecvtot.exists)
		design_add_tvar(d, f->irecvtot.name, VAR_IRECVTOT, f->irecvtot.window);
	if (f->isendtot.exists)
		design_add_tvar(d, f->isendtot.name, VAR_ISENDTOT, f->isendtot.window);
	if (f->nrecvtot.exists)
		design_add_tvar(d, f->nrecvtot.name, VAR_NRECVTOT, f->nrecvtot.intvls, f->nrecvtot.nintvl);
	if (f->nsendtot.exists)
		design_add_tvar(d, f->nsendtot.name, VAR_NSENDTOT, f->nsendtot.intvls, f->nsendtot.nintvl);
	assert(design_tvar_dim(d) == f->tvar_dim);
}

void design_test_teardown(struct design *d)
{
	design_deinit(d);
}
