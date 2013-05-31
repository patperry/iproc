#include "port.h"
#include <assert.h>
#include <stdlib.h>
#include <string.h>
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
	free((void *)f->trait_x);
}


void design_fixture_add_irecvtot(struct design_fixture *f, double window)
{
	f->irecvtot.exists = 1;
	f->irecvtot.name = "IRecvTot";
	f->irecvtot.window = window;
}

void design_test_setup(struct design *d, struct history *h, const struct design_fixture *f)
{
	design_init(d, h, f->count);
	design_add_traits(d, f->trait_names, f->trait_x, f->trait_dim);
	if (f->irecvtot.exists)
		design_add_tvar(d, f->irecvtot.name, VAR_IRECVTOT, f->irecvtot.window);
}

void design_test_teardown(struct design *d)
{
	design_deinit(d);
}
