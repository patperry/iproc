#include <stdlib.h>
#include "../../src/history.h"
#include "../../src/design.h"
#include "../../src/send_model.h"
#include "xalloc.h"
#include "history.h"
#include "design.h"
#include "send_model.h"


void send_model_fixture_setup(struct send_model_fixture	*f,
			      const struct design_fixture *d)
{
	size_t trait_dim = d->trait_dim;
	size_t tvar_dim = d->tvar_dim;
	
	f->params.coefs.traits = xcalloc(trait_dim, sizeof(double));
	f->params.coefs.tvars = xcalloc(tvar_dim, sizeof(double));
	f->trait_dim = trait_dim;
	f->tvar_dim = tvar_dim;
}

void send_model_fixture_set_rand(struct send_model_fixture *f)
{
	size_t i;
	int u;
	double x;

	for (i = 0; i < f->trait_dim; i++) {
		u = rand();
		x = ((double) (u % 21) - 10) / 10;
		f->params.coefs.traits[i] = x;
	}

	for (i = 0; i < f->tvar_dim; i++) {
		u = rand();
		x = ((double) (u % 21) - 10) / 10;
		f->params.coefs.tvars[i] = x;
	}
}

void send_model_fixture_teardown(struct send_model_fixture *f)
{
	free(f->params.coefs.tvars);
	free(f->params.coefs.traits);
}


void send_model_test_setup(struct send_model *m, struct design *d,
			   const struct send_model_fixture *f)
{
	send_model_init(m, &f->params, d);
}


void send_model_test_teardown(struct send_model *m)
{
	send_model_deinit(m);
}
