#include <stdlib.h>
#include "../../src/history.h"
#include "../../src/design.h"
#include "../../src/design2.h"
#include "../../src/recv_model.h"
#include "xalloc.h"
#include "history.h"
#include "design.h"
#include "design2.h"
#include "recv_model.h"


void recv_model_fixture_setup(struct recv_model_fixture	*f,
			      const struct design_fixture *r,
			      const struct design2_fixture *d)
{
	size_t recv_trait_dim = r->trait_dim;
	size_t recv_tvar_dim = r->tvar_dim;
	size_t dyad_trait_dim = d->trait_dim;
	size_t dyad_tvar_dim = d->tvar_dim;

	f->params.recv.traits = xcalloc(recv_trait_dim, sizeof(double));
	f->params.recv.tvars = xcalloc(recv_tvar_dim, sizeof(double));
	f->params.dyad.traits = xcalloc(dyad_trait_dim, sizeof(double));
	f->params.dyad.tvars = xcalloc(dyad_tvar_dim, sizeof(double));
	f->recv_trait_dim = recv_trait_dim;
	f->recv_tvar_dim = recv_tvar_dim;
	f->dyad_trait_dim = dyad_trait_dim;
	f->dyad_tvar_dim = dyad_tvar_dim;
}

void recv_model_fixture_set_rand(struct recv_model_fixture *f)
{
	size_t i;
	int u;
	double x;

	for (i = 0; i < f->recv_trait_dim; i++) {
		u = rand();
		x = ((double) (u % 21) - 10) / 10;
		f->params.recv.traits[i] = x;
	}

	for (i = 0; i < f->recv_tvar_dim; i++) {
		u = rand();
		x = ((double) (u % 21) - 10) / 10;
		f->params.recv.tvars[i] = x;
	}

	for (i = 0; i < f->dyad_trait_dim; i++) {
		u = rand();
		x = ((double) (u % 21) - 10) / 10;
		f->params.dyad.traits[i] = x;
	}

	for (i = 0; i < f->dyad_tvar_dim; i++) {
		u = rand();
		x = ((double) (u % 21) - 10) / 10;
		f->params.dyad.tvars[i] = x;
	}
}

void recv_model_fixture_teardown(struct recv_model_fixture *f)
{
	free(f->params.dyad.tvars);
	free(f->params.dyad.traits);
	free(f->params.recv.tvars);
	free(f->params.recv.traits);
}


void recv_model_test_setup(struct recv_model *m, struct design *r,
			   struct design2 *d,
			   const struct recv_model_fixture *f)
{
	recv_model_init(m, &f->params, r, d);
}


void recv_model_test_teardown(struct recv_model *m)
{
	recv_model_deinit(m);
}
