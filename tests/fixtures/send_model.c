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
	f->params.coefs.traits = xcalloc(d->trait_dim, sizeof(double));
	f->params.coefs.tvars = xcalloc(d->tvar_dim, sizeof(double));
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
