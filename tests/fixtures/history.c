#include "port.h"
#include <assert.h>
#include <stdlib.h>
#include "../../src/history.h"
#include "enron/actors.h"
#include "enron/messages.h"
#include "history.h"

void history_fixture_setup_enron(struct history_fixture *f)
{
	int ok;

	ok = enron_messages_init(0, (double **)&f->time, (size_t **)&f->from,
				 (size_t ***)&f->to, (size_t **)&f->nto,
				 (intptr_t **)&f->attr, &f->count);
	f->nsend = ENRON_ACTOR_COUNT;
	f->nrecv = ENRON_ACTOR_COUNT;
	assert(ok);
}


void history_fixture_teardown(struct history_fixture *f)
{
	size_t i, n = f->count;

	for (i = 0; i < n; i++) {
		free((void *)f->to[i]);
	}

	free((void *)f->attr);
	free((void *)f->nto);
	free((void *)f->to);
	free((void *)f->from);
	free((void *)f->time);
}


void history_test_setup(struct history *h, const struct history_fixture *f)
{
	size_t i, n = f->count;

	history_init(h, f->nsend, f->nrecv);

	for (i = 0; i < n; i++) {
		history_advance(h, f->time[i]);
		history_add(h, f->from[i], f->to[i], f->nto[i], f->attr[i]);
	}

	history_reset(h);
}


void history_test_teardown(struct history *h)
{
	history_deinit(h);
}
