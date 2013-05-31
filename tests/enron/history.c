#include "port.h"
#include <assert.h>
#include <stdlib.h>
#include "../../src/history.h"
#include "enron/messages.h"
#include "enron/history.h"

void history_fixture_setup_enron(struct history_fixture *f)
{
	int ok;

	ok = enron_messages_init(&f->nsend, &f->nrecv, &f->time, &f->from,
				 &f->to, &f->nto, &f->attr, &f->count, 0);
	assert(ok);
}


void history_fixture_teardown(struct history_fixture *f)
{
	size_t i, n = f->count;

	for (i = 0; i < n; i++) {
		free(f->to[i]);
	}

	free(f->attr);
	free(f->nto);
	free(f->to);
	free(f->from);
	free(f->time);
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


void history_test_teardown(struct history *h, const struct history_fixture *f)
{
	history_deinit(h);
}
