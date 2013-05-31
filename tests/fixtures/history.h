#ifndef FIXTURES_HISTORY_H
#define FIXTURES_HISTORY_H

#include <stddef.h>
#include <stdint.h>

struct history_fixture {
	size_t nsend;
	size_t nrecv;
	const double *time;
	const size_t *from;
	const size_t * const *to;
	const size_t *nto;
	const intptr_t *attr;
	size_t count;
};


void history_fixture_setup_enron(struct history_fixture *f);
void history_fixture_teardown(struct history_fixture *f);

void history_test_setup(struct history *h, const struct history_fixture *f);
void history_test_teardown(struct history *h);


#endif /* FIXTURES_HISTORY_H */
