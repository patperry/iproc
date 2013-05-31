#ifndef ENRON_HISTORY_H
#define ENRON_HISTORY_H

#include <stddef.h>
#include <stdint.h>

struct history_fixture {
	size_t nsend;
	size_t nrecv;
	double *time;
	size_t *from;
	size_t **to;
	size_t *nto;
	intptr_t *attr;
	size_t count;
};


void history_fixture_setup_enron(struct history_fixture *f);
void history_fixture_teardown(struct history_fixture *f);

void history_test_setup(struct history *h, const struct history_fixture *f);
void history_test_teardown(struct history *h, const struct history_fixture *f);


#endif /* ENRON_HISTORY_H */
