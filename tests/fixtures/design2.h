#ifndef FIXTURES_DESIGN2_H
#define FIXTURES_DESIGN2_H

#include <stddef.h>
#include <stdint.h>


struct design2_fixture {
	size_t count1, count2;

	const double *trait_x;
	const char * const *trait_names;
	size_t trait_dim;

	struct { int exists; const char *name; double window; } irecv;
	struct { int exists; const char *name; double window; } isend;
	struct { int exists; const char *name; const double *intvls; size_t nintvl; } nrecv;
	struct { int exists; const char *name; const double *intvls; size_t nintvl; } nsend;
	size_t tvar_dim;
};


void design2_fixture_setup(struct design2_fixture *f, size_t count1, size_t count2);
void design2_fixture_setup_enron(struct design2_fixture *f);
void design2_fixture_add_irecv(struct design2_fixture *f, double window);
void design2_fixture_add_isend(struct design2_fixture *f, double window);
void design2_fixture_add_nrecv(struct design2_fixture *f, const double *intvls, size_t nintvl);
void design2_fixture_add_nsend(struct design2_fixture *f, const double *intvls, size_t nintvl);
void design2_fixture_teardown(struct design2_fixture *f);

void design2_test_setup(struct design2 *d, struct history *h, const struct design2_fixture *f);
void design2_test_teardown(struct design2 *d);

#endif /* FIXTURES_DESIGN_H */