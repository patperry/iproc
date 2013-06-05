#ifndef FIXTURES_DESIGN_H
#define FIXTURES_DESIGN_H

#include <stddef.h>
#include <stdint.h>


struct design_fixture {
	size_t count;

	const double *trait_x;
	const char * const *trait_names;
	size_t trait_dim;

	struct { int exists; const char *name; double window; } irecvtot;
	struct { int exists; const char *name; double window; } isendtot;
	struct { int exists; const char *name; const double *intvls; size_t nintvl; } nrecvtot;
	struct { int exists; const char *name; const double *intvls; size_t nintvl; } nsendtot;
	size_t tvar_dim;
};


void design_fixture_setup(struct design_fixture *f, size_t count);
void design_fixture_setup_enron(struct design_fixture *f);
void design_fixture_add_traits_enron(struct design_fixture *f);
void design_fixture_add_irecvtot(struct design_fixture *f, double window);
void design_fixture_add_isendtot(struct design_fixture *f, double window);
void design_fixture_add_nrecvtot(struct design_fixture *f, const double *intvls, size_t nintvl);
void design_fixture_add_nsendtot(struct design_fixture *f, const double *intvls, size_t nintvl);
void design_fixture_teardown(struct design_fixture *f);

void design_test_setup(struct design *d, struct history *h, const struct design_fixture *f);
void design_test_teardown(struct design *d);

#endif /* FIXTURES_DESIGN_H */