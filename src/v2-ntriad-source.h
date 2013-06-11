#include "port.h"
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include "coreutil.h"
#include "ieee754.h"
#include "pqueue.h"
#include "xalloc.h"
#include "design2.h"


struct ntriad_thunk {
	struct design2 dyad1, dyad2;
	const struct var2 *v1, *v2;
};



static void ntriad_init(struct var_meta *meta, void **thunk, const char *name,
		       struct design2 *d, va_list ap)
{
	const double *intvls1 = va_arg(ap, double *);
	size_t nintvl1 = va_arg(ap, size_t);
	const double *intvls2 = va_arg(ap, double *);
	size_t nintvl2 = va_arg(ap, size_t);
	struct history *h = design2_history(d);
	struct ntriad_thunk *ntriad = xmalloc(sizeof(*ntriad));
	size_t dims[2] = { nintvl1, nintvl2 };
	var_meta_init(meta, VAR_TYPE_TVAR, name, dims, 2);
	size_t n = design2_count1(d);
	assert(design2_count2(d) == n);

	design2_init(&ntriad->dyad1, h, n, n);
	design2_init(&ntriad->dyad2, h, n, n);
	ntriad->v1 = design2_add_tvar(&ntriad->dyad1, "V1", TRIAD_V1, intvls1, nintvl1);
	ntriad->v2 = design2_add_tvar(&ntriad->dyad2, "V2", TRIAD_V2, intvls2, nintvl2);

	*thunk = ntriad;
}


static void ntriad_deinit(struct var_meta *meta, void *thunk, struct design2 *d)
{
	struct ntriad_thunk *ntriad = thunk;
	design2_deinit(&ntriad->dyad2);
	design2_deinit(&ntriad->dyad1);
	free(ntriad);
	var_meta_deinit(meta);
}


static int needs_update(struct tvar2 *tv, size_t i, double t0)
{
	struct ntriad_thunk *ntriad = tv->thunk;
	int upd = 0;

	if (!design2_count1(&ntriad->dyad1))
		return 0;

	const struct deltaset *ds1 = design2_changes(&ntriad->dyad1, i);
	if (deltaset_thead(ds1) >= t0) {
		upd = 1;
	} else {
		const double *x;
		const size_t *ind;
		size_t iz, nz;

		design2_get_tvar_matrix(&ntriad->dyad1, i, &x, &ind, &nz);
		for (iz = 0; iz < nz; iz++) {
			size_t h = ind[iz];
			const struct deltaset *ds2 = design2_changes(&ntriad->dyad2, h);
			if (deltaset_thead(ds2) >= t0) {
				upd = 1;
				break;
			}
		}
	}

	return upd;
}


static void clear(struct tvar2 *tv, size_t i)
{
	struct ntriad_thunk *ntriad = tv->thunk;
	struct design2 *d = tv->var.design;
	size_t size = tv->var.meta.size;
	const double *x1, *x2;
	const size_t *ind1, *ind2;
	size_t iz1, nz1, iz2, nz2;

	design2_get_tvar_matrix(&ntriad->dyad1, i, &x1, &ind1, &nz1);
	for (iz1 = 0; iz1 < nz1; iz1++) {
		size_t h = ind1[iz1];

		design2_get_tvar_matrix(&ntriad->dyad2, h, &x2, &ind2, &nz2);
		for (iz2 = 0; iz2 < nz2; iz2++) {
			size_t j = ind2[iz2];
			double *x = design2_make_active(d, tv, i, j);
			memset(x, 0, size * sizeof(double));
		}
	}
}

static double do_update(struct tvar2 *tv, size_t i)
{
	struct ntriad_thunk *ntriad = tv->thunk;
	struct design2 *d = tv->var.design;
	const size_t *dims = tv->var.meta.dims;
	size_t n1 = dims[0], n2 = dims[1];
	const double *x1, *x2;
	const size_t *ind1, *ind2;
	size_t iz1, nz1, iz2, nz2;
	double tnext = INFINITY;

	if (!n2)
		return INFINITY;

	clear(tv, i);

	const struct deltaset *ds1 = design2_changes(&ntriad->dyad1, i);
	double tnext1 = design2_next_time(&ntriad->dyad1, i);

	tnext = tnext1;

	design2_get_tvar_matrix(&ntriad->dyad1, i, &x1, &ind1, &nz1);
	for (iz1 = 0; iz1 < nz1; iz1++, x1 += n1) {
		size_t h = ind1[iz1];
		double tnext2 = design2_next_time(&ntriad->dyad2, h);
		double tlast1 = deltaset_tlast(ds1, h);
		const struct deltaset *ds2 = design2_changes(&ntriad->dyad2, h);

		tnext = MIN(tnext, tnext2);

		design2_get_tvar_matrix(&ntriad->dyad2, h, &x2, &ind2, &nz2);
		for (iz2 = 0; iz2 < nz2; iz2++, x2 += n2) {
			size_t j = ind2[iz2];
			double tlast2 = deltaset_tlast(ds2, j);
			double *x = design2_make_active(d, tv, i, j);

			blas_dger(n2, n1, 1.0, x2, 1, x1, 1, x, n2);

			double tlast = MAX(tlast1, tlast2);
			design2_update(d, tv, i, j, tlast);
		}
	}

	return tnext;
}


static double ntriad_update(struct tvar2 *tv, size_t i, double t0,
			    const struct history *h)
{
	double tnext = INFINITY;

	if (needs_update(tv, i, t0)) {
		tnext = do_update(tv, i);
	}

	return tnext;
}


static struct tvar2_type VAR2_NTRIAD_REP = {
	ntriad_init,
	ntriad_deinit,
	ntriad_update
};

