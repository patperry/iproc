#include "port.h"

#include <assert.h>
#include <stdlib.h>
#include "actors.h"

static void cohort_init(struct cohort *c)
{
	assert(c);
	array_init(&c->actors, sizeof(ssize_t));
}

static void cohort_init_copy(struct cohort *c, const struct cohort *src)
{
	assert(c);
	array_init_copy(&c->actors, &src->actors);
}

static void cohort_deinit(struct cohort *c)
{
	assert(c);
	array_deinit(&c->actors);
}

void actors_init(struct actors *a)
{
	assert(a);

	array_init(&a->actors, sizeof(struct actor));
	array_init(&a->cohorts, sizeof(struct cohort));
	refcount_init(&a->refcount);
}

void actors_init_copy(struct actors *a, const struct actors *src)
{
	assert(a);
	assert(src);
	array_init_copy(&a->actors, &src->actors);

	array_init(&a->cohorts, sizeof(struct cohort));
	struct cohort *csrc, *cdst;
	ARRAY_FOREACH(csrc, &src->cohorts) {
		cdst = array_add(&a->cohorts, NULL);
		cohort_init_copy(cdst, csrc);
	}

	refcount_init(&a->refcount);
}

void actors_clear(struct actors *a)
{
	array_clear(&a->actors);

	struct cohort *c;
	ARRAY_FOREACH(c, &a->cohorts) {
		cohort_deinit(c);
	}
	array_clear(&a->cohorts);
}

void actors_deinit(struct actors *a)
{
	assert(a);

	actors_clear(a);
	refcount_deinit(&a->refcount);
	array_deinit(&a->cohorts);
	array_deinit(&a->actors);
}

struct actors *actors_alloc(void)
{
	struct actors *actors = xcalloc(1, sizeof(*actors));
	actors_init(actors);
	return actors;
}

void actors_free(struct actors *a)
{
	if (!a)
		return;

	if (refcount_put(&a->refcount, NULL)) {
		refcount_get(&a->refcount);
		actors_deinit(a);
		xfree(a);
	}
}

struct actors *actors_ref(struct actors *actors)
{
	assert(actors);
	refcount_get(&actors->refcount);
	return actors;
}

ssize_t actors_add(struct actors *a, ssize_t cohort)
{
	assert(a);
	assert(0 <= cohort && cohort < actors_cohort_count(a));

	struct cohort *c = (struct cohort *)actors_cohorts(a);
	ssize_t id = array_count(&a->actors);
	array_add(&a->actors, &cohort);
	array_add(&c[cohort].actors, &id);
	return id;
}

ssize_t actors_add_cohort(struct actors *a)
{
	assert(a);

	ssize_t cid = array_count(&a->cohorts);
	struct cohort *c = array_add(&a->cohorts, NULL);
	cohort_init(c);
	return cid;
}

void actors_mul(double alpha, enum trans_op trans, const struct actors *a,
		const struct vector *x, double beta, struct vector *y)
{
	assert(a);
	assert(x);
	assert(y);
	assert(trans != TRANS_NOTRANS
	       || vector_dim(x) == actors_cohort_count(a));
	assert(trans != TRANS_NOTRANS || vector_dim(y) == actors_count(a));
	assert(trans == TRANS_NOTRANS || vector_dim(x) == actors_count(a));
	assert(trans == TRANS_NOTRANS
	       || vector_dim(y) == actors_cohort_count(a));

	if (beta == 0) {
		vector_fill(y, 0.0);
	} else if (beta != 1) {
		vector_scale(y, beta);
	}

	if (trans == TRANS_NOTRANS) {
		const struct actor *item = actors_items(a);
		ssize_t i, n = actors_count(a);

		for (i = 0; i < n; i++) {
			ssize_t c = item[i].cohort;
			double xc = vector_item(x, c);
			double *yi = vector_item_ptr(y, i);
			*yi += alpha * xc;
		}
	} else {
		const struct cohort *cohort = actors_cohorts(a);
		ssize_t c, p = actors_cohort_count(a);

		for (c = 0; c < p; c++) {
			double *yc = vector_item_ptr(y, c);
			const ssize_t *actors = array_to_ptr(&cohort[c].actors);
			ssize_t ic, nc = array_count(&cohort[c].actors);

			for (ic = 0; ic < nc; ic++) {
				ssize_t i = actors[ic];
				double xi = vector_item(x, i);
				*yc += alpha * xi;
			}
		}
	}
}

void
actors_muls(double alpha,
	    enum trans_op trans,
	    const struct actors *a,
	    const struct svector *x, double beta, struct svector *y)
{
	assert(a);
	assert(x);
	assert(y);
	assert(trans != TRANS_NOTRANS
	       || svector_dim(x) == actors_cohort_count(a));
	assert(trans != TRANS_NOTRANS || svector_dim(y) == actors_count(a));
	assert(trans == TRANS_NOTRANS || svector_dim(x) == actors_count(a));
	assert(trans == TRANS_NOTRANS
	       || svector_dim(y) == actors_cohort_count(a));

	if (beta == 0) {
		svector_clear(y);
	} else if (beta != 1) {
		svector_scale(y, beta);
	}

	if (trans == TRANS_NOTRANS) {
		const struct cohort *cohort = actors_cohorts(a);
		struct svector_iter itx;
		SVECTOR_FOREACH(itx, x) {
			ssize_t c = SVECTOR_IDX(itx);
			double xc = SVECTOR_VAL(itx);
			const ssize_t *actors = array_to_ptr(&cohort[c].actors);
			ssize_t ic, nc = array_count(&cohort[c].actors);

			for (ic = 0; ic < nc; ic++) {
				ssize_t i = actors[ic];
				double *yi = svector_item_ptr(y, i);
				*yi += alpha * xc;
			}
		}
	} else {
		const struct actor *item = actors_items(a);
		struct svector_iter itx;
		SVECTOR_FOREACH(itx, x) {
			ssize_t i = SVECTOR_IDX(itx);
			double xi = SVECTOR_VAL(itx);
			ssize_t c = item[i].cohort;
			double *yc = svector_item_ptr(y, c);
			*yc += alpha * xi;
		}
	}
}

void
actors_matmul(double alpha,
	      enum trans_op trans,
	      const struct actors *a,
	      const struct matrix *x, double beta, struct matrix *y)
{
	assert(a);
	assert(x);
	assert(y);
	assert(trans != TRANS_NOTRANS
	       || matrix_nrow(x) == actors_cohort_count(a));
	assert(trans != TRANS_NOTRANS || matrix_nrow(y) == actors_count(a));
	assert(trans == TRANS_NOTRANS || matrix_nrow(x) == actors_count(a));
	assert(trans == TRANS_NOTRANS
	       || matrix_nrow(y) == actors_cohort_count(a));
	assert(matrix_ncol(x) == matrix_ncol(y));

	ssize_t m = matrix_ncol(x);
	ssize_t j;

	if (beta == 0) {
		matrix_fill(y, 0.0);
	} else if (beta != 1) {
		matrix_scale(y, beta);
	}

	struct vector xcol;
	struct vector ycol;

	for (j = 0; j < m; j++) {
		xcol = matrix_col(x, j);
		ycol = matrix_col(y, j);
		actors_mul(alpha, trans, a, &xcol, 1.0, &ycol);
	}
}
