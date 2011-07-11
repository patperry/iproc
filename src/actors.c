#include "port.h"

#include <assert.h>
#include <stdlib.h>
#include "actors.h"

void actors_init(struct actors *actors, ssize_t dim)
{
	assert(actors);
	assert(dim >= 0);

	array_init(&actors->actors, sizeof(struct actor));
	array_init(&actors->cohorts, sizeof(struct cohort));
	intmap_init(&actors->trait_hashes, sizeof(struct intset),
		    alignof(struct intset));

	refcount_init(&actors->refcount);
	actors->dim = dim;
}

void actors_init_matrix(struct actors *actors, enum trans_op trans,
			const struct matrix *matrix)
{
	assert(actors);
	assert(matrix);

	ssize_t m = matrix_nrow(matrix);
	ssize_t n = matrix_ncol(matrix);
	ssize_t dim = trans == TRANS_NOTRANS ? n : m;

	actors_init(actors, dim);

	if (trans == TRANS_NOTRANS) {
		struct vector row;
		double *row_front;
		ssize_t i;

		vector_init(&row, n);
		row_front = vector_to_ptr(&row);

		for (i = 0; i < m; i++) {
			matrix_get_row(matrix, i, row_front);
			actors_add(actors, &row);
		}
		vector_deinit(&row);
	} else {
		struct vector col;
		ssize_t j;

		for (j = 0; j < n; j++) {
			col = matrix_col(matrix, j);
			actors_add(actors, &col);
		}
	}
}

void actors_init_copy(struct actors *actors, const struct actors *src)
{
	assert(actors);
	assert(src);

	actors_init(actors, src->dim);

	ssize_t i, n = actors_count(src);
	for (i = 0; i < n; i++) {
		actors_add(actors, actors_traits(src, i));
	}

	assert(actors_count(actors) == actors_count(src));
	assert(actors_cohort_count(actors) == actors_cohort_count(src));
}

static void cohorts_clear(struct array *cohorts)
{
	struct cohort *c;
	ARRAY_FOREACH(c, cohorts) {
		cohort_deinit(c);
	}
	array_clear(cohorts);
}

static void trait_hashes_clear(struct intmap *trait_hashes)
{
	struct intmap_iter it;
	INTMAP_FOREACH(it, trait_hashes) {
		struct intset *set = INTMAP_VAL(it);
		intset_deinit(set);
	}
	intmap_clear(trait_hashes);
}

void actors_clear(struct actors *a)
{
	array_clear(&a->actors);
	cohorts_clear(&a->cohorts);
	trait_hashes_clear(&a->trait_hashes);
}

void actors_deinit(struct actors *a)
{
	assert(a);

	actors_clear(a);
	refcount_deinit(&a->refcount);
	intmap_deinit(&a->trait_hashes);
	array_deinit(&a->cohorts);
	array_deinit(&a->actors);
}

struct actors *actors_alloc(ssize_t dim)
{
	assert(dim >= 0);
	struct actors *actors = xcalloc(1, sizeof(*actors));
	actors_init(actors, dim);
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

void actors_add(struct actors *actors, const struct vector *traits)
{
	assert(actors);
	assert(vector_dim(traits) == actors_dim(actors));

	ssize_t aid = array_count(&actors->actors);
	struct actor *a;
	ssize_t cid = -1;
	struct cohort *c;
	int32_t hash = vector_hash(traits);
	struct intmap_pos pos;
	struct intset *set = intmap_find(&actors->trait_hashes, hash, &pos);

	if (!set) {
		set = intmap_insert(&actors->trait_hashes, &pos, NULL);
		intset_init(set);
	}

	struct intset_iter it;
	INTSET_FOREACH(it, set) {
		cid = INTSET_KEY(it);
		c = array_item(&actors->cohorts, cid);
		const struct vector *x = cohort_traits(c);
		if (vector_equals(x, traits)) {
			goto found;
		}
	}

	/* cohort not found; create a new one */
	cid = array_count(&actors->cohorts);
	c = array_add(&actors->cohorts, NULL);
	cohort_init(c, traits);
	intset_add(set, cid);
found:
	cohort_add(c, aid);
	a = array_add(&actors->actors, NULL);
	a->cohort = cid;
}

void actors_mul(double alpha, enum trans_op trans, const struct actors *a,
		const struct vector *x, double beta, struct vector *y)
{
	assert(a);
	assert(x);
	assert(y);
	assert(trans != TRANS_NOTRANS || vector_dim(x) == actors_dim(a));
	assert(trans != TRANS_NOTRANS || vector_dim(y) == actors_count(a));
	assert(trans == TRANS_NOTRANS || vector_dim(x) == actors_count(a));
	assert(trans == TRANS_NOTRANS || vector_dim(y) == actors_dim(a));

	const struct vector *row;
	double alpha_dot, scale;
	struct cohort *c;
	struct cohort_iter c_it;
	ssize_t id;

	if (beta == 0) {
		vector_fill(y, 0.0);
	} else if (beta != 1) {
		vector_scale(y, beta);
	}

	if (trans == TRANS_NOTRANS) {
		ARRAY_FOREACH(c, &a->cohorts) {
			row = cohort_traits(c);
			alpha_dot = alpha * vector_dot(row, x);

			COHORT_FOREACH(c_it, c) {
				id = COHORT_KEY(c_it);
				*vector_item_ptr(y, id) += alpha_dot;
			}
		}
	} else {
		ARRAY_FOREACH(c, &a->cohorts) {
			row = cohort_traits(c);
			scale = 0.0;

			COHORT_FOREACH(c_it, c) {
				id = COHORT_KEY(c_it);
				scale += *vector_item_ptr(x, id);
			}

			vector_axpy(alpha * scale, row, y);
		}
	}
}

void
actors_muls(double alpha,
	    enum trans_op trans,
	    const struct actors *a,
	    const struct svector *x, double beta, struct vector *y)
{
	assert(a);
	assert(x);
	assert(y);
	assert(trans != TRANS_NOTRANS || svector_dim(x) == actors_dim(a));
	assert(trans != TRANS_NOTRANS || vector_dim(y) == actors_count(a));
	assert(trans == TRANS_NOTRANS || svector_dim(x) == actors_count(a));
	assert(trans == TRANS_NOTRANS || vector_dim(y) == actors_dim(a));

	const struct vector *row;
	double alpha_dot, entry;
	struct cohort *c;
	struct cohort_iter c_it;
	ssize_t id;

	if (beta == 0) {
		vector_fill(y, 0.0);
	} else if (beta != 1) {
		vector_scale(y, beta);
	}

	if (trans == TRANS_NOTRANS) {
		ARRAY_FOREACH(c, &a->cohorts) {
			row = cohort_traits(c);
			alpha_dot = alpha * svector_dot(x, row);

			COHORT_FOREACH(c_it, c) {
				id = COHORT_KEY(c_it);
				*vector_item_ptr(y, id) += alpha_dot;
			}
		}
	} else {
		/* NOTE: this could potentially be made more effecient by
		 * using the cohort structure.  Below, we assume that
		 * the sparsity in x is more important.
		 */
		struct svector_iter itx;
		SVECTOR_FOREACH(itx, x) {
			id = SVECTOR_IDX(itx);
			entry = SVECTOR_VAL(itx);
			row = actors_traits(a, id);
			vector_axpy(alpha * entry, row, y);

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
	assert(trans != TRANS_NOTRANS || matrix_nrow(x) == actors_dim(a));
	assert(trans != TRANS_NOTRANS || matrix_nrow(y) == actors_count(a));
	assert(trans == TRANS_NOTRANS || matrix_nrow(x) == actors_count(a));
	assert(trans == TRANS_NOTRANS || matrix_nrow(y) == actors_dim(a));
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
