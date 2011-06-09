#include "port.h"

#include <assert.h>
#include <stdlib.h>
#include "actors.h"

static uint32_t cohortp_traits_hash(const void *cohortp)
{
	return vector_hash(cohort_traits(*(struct cohort **)cohortp));
}

static bool cohortp_traits_equals(const void *cohortp1, const void *cohortp2)
{
	return vector_equals(cohort_traits(*(struct cohort **)cohortp1),
			     cohort_traits(*(struct cohort **)cohortp2));
}

void actors_init(struct actors *actors, ssize_t dim)
{
	assert(actors);
	assert(dim >= 0);

	array_init(&actors->actors, sizeof(struct actor));

	hashset_init(&actors->cohorts, cohortp_traits_hash,
		     cohortp_traits_equals, sizeof(struct cohort *));

	refcount_init(&actors->refcount);
	actors->dim = dim;
}

void actors_init_matrix(struct actors *actors, const struct matrix *matrix,
			enum trans_op trans)
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

	ssize_t i, n = actors_size(src);
	for (i = 0; i < n; i++) {
		actors_add(actors, actors_traits(src, i));
	}

	assert(actors_size(actors) == actors_cohorts_size(src));
	assert(actors_cohorts_size(actors) == actors_cohorts_size(src));
}

static void cohorts_clear(struct hashset *cohorts)
{
	struct hashset_iter it;
	HASHSET_FOREACH(it, cohorts) {
		struct cohort **cp = HASHSET_KEY(it);
		cohort_free(*cp);
	}
	hashset_clear(cohorts);
}

void actors_clear(struct actors *a)
{
	cohorts_clear(&a->cohorts);
	array_clear(&a->actors);
}

void actors_deinit(struct actors *a)
{
	assert(a);

	actors_clear(a);
	refcount_deinit(&a->refcount);
	hashset_deinit(&a->cohorts);
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
	if (refcount_get(&actors->refcount))
		return actors;
	return NULL;
}

ssize_t actors_size(const struct actors *a)
{
	assert(a);
	return array_count(&a->actors);
}

ssize_t actors_cohorts_size(const struct actors *a)
{
	assert(a);
	return hashset_count(&a->cohorts);
}

ssize_t actors_dim(const struct actors *actors)
{
	assert(actors);
	return actors->dim;
}

const struct actor *actors_at(const struct actors *a, ssize_t actor_id)
{
	assert(a);
	assert(0 <= actor_id && actor_id < actors_size(a));
	return array_item(&a->actors, actor_id);
}

const struct cohort *actors_cohort(const struct actors *a, ssize_t actor_id)
{
	return actors_at(a, actor_id)->cohort;
}

static struct cohort *actors_get_cohort(struct actors *actors,
					const struct vector *traits)
{
	struct cohort *key = container_of(traits, struct cohort, traits);
	struct hashset_pos pos;
	struct cohort *new_cohort;
	struct cohort **cohortp;

	if ((cohortp = hashset_find(&actors->cohorts, &key, &pos)))
		return *cohortp;

	if ((new_cohort = cohort_new(traits))) {
		if (hashset_insert(&actors->cohorts, &pos, &new_cohort)) {
			return new_cohort;
		}
		cohort_free(new_cohort);
	}
	return NULL;
}

double actors_get(const struct actors *a, ssize_t actor_id, ssize_t j)
{
	assert(a);
	assert(0 <= actor_id && actor_id < actors_size(a));
	assert(0 <= j && j < actors_dim(a));

	const struct vector *x = actors_traits(a, actor_id);
	double val = vector_item(x, j);
	return val;
}

bool actors_add(struct actors *actors, const struct vector *traits)
{
	assert(actors);
	assert(vector_dim(traits) == actors_dim(actors));

	struct actor a;
	ssize_t id;
	bool ok;

	a.cohort = actors_get_cohort(actors, traits);
	if (!a.cohort)
		goto fail_get_cohort;

	array_add(&actors->actors, &a);

	id = array_count(&actors->actors) - 1;
	cohort_add(a.cohort, id);
	ok = true;
	goto out;

fail_get_cohort:
	ok = false;
out:
	return ok;
}

const struct vector *actors_traits(const struct actors *a, ssize_t actor_id)
{
	assert(a);
	assert(0 <= actor_id && actor_id < actors_size(a));

	const struct actor *actor = actors_at(a, actor_id);
	return cohort_traits(actor->cohort);
}

void actors_mul(double alpha, enum trans_op trans, const struct actors *a,
		const struct vector *x, double beta, struct vector *y)
{
	assert(a);
	assert(x);
	assert(y);
	assert(trans != TRANS_NOTRANS || vector_dim(x) == actors_dim(a));
	assert(trans != TRANS_NOTRANS || vector_dim(y) == actors_size(a));
	assert(trans == TRANS_NOTRANS || vector_dim(x) == actors_size(a));
	assert(trans == TRANS_NOTRANS || vector_dim(y) == actors_dim(a));

	const struct vector *row;
	double alpha_dot, scale;
	struct hashset_iter it;
	struct cohort *c;
	struct cohort_iter c_it;
	ssize_t id;

	if (beta == 0) {
		vector_fill(y, 0.0);
	} else if (beta != 1) {
		vector_scale(y, beta);
	}

	if (trans == TRANS_NOTRANS) {
		HASHSET_FOREACH(it, &a->cohorts) {
			c = *(struct cohort **)HASHSET_KEY(it);
			row = cohort_traits(c);
			alpha_dot = alpha * vector_dot(row, x);

			cohort_iter_init(c, &c_it);
			while (cohort_iter_advance(c, &c_it)) {
				id = cohort_iter_current(c, &c_it);
				*vector_item_ptr(y, id) += alpha_dot;
			}
			cohort_iter_deinit(c, &c_it);
		}
	} else {
		HASHSET_FOREACH(it, &a->cohorts) {
			c = *(struct cohort **)HASHSET_KEY(it);
			row = cohort_traits(c);
			scale = 0.0;

			cohort_iter_init(c, &c_it);
			while (cohort_iter_advance(c, &c_it)) {
				id = cohort_iter_current(c, &c_it);
				scale += *vector_item_ptr(x, id);
			}
			cohort_iter_deinit(c, &c_it);

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
	assert(trans != TRANS_NOTRANS || vector_dim(y) == actors_size(a));
	assert(trans == TRANS_NOTRANS || svector_dim(x) == actors_size(a));
	assert(trans == TRANS_NOTRANS || vector_dim(y) == actors_dim(a));

	const struct vector *row;
	double alpha_dot, entry;
	struct hashset_iter it;
	struct cohort *c;
	struct cohort_iter c_it;
	ssize_t id;

	if (beta == 0) {
		vector_fill(y, 0.0);
	} else if (beta != 1) {
		vector_scale(y, beta);
	}

	if (trans == TRANS_NOTRANS) {
		HASHSET_FOREACH(it, &a->cohorts) {
			c = *(struct cohort **)HASHSET_KEY(it);
			row = cohort_traits(c);
			alpha_dot = alpha * svector_dot(x, row);

			cohort_iter_init(c, &c_it);
			while (cohort_iter_advance(c, &c_it)) {
				id = cohort_iter_current(c, &c_it);
				*vector_item_ptr(y, id) += alpha_dot;
			}
			cohort_iter_deinit(c, &c_it);
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
	assert(trans != TRANS_NOTRANS || matrix_nrow(y) == actors_size(a));
	assert(trans == TRANS_NOTRANS || matrix_nrow(x) == actors_size(a));
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
