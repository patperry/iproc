#ifndef _SPARSETABLE_H
#define _SPARSETABLE_H

/* An intmap maps the integers 0..n-1 to values.  The implementation on
 * Google's "sparsetable".
 *
 * Define SSIZE_MAX, ssize_t, bool and assert before including this file.
 */

#include <stddef.h>		// sizeof, size_t
#include "darray.h"

/* private */
#define SPARSETABLE_GROUP_SIZE 48

struct sparsegroup {
	void *group;				// (small) array of values
	uint16_t num_buckets;			// limits GROUP_SIZE to 64K
	uint8_t bitmap[(SPARSETABLE_GROUP_SIZE-1)/8 + 1];	// fancy math is so we round up
	uint8_t deleted[(SPARSETABLE_GROUP_SIZE-1)/8 + 1];	// indicates if a position was ever deleted
};

struct sparsegroup_pos {
	ssize_t index;
	ssize_t offset;
};

struct sparsegroup_iter {
	struct sparsegroup_pos pos;
	void *val;
};


/* public */
struct sparsetable {
	struct darray groups;			// our list of groups
	ssize_t table_size;			// how many buckets they want
	ssize_t num_buckets;			// number of non-empty buckets
	size_t elt_size;
};

struct sparsetable_pos {
	ssize_t index;
	struct sparsegroup *group;
	struct sparsegroup_pos group_pos;
};

struct sparsetable_iter {
	ssize_t index;
	struct sparsegroup *group;
	struct sparsegroup_iter group_it;
};



/* constructors */
bool sparsetable_init(struct sparsetable *t, ssize_t n, size_t elt_size);
bool sparsetable_init_copy(struct sparsetable *t, const struct sparsetable *src);
void sparsetable_deinit(struct sparsetable *t);

/* assign, copy, clear */
bool sparsetable_assign_copy(struct sparsetable *t,
			     const struct sparsetable *src);
void sparsetable_clear(struct sparsetable *t);

/* informative */
ssize_t sparsetable_empty(const struct sparsetable *t); // returns size == 0, NOT count == 0
ssize_t sparsetable_count(const struct sparsetable *t);
ssize_t sparsetable_size(const struct sparsetable *t);
ssize_t sparsetable_max_size(const struct sparsetable *t);
size_t sparsetable_elt_size(const struct sparsetable *t);

bool sparsetable_contains(const struct sparsetable *t, ssize_t index);
void *sparsetable_lookup(const struct sparsetable *t, ssize_t index);
void *sparsetable_lookup_with(const struct sparsetable *t, ssize_t index,
			      const void *val0);

/* modification */
bool sparsetable_add(struct sparsetable *t, ssize_t index, const void *val);
bool sparsetable_add_all(struct sparsetable *t, ssize_t *indexes,
			 const void *vals, ssize_t n);
void sparsetable_remove(struct sparsetable *t, ssize_t index);
void sparsetable_remove_all(struct sparsetable *t,
			    ssize_t *indexes,
			    ssize_t n);
void sparsetable_resize(struct sparsetable *t, ssize_t n);

/* position-based interface */
void *sparsetable_find(const struct sparsetable *t, ssize_t index,
			struct sparsetable_pos *pos);
bool sparsetable_insert(struct sparsetable *t,
			const struct sparsetable_pos *pos,
			const void *val);
void sparsetable_replace(struct sparsetable *t,
			 const struct sparsetable_pos *pos,
			 const void *val);
void sparsetable_erase(struct sparsetable *t,
		       const struct sparsetable_pos *pos);
bool sparsetable_deleted(const struct sparsetable *t,
			 const struct sparsetable_pos *pos);

/* iteration */
void sparsetable_iter_init(const struct sparsetable *t,
			   struct sparsetable_iter *it);
void sparsetable_iter_deinit(const struct sparsetable *t,
			     struct sparsetable_iter *it);
void sparsetable_iter_reset(const struct sparsetable *t,
			    struct sparsetable_iter *it);
bool sparsetable_iter_advance(const struct sparsetable *t,
			      struct sparsetable_iter *it);
ssize_t sparsetable_iter_skip(const struct sparsetable *t,
			      struct sparsetable_iter *it);
void *sparsetable_iter_current(const struct sparsetable *t,
			       const struct sparsetable_iter *it);
ssize_t sparsetable_iter_index(const struct sparsetable *t,
			       const struct sparsetable_iter *it);



#endif /* _SPARSETABLE_H */
