#ifndef _SPARSETABLE_H
#define _SPARSETABLE_H

/* An intmap maps the integers 0..n-1 to values.  The implementation on
 * Google's "sparsetable".
 *
 * Define SSIZE_MAX, ssize_t, bool and assert before including this file.
 */

#include <stddef.h>		// sizeof, size_t
#include "darray.h"


struct sparsetable {
	struct darray groups;			// our list of groups
	ssize_t table_size;			// how many buckets they want
	ssize_t num_buckets;			// number of non-empty buckets
	size_t elt_size;
};

bool sparsetable_init(struct sparsetable *t, size_t elt_size);
void sparsetable_deinit(struct sparsetable *t);

#endif /* _SPARSETABLE_H */
