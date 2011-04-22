#include "port.h"
#include <assert.h>
#include <stdlib.h>

#include "sparsetable.h"

static ssize_t charbit(ssize_t i)
{
	return i >> 3;
}

static ssize_t modbit(ssize_t i)
{
	return 1 << (i & 7);
}

// We need a small function that tells us how many set bits there are
// in positions 0..i-1 of the bitmap.  It uses a big table.
// We make it static so templates don't allocate lots of these tables.
// There are lots of ways to do this calculation (called 'popcount').
// The 8-bit table lookup is one of the fastest, though this
// implementation suffers from not doing any loop unrolling.  See, eg,
//   http://www.dalkescientific.com/writings/diary/archive/2008/07/03/hakmem_and_other_popcounts.html
//   http://gurmeetsingh.wordpress.com/2008/08/05/fast-bit-counting-routines/
static ssize_t index_to_offset(const uint8_t *bm, ssize_t index)
{
	// We could make these ints.  The tradeoff is size (eg does it overwhelm
	// the cache?) vs efficiency in referencing sub-word-sized array elements
	static const uint8_t bits_in[256] = {	// # of bits set in one char
		0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4,
		1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
		1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
		2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
		1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
		2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
		2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
		3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
		1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
		2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
		2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
		3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
		2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
		3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
		3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
		4, 5, 5, 6, 5, 6, 6, 7, 5, 6, 6, 7, 6, 7, 7, 8,
	};
	ssize_t retval = 0;

	// [Note: condition index > 8 is an optimization; convince yourself we
	// give exactly the same result as if we had index >= 8 here instead.]
	for (; index > 8; index -= 8)	// bm[0..index/8-1]
		retval += bits_in[*bm++];	// chars we want *all* bits in
	return retval + bits_in[*bm & ((1 << index) - 1)];	// the char that includes index
}

static ssize_t offset_to_index(const uint8_t *bm, ssize_t offset)
{
#ifndef NDEBUG
	const uint8_t *const bm0 = bm;
	const ssize_t offset0 = offset;
#endif

	static const uint8_t bits_in[256] = {	// # of bits set in one char
		0, 1, 1, 2, 1, 2, 2, 3, 1, 2, 2, 3, 2, 3, 3, 4,
		1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
		1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
		2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
		1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
		2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
		2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
		3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
		1, 2, 2, 3, 2, 3, 3, 4, 2, 3, 3, 4, 3, 4, 4, 5,
		2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
		2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
		3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
		2, 3, 3, 4, 3, 4, 4, 5, 3, 4, 4, 5, 4, 5, 5, 6,
		3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
		3, 4, 4, 5, 4, 5, 5, 6, 4, 5, 5, 6, 5, 6, 6, 7,
		4, 5, 5, 6, 5, 6, 6, 7, 5, 6, 6, 7, 6, 7, 7, 8,
	};
	ssize_t index = 0;
	uint8_t rem;

	while (bits_in[*bm] <= offset) {
		offset -= bits_in[*bm++];
		index += 8;
	}

	if (offset == 0)
		goto done;

	rem = *bm;
	if (bits_in[rem & 0x0F] <= offset) {
		offset -= bits_in[rem & 0x0F];
		index += 4;
		rem >>= 4;
	}
	if (bits_in[rem & 0x03] <= offset) {
		offset -= bits_in[rem & 0x03];
		index += 2;
		rem >>= 2;
	}
	if ((rem & 0x01) <= offset) {
		offset -= rem & 0x01;
		index += 1;
		rem >>= 1;
	}
done:
	assert(offset == 0);
	assert(index_to_offset(bm0, index) == offset0);

	return index;
}

/* private functions */

/* bitmap manipulation */
static int sparsegroup_bmtest(const struct sparsegroup *g, ssize_t i);
static void sparsegroup_bmset(struct sparsegroup *g, ssize_t i);
static void sparsegroup_bmclear(struct sparsegroup *g, ssize_t i);
static int sparsegroup_dtest(const struct sparsegroup *g, ssize_t i);
static void sparsegroup_dset(struct sparsegroup *g, ssize_t i);
static void sparsegroup_dclear(struct sparsegroup *g, ssize_t i);

/* group alloc/free */
static bool sparsegroup_allocate_group(struct sparsegroup *g, ssize_t n,
				       size_t elt_size);
static bool sparsegroup_realloc_group(struct sparsegroup *g, ssize_t n,
				      size_t elt_size);
static void sparsegroup_free_group(struct sparsegroup *g, size_t elt_size);

/* indexing */
ssize_t sparsegroup_index_to_offset(const struct sparsegroup *g, ssize_t index);
ssize_t sparsegroup_offset_to_index(const struct sparsegroup *g,
				    ssize_t offset);

/* public functions */

/* constructors */
static bool sparsegroup_init(struct sparsegroup *g);
static bool sparsegroup_init_copy(struct sparsegroup *g,
				  const struct sparsegroup *src,
				  size_t elt_size);
static void sparsegroup_deinit(struct sparsegroup *g, size_t elt_size);

/* assign, clear */
static bool sparsegroup_assign_copy(struct sparsegroup *g,
				    const struct sparsegroup *src,
				    size_t elt_size);
static void sparsegroup_clear(struct sparsegroup *g, size_t elt_size);

/* informative */
static bool sparsegroup_empty(const struct sparsegroup *g);
static ssize_t sparsegroup_count(const struct sparsegroup *g);
static ssize_t sparsegroup_size(const struct sparsegroup *g);
static ssize_t sparsegroup_max_size(const struct sparsegroup *g);
static bool sparsegroup_contains(const struct sparsegroup *g,
				 ssize_t index, size_t elt_size);
static void *sparsegroup_lookup(const struct sparsegroup *g, ssize_t index,
				size_t elt_size);
static void *sparsegroup_lookup_with(const struct sparsegroup *g, ssize_t index,
				     const void *val0, size_t elt_size);

/* modification */
static bool sparsegroup_add(struct sparsegroup *g, ssize_t index,
			    const void *val, size_t elt_size);
static ssize_t sparsegroup_add_all(struct sparsegroup *g,
				   ssize_t *indexes,
				   const void *vals, ssize_t n, size_t elt_size);
static void sparsegroup_remove(struct sparsegroup *g, ssize_t index,
			       size_t elt_size);
static void sparsegroup_remove_preserve(struct sparsegroup *g, ssize_t index,
					size_t elt_size);
static void sparsegroup_remove_all(struct sparsegroup *g,
				   ssize_t *indexes,
				   ssize_t n, size_t elt_size);
static void sparsegroup_remove_range(struct sparsegroup *g, ssize_t i,
				     ssize_t n, size_t elt_size);

/* position-based interface */
static void *sparsegroup_find(const struct sparsegroup *g, ssize_t index,
			      struct sparsegroup_pos *pos, size_t elt_size);
static bool sparsegroup_insert(struct sparsegroup *g,
			       const struct sparsegroup_pos *pos,
			       const void *val, size_t elt_size);
static void sparsegroup_replace(struct sparsegroup *g,
				const struct sparsegroup_pos *pos,
				const void *val, size_t elt_size);
static void sparsegroup_erase(struct sparsegroup *g,
			      const struct sparsegroup_pos *pos,
			      size_t elt_size);
static bool sparsegroup_deleted(const struct sparsegroup *g,
				const struct sparsegroup_pos *pos);

/* iteration */
static void sparsegroup_iter_init(const struct sparsegroup *g,
				  struct sparsegroup_iter *it);
static void sparsegroup_iter_deinit(const struct sparsegroup *g,
				    struct sparsegroup_iter *it);
static void sparsegroup_iter_reset(const struct sparsegroup *g,
				   struct sparsegroup_iter *it);
static bool sparsegroup_iter_advance(const struct sparsegroup *g,
				     struct sparsegroup_iter *it,
				     size_t elt_size);
static ssize_t sparsegroup_iter_skip(const struct sparsegroup *g,
				     struct sparsegroup_iter *it,
				     size_t elt_size);
static void *sparsegroup_iter_current(const struct sparsegroup *g,
				      const struct sparsegroup_iter *it);
static ssize_t sparsegroup_iter_index(const struct sparsegroup *g,
				      const struct sparsegroup_iter *it);

/* bitmap manipulation */
int sparsegroup_bmtest(const struct sparsegroup *g, ssize_t i)
{
	return g->bitmap[charbit(i)] & modbit(i);
}

void sparsegroup_bmset(struct sparsegroup *g, ssize_t i)
{
	g->bitmap[charbit(i)] |= modbit(i);
}

void sparsegroup_bmclear(struct sparsegroup *g, ssize_t i)
{
	g->bitmap[charbit(i)] &= ~modbit(i);
}

int sparsegroup_dtest(const struct sparsegroup *g, ssize_t i)
{
	return g->deleted[charbit(i)] & modbit(i);
}

void sparsegroup_dset(struct sparsegroup *g, ssize_t i)
{
	g->deleted[charbit(i)] |= modbit(i);
}

void sparsegroup_dclear(struct sparsegroup *g, ssize_t i)
{
	g->deleted[charbit(i)] &= ~modbit(i);
}

/* group alloc/free */
bool sparsegroup_allocate_group(struct sparsegroup *g, ssize_t n,
				size_t elt_size)
{
	g->group = malloc(n * elt_size);
	return g->group;
}

bool sparsegroup_realloc_group(struct sparsegroup *g, ssize_t n,
			       size_t elt_size)
{
	void *retval = realloc(g->group, n * elt_size);
	if (retval)
		g->group = retval;
	return retval;
}

void sparsegroup_free_group(struct sparsegroup *g, size_t elt_size)
{
	free(g->group);
	g->group = NULL;
}

/* indexing */
ssize_t sparsegroup_index_to_offset(const struct sparsegroup *g, ssize_t index)
{
	return index_to_offset(g->bitmap, index);
}

ssize_t sparsegroup_offset_to_index(const struct sparsegroup *g, ssize_t offset)
{
	return offset_to_index(g->bitmap, offset);
}

/* constructors */
bool sparsegroup_init(struct sparsegroup *g)
{
	g->group = NULL;
	g->num_buckets = 0;
	memset(g->bitmap, 0, sizeof(g->bitmap));
	memset(g->deleted, 0, sizeof(g->deleted));
	return true;
}

static bool sparsegroup_init_copy(struct sparsegroup *g,
				  const struct sparsegroup *src,
				  size_t elt_size)
{
	if (src->num_buckets) {
		g->num_buckets = src->num_buckets;
		if (!sparsegroup_allocate_group(g, g->num_buckets, elt_size))
			return false;
	} else {
		g->group = NULL;
		g->num_buckets = 0;
	}

	memcpy(g->bitmap, src->bitmap, sizeof(g->bitmap));
	memcpy(g->deleted, src->deleted, sizeof(g->deleted));
	return true;
}

void sparsegroup_deinit(struct sparsegroup *g, size_t elt_size)
{
	sparsegroup_free_group(g, elt_size);
}

/* assign, clear */
bool sparsegroup_assign_copy(struct sparsegroup *g,
			     const struct sparsegroup *src, size_t elt_size)
{
	if (g == src)
		return true;
	if (src->num_buckets == 0) {
		sparsegroup_free_group(g, elt_size);
	} else if (sparsegroup_realloc_group(g, src->num_buckets, elt_size)) {
		memcpy(g->group, src->group, sizeof(g->group));
	} else {
		return false;
	}

	memcpy(g->bitmap, src->bitmap, sizeof(g->bitmap));
	memcpy(g->deleted, src->deleted, sizeof(g->deleted));
	g->num_buckets = src->num_buckets;
	return true;
}

void sparsegroup_clear(struct sparsegroup *g, size_t elt_size)
{
	sparsegroup_free_group(g, elt_size);
	memset(g->group, 0, sizeof(g->group));
	memset(g->deleted, 0, sizeof(g->deleted));
	g->num_buckets = 0;
}

/* informative */
bool sparsegroup_empty(const struct sparsegroup *g)
{
	return sparsegroup_size(g) == 0;
}

ssize_t sparsegroup_count(const struct sparsegroup *g)
{
	return g->num_buckets;
}

ssize_t sparsegroup_size(const struct sparsegroup *g)
{
	return SPARSETABLE_GROUP_SIZE;
}

ssize_t sparsegroup_max_size(const struct sparsegroup *g)
{
	return SPARSETABLE_GROUP_SIZE;
}

bool sparsegroup_contains(const struct sparsegroup *g, ssize_t index,
			  size_t elt_size)
{
	return sparsegroup_lookup(g, index, elt_size);
}

void *sparsegroup_lookup(const struct sparsegroup *g,
			 ssize_t index, size_t elt_size)
{
	return sparsegroup_lookup_with(g, index, NULL, elt_size);
}

static void *sparsegroup_lookup_with(const struct sparsegroup *g, ssize_t index,
				     const void *val0, size_t elt_size)
{
	struct sparsegroup_pos pos;
	void *val = sparsegroup_find(g, index, &pos, elt_size);
	return val ? val : (void *)val0;
}

/* modification */
static bool sparsegroup_add(struct sparsegroup *g, ssize_t index,
			    const void *val, size_t elt_size)
{
	struct sparsegroup_pos pos;

	if (sparsegroup_find(g, index, &pos, elt_size)) {
		sparsegroup_replace(g, &pos, val, elt_size);
		return true;
	} else {
		return sparsegroup_insert(g, &pos, val, elt_size);
	}
}

static ssize_t sparsegroup_add_all(struct sparsegroup *g,
				   ssize_t *indexes,
				   const void *vals, ssize_t n, size_t elt_size)
{
	assert(n >= 0);

	struct sparsegroup_pos pos;
	ssize_t i;

	for (i = 0; i < n; i++) {
		if (sparsegroup_find(g, indexes[i], &pos, elt_size)) {
			sparsegroup_replace(g, &pos, (char *)vals + i * elt_size,
					    elt_size);

		} else {
			if (!sparsegroup_insert(g, &pos, (char *)vals + i * elt_size, elt_size))
				break;
		}
	}
	return i;
}

static void sparsegroup_remove(struct sparsegroup *g, ssize_t index,
			       size_t elt_size)
{
	struct sparsegroup_pos pos;

	if (sparsegroup_find(g, index, &pos, elt_size)) {
		sparsegroup_erase(g, &pos, elt_size);
	}
}

static void sparsegroup_remove_preserve(struct sparsegroup *g, ssize_t index,
					size_t elt_size)
{
	int dbit = sparsegroup_dtest(g, index);
	sparsegroup_remove(g, index, elt_size);
	if (!dbit)
		sparsegroup_dclear(g, index);
}

static void sparsegroup_remove_all(struct sparsegroup *g,
				   ssize_t *indexes, ssize_t n, size_t elt_size)
{
	ssize_t i;
	for (i = 0; i < n; i++) {
		sparsegroup_remove(g, indexes[i], elt_size);
	}
}

// This could be more efficient, but then we'd need to figure
// out if we spanned groups or not.  Doesn't seem worth it.
void sparsegroup_remove_range(struct sparsegroup *g, ssize_t i, ssize_t n,
			      size_t elt_size)
{
	assert(n >= 0);
	assert(0 <= i && i <= sparsegroup_size(g) - n);

	struct sparsegroup_pos pos;

	for (; i < n; i++) {
		if (sparsegroup_find(g, i, &pos, elt_size))
			sparsegroup_erase(g, &pos, elt_size);
	}
}

/* position-based interface */
void *sparsegroup_find(const struct sparsegroup *g, ssize_t index,
		       struct sparsegroup_pos *pos, size_t elt_size)
{
	pos->index = index;
	pos->offset = sparsegroup_index_to_offset(g, index);

	if (sparsegroup_bmtest(g, index)) {
		return (char *)g->group + pos->offset * elt_size;
	}
	return NULL;
}

bool sparsegroup_insert(struct sparsegroup *g,
			const struct sparsegroup_pos *pos,
			const void *val, size_t elt_size)
{
	assert(!sparsegroup_bmtest(g, pos->index));

	if (!sparsegroup_realloc_group(g, g->num_buckets + 1, elt_size))
		return false;

	memmove((char *)g->group + (pos->offset + 1) * elt_size,
		(char *)g->group + pos->offset * elt_size,
		(g->num_buckets - pos->offset) * elt_size);
	g->num_buckets++;
	sparsegroup_bmset(g, pos->index);
	memcpy((char *)g->group + pos->offset * elt_size, val, elt_size);
	assert(sparsegroup_contains(g, pos->index, elt_size));
	return true;
}

void sparsegroup_replace(struct sparsegroup *g,
			 const struct sparsegroup_pos *pos,
			 const void *val, size_t elt_size)
{
	assert(sparsegroup_bmtest(g, pos->index));
	memcpy((char *)g->group + pos->offset * elt_size, val, elt_size);
}

void sparsegroup_erase(struct sparsegroup *g,
		       const struct sparsegroup_pos *pos, size_t elt_size)
{
	assert(sparsegroup_bmtest(g, pos->index));

	if (g->num_buckets == 1) {
		sparsegroup_free_group(g, elt_size);
		g->group = NULL;
	} else {
		memmove((char *)g->group + pos->offset * elt_size,
			(char *)g->group + (pos->offset + 1) * elt_size,
			(g->num_buckets - pos->offset - 1) * elt_size);
		sparsegroup_realloc_group(g, g->num_buckets - 1, elt_size);
	}
	g->num_buckets--;
	sparsegroup_bmclear(g, pos->index);
	sparsegroup_dset(g, pos->index);
}

// true if the index has ever been deleted
bool sparsegroup_deleted(const struct sparsegroup *g,
			 const struct sparsegroup_pos *pos)
{
	return sparsegroup_dtest(g, pos->index);
}

/* iteration */
static void sparsegroup_iter_init(const struct sparsegroup *g,
				  struct sparsegroup_iter *it)
{
	sparsegroup_iter_reset(g, it);
}

static void sparsegroup_iter_deinit(const struct sparsegroup *g,
				    struct sparsegroup_iter *it)
{
}

static void sparsegroup_iter_reset(const struct sparsegroup *g,
				   struct sparsegroup_iter *it)
{
	it->pos.index = -1;
	it->pos.offset = -1;
	it->val = NULL;
}

static bool sparsegroup_iter_advance(const struct sparsegroup *g,
				     struct sparsegroup_iter *it,
				     size_t elt_size)
{
	if (it->pos.index < sparsegroup_size(g)) {
		it->val = sparsegroup_find(g, it->pos.index + 1, &it->pos,
					   elt_size);
		return true;
	} else {
		it->pos.offset = sparsegroup_count(g);
		it->val = NULL;
		return false;
	}
}

static ssize_t sparsegroup_iter_skip(const struct sparsegroup *g,
				     struct sparsegroup_iter *it,
				     size_t elt_size)
{
	it->pos.offset++;
	if (it->pos.offset < sparsegroup_count(g)) {
		ssize_t index0 = it->pos.index;
		it->pos.index = sparsegroup_offset_to_index(g, it->pos.offset);
		it->val = (char *)g->group + it->pos.offset * elt_size;
		return it->pos.index - index0;
	} else {
		it->pos.index = sparsegroup_size(g);
		it->val = NULL;
		return 0;
	}
}

static void *sparsegroup_iter_current(const struct sparsegroup *g,
				      const struct sparsegroup_iter *it)
{
	assert(it->pos.index >= 0);
	assert(it->pos.index < sparsegroup_size(g));
	return it->val;
}

static ssize_t sparsegroup_iter_index(const struct sparsegroup *g,
				      const struct sparsegroup_iter *it)
{
	return it->pos.index;
}

////////////////////////////////////////////////////////////////////////////////

// SPARSE-TABLE
// ------------
// The idea is that a table with (logically) t buckets is divided
// into t/M *groups* of M buckets each.  (M is a constant set in
// SPARSETABLE_GROUP_SIZE for efficiency.)  Each group is stored sparsely.
// Thus, inserting into the table causes some array to grow, which is
// slow but still constant time.  Lookup involves doing a
// logical-position-to-sparse-position lookup, which is also slow but
// constant time.  The larger M is, the slower these operations are
// but the less overhead (slightly).
//
// To store the sparse array, we store a bitmap B, where B[i] = 1 iff
// bucket i is non-empty.  Then to look up bucket i we really look up
// array[# of 1s before i in B].  This is constant time for fixed M.
//

/* How to deal with the proper group */
static ssize_t sparsetable_num_groups(ssize_t n);
static uint16_t sparsetable_pos_in_group(ssize_t i);
static ssize_t sparsetable_group_num(ssize_t i);
static struct sparsegroup *sparsetable_which_group(const struct sparsetable *t,
						   ssize_t i);

ssize_t sparsetable_num_groups(ssize_t n)
{				// how many to hold n buckets
	return n == 0 ? 0 : ((n - 1) / SPARSETABLE_GROUP_SIZE) + 1;
}

uint16_t sparsetable_pos_in_group(ssize_t i)
{
	return (uint16_t)(i % SPARSETABLE_GROUP_SIZE);
}

ssize_t sparsetable_group_num(ssize_t i)
{
	return i / SPARSETABLE_GROUP_SIZE;
}

struct sparsegroup *sparsetable_which_group(const struct sparsetable *t,
					    ssize_t i)
{
	struct sparsegroup *groups = darray_front(&t->groups);
	return &groups[sparsetable_group_num(i)];
}

bool sparsetable_init(struct sparsetable *t, ssize_t n, size_t elt_size)
{
	if (darray_init(&t->groups, sizeof(struct sparsegroup))) {
		t->table_size = 0;
		t->num_buckets = 0;
		t->elt_size = elt_size;
		if (sparsetable_resize(t, n)) {
			return true;
		}
		darray_deinit(&t->groups);
		return false;
	}
	return false;
}

bool sparsetable_init_copy(struct sparsetable *t, const struct sparsetable *src)
{
	size_t elt_size = src->elt_size;
	ssize_t i, n = src->table_size;
	struct sparsegroup *groups;
	const struct sparsegroup *src_groups;

	if (!sparsetable_init(t, n, elt_size))
		goto fail_init;

	groups = darray_front(&t->groups);
	src_groups = darray_front(&src->groups);

	for (i = 0; i < n; i++) {
		if (!sparsegroup_assign_copy(groups + i, src_groups + i,
					     elt_size))
			goto fail_assign;
	}

	t->num_buckets = src->num_buckets;
	return true;
fail_assign:
	sparsetable_deinit(t);
fail_init:
	return false;
}

void sparsetable_deinit(struct sparsetable *t)
{
	struct sparsegroup *groups = darray_front(&t->groups);
	ssize_t i, n = darray_size(&t->groups);
	size_t elt_size = t->elt_size;

	for (i = 0; i < n; i++) {
		sparsegroup_deinit(groups + i, elt_size);
	}
	darray_deinit(&t->groups);
}

bool sparsetable_assign_copy(struct sparsetable *t,
			     const struct sparsetable *src)
{
	struct sparsetable dst;

	if (sparsetable_init_copy(&dst, src)) {
		sparsetable_deinit(t);
		*t = dst;
		return true;
	}
	return false;
}

void sparsetable_clear(struct sparsetable *t)
{
	struct sparsegroup *groups = darray_front(&t->groups);
	ssize_t i, n = t->table_size;
	size_t elt_size = t->elt_size;

	for (i = 0; i < n; i++) {
		sparsegroup_clear(groups + i, elt_size);
	}
}

ssize_t sparsetable_empty(const struct sparsetable *t)
{
	return sparsetable_size(t) == 0;
}

ssize_t sparsetable_count(const struct sparsetable *t)
{
	return t->num_buckets;
}

ssize_t sparsetable_size(const struct sparsetable *t)
{
	return t->table_size;
}

ssize_t sparsetable_max_size(const struct sparsetable *t)
{
	ssize_t max_size = darray_max_size(&t->groups) * SPARSETABLE_GROUP_SIZE;
	return (max_size <= 0 ? SSIZE_MAX : max_size);
}

size_t sparsetable_elt_size(const struct sparsetable *t)
{
	return t->elt_size;
}

bool sparsetable_resize(struct sparsetable *t, ssize_t n)
{
	ssize_t num_groups0 = darray_size(&t->groups);
	ssize_t n0 = t->table_size;
	ssize_t num_groups = sparsetable_num_groups(n);
	ssize_t group = sparsetable_group_num(n);
	ssize_t index = sparsetable_pos_in_group(n);
	size_t elt_size = sparsetable_elt_size(t);
	ssize_t i;
	bool resize_ok;

	if (num_groups < num_groups0) {
		struct sparsegroup *groups = darray_front(&t->groups);

		for (i = group + 1; i < num_groups0; i++) {
			sparsegroup_deinit(&groups[i], elt_size);
		}
	}
	// the resize always succeeds when we shrink
	resize_ok = darray_resize(&t->groups, sparsetable_num_groups(n));

	if (n < n0) {
		// lower num_buckets, clear last group

		if (index > 0) {	// need to clear inside last group
			struct sparsegroup *g = darray_back(&t->groups);
			sparsegroup_remove_range(g, index,
						 sparsegroup_size(g) - index,
						 t->elt_size);
		}

		t->num_buckets = 0;	// refigure # of used buckets
		struct sparsegroup *groups = darray_front(&t->groups);
		ssize_t i, n = darray_size(&t->groups);

		for (i = 0; i < n; i++) {
			t->num_buckets += sparsegroup_count(&groups[i]);
		}
	}

	if (resize_ok) {
		t->table_size = n;
	}
	return resize_ok;
}

bool sparsetable_contains(const struct sparsetable *t, ssize_t index)
{
	return sparsetable_lookup(t, index);
}

void *sparsetable_lookup(const struct sparsetable *t, ssize_t index)
{
	return sparsetable_lookup_with(t, index, NULL);
}

void *sparsetable_lookup_with(const struct sparsetable *t, ssize_t index,
			      const void *val0)
{
	struct sparsetable_pos pos;
	void *val = sparsetable_find(t, index, &pos);
	return val ? val : (void *)val0;
}

bool sparsetable_add(struct sparsetable *t, ssize_t index, const void *val)
{
	struct sparsetable_pos pos;

	if (sparsetable_find(t, index, &pos)) {
		sparsetable_replace(t, &pos, val);
		return true;
	} else {
		return sparsetable_insert(t, &pos, val);
	}
}

static void sparsetable_remove_preserve(struct sparsetable *t, ssize_t index)
{
	struct sparsetable_pos pos;
	int dbit;

	if (sparsetable_find(t, index, &pos)) {
		dbit = sparsegroup_dtest(pos.group, pos.group_pos.index);
		sparsetable_erase(t, &pos);
		if (!dbit)
			sparsegroup_dclear(pos.group, pos.group_pos.index);
	}
}

ssize_t sparsetable_add_all(struct sparsetable *t, ssize_t *indexes,
			    const void *vals, ssize_t n)
{
	if (n <= 0)
		return true;

	size_t elt_size = sparsetable_elt_size(t);
	struct sparsetable_pos pos;
	ssize_t i;

	for (i = 0; i < n; i++, vals = (char *)vals + elt_size) {
		if (sparsetable_find(t, indexes[i], &pos)) {
			sparsetable_replace(t, &pos, vals);
		} else {
			if (!sparsetable_insert(t, &pos, vals))
				break;
		}
	}
	return i;
}

void sparsetable_remove(struct sparsetable *t, ssize_t index)
{
	struct sparsetable_pos pos;

	if (sparsetable_find(t, index, &pos)) {
		sparsetable_erase(t, &pos);
	}
}

void sparsetable_remove_all(struct sparsetable *t, ssize_t *indexes, ssize_t n)
{
	ssize_t i;
	for (i = 0; i < n; i++) {
		sparsetable_remove(t, indexes[i]);
	}
}

void *sparsetable_find(const struct sparsetable *t, ssize_t index,
		       struct sparsetable_pos *pos)
{
	ssize_t group_num = sparsetable_group_num(index);
	ssize_t index_in_group = sparsetable_pos_in_group(index);
	size_t elt_size = sparsetable_elt_size(t);

	pos->index = index;
	pos->group = darray_at(&t->groups, group_num);
	return sparsegroup_find(pos->group, index_in_group, &pos->group_pos,
				elt_size);
}

bool sparsetable_insert(struct sparsetable *t,
			const struct sparsetable_pos *pos, const void *val)
{
	if (sparsegroup_insert(pos->group, &pos->group_pos, val,
			       sparsetable_elt_size(t))) {
		t->num_buckets++;
		assert(sparsetable_contains(t, pos->index));
		return true;
	}
	return false;
}

void sparsetable_replace(struct sparsetable *t,
			 const struct sparsetable_pos *pos, const void *val)
{
	sparsegroup_replace(pos->group, &pos->group_pos, val,
			    sparsetable_elt_size(t));
}

void sparsetable_erase(struct sparsetable *t, const struct sparsetable_pos *pos)
{
	sparsegroup_erase(pos->group, &pos->group_pos, sparsetable_elt_size(t));
	t->num_buckets--;
}

bool sparsetable_deleted(const struct sparsetable *t,
			 const struct sparsetable_pos *pos)
{
	return sparsegroup_deleted(pos->group, &pos->group_pos);
}

/* iteration */
void sparsetable_iter_init(const struct sparsetable *t,
			   struct sparsetable_iter *it)
{
	it->index = -1;
	it->group = darray_front(&t->groups);
	sparsegroup_iter_init(it->group, &it->group_it);
}

void sparsetable_iter_deinit(const struct sparsetable *t,
			     struct sparsetable_iter *it)
{
	sparsegroup_iter_deinit(it->group, &it->group_it);
}

void sparsetable_iter_reset(const struct sparsetable *t,
			    struct sparsetable_iter *it)
{
	sparsetable_iter_deinit(t, it);
	sparsetable_iter_init(t, it);
}

bool sparsetable_iter_advance(const struct sparsetable *t,
			      struct sparsetable_iter *it)
{
	size_t elt_size = sparsetable_elt_size(t);

	it->index++;

	if (it->index < sparsetable_size(t)) {
		if (!sparsegroup_iter_advance
		    (it->group, &it->group_it, elt_size)) {
			sparsegroup_iter_deinit(it->group, &it->group_it);
			it->group++;
			sparsegroup_iter_init(it->group, &it->group_it);
			sparsegroup_iter_advance(it->group, &it->group_it,
						 elt_size);
		}
		return true;
	}
	return false;
}

ssize_t sparsetable_iter_skip(const struct sparsetable *t,
			      struct sparsetable_iter *it)
{
	size_t elt_size = sparsetable_elt_size(t);
	ssize_t group_skip = 0;
	ssize_t group_index0;
	ssize_t index0 = it->index;

	if (it->index == sparsetable_size(t))
		return 0;

	group_index0 = sparsegroup_iter_index(it->group, &it->group_it);
	group_skip = sparsegroup_iter_skip(it->group, &it->group_it, elt_size);

	while (group_skip == 0) {
		it->index += sparsegroup_size(it->group) - group_index0 - 1;

		if (it->index < sparsetable_size(t)) {
			sparsegroup_iter_deinit(it->group, &it->group_it);
			it->group++;
			sparsegroup_iter_init(it->group, &it->group_it);
			group_index0 = -1;
			group_skip = sparsegroup_iter_skip(it->group,
							   &it->group_it,
							   elt_size);
		} else {
			it->index = sparsetable_size(t);
			return 0;
		}
	}

	it->index += group_skip;
	return it->index - index0;
}

void *sparsetable_iter_current(const struct sparsetable *t,
			       const struct sparsetable_iter *it)
{
	return sparsegroup_iter_current(it->group, &it->group_it);
}

ssize_t sparsetable_iter_index(const struct sparsetable *t,
			       const struct sparsetable_iter *it)
{
	return it->index;
}
