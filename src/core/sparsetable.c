#include "port.h"
#include <assert.h>
#include <stdlib.h>

#include "util.h"
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

/* group alloc/free */
static void sparsegroup_realloc_group(struct sparsegroup *g, ssize_t n,
				      size_t elt_size);
static void sparsegroup_free_group(struct sparsegroup *g, size_t elt_size);

/* indexing */
static ssize_t sparsegroup_index_to_offset(const struct sparsegroup *g,
					   ssize_t index);
static ssize_t sparsegroup_offset_to_index(const struct sparsegroup *g,
					   ssize_t offset);

/* public functions */

/* constructors */
static void sparsegroup_init(struct sparsegroup *g);
static void sparsegroup_deinit(struct sparsegroup *g, size_t elt_size);

/* assign, clear */
static void sparsegroup_assign_copy(struct sparsegroup *g,
				    const struct sparsegroup *src,
				    size_t elt_size);
static void sparsegroup_clear(struct sparsegroup *g, size_t elt_size);

/* informative */
static ssize_t sparsegroup_count(const struct sparsegroup *g);
static ssize_t sparsegroup_size(const struct sparsegroup *g);

/* modification */
static ssize_t sparsegroup_remove_range(struct sparsegroup *g, ssize_t i,
					ssize_t n, size_t elt_size);

/* position-based interface */
static void *sparsegroup_find(const struct sparsegroup *g, ssize_t index,
			      struct sparsegroup_pos *pos, size_t elt_size);
static void *sparsegroup_insert(struct sparsegroup *g,
				const struct sparsegroup_pos *pos,
				const void *val, size_t elt_size);
static void *sparsegroup_replace(struct sparsegroup *g,
				 const struct sparsegroup_pos *pos,
				 const void *val, size_t elt_size);
static void sparsegroup_erase(struct sparsegroup *g,
			      const struct sparsegroup_pos *pos,
			      size_t elt_size);
static bool sparsegroup_deleted(const struct sparsegroup *g,
				const struct sparsegroup_pos *pos);

/* iteration */
static struct sparsegroup_iter sparsegroup_iter_make(const struct sparsegroup
						     *g, size_t elt_size);
static void sparsegroup_iter_reset(const struct sparsegroup *g,
				   struct sparsegroup_iter *it,
				   size_t elt_size);
static bool sparsegroup_iter_advance(const struct sparsegroup *g,
				     struct sparsegroup_iter *it,
				     size_t elt_size);

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

/* group alloc/free */
void sparsegroup_realloc_group(struct sparsegroup *g, ssize_t n,
			       size_t elt_size)
{
	g->group = xrealloc(g->group, n * elt_size);
}

void sparsegroup_free_group(struct sparsegroup *g, size_t elt_size)
{
	xfree(g->group);
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
void sparsegroup_init(struct sparsegroup *g)
{
	g->group = NULL;
	g->num_buckets = 0;
	memset(g->bitmap, 0, sizeof(g->bitmap));
	memset(g->deleted, 0, sizeof(g->deleted));
}

void sparsegroup_deinit(struct sparsegroup *g, size_t elt_size)
{
	sparsegroup_free_group(g, elt_size);
}

/* assign, clear */
void sparsegroup_assign_copy(struct sparsegroup *g,
			     const struct sparsegroup *src, size_t elt_size)
{
	if (g == src)
		return;
	if (src->num_buckets == 0) {
		sparsegroup_free_group(g, elt_size);
	} else {
		sparsegroup_realloc_group(g, src->num_buckets, elt_size);
		memcpy(g->group, src->group, sizeof(g->group));
	}

	memcpy(g->bitmap, src->bitmap, sizeof(g->bitmap));
	memcpy(g->deleted, src->deleted, sizeof(g->deleted));
	g->num_buckets = src->num_buckets;
}

void sparsegroup_clear(struct sparsegroup *g, size_t elt_size)
{
	sparsegroup_free_group(g, elt_size);
	memset(g->bitmap, 0, sizeof(g->bitmap));
	memset(g->deleted, 0, sizeof(g->deleted));
	g->num_buckets = 0;
}

/* informative */
ssize_t sparsegroup_count(const struct sparsegroup *g)
{
	return g->num_buckets;
}

ssize_t sparsegroup_size(const struct sparsegroup *g)
{
	return SPARSETABLE_GROUP_SIZE;
}

// This could be more efficient, but then we'd need to figure
// out if we spanned groups or not.  Doesn't seem worth it.
ssize_t sparsegroup_remove_range(struct sparsegroup *g, ssize_t i, ssize_t n,
				 size_t elt_size)
{
	assert(n >= 0);
	assert(0 <= i && i <= sparsegroup_size(g) - n);

	ssize_t res = 0;
	struct sparsegroup_pos pos;

	for (; i < n; i++) {
		if (sparsegroup_find(g, i, &pos, elt_size)) {
			sparsegroup_erase(g, &pos, elt_size);
			res++;
		}
	}
	return res;
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

void *sparsegroup_insert(struct sparsegroup *g,
			 const struct sparsegroup_pos *pos,
			 const void *val, size_t elt_size)
{
	assert(!sparsegroup_bmtest(g, pos->index));

	sparsegroup_realloc_group(g, g->num_buckets + 1, elt_size);

	memmove((char *)g->group + (pos->offset + 1) * elt_size,
		(char *)g->group + pos->offset * elt_size,
		(g->num_buckets - pos->offset) * elt_size);
	g->num_buckets++;
	sparsegroup_bmset(g, pos->index);

	void *res = (char *)g->group + pos->offset * elt_size;

	if (val) {
		memcpy(res, val, elt_size);
	}

	return res;
}

void *sparsegroup_replace(struct sparsegroup *g,
			  const struct sparsegroup_pos *pos,
			  const void *val, size_t elt_size)
{
	assert(sparsegroup_bmtest(g, pos->index));

	void *res = (char *)g->group + pos->offset * elt_size;

	if (val)
		memcpy(res, val, elt_size);

	return res;
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
static struct sparsegroup_iter sparsegroup_iter_make(const struct sparsegroup
						     *g, size_t elt_size)
{
	struct sparsegroup_iter it;
	sparsegroup_iter_reset(g, &it, elt_size);
	return it;
}

static void sparsegroup_iter_reset(const struct sparsegroup *g,
				   struct sparsegroup_iter *it, size_t elt_size)
{
	it->val = (char *)g->group - elt_size;
	it->pos.index = -1;
	it->pos.offset = -1;
}

static bool sparsegroup_iter_advance(const struct sparsegroup *g,
				     struct sparsegroup_iter *it,
				     size_t elt_size)
{
	it->pos.offset++;
	if (it->pos.offset < sparsegroup_count(g)) {
		it->pos.index = sparsegroup_offset_to_index(g, it->pos.offset);
		it->val = (char *)it->val + elt_size;
		assert(it->val == (char *)g->group + it->pos.offset * elt_size);
		return true;
	} else {
		it->pos.index = sparsegroup_size(g);
		it->val = NULL;
		return false;
	}
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
	struct sparsegroup *groups = array_item(&t->groups, 0);
	return &groups[sparsetable_group_num(i)];
}

void sparsetable_init(struct sparsetable *t, ssize_t n, size_t elt_size)
{
	array_init(&t->groups, sizeof(struct sparsegroup));
	t->table_size = 0;
	t->num_buckets = 0;
	t->elt_size = elt_size;
	sparsetable_set_size(t, n);
}

void sparsetable_init_copy(struct sparsetable *t, const struct sparsetable *src)
{
	size_t elt_size = src->elt_size;
	ssize_t i, n = src->table_size;
	struct sparsegroup *groups;
	const struct sparsegroup *src_groups;

	sparsetable_init(t, n, elt_size);

	if (n > 0) {
		groups = array_item(&t->groups, 0);
		src_groups = array_item(&src->groups, 0);

		for (i = 0; i < n; i++) {
			sparsegroup_assign_copy(groups + i, src_groups + i,
						elt_size);
		}
	}

	t->num_buckets = src->num_buckets;
}

void sparsetable_deinit(struct sparsetable *t)
{
	struct sparsegroup *groups = array_item(&t->groups, 0);
	ssize_t i, n = array_count(&t->groups);
	size_t elt_size = t->elt_size;

	for (i = 0; i < n; i++) {
		sparsegroup_deinit(groups + i, elt_size);
	}
	array_deinit(&t->groups);
}

void sparsetable_assign_copy(struct sparsetable *t,
			     const struct sparsetable *src)
{
	sparsetable_deinit(t);
	sparsetable_init_copy(t, src);
}

void sparsetable_clear(struct sparsetable *t)
{
	struct sparsegroup *groups = array_item(&t->groups, 0);
	ssize_t i, n = array_count(&t->groups);
	size_t elt_size = t->elt_size;

	for (i = 0; i < n; i++) {
		sparsegroup_clear(groups + i, elt_size);
	}

	t->num_buckets = 0;
}

ssize_t sparsetable_count(const struct sparsetable *t)
{
	return t->num_buckets;
}

ssize_t sparsetable_size(const struct sparsetable *t)
{
	return t->table_size;
}

size_t sparsetable_elt_size(const struct sparsetable *t)
{
	return t->elt_size;
}

void sparsetable_set_size(struct sparsetable *t, ssize_t n)
{
	ssize_t num_groups0 = array_count(&t->groups);
	ssize_t n0 = t->table_size;
	ssize_t num_groups = sparsetable_num_groups(n);
	ssize_t index = sparsetable_pos_in_group(n);
	size_t elt_size = sparsetable_elt_size(t);
	ssize_t i;

	if (num_groups <= num_groups0) {
		struct sparsegroup *g;

		for (i = num_groups; i < num_groups0; i++) {
			g = array_item(&t->groups, i);
			assert(!sparsegroup_count(g));
			sparsegroup_deinit(g, elt_size);
		}

		array_remove_range(&t->groups, num_groups,
				   num_groups0 - num_groups);
	} else {
		array_add_range(&t->groups, NULL, num_groups - num_groups0);
		// no need to init groups since add_range zeros space
	}

	if (n < n0) {
		// lower num_buckets, clear last group

		if (index > 0) {	// make sure nothing remains in last group
			struct sparsegroup *g =
			    array_item(&t->groups, array_count(&t->groups) - 1);
			assert(!sparsegroup_remove_range(g, index,
							 sparsegroup_size(g) -
							 index, t->elt_size));
		}
	}

	t->table_size = n;
}

void *sparsetable_find(const struct sparsetable *t, ssize_t index,
		       struct sparsetable_pos *pos)
{
	ssize_t group_num = sparsetable_group_num(index);
	ssize_t index_in_group = sparsetable_pos_in_group(index);
	size_t elt_size = sparsetable_elt_size(t);

	pos->index = index;
	pos->group = array_item(&t->groups, group_num);
	return sparsegroup_find(pos->group, index_in_group, &pos->group_pos,
				elt_size);
}

void *sparsetable_insert(struct sparsetable *t,
			 const struct sparsetable_pos *pos, const void *val)
{
	void *res;
	if ((res = sparsegroup_insert(pos->group, &pos->group_pos, val,
				      sparsetable_elt_size(t)))) {
		t->num_buckets++;
	}
	return res;
}

void *sparsetable_replace(struct sparsetable *t,
			  const struct sparsetable_pos *pos, const void *val)
{
	return sparsegroup_replace(pos->group, &pos->group_pos, val,
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
struct sparsetable_iter sparsetable_iter_make(const struct sparsetable *t)
{
	assert(t);

	struct sparsetable_iter it;
	it.table = t;
	sparsetable_iter_reset(&it);
	return it;
}

void sparsetable_iter_reset(struct sparsetable_iter *it)
{
	size_t elt_size = sparsetable_elt_size(it->table);

	it->index = -1;
	it->group = array_item(&it->table->groups, 0);
	it->group_it = sparsegroup_iter_make(it->group, elt_size);
}

bool sparsetable_iter_advance(struct sparsetable_iter *it)
{
	size_t elt_size = sparsetable_elt_size(it->table);
	bool group_adv;
	ssize_t group_idx0, skip;

	if (it->index == sparsetable_size(it->table))
		return false;

	group_idx0 = SPARSEGROUP_IDX(it->group_it);
	group_adv =
	    sparsegroup_iter_advance(it->group, &it->group_it, elt_size);

	while (!group_adv) {
		it->index += sparsegroup_size(it->group) - group_idx0 - 1;

		if (it->index < sparsetable_size(it->table)) {
			it->group++;
			it->group_it =
			    sparsegroup_iter_make(it->group, elt_size);
			group_idx0 = -1;
			group_adv = sparsegroup_iter_advance(it->group,
							     &it->group_it,
							     elt_size);
		} else {
			it->index = sparsetable_size(it->table);
			return false;
		}
	}

	skip = SPARSEGROUP_IDX(it->group_it) - group_idx0;
	it->index += skip;
	assert(skip > 0);
	return true;
}
