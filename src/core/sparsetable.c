#include "port.h"
#include <assert.h>
#include <stdlib.h>

#include "sparsetable.h"

#define GROUP_SIZE 48

struct sparsegroup_item {
	ssize_t index;
	void *val;
};

struct sparsegroup {
	void *group;				// (small) array of values
	uint16_t num_buckets;			// limits GROUP_SIZE to 64K
	uint8_t bitmap[(GROUP_SIZE-1)/8 + 1];	// fancy math is so we round up
};

struct sparsegroup_pos {
	ssize_t index;
	ssize_t offset;
};

struct sparsegroup_iter {
	ssize_t offset;
};

static ssize_t charbit(ssize_t i)  { return i >> 3; }
static ssize_t modbit(ssize_t i)   { return 1 << (i&7); }

// We need a small function that tells us how many set bits there are
// in positions 0..i-1 of the bitmap.  It uses a big table.
// We make it static so templates don't allocate lots of these tables.
// There are lots of ways to do this calculation (called 'popcount').
// The 8-bit table lookup is one of the fastest, though this
// implementation suffers from not doing any loop unrolling.  See, eg,
//   http://www.dalkescientific.com/writings/diary/archive/2008/07/03/hakmem_and_other_popcounts.html
//   http://gurmeetsingh.wordpress.com/2008/08/05/fast-bit-counting-routines/
static ssize_t index_to_offset(const uint8_t *bm, ssize_t index) {
	// We could make these ints.  The tradeoff is size (eg does it overwhelm
	// the cache?) vs efficiency in referencing sub-word-sized array elements
	static const uint8_t bits_in[256] = {      // # of bits set in one char
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
	for ( ; index > 8; index -= 8 )                    // bm[0..index/8-1]
		retval += bits_in[*bm++];              // chars we want *all* bits in
	return retval + bits_in[*bm & ((1 << index)-1)]; // the char that includes index
}

static ssize_t offset_to_index(const uint8_t *bm, ssize_t offset) {
#ifndef NDEBUG
	const uint8_t * const bm0 = bm;
	const ssize_t offset0 = offset;
#endif
	
	static const uint8_t bits_in[256] = {      // # of bits set in one char
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

/* group alloc/free */
static void *sparsegroup_allocate_group(const struct sparsegroup *g, ssize_t n,
					size_t elt_size);
static void *sparsegroup_realloc_group(struct sparsegroup *g, ssize_t n,
				       size_t elt_size);
static void sparsegroup_free_group(struct sparsegroup *g, size_t elt_size);

/* indexing */
ssize_t sparsegroup_index_to_offset(const struct sparsegroup *g, ssize_t index);
ssize_t sparsegroup_offset_to_index(const struct sparsegroup *g, ssize_t offset);


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
static ssize_t sparsegroup_size(const struct sparsegroup *g);
static ssize_t sparsegroup_max_size(const struct sparsegroup *g);
static bool sparsegroup_contains(const struct sparsegroup *g,
				 ssize_t index,
				 size_t elt_size);
static const void *sparsegroup_lookup(const struct sparsegroup *g,
				      ssize_t index,
				      size_t elt_size);

/* modification */
static bool sparsegroup_add(struct sparsegroup *g, ssize_t index,
			    const void *val, size_t elt_size);
static bool sparsegroup_add_all(struct sparsegroup *g,
				ssize_t *indexes,
				const void *vals,
				ssize_t n,
				size_t elt_size);
static void sparsegroup_remove(struct sparsegroup *g, ssize_t index,
			       size_t elt_size);
static void sparsegroup_remove_all(struct sparsegroup *g,
				   ssize_t *indexes,
				   ssize_t n,
				   size_t elt_size);

/* position-based interface */
static const void *sparsegroup_find(const struct sparsegroup *g, ssize_t index,
				    struct sparsegroup_pos *pos,
				    size_t elt_size);
static bool sparsegroup_insert(struct sparsegroup *g,
			       const struct sparsegroup_pos *pos,
			       const void *val,
			       size_t elt_size);
static void sparsegroup_replace(struct sparsegroup *g,
				const struct sparsegroup_pos *pos,
				const void *val,
				size_t elt_size);
static void sparsegroup_erase(struct sparsegroup *g,
			      const struct sparsegroup_pos *pos,
			      size_t elt_size);

/* iteration */
static void sparsegroup_iter_init(const struct sparsegroup *g,
				  struct sparsegroup_iter *it);
static void sparsegroup_iter_deinit(const struct sparsegroup *g,
				    struct sparsegroup_iter *it);
static void sparsegroup_iter_reset(const struct sparsegroup *g,
				   struct sparsegroup_iter *it);
static bool sparsegroup_iter_advance(const struct sparsegroup *g,
				     struct sparsegroup_iter *it);
static struct sparsegroup_item sparsegroup_iter_current(const struct sparsegroup *g,
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

/* group alloc/free */
void *sparsegroup_allocate_group(const struct sparsegroup *g, ssize_t n,
				 size_t elt_size)
{
	return malloc(n * elt_size);
}

void *sparsegroup_realloc_group(struct sparsegroup *g, ssize_t n,
				size_t elt_size)
{
	return realloc(g->group, n * elt_size);
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
	memset(&g->bitmap, 0, sizeof(g->bitmap));
	return true;
}

static bool sparsegroup_init_copy(struct sparsegroup *g,
				  const struct sparsegroup *src,
				  size_t elt_size)
{
	if (src->num_buckets) {
		g->num_buckets = src->num_buckets;		
		g->group = sparsegroup_allocate_group(g, g->num_buckets,
						      elt_size);
		if (!g->group)
			return false;
	} else {
		g->group = NULL;
		g->num_buckets = 0;
	}
	
	memcpy(g->bitmap, src->bitmap, sizeof(g->bitmap));
	return true;
}


void sparsegroup_deinit(struct sparsegroup *g, size_t elt_size)
{
	sparsegroup_free_group(g, elt_size);
}

/* assign, clear */
bool sparsegroup_assign_copy(struct sparsegroup *g,
			     const struct sparsegroup *src,
			     size_t elt_size)
{
	if (g == src) return true;
	if (src->num_buckets == 0) {
		sparsegroup_free_group(g, elt_size);
	} else if (sparsegroup_realloc_group(g, src->num_buckets, elt_size)) {
		memcpy(g->group, src->group, sizeof(g->group));
	} else {
		return false;
	}
	
	memcpy(g->bitmap, src->bitmap, sizeof(g->bitmap));
	g->num_buckets = src->num_buckets;
	return true;
}

void sparsegroup_clear(struct sparsegroup *g, size_t elt_size)
{
	sparsegroup_free_group(g, elt_size);
	memset(g->group, 0, sizeof(g->group));
	g->num_buckets = 0;
}


/* informative */
bool sparsegroup_empty(const struct sparsegroup *g)
{
	return g->num_buckets == 0;
}

ssize_t sparsegroup_size(const struct sparsegroup *g)
{
	return g->num_buckets;
}

ssize_t sparsegroup_max_size(const struct sparsegroup *g)
{
	return GROUP_SIZE;
}

bool sparsegroup_contains(const struct sparsegroup *g, ssize_t index,
			  size_t elt_size)
{
	return sparsegroup_lookup(g, index, elt_size);
}

const void *sparsegroup_lookup(const struct sparsegroup *g,
			       ssize_t index,
			       size_t elt_size)
{
	struct sparsegroup_pos pos;
	return sparsegroup_find(g, index, &pos, elt_size);
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

static bool sparsegroup_add_all(struct sparsegroup *g,
				ssize_t *indexes,
				const void *vals,
				ssize_t n,
				size_t elt_size)
{
	if (n <= 0)
		return true;
	
	char oldvals[n * elt_size];
	const void *oldval;
	bool exists[n];
	struct sparsegroup_pos pos;
	ssize_t i;
	
	for (i = 0; i < n; i++) {
		if ((oldval = sparsegroup_find(g, indexes[i], &pos, elt_size))) {
			exists[i] = true;
			memcpy(oldvals + i * elt_size, oldval, elt_size);
			sparsegroup_replace(g, &pos, vals + i * elt_size, elt_size);
			
		} else {
			exists[i] = false;
			if (!sparsegroup_insert(g, &pos, vals + i * elt_size, elt_size))
				goto rollback;
		}
	}
	return true;
rollback:
	for(; i > 0; i--) {
		if (exists[i-1]) {
			sparsegroup_add(g, indexes[i-1], oldvals + (i-1) * elt_size, elt_size);
		} else {
			sparsegroup_remove(g, indexes[i-1], elt_size);
		}
	}
	return false;
}


static void sparsegroup_remove(struct sparsegroup *g, ssize_t index,
			       size_t elt_size)
{
	struct sparsegroup_pos pos;

	if (sparsegroup_find(g, index, &pos, elt_size)) {
		sparsegroup_erase(g, &pos, elt_size);
	}
}

static void sparsegroup_remove_all(struct sparsegroup *g,
				   ssize_t *indexes,
				   ssize_t n,
				   size_t elt_size)
{
	ssize_t i;
	for (i = 0; i < n; i++) {
		sparsegroup_remove(g, indexes[i], elt_size);
	}
}


/* position-based interface */
const void *sparsegroup_find(const struct sparsegroup *g, ssize_t index,
			     struct sparsegroup_pos *pos, size_t elt_size)
{
	pos->index = index;
	pos->offset = sparsegroup_index_to_offset(g, index);
	
	if (sparsegroup_bmtest(g, index)) {
		return NULL;
	} else {
		return g->group + pos->offset * elt_size;
	}
}

bool sparsegroup_insert(struct sparsegroup *g,
			const struct sparsegroup_pos *pos,
			const void *val,
			size_t elt_size)
{
	assert(!sparsegroup_bmtest(g, pos->index));
	
	if (!sparsegroup_realloc_group(g, g->num_buckets + 1, elt_size))
		return false;
	
	memmove(g->group + (pos->offset + 1) * elt_size,
		g->group + pos->offset * elt_size,
		(g->num_buckets - pos->offset) * elt_size);
	g->num_buckets++;
	sparsegroup_bmset(g, pos->index);
	memcpy(g->group + pos->offset * elt_size, val, elt_size);
	return true;
}


void sparsegroup_replace(struct sparsegroup *g,
			 const struct sparsegroup_pos *pos,
			 const void *val,
			 size_t elt_size)
{
	assert(sparsegroup_bmtest(g, pos->index));
	memcpy(g->group + pos->offset * elt_size, val, elt_size);
}

void sparsegroup_erase(struct sparsegroup *g,
		       const struct sparsegroup_pos *pos,
		       size_t elt_size)
{
	assert(sparsegroup_bmtest(g, pos->index));
	
	if (g->num_buckets == 1) {
		sparsegroup_free_group(g, elt_size);
		g->group = NULL;
	} else {
		memmove(g->group + pos->offset * elt_size,
			g->group + (pos->offset + 1) * elt_size,
			(g->num_buckets - pos->offset - 1) * elt_size);
		sparsegroup_realloc_group(g, g->num_buckets - 1,
					  elt_size);
	}
	g->num_buckets--;
	sparsegroup_bmclear(g, pos->index);
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
	it->offset = -1;
}

static bool sparsegroup_iter_advance(const struct sparsegroup *g,
				     struct sparsegroup_iter *it)
{
	it->offset++;
	return it->offset < g->num_buckets;
}

static struct sparsegroup_item sparsegroup_iter_current(const struct sparsegroup *g,
							struct sparsegroup_iter *it,
							size_t elt_size)
{
	assert(it->offset >= 0);
	assert(it->offset < sparsegroup_size(g));
	
	struct sparsegroup_item item;
	
	item.index = sparsegroup_offset_to_index(g, it->offset);
	item.val = g->group + item.index * elt_size;
	return item;
}
