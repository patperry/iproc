#ifndef _SARRAY_H
#define _SARRAY_H

/* Sparse Array type
 *
 * Supports O(log n) insert and delete
 *
 * Define SSIZE_MAX, ssize_t, bool and assert before including this file.
 */

#include <stddef.h>    // sizeof, size_t
#include <string.h>    // memcpy
#include "darray.h"


struct intmap {
    struct darray keys;    
    struct darray values;
};

struct intmap {
    intptr_t key;
    ssize_t  index;
};

void * sarray_lookup (const struct sarray *a, ssize_t key, struct sarray_pos *pos);

static inline bool sarray_pos_exists (const struct sarray_pos *pos)
{
    assert(pos);
    return (pos->index >= 0);
}

static inline bool sarray_pos_key (const struct sarray_pos *pos)
{
    assert(pos);
    return pos->key;
}

void * sarray_insert (struct sarray *a, struct sarray_pos *pos, const void *val);


#ifdef NDEBUG
# define sarray_index(a,t,pos) \
         darray_index(&(a)->values, t, (pos)->index)
#else
# define sarray_index(a,t,pos) \
         (*(t *)_sarray_index(a, pos))
#endif


static inline void * _sarray_index (struct sarray *a, const struct sarray_pos *pos)
{
    assert(sarray_pos_exists(pos));
    return darray_ptr(&a->values, pos->index);
}


/* create, destroy */
#define               sarray_init(a,t)    (_sarray_init(a,sizeof(t)))
struct sarray *       sarray_init_copy    (struct sarray       *a,
                                           const struct sarray *src);
void                  sarray_deinit       (struct sarray *a);


/* assignment, copy */
struct sarray *       sarray_assign_copy     (struct sarray       *a,
                                              const struct sarray *src);
void *                sarray_copy_to         (const struct sarray *a,
                                              void                *dst);
ssize_t *             sarray_copy_pattern_to (const struct sarray *a,
                                              ssize_t             *dst);

/* index */
#define sarray_index(a,t,i)           (*((t *)_sarray_index(a, i)))
#define sarray_index_with(a,t,i,val0) (*((t *)_sarray_index_with(a, i, val0)))

/* informative */
static inline ssize_t sarray_size      (const struct sarray *a);
static inline bool    sarray_empty     (const struct sarray *a);
static inline ssize_t sarray_max_size  (const struct sarray *a);
static inline size_t  sarray_elt_size  (const struct sarray *a);

/* operations */
void * sarray_find      (const struct sarray *a, ssize_t i);
void * sarray_find_with (const struct sarray *a, ssize_t i, const void *val0);
void * sarray_replace   (struct sarray *a, ssize_t i, void *src);
void   sarray_remove    (struct sarray *a, ssize_t i);

void * sarray_scatter (const struct sarray *a, void *dst);
void * sarray_gather  (struct sarray *a, const void *src);


/* iteration */
static inline void * sarray_begin (const struct sarray *a);
static inline void * sarray_end   (const struct sarray *a);
static inline void * sarray_ptr   (const struct sarray *a, ssize_t i);

static inline const ssize_t * sarray_keys_begin (const struct sarray *a);
static inline const ssize_t * sarray_keys_end   (const struct sarray *a);
static inline const ssize_t * sarray_keys_ptr   (const struct sarray *a, ssize_t i);


/* private functions */
struct sarray * _sarray_init (struct sarray *a, size_t elt_size);
void * _sarray_index      (struct sarray *a, ssize_t i);
void * _sarray_index_with (struct sarray *a, ssize_t i, const void *val0);


/* inline function definitions */
ssize_t sarray_size     (const struct sarray *a) { return darray_size(&a->values); }
bool    sarray_empty    (const struct sarray *a) { return darray_empty(&a->values); }
ssize_t sarray_max_size (const struct sarray *a) { return darray_max_size(&a->values); }
size_t  sarray_elt_size (const struct sarray *a) { return darray_elt_size(&a->values); }

void * sarray_begin (const struct sarray *a) { return darray_begin(&a->values); }
void * sarray_end   (const struct sarray *a) { return darray_end(&a->values); }
void * sarray_ptr   (const struct sarray *a, ssize_t i) { return darray_ptr(&a->values, i); }

const ssize_t * sarray_keys_begin (const struct sarray *a) { return darray_begin(&a->keys); }
const ssize_t * sarray_keys_end   (const struct sarray *a) { return darray_end(&a->keys); }
const ssize_t * sarray_keys_ptr   (const struct sarray *a, ssize_t i) { return darray_ptr(&a->keys, i); }



#endif /* _SARRAY_H */
