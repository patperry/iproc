#ifndef _HASHSET_H
#define _HASHSET_H

#include <stddef.h>

#include "compare.h"
#include "darray.h"
#include "hash.h"
#include "sarray.h"
#include "util.h"

struct hashset_bucket {
    struct darray values;
};

struct hashset {
    struct sarray         buckets;
    struct hashset_bucket empty_bucket;
    hash_fn               hash;
    equal_fn              equal;
    destroy_fn            destroy;
};


#define hashset_init(s, type, hash, equal) \
        _hashset_init(s, hash, equal, sizeof(type))

#define hashset_init_with_destroy(s, type, hash, equal, destroy) \
        _hashset_init_with_destroy(s, hash, equal, destroy, sizeof(type))


void hashset_deinit (struct hashset *s);



static inline bool    hashset_empty    (const struct hashset *s);
static inline ssize_t hashset_max_size (const struct hashset *s);
static inline size_t  hashset_elt_size (const struct hashset *s);

ssize_t hashset_size (const struct hashset *s);


const void * hashset_find        (const struct hashset *s, const void *key);
const void * hashset_find_with   (const struct hashset *s, const void *key, const void *val0);
void *       hashset_add         (struct hashset *s, const void *val);
void *       hashset_add_all     (struct hashset *s, const void *ptr, ssize_t n);
void         hashset_remove      (struct hashset *s, const void *key);
void         hashset_remove_all  (struct hashset *s, const void *ptr, ssize_t n);
void         hashset_clear       (struct hashset *s);



/* private functions */
struct hashset * _hashset_init              (struct hashset *s,
                                             hash_fn         hash,
                                             equal_fn        equal,
                                             size_t          elt_size);

struct hashset * _hashset_init_with_destroy (struct hashset *s,
                                             hash_fn         hash,
                                             equal_fn        equal,
                                             destroy_fn      destroy,
                                             size_t          elt_size);

/* inline function definitions */
bool    hashset_empty    (const struct hashset *s) { return sarray_empty(&s->buckets); }
size_t  hashset_elt_size (const struct hashset *s) { return darray_elt_size(&(s->empty_bucket).values); }

ssize_t hashset_max_size (const struct hashset *s)
{
    return MIN(sarray_max_size(&s->buckets),
               darray_max_size(&(s->empty_bucket).values));
}


#endif /* _HASHSET_H */
