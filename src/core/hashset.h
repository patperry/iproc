#ifndef _HASHSET_H
#define _HASHSET_H

#include <stddef.h>

#include "compare.h"
#include "darray.h"
#include "hash.h"
#include "intmap.h"
#include "trashstack.h"
#include "util.h"


struct hashset {
    struct intmap buckets;
    hash_fn       hash;
    equals_fn     equals;
    size_t        elt_size;
    size_t        elt_offset;
    struct        trashstack nodes;
};

// a position in a hash set
struct hashset_pos {
    struct intmap_pos      key;
    struct hashset_bucket *bucket;
    union { struct hashset_node   *n;
            struct hashset_bucket *b; } node;
};

// a hash set iterator
struct hashset_it {
    struct hashset_bucket *bucket, *end;
    struct hashset_node   *node;
};


/* create, destroy */
#define hashset_init(set, type, hash, equal) \
        _hashset_init(set, hash, equal, sizeof(type), \
            offsetof(struct { struct node *next; type t; }, t)

// TODO: bool hashset_init_copy (struct hashset *s, const struct hashset *src);
void hashset_deinit (struct hashset *s);


/* assign, copy, clear */
// bool                  hashset_assign_copy (struct hashset       *s,
//                                           const struct hashset *src);
// void *                hashset_copy_to     (const struct hashset *s,
//                                            void                 *dst);
void                  hashset_clear       (struct hashset *s);
bool                  hashset_reserve     (struct hashset *s, ssize_t n);


/* informative */
              size_t   hashset_size     (const struct hashset *s);
static inline bool     hashset_empty    (const struct hashset *s);
static inline ssize_t  hashset_max_size (const struct hashset *s);
static inline size_t   hashset_elt_size (const struct hashset *s);
static inline uint32_t hashset_hash     (const struct hashset *s, const void *val);
static inline int      hashset_equals   (const struct hashset *s, const void *val1, const void *val2);

/* modification */
bool    hashset_contains    (const struct hashset *s, const void *key);
void *  hashset_lookup      (const struct intmap *m,  const void *key);
void *  hashset_lookup_with (const struct intmap *m,  const void *key, const void *val0);


/* position-based operations */
void * hashset_add        (struct hashset *s, const void *val);
bool   hashset_add_all    (struct hashset *s, const void *vals, ssize_t n);
void   hashset_remove     (struct hashset *s, const void *key);
void   hashset_remove_all (struct hashset *s, const void *keys, ssize_t n);


// return pointer to value if it exists
void * hashset_find (const struct hashset *s, const void *key, struct hashset_pos *pos);

// create storage location and return it; valid when hashset_find returns NULL
void * hashset_insert  (struct hashset *s, struct hashset_pos *pos, const void *val);

// if val is non-null, replace value at storage location and return it;
// otherwise, remove the storage location and return NULL.
// valid when hashset_find returns non-NULL
void * hashset_replace (struct hashset *s, struct hashset_pos *pos, const void *val);

// remove storage location; valid when hashset_find returns non-NULL
void   hashset_erase  (struct hashset *s, struct hashset_pos *pos);


/* iteration */
bool   hashset_it_init    (const struct hashset *s, struct hashset_it *it);
void   hashset_it_deinit  (const struct hashset *s, struct hashset_it *it);
void   hashset_it_reset   (const struct hashset *s, struct hashset_it *it);
bool   hashset_it_advance (const struct hashset *s, struct hashset_it *it);
void * hashset_it_current (const struct hashset *s, const struct hashset_it *it);



/* private functions */
bool _hashset_init (struct hashset *s, hash_fn hash, equals_fn equals,
                    size_t elt_size, size_t elt_offset);



/* inline function definitions */
bool    hashset_empty    (const struct hashset *s) { return intmap_empty(&s->buckets); }
size_t  hashset_elt_size (const struct hashset *s) { return s->elt_size; }
ssize_t hashset_max_size (const struct hashset *s) { return intmap_max_size(&s->buckets); }


static inline uint32_t hashset_hash (const struct hashset *s, const void *val)
{
    if (s->hash) {
        return s->hash(val);
    } else {
        return memory_hash(val, hashset_elt_size(s));
    }
}


static inline int hashset_equals (const struct hashset *s, const void *val1, const void *val2)
{
    if (s->equals) {
        return s->equals(val1, val2);
    } else {
        return (memcmp(val1, val2, hashset_elt_size(s)) == 0);
    }
}



#endif /* _HASHSET_H */
