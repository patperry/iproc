#ifndef _intmap_H
#define _intmap_H

/* 
 * Define SSIZE_MAX, ssize_t, bool and assert before including this file.
 */

#include <stddef.h>    // sizeof, size_t
#include "intset.h"
#include "darray.h"
#include "util.h"


struct intmap {
    struct intset keys;
    struct darray vals;
};

struct intmap_pos {
    struct intset_pos key;
};

/* create, destroy */
#define intmap_init(m,t) \
        _intmap_init(m,sizeof(t))

bool intmap_init_copy (struct intmap *m, const struct intmap *src);
void intmap_deinit    (struct intmap *m);


/* assign, copy, clear */
bool                  intmap_assign_copy  (struct intmap       *a,
                                           const struct intmap *src);
void *                intmap_copy_vals_to (const struct intmap *m,
                                           void                *dst);
intptr_t *            intmap_copy_keys_to (const struct intmap *m,
                                           intptr_t            *dst);
void                  intmap_clear        (struct intmap *m);
bool                  intmap_reserve      (struct intmap *m, ssize_t n);


/* informative */
static inline ssize_t  intmap_size      (const struct intmap *m);
static inline bool     intmap_empty     (const struct intmap *m);
static inline ssize_t  intmap_max_size  (const struct intmap *m);
static inline size_t   intmap_elt_size  (const struct intmap *m);
static inline intptr_t intmap_min_key   (const struct intmap *m);
static inline intptr_t intmap_max_key   (const struct intmap *m);
static inline intptr_t intmap_key_at    (const struct intmap *m, ssize_t i);
#define intmap_value_at(m,t,i) \
        darray_index((m)->values, t, i)

bool    intmap_contains    (const struct intmap *m, intptr_t key);
ssize_t intmap_index       (const struct intmap *m, intptr_t key);
void *  intmap_lookup      (const struct intmap *m, intptr_t key);
void *  intmap_lookup_with (const struct intmap *m, intptr_t key, const void *val0);


/* modification */
bool intmap_add        (struct intmap *m, intptr_t key, const void *val);
bool intmap_add_all    (struct intmap *m, const intptr_t *keys, const void *vals, ssize_t n);
void intmap_remove     (struct intmap *m, intptr_t key);
void intmap_remove_all (struct intmap *m, const intptr_t *keys, ssize_t n);

/* iteration */
static inline void * intset_vals_begin (const struct intmap *m);
static inline void * intset_vals_end   (const struct intmap *m);
static inline void * intset_vals_ptr   (const struct intmap *m, ssize_t i);

static inline const intptr_t * intset_keys_begin (const struct intmap *m);
static inline const intptr_t * intset_keys_end   (const struct intmap *m);
static inline const intptr_t * intset_keys_ptr   (const struct intmap *m, ssize_t i);


/* position-based operations */
void * intmap_find (const struct intmap *m, intptr_t key, struct intmap_pos *pos);
bool intmap_insert (      struct intmap *m, const struct intmap_pos *pos, const void *val);
void intmap_erase  (      struct intmap *m, const struct intmap_pos *pos);


/* private functions */
bool _intmap_init (struct intmap *m, size_t elt_size);


/* inline function definitions */
ssize_t intmap_size     (const struct intmap *m) { return intset_size(&m->keys); }
bool    intmap_empty    (const struct intmap *m) { return intset_empty(&m->keys); }
ssize_t intmap_max_size (const struct intmap *m) { return MIN(intset_max_size(&m->keys),
                                                              darray_max_size(&m->vals)); }
size_t  intmap_elt_size (const struct intmap *m) { return darray_elt_size(&m->vals); }
intptr_t intmap_min_key  (const struct intmap *m) { return intset_min(&m->keys); }
intptr_t intmap_max_key  (const struct intmap *m) { return intset_max(&m->keys); }
intptr_t intmap_key_at  (const struct intmap *m, ssize_t i) { return intset_at(&m->keys, i); }

void * intmap_vals_begin (const struct intmap *m) { return darray_begin(&m->vals); }
void * intmap_vals_end   (const struct intmap *m) { return darray_end(&m->vals); }
void * intmap_vals_ptr   (const struct intmap *m, ssize_t i) { return darray_ptr(&m->vals, i); }

const intptr_t * intmap_keys_begin (const struct intmap *m) { return intset_begin(&m->keys); }
const intptr_t * intmap_keys_end   (const struct intmap *m) { return intset_end(&m->keys); }
const intptr_t * intmap_keys_ptr   (const struct intmap *m, ssize_t i) { return intset_ptr(&m->keys, i); }


#endif /* _intmap_H */
