#ifndef _INTSET_H
#define _INTSET_H

#include <stddef.h>
#include "darray.h"

struct intset {
    struct darray values;
};

struct intset_pos {
    intptr_t value;
    ssize_t  index;
};


/* create, destroy */
struct intset * intset_init      (struct intset *s);
struct intset * intset_init_copy (struct intset *s, const struct intset *src);
void            intset_deinit    (struct intset *s);


/* assign, copy, clear */
struct intset * intset_assign        (struct intset *s, const intptr_t *ptr, ssize_t n);
struct intset * intset_assign_sorted (struct intset *s, const intptr_t *ptr, ssize_t n);
struct intset * intset_assign_copy   (struct intset *s, const struct intset *src);
intptr_t *      intset_copy_to       (const struct intset *s, intptr_t *dst);
void            intset_clear         (struct intset *s);


/* informative */
static inline bool     intset_empty    (const struct intset *s);
static inline ssize_t  intset_size     (const struct intset *s);
static inline ssize_t  intset_max_size (const struct intset *s);
static inline intptr_t intset_min      (const struct intset *s);
static inline intptr_t intset_max      (const struct intset *s);
static inline intptr_t intset_at       (const struct intset *s, ssize_t i);

bool     intset_contains   (const struct intset *s, intptr_t val);
ssize_t  intset_find_index (const struct intset *s, intptr_t val);


/* modification */
bool       intset_add        (struct intset *s, intptr_t val);
intptr_t * intset_add_all    (struct intset *s, intptr_t *ptr, ssize_t n);
void       intset_remove     (struct intset *s, intptr_t val);
void       intset_remove_all (struct intset *s, intptr_t *ptr, ssize_t n);


/* iteration */
static inline const intptr_t * intset_begin (const struct intset *s);
static inline const intptr_t * intset_end   (const struct intset *s);
static inline const intptr_t * intset_ptr   (const struct intset *s, ssize_t i);



/* position-based operations */
bool intset_find   (const struct intset *s, intptr_t val, struct intset_pos *pos);
bool intset_insert (      struct intset *s, const struct intset_pos *pos);
void intset_erase  (      struct intset *s, const struct intset_pos *pos);


/* inline function definitions */
bool    intset_empty    (const struct intset *s) { return darray_empty(&s->values); }
ssize_t intset_size     (const struct intset *s) { return darray_size(&s->values); }
ssize_t intset_max_size (const struct intset *s) { return darray_max_size(&s->values); }
intptr_t intset_min     (const struct intset *s) { return darray_front(&s->values, intptr_t); }
intptr_t intset_max     (const struct intset *s) { return darray_back(&s->values, intptr_t); }
intptr_t intset_at      (const struct intset *s, ssize_t i) { return darray_index(&s->values, intptr_t, i); }
const intptr_t * intset_begin (const struct intset *s)            { return darray_begin(&s->values); }
const intptr_t * intset_end   (const struct intset *s)            { return darray_end(&s->values); }
const intptr_t * intset_ptr   (const struct intset *s, ssize_t i) { return darray_ptr(&s->values, i); }



#endif /* _INTSET_H */