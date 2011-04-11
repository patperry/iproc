#ifndef _ARRAY_H
#define _ARRAY_H

/* Generic Array type.
 *
 * Define SSIZE_MAX, ssize_t, bool and assert before including this file.
 */

#include <stddef.h>    // sizeof, size_t
#include <string.h>    // memcpy
#include "compare.h"


struct array {
    void    *data;
    ssize_t  size;
    size_t   elt_size;
    bool     owner;
};

/* create, destroy */
#define               array_init(a,t,n)           _array_init(a, n, sizeof(t))
#define               array_init_view(a,t,ptr,n)  _array_init_view(a,ptr,n,sizeof(t))
struct array *        array_init_slice        (struct array       *a,
                                               const struct array *parent,
                                               ssize_t             i,
                                               ssize_t             n);
struct array *        array_init_copy         (struct array *a,
                                               const struct array *src);
void                  array_deinit            (struct array *a);


/* assign, copy, fill */
struct array *        array_assign_copy (struct array       *a,
                                         const struct array *src);
void *                array_copy_to     (const struct array *a,
                                         void               *dst);
void                  array_fill        (struct array *a,
                                         const void   *val);
void                  array_fill_range  (struct array *a,
                                         ssize_t       i,
                                         ssize_t       n,
                                         const void   *src);

/* index */
#ifdef NDEBUG
# define              array_index(a,t,i)  (((t *)((a)->data))[(i)])
#else
# define              array_index(a,t,i)  (*((t *)array_ptr(a, i)))
#endif
static inline void    array_get        (const struct array *a,
                                        ssize_t             i,
                                        void               *dst);
static inline void *  array_get_range  (const struct array *a,
                                        ssize_t             i,
                                        ssize_t             n,
                                        void               *dst);
static inline void    array_set        (struct array *a,
                                        ssize_t       i,
                                        const void   *src);
static inline void *  array_set_range  (struct array *a,
                                        ssize_t       i,
                                        ssize_t       n,
                                        const void   *src);


/* informative */
#define               array_front(a,t)       (array_index(a, t, 0))
#define               array_get_front(a,dst) (array_get(a, 0, dst))
#define               array_set_front(a,src) (array_set(a, 0, src))

#define               array_back(a,t)        (array_index(a, t, (a)->size - 1))
#define               array_get_back(a,dst)  (array_get(a, (a)->size - 1, dst))
#define               array_set_back(a,src)  (array_set(a, (a)->size - 1, src))

static inline ssize_t array_size     (const struct array *a);
static inline bool    array_empty    (const struct array *a);
static inline size_t  array_elt_size (const struct array *a);
static inline ssize_t array_max_size (const struct array *a);
static inline bool    array_owner    (const struct array *a);
static inline bool    array_overlaps (const struct array *a,
                                      ssize_t             i,
                                      ssize_t             n,
                                      const void         *ptr,
                                      ssize_t             nel);

/* operations */
void                  array_swap    (struct array *a, ssize_t i, ssize_t j);
void                  array_reverse (struct array *a);


/* iteration */
static inline void *  array_begin (const struct array *a);
static inline void *  array_end   (const struct array *a);
static inline void *  array_ptr   (const struct array *a,
                                   ssize_t             i);


/* searching, sorting */
bool                  array_contains        (const struct array *a,
                                             const void         *key,
                                             equal_fn            equal);
void *                array_find            (const struct array *a,
                                             const void         *key,
                                             equal_fn            equal);
ssize_t               array_find_index      (const struct array *a,
                                             const void         *key,
                                             equal_fn            equal);
void *                array_find_last       (const struct array *a,
                                             const void         *key,
                                             equal_fn            equal);
ssize_t               array_find_last_index (const struct array *a,
                                             const void         *key,
                                             equal_fn            equal);
ssize_t               array_binary_search   (const struct array *a,
                                             const void         *key,
                                             compare_fn          compar);


/* private functions */
struct array *        _array_init      (struct array *a,
                                        ssize_t       size,
                                        size_t        elt_size);
struct array *        _array_init_view (struct array *a,
                                        const void   *data,
                                        ssize_t       size,
                                        size_t        elt_size);
struct array *        _array_reinit    (struct array *a,
                                        ssize_t       n,
                                        size_t        elt_size);




/* inline function defs */
ssize_t array_size      (const struct array *a) { return a->size; }
bool    array_empty     (const struct array *a) { return a->size == 0; }
size_t  array_elt_size  (const struct array *a) { return a->elt_size; }
ssize_t array_max_size  (const struct array *a) { return SSIZE_MAX / a->elt_size; }
bool    array_owner     (const struct array *a) { return a->owner; }


bool array_overlaps (const struct array *a, ssize_t i, ssize_t n,
                     const void *ptr, ssize_t nel)
{
    assert(a);
    assert(n >= 0);
    assert(0 <= i && i <= array_size(a) - n);
    assert(0 <= nel && nel < array_max_size(a));

    const void *begin1 = array_ptr(a, i);
    const void *end1 = array_ptr(a, i + n);
    const void *begin2 = ptr;
    const void *end2 = begin2 + nel * array_elt_size(a);
    
    return ((begin1 <= begin2 && begin2 < end1)
            || (begin2 <= begin1 && begin1 < end2));
}

void array_get (const struct array *a, ssize_t i, void *dst)
{
    assert(0 <= i && i < array_size(a));
    assert(!array_overlaps(a, i, 1, dst, 1));

    array_get_range(a, i, 1, dst);
}


void * array_get_range (const struct array *a, ssize_t i, ssize_t n, void *dst)
{
    assert(n >= 0);
    assert(0 <= i && i <= array_size(a) - n);
    assert(!array_overlaps(a, i, n, dst, n));

    size_t nbytes = n * array_elt_size(a);
    memcpy(dst, array_ptr(a, i), nbytes);
    return (void *)dst + nbytes;
}

void array_set (struct array *a, ssize_t i, const void *src)
{
    assert(0 <= i && i < array_size(a));
    assert(!array_overlaps(a, i, 1, src, 1));

    array_set_range(a, i, 1, src);
}


void * array_set_range (struct array *a, ssize_t i, ssize_t n, const void *src)
{
    assert(n >= 0);
    assert(0 <= i && i <= array_size(a) - n);
    assert(!array_overlaps(a, i, n, src, n));

    size_t nbytes = n * array_elt_size(a);
    memcpy(array_ptr(a, i), src, nbytes);
    return (void *)src + nbytes;
}


void * array_begin (const struct array *a)
{
    return a->data;
}


void * array_end (const struct array *a)
{
    return array_ptr(a, a->size);
}


void * array_ptr (const struct array *a, ssize_t i)
{
    assert(0 <= i && i <= array_size(a));
    return a->data + i * a->elt_size;
    
}


#endif /* _ARRAY_H */
