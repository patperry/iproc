#ifndef _ARRAY_H
#define _ARRAY_H

#include <stddef.h>    // sizeof, size_t
#include <string.h>    // memcpy
#include "compare.h"


struct array {
    void    *data;
    ssize_t  size;
    ssize_t  elt_size;
    bool     owner;
};

/* create, destroy */
#define               array_init(a,t,n)           _array_init(a, n, sizeof(t))
#define               array_init_view(a,t,ptr,n)  _array_init_view(a,ptr,n,sizeof(t))
struct array *        array_init_copy         (struct array *a, const struct array *src);
struct array *        array_init_slice        (struct array *a,
                                               const struct array *parent,
                                               ssize_t             i,
                                               ssize_t             n);
void                  array_deinit            (struct array *a);
bool                  array_realloc           (struct array *a,
                                               ssize_t       n);




/* assignment, swap */
void *                array_assign            (struct array *a,
                                               ssize_t       n,
                                               const void   *val);
void *                array_assign_with       (struct array *a,
                                               ssize_t       n,
                                               const void   *val,
                                               const void   *val0);
void *                array_assign_array      (struct array *a,
                                               const void   *ptr,
                                               ssize_t       n);
void *                array_assign_array_with (struct array *a,
                                               const void   *ptr,
                                               ssize_t       n,
                                               const void   *val0);

void                  array_copy              (const struct array *a,
                                               struct array       *dst);
void                  array_copy_with         (const struct array *a,
                                               struct array       *dst,
                                               const void         *val0);

void                  array_swap              (struct array       *a,
                                               struct array       *b);

/* index */
#define               array_index(a,t,i)  (((t *)((a)->data))[(i)])
static inline void *  array_get (const struct array *a,
                                 ssize_t             i,
                                 void               *dst);
static inline void *  array_set (struct array *a,
                                 ssize_t       i,
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
static inline ssize_t array_max_size (const struct array *a);
static inline bool    array_owner    (const struct array *a);


/* iteration */
static inline void *  array_begin (const struct array *a);
static inline void *  array_end   (const struct array *a);
static inline void *  array_ptr   (const struct array *a,
                                   ssize_t             i);


/* searching */
ssize_t               array_lfind   (const struct array *a,
                                     const void         *key,
                                     compare_fn          compar);
ssize_t               array_bsearch (const struct array *a,
                                     const void         *key,
                                     compare_fn          compar);


/* private functions */
struct array *        _array_new       (ssize_t size,
                                        ssize_t elt_size);
struct array *        _array_init      (struct array *a,
                                        ssize_t       size,
                                        ssize_t       elt_size);
struct array *        _array_init_view (struct array *a,
                                        const void   *data,
                                        ssize_t       size,
                                        ssize_t       elt_size);

static inline ssize_t _array_elt_size (const struct array *a);


/* inline function defs */
ssize_t _array_elt_size (const struct array *a) { return a->elt_size; }
ssize_t array_size      (const struct array *a) { return a->size; }
bool    array_empty     (const struct array *a) { return a->size == 0; }
ssize_t array_max_size  (const struct array *a) { return SSIZE_MAX / a->elt_size; }
bool    array_owner     (const struct array *a) { return a->owner; }



void * array_get (const struct array *a, ssize_t i, void *dst)
{
    memcpy(dst, array_ptr(a, i), _array_elt_size(a));
    return (void *)dst + _array_elt_size(a);
}


void * array_set (struct array *a, ssize_t i, const void *src)
{
    memcpy(array_ptr(a, i), src, _array_elt_size(a));
    return (void *)src + _array_elt_size(a);
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
    return a->data + i * a->elt_size;
    
}


#endif /* _ARRAY_H */
