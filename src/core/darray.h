#ifndef _DARRAY_H
#define _DARRAY_H

/* define SSIZE_MAX, ssize_t, and bool before including this file */

#include <stddef.h>    // sizeof, size_t
#include <string.h>    // memcpy
#include "array.h"


struct darray {
    struct array array;
    ssize_t      size;
};

/* create, destroy */
#define               darray_init(a,t)    (_darray_init(a,sizeof(t)))
struct darray *       darray_init_copy    (struct darray       *a,
                                           const struct darray *src);
void                  darray_deinit       (struct darray *a);


/* assignment, swap */
void *                darray_assign       (struct darray       *a,
                                           ssize_t              n,
                                           const void          *val);
void *                darray_assign_array (struct darray       *a,
                                           const void          *ptr,
                                           ssize_t              n);
void                  darray_copy         (const struct darray *a,
                                           struct darray       *dst);
void                  darray_swap         (struct darray       *a,
                                           struct darray       *b);
                                           


/* index */
#define               darray_index(a,t,i) array_index((&(a)->array),t,i)
static        void *  darray_get (const struct darray *a,
                                  ssize_t              i,
                                  void                *dst);
static inline void *  darray_set (struct darray       *a,
                                  ssize_t              i,
                                  const void          *src);

/* informative */
#define               darray_front(a,t)                        (darray_index(a,t,0))
#define               darray_back(a,t)                         (darray_index(a,t,(a)->size - 1))

static inline ssize_t darray_size     (const struct darray *a);
static inline bool    darray_empty    (const struct darray *a);
static inline ssize_t darray_capacity (const struct darray *a);
static inline ssize_t darray_max_size (const struct darray *a);

/* standard operations */
void *                darray_insert       (struct darray *a,
                                           ssize_t        i,
                                           const void    *val);
void *                darray_insert_many  (struct darray *a,
                                           ssize_t        i,
                                           ssize_t        n,
                                           const void    *val);
void *                darray_insert_array (struct darray *a,
                                           ssize_t        i,
                                           const void    *ptr,
                                           ssize_t        n);
void *                darray_push_back    (struct darray *a,
                                           const void    *val);
void                  darray_erase        (struct darray *a,
                                           ssize_t        i);
void                  darray_erase_range  (struct darray *a,
                                           ssize_t        i,
                                           ssize_t        n);
void                  darray_pop_back     (struct darray *a);
void                  darray_clear        (struct darray *a);


/* allocation/size modification */
void                  darray_reserve     (struct darray *a,
                                          ssize_t        n);
void                  darray_resize      (struct darray *a,
                                          ssize_t        n);
void                  darray_resize_with (struct darray *a,
                                          ssize_t        n,
                                          const void    *val);

/* iteration */
static inline void *  darray_begin (const struct darray *a);
static inline void *  darray_end   (const struct darray *a);
static inline void *  darray_ptr   (const struct darray *a,
                                    ssize_t              i);

/* searching */
ssize_t               darray_lfind   (const struct darray *a,
                                      const void          *key,
                                      int        (*compar) (const void *, const void *));
ssize_t               darray_bsearch (const struct darray *a,
                                      const void           *key,
                                      int        (*compar) (const void *, const void *));



/* private functions */
struct darray *       _darray_init     (struct darray *a,
                                        size_t         elt_size);
static inline ssize_t _darray_elt_size (const struct darray *a);


/* inline function definitions */
ssize_t darray_size      (const struct darray *a) { return a->size; }
bool    darray_empty     (const struct darray *a) { return a->size == 0; }
ssize_t darray_capacity  (const struct darray *a) { return array_size(&a->array); }
ssize_t darray_max_size  (const struct darray *a) { return array_max_size(&a->array); }
ssize_t _darray_elt_size (const struct darray *a) { return _array_elt_size(&a->array); }


void * darray_get (const struct darray *a, ssize_t i, void *dst) { return array_get(&a->array, i, dst); }
void * darray_set (struct darray *a, ssize_t i, const void *src) { return array_set(&a->array, i, src); }


void * darray_begin (const struct darray *a)            { return array_begin(&a->array); }
void * darray_end   (const struct darray *a)            { return array_ptr(&a->array, a->size); }
void * darray_ptr   (const struct darray *a, ssize_t i) { return array_ptr(&a->array, i); }


#endif /* _DARRAY_H */
