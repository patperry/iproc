#ifndef _DARRAY_H
#define _DARRAY_H

/* define SSIZE_MAX, ssize_t, and bool before including this file */

#include <stddef.h>    // sizeof, size_t
#include <string.h>    // memcpy


struct darray {
    void    *data;
    size_t   elt_size;    
    ssize_t  size;
    ssize_t  capacity;
};

/* create, destroy */
#define               darray_new(t)      (_darray_new(sizeof(t)))
struct darray *       darray_new_copy    (const struct darray *a);
void                  darray_free        (struct darray *a);

/* assignment, swap */
void                  darray_assign       (struct darray       *a,
                                           ssize_t              n,
                                           const void          *val);
void                  darray_assign_array (struct darray       *a,
                                           const void          *ptr,
                                           ssize_t              n);
void                  darray_copy         (const struct darray *a,
                                           struct darray       *dst);
void                  darray_swap         (struct darray       *a,
                                           struct darray       *b);
                                           


/* index */
#define               darray_index(a,t,i) (((t *)((a)->data))[(i)])
static inline void    darray_set (struct darray *a,
                                  ssize_t        i,
                                  const void    *val);

/* informative */
#define               darray_front(a,t)                        (darray_index(a,t,0))
#define               darray_back(a,t)                         (darray_index(a,t,(a)->size - 1))
static inline size_t  darray_elt_size (const struct darray *a);
static inline ssize_t darray_size     (const struct darray *a);
static inline bool    darray_empty    (const struct darray *a);
static inline ssize_t darray_capacity (const struct darray *a);
static inline ssize_t darray_max_size (const struct darray *a);

/* standard operations */
void                  darray_insert       (struct darray *a,
                                           ssize_t        i,
                                           const void    *val);
void                  darray_insert_many  (struct darray *a,
                                           ssize_t        i,
                                           ssize_t        n,
                                           const void    *val);
void                  darray_insert_array (struct darray *a,
                                           ssize_t        i,
                                           const void    *ptr,
                                           ssize_t        n);
void                  darray_push_back    (struct darray *a,
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
struct darray *       _darray_new (size_t elt_size);


/* inline function definitions */
size_t  darray_elt_size (const struct darray *a) { return a->elt_size; }
ssize_t darray_size     (const struct darray *a) { return a->size; }
bool    darray_empty    (const struct darray *a) { return a->size == 0; }
ssize_t darray_capacity (const struct darray *a) { return a->capacity; }
ssize_t darray_max_size (const struct darray *a) { return SSIZE_MAX / a->elt_size; }


void darray_set (struct darray *a, ssize_t i, const void *val)
{
    memcpy(darray_ptr(a, i), val, darray_elt_size(a));
}


void * darray_begin (const struct darray *a)
{
    return a->data;
}


void * darray_end (const struct darray *a)
{
    return darray_ptr(a, a->size);
}


void * darray_ptr (const struct darray *a, ssize_t i)
{
    return a->data + i * a->elt_size;

}


#endif /* _DARRAY_H */
