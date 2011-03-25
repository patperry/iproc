#ifndef _IPROC_COMPARE_H
#define _IPROC_COMPARE_H

#include <stddef.h>
#include <stdint.h>

typedef int (*iproc_compare_fn) (void *px, void *py);


static inline int
iproc_int64_compare (void *px, void *py)
{
    int64_t x = *(int64_t *)px;
    int64_t y = *(int64_t *)py;
    
    if (x < y) {
        return -1;
    } else if (x > y) {
        return +1;
    } else {
        return 0;
    }
}

static inline int
iproc_uint64_compare (void *px, void *py)
{
    uint64_t x = *(int64_t *)px;
    uint64_t y = *(int64_t *)py;
    
    if (x < y) {
        return -1;
    } else if (x > y) {
        return +1;
    } else {
        return 0;
    }
}

static inline int
iproc_size_compare (void *px, void *py)
{
    size_t x = *(size_t *)px;
    size_t y = *(size_t *)py;
    
    if (x < y) {
        return -1;
    } else if (x > y) {
        return +1;
    } else {
        return 0;
    }
}

/* compare using uint64_t instead of double to avoid dealing with NaN
 * comparisons; this relies on IEEE doubles being 8 bytes and lexicographically
 * ordered, and uint64_t having the same endianness as double */
#define iproc_double_compare iproc_uint64_compare


#endif /* _IPROC_COMPARE_H */
