#ifndef _IPROC_COMPARE_H
#define _IPROC_COMPARE_H

#include <stddef.h>
#include <stdint.h>
#include <sys/types.h>

typedef int (*iproc_compare_fn) (const void *px, const void *py);


static inline int
iproc_int64_compare (const void *px, const void *py)
{
    int64_t x = *(const int64_t *)px;
    int64_t y = *(const int64_t *)py;
    
    if (x < y) {
        return -1;
    } else if (x > y) {
        return +1;
    } else {
        return 0;
    }
}

static inline int
iproc_uint64_compare (const void *px, const void *py)
{
    uint64_t x = *(const int64_t *)px;
    uint64_t y = *(const int64_t *)py;
    
    if (x < y) {
        return -1;
    } else if (x > y) {
        return +1;
    } else {
        return 0;
    }
}

static inline int
iproc_size_compare (const void *px, const void *py)
{
    size_t x = *(const size_t *)px;
    size_t y = *(const size_t *)py;
    
    if (x < y) {
        return -1;
    } else if (x > y) {
        return +1;
    } else {
        return 0;
    }
}

static inline int
iproc_ssize_compare (const void *px, const void *py)
{
    ssize_t x = *(const ssize_t *)px;
    ssize_t y = *(const ssize_t *)py;
    
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
