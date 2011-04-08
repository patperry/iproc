#ifndef _UTIL_H
#define _UTIL_H

#include <stddef.h>
#include "compare.h"

#define MAX(x,y) ((y) > (x) ? (y) : (x))
#define MIN(x,y) ((y) < (x) ? (y) : (x))

#define container_of(ptr, type, member) ({ \
            const typeof( ((type *)0)->member ) *__mptr = (ptr);    \
            (type *)( (char *)__mptr - offsetof(type,member) );})


ssize_t               find_index      (const void *begin,
                                       ssize_t     size,
                                       const void *key,
                                       compare_fn  compar,                                       
                                       size_t      elt_size);
ssize_t               find_last_index (const void *begin,
                                       ssize_t     size,
                                       const void *key,
                                       compare_fn  compar,                                       
                                       size_t      elt_size);
ssize_t               binary_search   (const void *begin,
                                       ssize_t     size,
                                       const void *key,
                                       compare_fn  compar,
                                       size_t      elt_size);                                     


#endif /* _UTIL_H */
