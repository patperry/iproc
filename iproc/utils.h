#ifndef _IPROC_UTILS_H
#define _IPROC_UTILS_H

#include <stddef.h>

#define IPROC_MAX(x,y) ((y) > (x) ? (y) : (x))
#define IPROC_MIN(x,y) ((y) < (x) ? (y) : (x))

#define container_of(ptr, type, member) ({ \
            const typeof( ((type *)0)->member ) *__mptr = (ptr);    \
            (type *)( (char *)__mptr - offsetof(type,member) );})

#endif /* _IPROC_UTILS_H */
