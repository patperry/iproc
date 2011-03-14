#ifndef _IPROC_HASH_H
#define _IPROC_HASH_H

#include <stddef.h>

size_t iproc_hash_double  (double v);

size_t iproc_hash_combine (size_t seed,
                           size_t hash_value);


#endif /* _IPROC_HASH_H */
