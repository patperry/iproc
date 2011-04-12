#ifndef _REFCOUNT_H
#define _REFCOUNT_H

#include <limits.h>
#include "util.h"

#define REFCOUNT_MAX INT_MAX

struct refcount {
    int count;
};

bool refcount_init (struct refcount *refcount);
void refcount_set  (struct refcount *refcount,
                    int              count);
bool refcount_get  (struct refcount *refcount);
bool refcount_put  (struct refcount *refcount,
                    void (*release) (struct refcount *refcount));


#endif /* _REFCOUNT_H */
