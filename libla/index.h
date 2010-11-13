#ifndef _LA_INDEX_H
#define _LA_INDEX_H

#include <limits.h>

typedef int          la_index;
#define LA_INDEX_MAX INT_MAX
#define LA_INDEX_MIN INT_MIN
#define LA_INDEX_FMT "%d"

typedef la_index     la_size;
#define LA_SIZE_MAX  LA_INDEX_MAX
#define LA_SIZE_MIN  LA_INDEX_MIN
#define LA_SIZE_FMT  LA_INDEX_FMT

#endif /* _LA_INDEX_H */
