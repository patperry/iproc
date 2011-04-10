#include "config.h"

#ifndef HAVE_IMAXABS

#ifdef HAVE_INTTYPES_H
# include <inttypes.h>
#endif

#include <stdlib.h>

intmax_t imaxabs (intmax_t i)
{
#if SIZEOF_INTMAX_T == SIZEOF_INT
    return abs(i);
#elif SIZEOF_INTMAX_T == SIZEOF_LONG
    return labs(i);
#elif SIZEOF_INTMAX_T == SIZEOF_LONG_LONG
    return llabs(i);
#else
# error "Need definition for imaxabs; no support for sizeof(intmax_t) > sizeof(long long)"
#endif
}

#endif /* HAVE_IMAXABS */
