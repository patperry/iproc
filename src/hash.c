
#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <assert.h>
#include <errno.h>
#include <float.h>
#include <limits.h>
#include <math.h>

#include "hash.h"


/* C version of the Boost hash_combine function from
 * <boost/functional/hash.hpp>
 */
size_t
iproc_hash_combine (size_t seed, size_t hash_value)
{
    seed ^= hash_value + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    return seed;
}

/* Ported from boost/functional/detail/hash_float.hpp
 *   C++ version written by Daniel James (Boost License 1.0)
 */

static size_t
hash_float_combine (size_t seed, size_t value)
{
    seed ^= value + (seed<<6) + (seed>>2);
    return seed;
}

static size_t
double_hash_impl (double v)
{
    int exp = 0;
    errno = 0;
    v = frexp(v, &exp);
    if (errno) return 0;

    size_t seed = 0;

    int int_digits = sizeof(int) * CHAR_BIT - 1;
    size_t length = (DBL_MANT_DIG + int_digits - 1) / int_digits;
    size_t i;

    for (i = 0; i < length; i++) {
        v = ldexp(v, int_digits);
        int part = (int)v;
        v -= part;
        seed = hash_float_combine(seed, part);
    }
 
    seed = hash_float_combine(seed, exp);

    return seed;
}

size_t
iproc_hash_double (double v)
{
    switch (fpclassify(v)) {
    case FP_ZERO:
        return 0;
    case FP_INFINITE:
        return (size_t)(v > 0 ? -1 : -2);
    case FP_NAN:
        return (size_t)(-3);
    case FP_NORMAL:
    case FP_SUBNORMAL:
        return double_hash_impl(v);
    default:
        assert(0);
        return 0;
    }
}
