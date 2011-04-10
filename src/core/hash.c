#include "port.h"

#include <assert.h>
#include <errno.h>
#include <float.h>
#include <limits.h>
#include <math.h>

#include "hash.h"

/* From http://www.azillionmonkeys.com/qed/hash.html */
#undef get16bits
#if (defined(__GNUC__) && defined(__i386__)) || defined(__WATCOMC__) \
|| defined(_MSC_VER) || defined (__BORLANDC__) || defined (__TURBOC__)
#define get16bits(d) (*((const uint16_t *) (d)))
#endif

#if !defined (get16bits)
#define get16bits(d) ((((uint32_t)(((const uint8_t *)(d))[1])) << 8)\
+(uint32_t)(((const uint8_t *)(d))[0]) )
#endif

uint32_t memory_hash (const void *ptr, ssize_t len)
{
    const char *data = ptr;
    uint32_t hash = (uint32_t)len;
    uint32_t tmp;
    int rem;
    
    if (len <= 0 || data == NULL) return 0;
    
    rem = len & 3;
    len >>= 2;
    
    /* Main loop */
    for (;len > 0; len--) {
        hash  += get16bits (data);
        tmp    = (get16bits (data+2) << 11) ^ hash;
        hash   = (hash << 16) ^ tmp;
        data  += 2*sizeof (uint16_t);
        hash  += hash >> 11;
    }
    
    /* Handle end cases */
    switch (rem) {
        case 3: hash += get16bits (data);
            hash ^= hash << 16;
            hash ^= data[sizeof (uint16_t)] << 18;
            hash += hash >> 11;
            break;
        case 2: hash += get16bits (data);
            hash ^= hash << 11;
            hash += hash >> 17;
            break;
        case 1: hash += *data;
            hash ^= hash << 10;
            hash += hash >> 1;
    }
    
    /* Force "avalanching" of final 127 bits */
    hash_finalize(&hash);
    
    return hash;
}


/* Ported from boost/functional/detail/hash_float.hpp
 *   C++ version written by Daniel James (Boost License 1.0)
 */
static void hash_float_combine (uint32_t *seedp, uint32_t value)
{
    uint32_t seed = *seedp;
    seed ^= value + (seed<<6) + (seed>>2);
    *seedp = seed;
}


static uint32_t double_hash_impl (double v)
{
    int exp = 0;
    errno = 0;
    v = frexp(v, &exp);
    if (errno) return 0;

    uint32_t seed = 0;

    int int_digits = sizeof(int) * CHAR_BIT - 1;
    size_t length = (DBL_MANT_DIG + int_digits - 1) / int_digits;
    size_t i;

    for (i = 0; i < length; i++) {
        v = ldexp(v, int_digits);
        int part = (int)v;
        v -= part;
        hash_float_combine(&seed, part);
    }
 
    hash_float_combine(&seed, exp);

    return seed;
}


uint32_t double_hash (const void *val)
{
    double v = *(double *)val;
    
    switch (fpclassify(v)) {
    case FP_ZERO:
        return 0;
    case FP_INFINITE:
        return (uint32_t)(v > 0 ? -1 : -2);
    case FP_NAN:
        return (uint32_t)(-3);
    case FP_NORMAL:
    case FP_SUBNORMAL:
        return double_hash_impl(v);
    default:
        assert(0);
        return 0;
    }
}
