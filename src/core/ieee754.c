#ifdef HAVE_CONFIG_H
# include <config.h>
#endif

#include <float.h>
#include <math.h>
#include <sys/types.h>
#include <sys/param.h>
#include "ieee754.h"


#define DBL_EXPMASK            ((uint16_t) 0x7FF0)
#define DBL_SIGNMASK           ((uint16_t) 0x8000)
#define DBL_EXPBIAS            ((uint16_t) 0x3FE0)
#define DBL_EXPBIAS_INT32      ((uint32_t) 0x7FF00000)
#define DBL_MANTISSAMASK_INT32 ((uint32_t) 0x000FFFFF) /* for the MSB only */
#if BYTE_ORDER == BIG_ENDIAN
# define DBL_EXPPOS_INT16 0
# define DBL_SIGNPOS_BYTE 0
#else
# define DBL_EXPPOS_INT16 3
# define DBL_SIGNPOS_BYTE 7
#endif

union double_uint64 {
    double d;
    uint64_t w;
};

int
iproc_identical (double x,
                 double y)
{
    union double_uint64 ux = { x };
    union double_uint64 uy = { y };
    return ux.w == uy.w;
}

/* ported from tango/math/IEEE.d */
double
iproc_nextup (double x)
{
    union double_uint64 ps = { x };

    if ((ps.w & UINT64_C(0x7FF0000000000000)) == UINT64_C(0x7FF0000000000000)) {
        /* First, deal with NANs and infinity */
        if (x == -INFINITY) return -DBL_MAX;
        return x; // +INF and NAN are unchanged.
    }
    if (ps.w & UINT64_C(0x8000000000000000))  { /* Negative number */
        if (ps.w == UINT64_C(0x8000000000000000)) { /* it was negative zero */
            ps.w = UINT64_C(0x0000000000000001); /* change to smallest subnormal */
            return ps.d;
        }
        --ps.w;
    } else { /* Positive number */
        ++ps.w;
    }
    return ps.d;
}

/* ported from tango/math/IEEE.d */
double
iproc_nextdown (double x)
{
    return -iproc_nextup(-x);
}

/* ported from tango/math/IEEE.d */
double
iproc_ieeemean (double x,
                double y)
{
    if (!((x>=0 && y>=0) || (x<=0 && y<=0))) return NAN;

    union double_uint64 ul;
    union double_uint64 xl = { x };
    union double_uint64 yl = { y };
    uint64_t m = ( (xl.w & UINT64_C(0x7FFFFFFFFFFFFFFF))
                 + (yl.w & UINT64_C(0x7FFFFFFFFFFFFFFF)) ) >> 1;
    m |= (xl.w & UINT64_C(0x8000000000000000));
    ul.w = m;

    return ul.d;
}

double
iproc_mknan (uint64_t payload)
{
    union double_uint64 x = { NAN };

    /* get sign, exponent, and quiet bit from NAN */
    x.w &= UINT64_C(0xFFF8000000000000);

    /* ignore sign, exponent, and quiet bit in payload */
    payload &= UINT64_C(0x0007FFFFFFFFFFFF);
    x.w |= payload;

    return x.d;
}

uint64_t
iproc_getnan (double x)
{
    union double_uint64 payload = { x };

    /* clear sign, exponent, and quiet bit */
    payload.w &= UINT64_C(0x0007FFFFFFFFFFFF);
    return payload.w;
}

int
iproc_feqrel (double x,
              double y)
{
    /* Public Domain. Original Author: Don Clugston, 18 Aug 2005.
     * Ported to C by Patrick Perry, 26 Feb 2010.
     */
    if (x == y) return DBL_MANT_DIG; /* ensure diff!= 0, cope with INF. */

    double diff = fabs(x - y);

    union { double r; uint16_t w[sizeof(double)/2]; } pa = { x };
    union { double r; uint16_t w[sizeof(double)/2]; } pb = { y };
    union { double r; uint16_t w[sizeof(double)/2]; } pd = { diff };

    /* The difference in abs(exponent) between x or y and abs(x-y)
     * is equal to the number of significand bits of x which are
     * equal to y. If negative, x and y have different exponents.
     * If positive, x and y are equal to 'bitsdiff' bits.
     * AND with 0x7FFF to form the absolute value.
     * To avoid out-by-1 errors, we subtract 1 so it rounds down
     * if the exponents were different. This means 'bitsdiff' is
     * always 1 lower than we want, except that if bitsdiff==0,
     * they could have 0 or 1 bits in common.
     */
    int bitsdiff = (( ((pa.w[DBL_EXPPOS_INT16] & DBL_EXPMASK)
                     + (pb.w[DBL_EXPPOS_INT16] & DBL_EXPMASK)
                     - ((uint16_t) 0x8000 - DBL_EXPMASK)) >> 1)
                   - (pd.w[DBL_EXPPOS_INT16] & DBL_EXPMASK)) >> 4;

    if ((pd.w[DBL_EXPPOS_INT16] & DBL_EXPMASK) == 0) {
        /* Difference is denormal
         * For denormals, we need to add the number of zeros that
         * lie at the start of diff's significand.
         * We do this by multiplying by 2^DBL_MANT_DIG
         */
        pd.r *= (1 / DBL_EPSILON);

        return (bitsdiff + DBL_MANT_DIG
                - (pd.w[DBL_EXPPOS_INT16] >> 4));
    }

    if (bitsdiff > 0)
        return bitsdiff + 1; /* add the 1 we subtracted before */

    /* Avoid out-by-1 errors when factor is almost 2. */
    return (bitsdiff == 0
            && !((pa.w[DBL_EXPPOS_INT16]
                  ^ pb.w[DBL_EXPPOS_INT16]) & DBL_EXPMASK)) ? 1 : 0;
}
